function T = Step4_Parameter_uncertainty(T)
%%
%  File: mf_epid_tracking_SNMPC_init.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2021. November 16. (2021b)
% 
% Forked from:
% 
%  File: epid_tracking_SNMPC_LTV.m
%  Author: Peter Polcz (ppolcz@gmail.com)
%  Reviewed: 2021. October 12. (2021b)
%

%%

import casadi.*

xls = "Parameters/Par_HUN_Reconstruction.xlsx";

% Uncertainty of each parameter
U = readtable(xls,"Sheet","Uncertainty");
u = U.Variables;
T.ParamStd = u .* T.Param / 200;


He = @(A) A + A';

%%%
% Dynamical model

[f,h,hvar,J] = epid_ode_model_8comp;

[x,p,beta,nu,omega] = deal(f.args{:});

% Symbolic (to test)
%{
    [x,p,beta,nu,omega] = deal(f.Input_1,f.Input_2,f.Input_3,f.Input_4,f.Input_5);
    f.val = f.fh(x,p,beta,nu,omega);
%}

% New arguments
u = [beta;omega];
v = nu;

J.nu = numel(u);

A = jacobian(f.val,x);
B = jacobian(f.val,u);
E = jacobian(f.val,p);

Fn_A = Function('Fn_A',{x,u,p,v},{A},...
    {'x','u','p','V'},'A');
Fn_B = Function('Fn_B',{x,u,p,v},{B},...
    {'x','u','p','V'},'Bk');
Fn_E = Function('Fn_E',{x,u,p,v},{E},...
    {'x','u','p','V'},'E');

Sx0 = blkdiag(7,eye(7));

% Symbolic placeholder variables: mu_u(k), mu_x(k), Sigma_x(k), Sigma_xp(k)
Mu = SX.sym('Mu_free',size(u));
Mx = SX.sym('Mx_free',size(x));
[Sx,Sx_half,Fn_Sx_vec] = Pcz_CasADi_Helper.create_sym('Sx',J.nx,'str','sym','r3','f_vec');
[Sxp,Fn_Sxp_vec] = Pcz_CasADi_Helper.create_sym('Sxp',[J.nx J.np],'str','full','r2','f_vec');

% Symbolic placeholder variables: mu_u(k+1), Sigma_xp(k+1)
Mu_pp = SX.sym('mu_u_kp1',size(u));
[Sx_pp,Sx_pp_half] = Pcz_CasADi_Helper.create_sym('Sigma_x_pp',J.nx,'str','sym');

% Symbolic placeholder variables for K(k) and K(k+1)
K = SX.sym('K_k',size(B'));
K_pp = SX.sym('K_kp1',size(B'));

%%%
% Recursive equation for the next expected state vector.
Fn_Mx_pp = Function('Fn_Mx_pp',...
    {x,u,p,v},...
    {f.val},...
    {'Mx','Mu','Mp','v'},...
    {'Mx_pp'});

%%%
% Recursive equation for the covariance between the next state and parameter vectors.
Sp = diag(SX.sym('Sp',J.np,1));
Fn_Sxp_pp = Function('Fn_Sxp_pp',...
    {x,u,p,v,Sxp,Sp,K},...
    {(A-B*K)*Sxp + E*Sp},...
    {'Mx','Mu','Mp','v','Sxp','Sp','K'},...
    {'Sxp_pp'});

%%%
% Recursive equation for the variance of the next state vector.
Fn_Sx_pp = Function('Fn_Sx_pp',...
    {x,u,p,v,Sxp,Sx_half,Sp,K},...
    {(A-B*K)*Sx*(A-B*K)' + E*Sp*E' + He( (A-B*K)*Sxp*E' )},...
    {'Mx','Mu','Mp','v','Sxp','Sx','Sp','K'},...
    {'Sx_pp'});

% 2021.12.14. (december 14, kedd), 16:08
Sxup_fh = @(Sx,K,Sxp,Sp) ...
    blkdiag([eye(J.nx) -K'],eye(J.np))' * [ 
        Sx   Sxp
        Sxp' Sp
    ] * blkdiag([eye(J.nx) -K'],eye(J.np));
% Sxup_fh = @(Sx,K,Sxp) [
%      Sx   -Sx*K'    Sxp
%     -K*Sx  K*Sx*K' -K*Sxp
%     Sxp'  -Sxp*K    Sp
%     ];

%%%
% Compute variances

% Controllable_modes
Ctrl_idx = [J.S,J.L,J.P,J.I,J.A,J.H,J.R];

N = height(T);

Mx_pp_error = zeros(N,1);
K_val   = cell(1,N);
Sx_val  = [ Sx0              cell(1,N-1) ];
Sxp_val = [ zeros(J.nx,J.np) cell(1,N-1) ];
Su_val  = cell(1,N);
Sy_val  = cell(1,N);

T.TrRate(end) = T.TrRate(end-1);
T.ImLossRate(end) = T.ImLossRate(end-1);
T.Rc(end) = T.Rc(end-1);
T.New_Cases(end) = T.New_Cases(end-1);
T.Infected(end) = T.Infected(end-1);
T.Infectious(end) = T.Infectious(end-1);

State = T(:,["S","L","P","I","A","H","D","R"]);
T.x = State.Variables;
T.u = [T.TrRate , T.ImLossRate];

R = eye(J.nu);

status = PStatus(N,'Compute variances for initial guess');
for k = 1:N-1
    Mx_val = T.x(k,:)';
    Mu_val = T.u(k,:)';
    Mp_val = T.Param(k,:)';
    Sp_val = diag(T.ParamStd(k,:)'.^2);
    v_val  = T.ImGainRate(k);

    % 2023.05.17. (május 17, szerda), 13:19
    % A small correction
    Mx_val(J.R) = max(Mx_val(J.R),1);

    % Prediction error.
    mu_x_kp1 = full(Fn_Mx_pp(Mx_val,Mu_val,Mp_val,v_val));
    Mx_pp_error(k) = norm(mu_x_kp1 - T.x(k+1,:)')^2;
    
    % Compute feedback gain
    A_num = full(Fn_A(Mx_val,Mu_val,Mp_val,v_val));
    B_num = full(Fn_B(Mx_val,Mu_val,Mp_val,v_val));    
    K_val{k} = zeros(J.nu,J.nx);
    K_val{k}(:,Ctrl_idx) = dlqr(A_num(Ctrl_idx,Ctrl_idx),B_num(Ctrl_idx,:),eye(numel(Ctrl_idx)),R);
    
    % Joint variance of the actual state, input, and parameter    
    Sxup_val = Sxup_fh(Sx_val{k},K_val{k},Sxp_val{k},Sp_val);
    
    Sy_val{k} = full(hvar.Fn(Mx_val,Mp_val,Mu_val(1),Mu_val(2),Sxup_val));

    % Compute the input (transmission rate, immunity loss rate) covariance 
    Su_val{k} = K_val{k} * Sx_val{k} * K_val{k}';

    % Compute state-parameter covariance
    Sxp_val{k+1} = full(Fn_Sxp_pp(Mx_val,Mu_val,Mp_val,v_val,Sxp_val{k},Sp_val,K_val{k}));
    
    % Compute state variance
    Sx_val{k+1} = full(Fn_Sx_pp(Mx_val,Mu_val,Mp_val,v_val,Sxp_val{k},Sx_val{k},Sp_val,K_val{k}));

    status.progress(k);
end

Sx_Std = cellfun(@(S) {sqrt(diag(S)')}, Sx_val);
T.StateStd = vertcat(Sx_Std{:});

Su_val{end} = Su_val{end-1};
Su_Std = cellfun(@(S) {sqrt(diag(S)')}, Su_val);
Su_Std = vertcat(Su_Std{:});
T.TrRateStd = Su_Std(:,1);
T.ImLossRateStd = Su_Std(:,2);

Sy_val{end} = Sy_val{end-1};
T.OutputStd = sqrt([ Sy_val{:} ]');

% fig = figure(5);
% delete(fig.Children)
% ax = subplot(511); hold on, plot_mean_var(T.Date,T.ImLossRate,T.ImLossRateStd), ylim([0,0.02])
% ax = subplot(512); hold on, plot_mean_var(T.Date,T.TrRate,T.TrRateStd), ylim([0,3])
% ax = subplot(513); hold on, plot_mean_var(T.Date,T.L,T.StateStd(:,2)), ax.YLim(1) = 0;
% ax = subplot(514); hold on, plot_mean_var(T.Date,T.New_Cases,T.OutputStd(:,2)), ax.YLim(1) = 0;
% ax = subplot(515); hold on, plot_mean_var(T.Date,T.H,real(T.StateStd(:,J.H))), ax.YLim(1) = 0;



end
