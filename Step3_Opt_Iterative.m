function T = Step3_Opt_Iterative(T)
%%
%  File: Step23_Convergence.m
%  Directory: /home/ppolcz/T/_Epid/RecPred_UIO_then_Opt/Ver_2022_09_09_Opt_SzennyvizRef
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. September 19. (2022a)
%

Timer_m4mu = pcz_dispFunctionName;

%%
T.New_Cases = max(0,fillmissing(T.New_Cases,"previous"));

tau_12 = 21;
tau_12_pm = 7;
tau_e = 14;
tau_e_pm = 4;
Verbose = true;
nr_sigma = 3;

sigma12e = round(sqrt(tau_12_pm^2 + tau_e_pm^2));

nr = 2;
Phi = @(x) normcdf(x,0,1);
Pr = @(k,sigma) ( Phi((2*k+1)/(2*sigma)) - Phi((2*k-1)/(2*sigma)) ) ./ ( 2*Phi((2*nr*sigma+1)/(2*sigma)) - 1 );
w12e = Pr(-nr*sigma12e:sigma12e*nr,sigma12e);
we = Pr(-nr*tau_e_pm:tau_e_pm*nr,tau_e_pm);

%% Vaccination model
% 2022.05.05. (május  5, csütörtök), 15:20

Vac_Delay = tau_12 + tau_e + tau_12_pm + tau_e_pm;
if T.Date(end) < T.Properties.UserData.Date_Last_Available_H + Vac_Delay
    T = tb_expand_timerange(T,[T.Date(1) T.Properties.UserData.Date_Last_Available_H+Vac_Delay]);
    T = gen_future_nans(T);
end

% Shift vaccination data
T_V = lag(T(:,"V_first"),tau_12 + tau_e);
T_W = lag(T(:,"V_boosted"),tau_e);

T_V = fillmissing(T_V,"constant",0);
T_W = fillmissing(T_W,"constant",0);

% Multiply first vaccination numbers with the probability to receive the
% second dose
p_v2 = T.V_second(end) / T.V_first(end); if isnan(p_v2), p_v2 = 1; end;
T_V.V_first = T_V.V_first * p_v2;

% Append variables to data table
T.Vac0_det = T.Np - T_V.V_first; % Not fully vaccinated
T.VacF_det = T_V.V_first;        % Fully vaccinated
T.VacB_det = T_W.V_boosted;      % Boosted

T.New_VacF = conv([0 ; diff(T.VacF_det)],w12e,'same');
T.New_VacB = conv([0 ; diff(T.VacB_det)],we,'same');

T.VacB = cumsum(T.New_VacB);
T.VacF = cumsum(T.New_VacF);
T.Vac0 = T.Np - T.VacF;

% Initial immunity loss rate value 
nu0 = T.New_VacF ./ (T.Np - T.VacF); % <- this should be updated iteratively 

%% Further data available from OWID

% A theoretical transmission rate computed from the stringency index (OWID)
beta_Is = T.beta0 .* T.Stringency;

% Filtering steps (diffusion) to obtain a smooth transmission rate
rwin = 15;
win = hamming(2*rwin+1);
win = win / sum(win);
beta_Is_filt = [ beta_Is(1)*ones(rwin,1) ; beta_Is ; beta_Is(end)*ones(rwin,1) ];
beta_Is_filt = conv(beta_Is_filt,win,'valid');

% Lower and upper bounds for beta
beta_lb_thry = 1/3 * (1 - 0.82);
beta_lb = max(beta_Is_filt*(1 - 0.2) - 0.05 , beta_lb_thry);
beta_ub = beta_Is_filt*(1 + 0.2) + 0.15;

% Allow delays in both directions
for i = 2:5
    beta_lb(i+1:end) = min(beta_lb(1:end-i),beta_lb(i+1:end));
    beta_lb(1:end-i) = min(beta_lb(1:end-i),beta_lb(i+1:end));

    beta_ub(i+1:end) = max(beta_ub(1:end-i),beta_ub(i+1:end));
    beta_ub(1:end-i) = max(beta_ub(1:end-i),beta_ub(i+1:end));
end

% Remove transient
Idx_Transient = 1:50;
beta_lb(Idx_Transient) = beta_lb_thry;
beta_ub(Idx_Transient) = Inf;

%% Construct optimization problem -- 1. Initial guess for optimization
% Egyébként ezt egy jó stratégiának gondolom, hogy még mielőtt definiálnám
% a szabad változókat, már akkor van egy lehetséges sejtés azok értékére.

% Initial guess for susceptibles
S0 = T.Np(1) - T.K(1);
S_fh = @(S) [ S0 ; S ];
S_var_fh = @(S) S(2:end);

% Immunity loss rate
w_fh = @(w) [ 0 ; w ; NaN ];
w_var_fh = @(w) w(2:end-1);

% Transmission rate
beta0 = 1/3;
beta_fh = @(beta) [ beta0 ; beta ; NaN ];
beta_var_fh = @(beta) beta(2:end-1);

% New cases as an independent variable
z1_fh = @(z) [ z ; NaN ];
z1_var_fh = @(z) z(1:end-1);

% Initial guess for susceptibles
S = S_fh(zeros(height(T)-1,1)+S0);
w0 = 0.01;
for k = 1:height(T)-1
    S(k+1) = (1 - nu0(k) - w0) * S(k) + w0*(T.Np(1) - T.K(k)) - T.New_Cases(k);
end
S_guess = S_var_fh(S);

% Initial guess for immunity loss rate
w_guess = zeros(height(T)-2,1) + w0;

% Initial guess for beta
beta_guess = beta_var_fh(beta_Is_filt);

% Lower and upper bounds for beta
beta_lb_var = beta_var_fh(beta_lb);
beta_ub_var = beta_var_fh(beta_ub);

% Initial guess for the new cases
z1_guess = z1_var_fh(T.New_Cases);

% These variables will be reassigned
clear S w

%% Construct optimization problem -- 2. Free decision variables
% (Column-wise along the horizon)

import casadi.*

helper = Pcz_CasADi_Helper('SX');

% Susceptibles
S_var = helper.new_var('S',size(S_guess),1,'str','full','lb',0,'ub',T.Np(1));
S = S_fh(S_var);

% Imminity loss rate
w_var = helper.new_var('w',size(w_guess),1,'str','full','lb',0,'ub',1);
w = w_fh(w_var);

% Transmission rate
beta_var = helper.new_var('beta',size(beta_guess),1,'str','full','lb',beta_lb_var,'ub',beta_ub_var);
beta = beta_fh(beta_var);

% New cases rate
z1_var = helper.new_var('z1',size(z1_guess),1,'str','full','lb',0);
z1 = z1_fh(z1_var);

Fn_var = helper.gen_var_mfun;

%% Construct optimization problem -- 3. Parameters

% Immunity gain rate as a parameter of the optimization problem
nu = helper.new_par('nu',height(T),1,'str','full');

%% Objective function

% Tracking error
w_beta = 1e1;
w_nce = 1e-7;

% Immunity loss cost
w_omega = 1e3;

% Smoothness
w_domega = 1e7;
w_dbeta = 1e3;

% Tracking error: the new cases (computed by LTI inversion)
Tracking_error_z1 = z1_var_fh(T.New_Cases - z1);
helper.add_obj('Tracking_error_z1',sumsqr(Tracking_error_z1),w_nce)

% Tracking error: beta (inferred from the stringency index)
Tracking_error_beta = beta_var_fh(beta - beta_Is_filt);
helper.add_obj('Tracking_error_beta',sumsqr(Tracking_error_beta),w_beta)

% Minimize the immunity loss rate (w)
helper.add_obj('omega',sumsqr(w_var),w_omega)

% Minimize the slope of the immunity loss rate (w)
dw = w(2:end-1) - w(1:end-2);
helper.add_obj('omega_rate',sumsqr(dw),w_domega)

% Minimize the slope of the transmission rate (beta)
dbeta = beta(2:end-1) - beta(1:end-2);
helper.add_obj('beta_rate',sumsqr(dbeta),w_dbeta)

%% Inequality constraints

% Hard constraint on the slope of the transmission rate
% helper.add_ineq_con('D_TrRate',dbeta,"lb",-0.005,"ub",0.005);

%% Equality constraints

% Dynamic equation for S
S_kp1_expr = (1 - nu - w) .* S + w.*(T.Np(1)-T.K) - z1;
Eq_con_S = S_kp1_expr(1:end-1) - S(2:end);
helper.add_eq_con(Eq_con_S);

% Algebraic equation for z1
z1_expr = beta .* S .* T.Infectious / T.Np(1);
Eq_con_z1 = z1_var_fh(z1_expr) - z1_var;
helper.add_eq_con(Eq_con_z1);

% 2023.01.03. (január  3, kedd), 08:36
% The sum of new cases should be preserved
helper.add_eq_con(sum(Tracking_error_z1)); % (Ez nagyon fontos!)

% Check whether the initial guess is a perfect solution or not.
%{

Fn_h = Function('h',{helper.x},{helper.h});
full(Fn_h(sol_guess))

%}

%% Generate NLP solver

NL_solver = helper.get_nl_solver("Verbose",Verbose);

%% Iterate

K = Epid_Par.GetK;

Iterated_nu = nu0;
for Iteration = 1:4    
    %% Solve the NLP and retrieve solution
    
    sol_guess = full(Fn_var(S_guess,w_guess,beta_guess,z1_guess));
    ret = NL_solver.solve(Iterated_nu(:,end),sol_guess);
    
    % Retrieve solution
    w = w_fh(helper.get_value('w'));
    S = S_fh(helper.get_value('S'));
    beta = beta_fh(helper.get_value('beta'));
    z1 = z1_fh(helper.get_value('z1'));
    
    % Recovered people
    R = T.Np(1) - S - T.K;
    
    % Compute achieved cost
    [f,F_Opt] = helper.get_obj;
    
    %% Further quantities and store data
    
    T.S = S;
    T.R = R;
    
    T.ImGainRate = Iterated_nu(:,end);
    T.ImLossRate = w;
    
    T.TrRate = beta;
    T.TrRate_Ref = beta_Is_filt;
    T.TrRate_bounds = [beta_lb , beta_ub];
    
    T.New_Cases = z1;
    
    % 2023.04.24. (április 24, hétfő), 11:55  
    T.Rc = T.TrRate .* ( ...
        1 ./ T.Param(:,K.P_iPeriod) ...
        + T.Param(:,K.Pr_I) ./ T.Param(:,K.I_iPeriod) ...
        + T.Param(:,K.Rel_beta_A).*(1-T.Param(:,K.Pr_I)) ./ T.Param(:,K.A_iPeriod) ...
        ) .* T.S / T.Np(1);

    T.New_ImGain = T.ImGainRate .* S;
    T.New_ImLoss = w .* R;
    
    T.Cum_ImGain = cumsum(T.New_ImGain);
    T.Cum_ImLoss = cumsum(T.New_ImLoss);
    T.Cum_Cases = cumsum(T.New_Cases);
    
%{
    figure
    plot(T.Date, [T.TrRate , T.New_Cases * T.Np(1) ./ T.Infectious ./ T.S ]);

    figure
    ncL = diff(T.L) + T.L_iPeriod(1:end-1) .* T.L(1:end-1);
    ncL = [ncL ; ncL(end)];
    plot(T.Date,[T.New_Cases ncL])
    keyboard
%}

% Solve recursion -- the simplest way of prediction
%{
    StateNames = ["L" "P" "I" "A" "H" "D"];
    for It = 1:100
        Tx = T(:,StateNames);        
        x = Tx.Variables;
        for i = 1:height(Tx)-1
            T.S(i+1) = T.S(i) - T.New_Cases(i) - T.ImGainRate(i)*T.S(i) + T.ImLossRate(i)*T.R(i);
            x(i+1,:) = x(i,:) + Fn_LPIAHD_ode(x(i,:)',T.Param(i,:)',T.New_Cases(i))';
            T.R(i+1) = T.Np(i+1) - sum(x(i+1,:)) - T.S(i+1);
        end
        Tx.Variables = x;
    
        T(:,StateNames) = Tx;
        T = gen_depvars_from_LPIAHD(T);
    
        plot([T.New_Cases T.TrRate .* T.Infectious .* T.S ./ T.Np])
        drawnow
    
        T.New_Cases = T.TrRate .* T.Infectious .* T.S ./ T.Np;        
    end
%}

% Solve recursion -- the simplest way of prediction
%{        
    StateNames = ["S" "L" "P" "I" "A" "H" "D" "R"];
    Tx = T(:,StateNames);
    
    x = Tx.Variables;
    for i = 1:height(T)-1
        x(i+1,:) - ( x(i,:) + Fn_SLPIAHDR_ode(x(i,:)',T.Param(i,:)',T.TrRate(i),T.ImGainRate(i),T.ImLossRate(i))' )
    end
    Tx.Variables = x;
    
    T(:,StateNames) = Tx;
    T = gen_depvars_from_LPIAHD(T);
%}
    
    %%
    
    S = T.S;
    S0 = T.S;
    Sf = T.S * 0;
    Sb = T.S * 0;
    
    R = T.R;
    R0 = T.R;
    Rf = T.R * 0;
    Rb = T.R * 0;
    
    Vf = T.New_VacF;
    Vb = T.New_VacB;
    % Vf = [0 ; diff(T.VacF_det)];
    % Vb = [0 ; diff(T.VacB_det)];
    
    z1 = T.New_Cases;
    z2 = T.New_Recoveries;
    
    w = T.ImLossRate;
    
    % 2023.04.24. (április 24, hétfő), 11:55
    Length_of_illness = ...
        1 ./ T.Param(:,K.L_iPeriod) ...
        + 1 ./ T.Param(:,K.P_iPeriod) ...
        + T.Param(:,K.Pr_I) ./ T.Param(:,K.I_iPeriod) ...
        + (1 - T.Param(:,K.Pr_I)) ./ T.Param(:,K.A_iPeriod) ...
        + T.Param(:,K.Pr_H) .* T.Param(:,K.Pr_I) ./ T.Param(:,K.H_iPeriod);
    
    Successfull_Vaccination = S0 * 0;
    nu = S0 * 0;
    
    for i = 1:height(T)-1
    
        Si = S0(i)+Sf(i)+Sb(i);
        Ri = R0(i)+Rf(i)+Rb(i);
        
        % Correction (gyorsitja, de lenyegeben nem befolyasolja a konvergenciat)
        if Si > 1
            S0(i) = S0(i) * S(i) / Si;
            Sf(i) = Sf(i) * S(i) / Si;
            Sb(i) = Sb(i) * S(i) / Si;
        end
        if Ri > 1
            R0(i) = R0(i) * R(i) / Ri;
            Rf(i) = Rf(i) * R(i) / Ri;
            Rb(i) = Rb(i) * R(i) / Ri;
        end
        
        Si = S0(i)+Sf(i)+Sb(i);
        Ri = R0(i)+Rf(i)+Rb(i);
        Vaccinable_f = 1/(S0(i)+R0(i));
        Vaccinable_b = 1/(Sf(i)+Rf(i)+Sb(i));
    
        j = max(round(i - Length_of_illness(i)),1);
        Sj = S0(j)+Sf(j)+Sb(j);
    
        if isinf(Vaccinable_f)
            Vaccinable_f = 0;
        end
    
        if isinf(Vaccinable_b)
            Vaccinable_b = 0;
        end
    
        S0(i+1) = S0(i) - z1(i) * S0(i)/Si ...              infection
                        - Vf(i) * S0(i)*Vaccinable_f ...    full vaccination
                        + w(i) * R0(i); ...                 immunity loss
    
        Sf(i+1) = Sf(i) - z1(i) * Sf(i)/Si ...              infection
                        - Vb(i) * Sf(i)*Vaccinable_b ...    booster
                        + w(i) * Rf(i); ...                 immunity loss
    
        Sb(i+1) = Sb(i) - z1(i) * Sb(i)/Si ...              infection
                        - Vb(i) * Sb(i)*Vaccinable_b ...    booster
                        + w(i) * Rb(i); ...                 immunity loss
        
        R0(i+1) = R0(i) + z2(i) * S0(j)/Sj ...              recovery
                        - Vf(i) * R0(i)*Vaccinable_f ...    full vaccination
                        - w(i) * R0(i); ...                 immunity loss
    
        Rf(i+1) = Rf(i) + z2(i) * Sf(j)/Sj ...              recovery
                        - Vb(i) * Rf(i)*Vaccinable_b ...    booster dose
                        + Vf(i) ...                         full vaccination
                        - w(i) * Rf(i); ...                 immunity loss
        
        Rb(i+1) = Rb(i) + z2(i) * Sb(j)/Sj ...              recovery
                        + Vb(i) ...                         booster dose
                        - w(i) * Rb(i); ...                 immunity loss
        
        Successfull_Vaccination(i) = 0 ...
            + Vf(i) * S0(i)*Vaccinable_f ...          fully vaccinated
            + Vb(i) * (Sf(i)+Sb(i))*Vaccinable_b; ... boosted
        nu(i) = Successfull_Vaccination(i) / Si;
    
    end
    Successfull_Vaccination(end) = Successfull_Vaccination(end-1);
    
    %%
    
    % nu_filt = conv(nu,hamming(31),'same');
    % nu_filt = nu_filt / sum(nu_filt) * sum(nu);
    % nu = nu_filt;
    
    T.ImGainRate = nu;

end
%%

%%

pcz_dispFunctionEnd(Timer_m4mu);
end

