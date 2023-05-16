function [T,Q] = Step2_LTI_inversion_Opt(T,Q,args)
arguments
    T,Q
    args.UseWasteWater = false;
end
%%
%  Author: Peter Polcz (ppolcz@gmail.com)
%  Revised on 2022. September 12. (2022a)
% 
% 

Timer_5bix = pcz_dispFunctionName;

LPIAHD_Str = ["L","P","I","A","H","D"];

%% 

% Q_Szv = Epid_Par.Read("HUN_22Szept_Szv_Mtp");

% 2022.12.14. (december 14, szerda), 11:27
Q_Szv = Epid_Par.Read("HUN_Wastewater_Mtp");

%%

import casadi.*

Verbose = true;

% Idorendi sorrendben a kovetkezo datumok a fontosak:
Date_Start_MPC = max(datetime(2020,03,03),T.Date(1));
Date_Last_H = T.Properties.UserData.Date_Last_Available_H;
Date_Last_Szv = T.Properties.UserData.Date_Last_Szennyviz;
Date_Last_Relevant_New_Cases = max(Date_Last_H - 4,Date_Last_Szv - 2);
Date_End_REC = max(Date_Last_H,Date_Last_Szv);
Date_End_PRED = Date_End_REC + 10;

%% Adjust and expand time range of data table (T) if necessary
% Adjust end date of reconstruction

% Expand time range of data table
T = tb_expand_timerange(T,[T.Date(1) Date_End_PRED]);
T = gen_future_nans(T);

%% Match data table (T) and MPC dates

T.Day_MPC = T.Day - days(Date_Start_MPC - T.Date(1));
T.Label = categorical(zeros(size(T.Day)));
T = movevars(T,["Day_MPC","Label"],"After","Day");

%% Dynamical model

[f,h,~,J] = epid_LPIAHD_submodel;
K = Epid_Par.GetK;

%% Categorize time labels
% NREL: not relevant (not included in the MPC)
% REC:  reconstruct (the past)
% PRED: predict (the future)

T.Label(T.Day_MPC < 0) = "NREL";
T.Label(T.Day_MPC >= 0) = "REC";
T.Label(T.Date > Date_End_REC) = "PRED";
T.Label(T.Date > Date_End_PRED) = "NREL";

% Find relevant rows of the data table (which actually values take part in the MPC).
ldx_MPC = T.Label ~= "NREL";
ldx_Rec = T.Label == "REC";

%% Prediction horizon.

Idx_last_rec_z1 = days(Date_Last_Relevant_New_Cases - Date_Start_MPC) + 1;
Idx_last_pred_z1 = days(Date_End_PRED - Date_Start_MPC);
N_rec = days(Date_End_REC - Date_Start_MPC);
N_pred = days(Date_End_PRED - Date_End_REC);
N_recpred = N_rec + N_pred;

%% Log

if Date_Start_MPC < Date_End_REC
    pcz_dispFunction('Reconstruction from %s to %s (N_rec = %d days)', Date_Start_MPC, Date_End_REC, N_rec);
end

if Date_End_REC < Date_End_PRED
    pcz_dispFunction('Prediction to %s (N_pred = %d days)', Date_End_PRED, N_pred);
end

pcz_dispFunction('Horizon: 1+N_rec+N_pred = %d',N_recpred+1);

%% Initial guess for optimization
% Egyébként ezt egy jó stratégiának gondolom, hogy még mielőtt definiálnám
% a szabad változókat, már akkor van egy lehetséges sejtés azok értékére.

% Initial values
T0 = table2struct(timetable2table(T(T.Date == Date_Start_MPC,:)));

% Precomputed state values relavant in the MPC
T_MPC_x = T(ldx_MPC,LPIAHD_Str);
T = T(ldx_MPC,:);

% ---

% [GUESS] Initial guess for the state evaluation
x_guess = fillmissing(T_MPC_x.Variables,"previous");

% [Prepare variable]
x0 = x_guess(1,:)';
x_fh = @(x_var) [ x0' ; x_var ];
x_var_fh = @(x) x(2:end,:);

% [GUESS] cont...
x_var_guess = x_var_fh(x_guess);

% ---

% [Prepare variable]
z1_fh = @(z1_var) [ T0.New_Cases ; z1_var ; NaN ];
z1_var_fh = @(z1) z1(2:end-1);

% [GUESS] Initial guess for the unknown input
z1_guess = fillmissing(T.New_Cases,"previous");
z1_var_guess = z1_var_fh(z1_guess);

% ---

% [Prepare variable]
z2_fh = @(z2_var) [ T0.Infectious ; z2_var ; NaN ];
z2_var_fh = @(z2) z2(2:end-1);

% [GUESS] Initial guess for the instrumental variable infectious people
z2_guess = T.Infectious + randn(N_recpred+1,1);
z2_var_guess = z2_var_fh(z2_guess);

% ---

% Relevant section of the parameter table
ldx_Q = Date_Start_MPC < Q.Date + Q.TransitionLength;
ldx_Q(1:end-1) = Date_Start_MPC < Q.Date(2:end) + Q.TransitionLength(2:end);
Q_MPC = Q(ldx_Q,:);

% Relevant section of the parameter table
ldx_Q_Szv = Date_Start_MPC < Q_Szv.Date + Q_Szv.TransitionLength;
ldx_Q_Szv(1:end-1) = Date_Start_MPC < Q_Szv.Date(2:end) + Q_Szv.TransitionLength(2:end);
Q_Szv_MPC = Q_Szv(ldx_Q_Szv,:);

% Variant dominance pattern (+visualization)
P = epid_get_variant_dominance_pattern(Q_MPC,T.Date,"EnableOverlapping",true);
P_Szv = epid_get_variant_dominance_pattern(Q_Szv_MPC,T.Date,"EnableOverlapping",true);

% Check dimension of the relevant parameter table
nPattern = width(P.Pattern);
nQ_MPC = height(Q_MPC);
assert(nPattern == nQ_MPC,"Relevant parameter table width (%d) != (%d) No. patterns",nQ_MPC,nPattern)

% [GUESS] Initial guess for pI
pI_var_guess = Q_MPC.Pr_I;
pI_var_bounds = pI_var_guess .* (1 + [-1 1]*0.25);
pI_guess = P.Pattern * pI_var_guess;

% [GUESS] Initial guess for pH
pH_var_guess = Q_MPC.Pr_H;
pH_var_bounds = pH_var_guess .* (1 + [-1 1]*0.5);
pH_guess = P.Pattern * pH_var_guess;

% [GUESS] Initial guess for pD
pD_var_guess = Q_MPC.Pr_D;
pD_var_bounds = pD_var_guess .* (1 + [-1 1]*0.5);
pD_guess = P.Pattern * pD_var_guess;

% [GUESS] Initial guess for the szennyviz multiplier
wconst_var_guess = Q_Szv_MPC.Szennyviz_Mtp*0 + 1;
wconst_var_bounds = wconst_var_guess .* (1 + [-1 1]*1);

%% Free decision variables along the horizon (column-wise).

helper = Pcz_CasADi_Helper('SX');

x_var = helper.new_var('x',size(x_var_guess),1,'str','full','lb',0);
x = x_fh(x_var);

z1_var = helper.new_var('z1',size(z1_var_guess),1,'str','full','lb',0);
z1 = z1_fh(z1_var);

z2_var = helper.new_var('z2',size(z2_var_guess),1,'str','full','lb',0);
z2 = z2_fh(z2_var);

% 2022.12.14. (december 14, szerda), 15:27
SelectorI  = [
    0 % Transient
    0 % Wild type
    0 % Alpha
    0 % Delta
    0 % BA.1
    0 % BA.2 
    0 % BA.4
    0 % BA.5
    0 % School
    0 % Future
    ]*args.UseWasteWater*0;
pI_var = helper.new_var('pI',size(pI_var_guess),1,'str','full', ...
    'lb',pI_var_bounds(:,1),'ub',pI_var_bounds(:,2));
pI = P.Pattern * (pI_var.*SelectorI + pI_var_guess.*(1-SelectorI));

% 2023.04.24. (április 24, hétfő), 11:34
SelectorH  = [
    0 % Transient
    0 % Wild type
    1 % Alpha
    1 % Delta
    1 % BA.1
    1 % BA.2 
    1 % BA.5
    1 % Turn point 1
    1 % Turn point 2
    1 % Future
    ]*args.UseWasteWater;
pH_var = helper.new_var('pH',size(pH_var_guess),1,'str','full', ...
    'lb',pH_var_bounds(:,1),'ub',pH_var_bounds(:,2));
pH = P.Pattern * (pH_var.*SelectorH + pH_var_guess.*(1-SelectorH));

pD_var = helper.new_var('pD',size(pD_var_guess),1,'str','full', ...
    'lb',pD_var_bounds(:,1),'ub',pD_var_bounds(:,2));
pD = P.Pattern * pD_var;

wconst_var = helper.new_var('wconst',size(wconst_var_guess),1,'str','full', ...
    'lb',wconst_var_bounds(:,1),'ub',wconst_var_bounds(:,2));
wconst = P_Szv.Pattern * wconst_var;

Idx_Alpha = find(strcmp(Q_Szv_MPC.Properties.RowNames,'Alpha'));
Idx_Delta = find(strcmp(Q_Szv_MPC.Properties.RowNames,'Delta'));
Idx_OmicronBA1 = Idx_Delta + 1;
Idx_OmicronBA5 = Idx_OmicronBA1 + 1;

% helper.add_ineq_con('wconst',...
%     wconst_var(Idx_Alpha)*1.2 - wconst_var(Idx_Delta), ...
%     "ub",0)
% 
% helper.add_ineq_con('wconst',...
%     wconst_var(Idx_Delta)*1.5 - wconst_var(Idx_OmicronBA1), ...
%     "ub",0)
% 
% helper.add_ineq_con('wconst',...
%     wconst_var(Idx_Delta)*1.5 - wconst_var(Idx_OmicronBA5), ...
%     "ub",0)

w_BA1 = wconst_var(Idx_OmicronBA1);
w_BA5 = wconst_var(Idx_OmicronBA5);

% 2023.04.24. (április 24, hétfő), 12:50
% Fraction = 5;
% helper.add_ineq_con('wconst',...
%     (Fraction - 1) / Fraction * w_BA1 - w_BA5, ...
%     "ub",0)
% helper.add_ineq_con('wconst',...
%     (Fraction + 1) / Fraction * w_BA1 - w_BA5, ...
%     "lb",0)

% wconst_var = helper.new_var('wconst',1,1,'lb',-0.6,'ub',0.6);
% wconst = P.Pattern * wconst_var_guess + wconst_var;

%% Objective function

% Eddig ez mukodott es nagyon jol
wrel_SzvO = [
    1 % Transient
    1 % Wild type
    1 % Alpha
    1 % Delta
    5 % BA.1
    5 % BA.2 
    10 % BA.4
    10 % BA.5
    10 % Future
    ]*args.UseWasteWater; 

% 2023.01.14. (január 14, szombat), 10:57
wrel_SzvO = [
    5 % Transient
    15 % Wild type
    15 % Alpha
    15 % Delta
    5 % BA.1
    5 % BA.2 
    5 % BA.4
    5 % BA.5
    5 % Future
    ]*args.UseWasteWater; 

wrel_H = [
    1 % Transient
    1 % Wild type
    1 % Alpha
    1 % Delta
    1 % BA.1
    1 % BA.2 
    1 % BA.4
    1 % BA.5
    1 % Future
    ]*args.UseWasteWater + (1-args.UseWasteWater);

wRefH   = 1 * P.Pattern(:,1:end-1) * wrel_H(ldx_Q(1:end-1));
wRefD   = 1e-5;
wRefSzv_O = 1e-5 * P.Pattern(:,1:end-1) * wrel_SzvO(ldx_Q(1:end-1));
wRefSzv_Nv = wRefSzv_O * 0;
wsthpH  = 1e10;
wsthz1  = 1e-3;
wdevpH  = [
    1 % Transient
    1 % Wild type
    1 % Alpha
    1 % Delta
    1 % BA.1
    1 % BA.2 
    1 % BA.4
    1 % BA.5
    1 % School
    1 % Future
    ]*args.UseWasteWater*1e9;
wdevpI  = 1e9;

% Q_MPC

Idx_Href = find(T.Day_MPC > 0 & T.Date <= Date_Last_H);
H_ref = T.H_off_ma(Idx_Href);
helper.add_obj('H_error',(x(Idx_Href,J.H) - H_ref).^2,wRefH(Idx_Href));

D_ref = T.D_off_NNK(Idx_Href);
helper.add_obj('D_error',(x(Idx_Href,J.D) - D_ref).^2,wRefD);

Idx_SzvRef = find(T.Day_MPC > 0 & T.Date <= Date_Last_Szv & ~isnan(T.Szennyviz));
helper.add_obj('Szv_error',(z2(Idx_SzvRef) - wconst(Idx_SzvRef).*T.Szennyviz(Idx_SzvRef)).^2,wRefSzv_O(Idx_SzvRef));

% Idx_SzvRef_Nv = setdiff(Idx_SzvRef,find(isnan(T.Szennyviz_Nv)));
% helper.add_obj('Szv_error_Nv',(z2(Idx_SzvRef_Nv) - wconst(Idx_SzvRef_Nv).*T.Szennyviz_Nv(Idx_SzvRef_Nv)).^2,wRefSzv_Nv(Idx_SzvRef_Nv));

% Smooth daily new cases
Delta_z1 = (z1(1:end-2) - z1(2:end-1)).^2;
helper.add_obj('Smooth_z1',Delta_z1,wsthz1);

% Delta_pH = sumsqr(pH(1:end-2) - pH(2:end-1));
% helper.add_obj('Smooth_pH',Delta_pH,wsthpH);

% Delta_pI = sumsqr(pI(1:end-2) - pI(2:end-1));
% helper.add_obj('Smooth_pI',Delta_pI,wsthpH);

% Deviation of pH
Deviation_pH = (pH_var - pH_var_guess).^2;
helper.add_obj('Deviation_pH',Deviation_pH,wdevpH)

% Deviation of pI
Deviation_pI = sumsqr(pI_var - pI_var_guess);
helper.add_obj('Deviation_pI',Deviation_pI,wdevpI)

% 2023.04.24. (április 24, hétfő), 12:53
helper.add_obj('Deviation_wconst_Omicron',(w_BA1-w_BA5).^2,1e9);

%% Equality constraints (initialize)

% Enforce the output equation for the instrumental variable, the number of
% infectious people
Infectious = x(:,J.P) + x(:,J.I) + T.Param(:,K.Rel_beta_A) .* x(:,J.A);
Eq_z2 = Infectious - z2;
helper.add_eq_con( Eq_z2(Idx_SzvRef) )

% Update the parameter trajectory with the free decision parameter
% variables
Param = SX(T.Param);
Param(:,K.pH) = pH;
Param(:,K.pI) = pI;
Param(:,K.pD) = pD;

% Enforce the state equations
for k = 1:N_recpred
    x_kp1 = f.Fn(x(k,:),Param(k,:),z1(k));
    helper.add_eq_con( x_kp1 - x(k+1,:)' );
end

% Enforce the auxiliary constraint for the new cases, namely, the future
% new cases should follow the same curve that it showed in the near past.
for k = Idx_last_rec_z1+1:Idx_last_pred_z1
    helper.add_eq_con( z1(2*Idx_last_rec_z1 - k) + z1(k) - 2*z1(Idx_last_rec_z1) );
end

% Construct the nonlinear solver object
NL_solver = helper.get_nl_solver("Verbose",Verbose);

% Retain the mapping for the free variables, which allows to construct an
% initial guess vector for the nonlinear solver.
Fn_var = NL_solver.helper.gen_var_mfun;

%% Solve MPC

Solution_found = false;
try                       % x[383x5],   z1[382],     z2[382],     pI[5],       pH[5],       pD[5],       wconst[3]
    sol_guess = full(Fn_var(x_var_guess,z1_var_guess,z2_var_guess,pI_var_guess,pH_var_guess,pD_var_guess,wconst_var_guess));

    Timer_2xw3 = pcz_dispFunctionName('Calling the solver (with initial guess)');
    ret = NL_solver.solve([],sol_guess);
    ret.Elapsed_time = toc(Timer_2xw3);
    ret.Elapsed_time__desc = 'Solver initialized';
    pcz_dispFunctionEnd(Timer_2xw3);

    Solution_found = true;
catch e
    getReport(e)
end

if ~Solution_found
    Timer_2xw3 = pcz_dispFunctionName('Calling the solver (without initial guess)'); 
    ret = NL_solver.solve();
    ret.Elapsed_time = toc(Timer_2xw3);
    ret.Elapsed_time__desc = 'Solution from scratch';
    pcz_dispFunctionEnd(Timer_2xw3);
end

[~,F_MPC] = helper.get_obj;
display(F_MPC,'Achieved cost')
    
%% Retain and store solution

pI = helper.get_value('pI');
Q.Pr_I(ldx_Q) = pI;

pH = helper.get_value('pH');
Q.Pr_H(ldx_Q) = pH;

pD = helper.get_value('pD');
Q.Pr_D(ldx_Q) = pD;

wconst = helper.get_value('wconst');
Q.Szennyviz_Mtp(ldx_Q) = wconst([(1:end-2) end-2 end-1 end-1 end-1 end]);

x = x_fh(helper.get_value('x'));
z1 = z1_fh(helper.get_value('z1'));
z2 = z2_fh(helper.get_value('z2'));
pH = P.Pattern * pH;
pI = P.Pattern * pI;
pD = P.Pattern * pD;
Szennyviz_Mtp = P_Szv.Pattern * wconst;

for C = LPIAHD_Str
    T.(C) = x(:,J.(C));
end
T = gen_depvars_from_LPIAHD(T);

T.New_Cases = z1;
% T.Pr_H = pH;
% T.Pr_I = pI;
% T.Pr_D = pD;
T.Param(:,[K.pI,K.pH,K.pD]) = [pI,pH,pD];

T.Szennyviz_Mtp = Szennyviz_Mtp;

%%

pcz_dispFunctionEnd(Timer_5bix);

end
