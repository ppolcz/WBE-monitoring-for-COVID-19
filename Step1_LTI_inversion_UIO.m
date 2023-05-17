function [T] = Step1_LTI_inversion_UIO(T,s)
%%
%  File: Main.m
%  Directory: /home/ppolcz/T/_Epid/Approach_LPIH_inv/Ver_2022_05_31
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. May 31. (2022a)
%

Timer_eoxg = pcz_dispFunctionName;

%% Load OWID and NNK data
% Generates
% # data table T
% # parameter table P
% # variant table Q
% # a struct Jp, with the important compartment and parameter order
% # a few string constants like cA = '\mathbf{A}'

%% State space model and filter design
%
%  u: L (latent)
%  x1: P (pre-symptomatic)
%  x2: I (symptomatic infected)
%  x3: H (hospitalized)

Jp = Epid_Par.GetK;
alpha = Jp.L_iPeriod;
zeta = Jp.P_iPeriod;
rhoI = Jp.I_iPeriod;
rhoA = Jp.A_iPeriod;
lambda = Jp.H_iPeriod;
gamma = Jp.Pr_I;
eta = Jp.Pr_H;
mu = Jp.Pr_D;

A = @(p) [
   -p(zeta)           0               0         % P
    p(gamma)*p(zeta) -p(rhoI)         0         % I
    0                 p(rhoI)*p(eta) -p(lambda) % H
    ];

B = @(p) [
    p(alpha)
    0
    0
    ];

C = [ 0 0 1 ];

% Initial states (L0 is not necessary)
PIH0 = [
    10
    10
    0
    ];
A0 = 10;
D0 = 0;

%%%
% Unknown Input Observer (UIO) design (Chen.Patton_1999)

Poles_Obsv = [-3.5; -3.51; -3.55]-4;

% Observability matrix
On = @(p) obsv(A(p),C);

% Eqs. (3.5)-(3.8)
H = @(p) B(p) / (On(p)*B(p));
K1_pardep = @(p) place((A(p) - H(p)*On(p)*A(p))',On(p)',Poles_Obsv)';
% K1_pardep = @(p) lqr((A(p) - H(p)*On(p)*A(p))',On(p)',eye(3),10);
K1 = K1_pardep(T.Param(1,:)'); % Gain K1 is computed for the original variant
F = @(p) A(p) - H(p)*On(p)*A(p) - K1*On(p);
K2 = @(p) F(p)*H(p);
K = @(p) K1 + K2(p);

Eig = [ eig(F(T.Param(1,:)')) eig(F(T.Param(end,:)'))];
display(Eig, ...
    'Poles of the error dynamics (for the original and the last variant)')
assert(all(Eig(:) < 0),'Observer not stable for all parameter configuration')

% Parameter function
p_fh = @(t) interp1(T.Day,T.Param,t)';

% Hospitalization function and its first two derivatives
Y_fh = @(t) [
    fnval(s.H_sp,t)
    fnval(s.d1_H_sp,t)
    fnval(s.d2_H_sp,t)
    ];

% State equation of the Unknown Input Observer
f_uio = @(t,x) F(p_fh(t)) * x + K(p_fh(t)) * Y_fh(t);

% Output function of the Unknown Input Observer
h_uio = @(t,x) x + H(p_fh(t)) * Y_fh(t);

% Unknown input expressed from the third derivative and the states (PIH)
g_uio = @(t,x) (C*A(p_fh(t))^2*B(p_fh(t))) \ ( fnval(s.d3_H_sp,t) - C*A(p_fh(t))^3 * x );

%% Simulation (continuous-time)

tspan = T.Day([1,end]);

% Simulate UIO dynamics with the computed output function (Y(t))
[tt,State_uio] = ode45(f_uio,tspan,PIH0 - H(p_fh(0))*Y_fh(0));

% Compute the output (i.e., [P,I,H]) of the UIO 
PIH_cell = cellfun(h_uio,num2cell(tt'),num2cell(State_uio',1), ...
    "UniformOutput",false);
PIH = [PIH_cell{:}]';

% As the state PIH is available, compute the input L from the third deriv.
L = cellfun(g_uio,num2cell(tt'),PIH_cell)';

% Simulate A
dA = @(p,A,P) (1-p(gamma))*p(zeta)*P - p(rhoA)*A;
dA_ode = @(t,A) dA(p_fh(t),A,interp1(tt,PIH(:,1),t));
[ttA,AA] = ode45(dA_ode,tspan,A0);
A = interp1(ttA,AA,tt);

% Simulate D
dD = @(p,D,H) p(mu)*p(lambda)*H;
dD_ode = @(t,D) dD(p_fh(t),D,interp1(tt,PIH(:,3),t));
[ttD,DD] = ode45(dD_ode,tspan,D0);
D = interp1(ttD,DD,tt);

% Append the state variables
LPIAHD = [L PIH(:,1:2) A PIH(:,3) D];

% A small correction:
LPIAHD(LPIAHD < 0) = 0;

%% Further quantities (discrete-time)
% which are append to the subtable T.UIO

LPIAHD = interp1(tt,LPIAHD,T.Day);

UIO = array2timetable(LPIAHD,"RowTimes",T.Date,"VariableNames",["L","P","I","A","H","D"]);
T = synchronize(T,UIO);
T = gen_depvars_from_LPIAHD(T);

%%

Date_Last_H = T.Properties.UserData.Date_Last_Available_H;
Date_Last_WW = T.Properties.UserData.Date_Last_WW;
Date_Last_Relevant_New_Cases = max(Date_Last_H - 4,Date_Last_WW - 2);

T.H(T.Date > T.Properties.UserData.Date_Last_Available_H) = NaN;
T.Infected(T.Date > Date_Last_Relevant_New_Cases) = NaN;
T.New_Cases(T.Date > Date_Last_Relevant_New_Cases) = NaN;

pcz_dispFunctionEnd(Timer_eoxg);
end
