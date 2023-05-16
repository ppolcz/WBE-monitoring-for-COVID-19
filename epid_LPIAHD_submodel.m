function [f,h,hvar,J] = mf_epid_LPIAHD_submodel
%%
%  File: mf_epid_ode_model.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. October 12. (2021b)
%

import casadi.*

[~,np] = Epid_Par.GetK;
[s,p,~,p_Cas] = Epid_Par.Symbolic_New;

Cnt(0);
% State variables:
L = sym('L'); J.L = Cnt;
P = sym('P'); J.P = Cnt;
I = sym('I'); J.I = Cnt;
A = sym('A'); J.A = Cnt;
H = sym('H'); J.H = Cnt;
D = sym('D'); J.D = Cnt;
nx = Cnt - 1;
% -----
x = [L;P;I;A;H;D];

x_Cas = SX(nx,1);
for xi = x.'
    name = char(xi);
    x_Cas(J.(name)) = SX.sym(name);
end

% Unknown time-dependent parameters:
z1 = sym('z1');   % new_cases
z1_Cas = SX.sym('z1');

dL = z1 - s.tauL*L;
dP = s.tauL*L - s.tauP*P;
dI = s.pI*s.tauP*P - s.tauI*I;
dA = (1-s.pI)*s.tauP*P - s.tauA*A;
dH = s.pH*s.tauI*I - s.tauH*H;
dD = s.pD*s.tauH*H;

f_sym = [dL dP dI dA dH dD].';

Cnt(0);
J.Daily_All = Cnt;
J.Infectious = Cnt;
J.Rt = Cnt;

h_sym = [
    L + P + I + A + H
    P + I + s.qA*A
    z1/(P+I+s.qA*A) * (1/s.tauI + s.pI/s.tauI + s.qA*(1-s.pI)/s.tauA)
    ];

matlabFunction(f_sym,'File','Fn_LPIAHD_ode','Vars',{x,p,z1});
matlabFunction(h_sym,'File','Fn_LPIAHD_out','Vars',{x,p,z1});

Ts = 1;

f = {};
f.desc = 'f(x=LPIAHD,p,u=z1)';
f.val = x_Cas + Ts*Fn_LPIAHD_ode(x_Cas,p_Cas,z1_Cas);
f.Fn = Function('f',{x_Cas,p_Cas,z1_Cas},{f.val},{'x','p','u'},{f.desc});

f.      Input_1 = x_Cas;
f.desc__Input_1 = 'State vector (x = [L,P,I,A,H,D])';
f.      Input_2 = p_Cas;
f.desc__Input_2 = 'Parameter vector (p)';
f.      Input_3 = z1_Cas;
f.desc__Input_3 = 'Input (u = z1 = new cases)';

h = {};
h.desc = '[all_inf,intectious,Rt]';
h.val = Fn_LPIAHD_out(x_Cas,p_Cas,z1_Cas);
h.Fn = Function('h',{x_Cas,p_Cas,z1_Cas},{h.val},{'x','p','u'},{h.desc});

h.      Input_1 = x_Cas;
h.desc__Input_1 = 'State vector (x = [L,P,I,A,H,D])';
h.      Input_2 = p_Cas;
h.desc__Input_2 = 'Parameter vector (p)';
h.      Input_3 = z1_Cas;
h.desc__Input_3 = 'Input (u = z1 = new cases)';

h.desc__Output_dim1 = 'All infected';
h.desc__Output_dim2 = 'Infectious people';
h.desc__Output_dim3 = 'Time-dependent reproduction number (Rt)';

hvar = {};
hvar.val = zeros(size(h.val));
hvar.Fn = [];
vars = [x_Cas;p_Cas;z1_Cas];
Jh = jacobian(h.val,vars);

[Sigma,Sigma_half] = Pcz_CasADi_Helper.create_sym('Sigma',numel(vars),1,'str','sym');

hvar.val = diag(Jh * Sigma * Jh');
hvar.Fn = Function('hvar',{x_Cas,p_Cas,z1_Cas,Sigma_half},{hvar.val},{'x','p','z1','Sigma'},{'[Var_all_inf,Var_intectious,Var_Rt]'});

J.nx = numel(x);
J.ny = numel(h.val);
J.np = np;

end
