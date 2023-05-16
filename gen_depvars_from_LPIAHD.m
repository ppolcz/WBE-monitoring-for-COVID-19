function T = gen_depvars_from_LPIAHD(T)
%%
%  File: gen_depvars_from_LPIAHD.m
%  Directory: /home/ppolcz/T/_Epid/RecPred_UIO_then_Opt/Ver_2022_07_12_UIO_then_Opt
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 12. (2022a)
%

K = Epid_Par.GetK;

% These notation ease further computations
L = T.L;
P = T.P;
I = T.I;
A = T.A;
H = T.H;
D = T.D;

% All infected + deceased
T.K = L + P + I + A + H + D;

% All infected
T.Infected = L + P + I + A + H;

% Daily new cases inferred by the UIO
T.New_Cases = [ L(2:end) + (T.Param(1:end-1,K.L_iPeriod) - 1) .* L(1:end-1) ; 0 ];
T.New_Recoveries = I .* T.Param(:,K.I_iPeriod) .* (1 - T.Param(:,K.Pr_H)) + A .* T.Param(:,K.A_iPeriod) + H .* T.Param(:,K.H_iPeriod) .* (1 - T.Param(:,K.Pr_D));

T.Cum_Cases = cumsum(T.New_Cases);
T.Cum_Recoveries = cumsum(T.New_Recoveries);

% These notation ease further computations
T.Infectious = P + I + T.Param(:,K.Rel_beta_A).*A;
betaSpN = T.New_Cases ./ T.Infectious;

% Effective reproduction rate inferred by the UIO
T.Rc = betaSpN .* ( ...
    1 ./ T.Param(:,K.P_iPeriod) ...
    + T.Param(:,K.Pr_I) ./ T.Param(:,K.I_iPeriod) ...
    + T.Param(:,K.Rel_beta_A).*(1-T.Param(:,K.Pr_I)) ./ T.Param(:,K.A_iPeriod));

end