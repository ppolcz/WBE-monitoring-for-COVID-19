function f_sym = Fn_SLPIAHDR_ode(in1,in2,TrRate,ImGainRate,ImLossRate)
%Fn_SLPIAHDR_ode
%    F_SYM = Fn_SLPIAHDR_ode(IN1,IN2,TrRate,ImGainRate,ImLossRate)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    17-May-2023 16:25:34

A = in1(5,:);
H = in1(6,:);
I = in1(4,:);
L = in1(2,:);
P = in1(3,:);
R = in1(8,:);
S = in1(1,:);
alpha = in2(1,:);
delta = in2(6,:);
eta = in2(8,:);
gamma = in2(7,:);
lambda = in2(5,:);
mu = in2(9,:);
rhoA = in2(3,:);
rhoI = in2(4,:);
zeta = in2(2,:);
t2 = A.*delta;
t3 = L.*alpha;
t4 = A.*rhoA;
t5 = ImLossRate.*R;
t6 = ImGainRate.*S;
t7 = I+P+t2;
t8 = (S.*TrRate.*t7)./9.967304e+6;
f_sym = [t5-t6-t8;-t3+t8;t3-P.*zeta;-I.*rhoI+P.*gamma.*zeta;-t4-P.*zeta.*(gamma-1.0);-H.*lambda+I.*eta.*rhoI;H.*lambda.*mu;t4-t5+t6-I.*rhoI.*(eta-1.0)-H.*lambda.*(mu-1.0)];
