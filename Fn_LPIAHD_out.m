function h_sym = Fn_LPIAHD_out(in1,in2,z1)
%Fn_LPIAHD_out
%    H_SYM = Fn_LPIAHD_out(IN1,IN2,Z1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    17-May-2023 16:11:13

A = in1(4,:);
H = in1(5,:);
I = in1(3,:);
L = in1(1,:);
P = in1(2,:);
pI = in2(7,:);
qA = in2(6,:);
tauA = in2(3,:);
tauI = in2(4,:);
t2 = A.*qA;
t3 = 1.0./tauI;
t4 = I+P+t2;
h_sym = [A+H+I+L+P;t4;(z1.*(t3+pI.*t3-(qA.*(pI-1.0))./tauA))./t4];
