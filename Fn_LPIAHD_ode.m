function f_sym = Fn_LPIAHD_ode(in1,in2,z1)
%Fn_LPIAHD_ode
%    F_SYM = Fn_LPIAHD_ode(IN1,IN2,Z1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    03-May-2023 09:37:15

A = in1(4,:);
H = in1(5,:);
I = in1(3,:);
L = in1(1,:);
P = in1(2,:);
pD = in2(9,:);
pH = in2(8,:);
pI = in2(7,:);
tauA = in2(3,:);
tauH = in2(5,:);
tauI = in2(4,:);
tauL = in2(1,:);
tauP = in2(2,:);
t2 = L.*tauL;
f_sym = [-t2+z1;t2-P.*tauP;-I.*tauI+P.*pI.*tauP;-A.*tauA-P.*tauP.*(pI-1.0);-H.*tauH+I.*pH.*tauI;H.*pD.*tauH];
