% Test for lsq_lut_piecewise function

rng(1)

Nr_intervals = 15;
Nr_samples = 1000;

% generate noisy x and y data
t = [-10:0.1:10]';
t = t+2*(rand(size(t))-0.5);
t = max(min(t,10),-10);
x = t+t.^2+10*sin(t)+5*(rand(size(t))-0.5);

% vector of 1-D look-up table "x" points
t_int = linspace(min(t),max(t),Nr_intervals);

% obtain vector of 1-D look-up table "y" points
x_int = lsq_lut_piecewise( t, x, t_int );


t_dense = linspace(min(t),max(t),Nr_samples);
x_dense = interp1(t_int,x_int,t_dense,'spline');

% plot fit
plot(t,x,'.',t_int,x_int,'+-',t_dense,x_dense)
legend('experimental data (x,y(x))','LUT points (XI,YI)')
title('Piecewise 1-D look-up table least square estimation')
