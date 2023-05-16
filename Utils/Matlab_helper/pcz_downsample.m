function [t,x,ppx,x_ZoH] = pcz_downsample(tt,xx,Ts,args)
arguments
    tt,xx,Ts
    args.MovMeanN = 5;
    args.MatchIntegral = false;
end
%%
%  File: pcz_downsample.m
%  Directory: /home/ppolcz/_SZTAKI_Doc04/sources/Main/11_QArm/Identification
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. April 25. (2022a)
%

t0 = tt(1);
T = tt(end);
t = (t0:Ts:T)';

% Moving average filter on the original dense data
xx_smooth = movmean(xx,args.MovMeanN);

% Downsample the smoothed dense data
x = interp1(tt,xx_smooth,t);

if args.MatchIntegral

    Ts_dense = round(mean(diff(tt)),4,'significant');

    int_xx = cumsum(xx).*Ts_dense;
    int_x = interp1(tt,int_xx,t);
    x_ZoH = [ diff(int_x)./Ts ; xx(end,:) ];

else
    x_ZoH = x;
end

% Fit a spline to the downsampled data
ppx = spline(t,x');

return
%%

% Evaluate the piecewise polynomial fit
xx_spline = ppval(ppx,tt)';

fig = figure(1);
delete(fig.Children)
plot(tt,xx)
hold on
plot(t,x);
stairs(t,x_ZoH);
plot(tt,xx_spline)

end