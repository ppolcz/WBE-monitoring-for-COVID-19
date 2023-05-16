function w_Mon = pf_week2date(y,w)
%%
%  File: pf_week2date.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study31_Szennyviz_analizis
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. April 26. (2022a)
%
% 
% Sunday:   weekday = 1
% Monday:   weekday = 2
% ...
% Saturday: weekday = 7
% 

y1 = datetime(y,1,1);
wd1 = weekday(y1);
wd1_HUF = wd1-1 + 7*(1-sign(wd1-1));

w1_Mon = y1 - wd1_HUF + 1;

Idx = y == 2021 | y == 2022;
w1_Mon(Idx) = w1_Mon(Idx) + 7;

w_Mon = w1_Mon + (w - 1)*7;

end