function T = gen_future_nans(T)
%%
%  File: gen_future_nans.m
%  Directory: /home/ppolcz/T/_Epid/RecPred_UIO_then_Opt/Ver_2022_07_12_UIO_then_Opt
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 16. (2022a)
%

T.H_off_ma(T.Date > T.Properties.UserData.Date_Last_Available_H) = NaN;
T.H_off(T.Date > T.Properties.UserData.Date_Last_Available_H) = NaN;
T.D_off(T.Date > T.Properties.UserData.Date_Last_OWID) = NaN;
T.Szennyviz(T.Date > T.Properties.UserData.Date_Last_Szennyviz) = NaN;
% T.Szennyviz_Nv(T.Date > T.Properties.UserData.Date_Last_Szennyviz_Nv) = NaN;

end