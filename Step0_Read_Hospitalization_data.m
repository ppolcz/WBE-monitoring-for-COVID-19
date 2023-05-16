function T = Step0_Read_Hospitalization_data(xls_hosp,xls_new_adm,T,r)
arguments
    xls_hosp,xls_new_adm,T
    r.Clip = false;
    r.Date_End = datetime("today");
end
%%
%  File: Tkn_Read_NNK.m
%  Directory: /home/ppolcz/T/_Epid/RecPred_UIO_then_Opt/Ver_2022_05_31_UIO_then_Opt
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 07. (2022a)
%
%  This is a possible and optional TUNING KNOB (Tkn).
% 

K = Epid_Par.GetK;

T_hosp = readtimetable(xls_hosp);
T_new_adm = readtimetable(xls_new_adm);

Date = T.Date;
Th_OWID = timetable(Date);
Th_OWID.H_OWID = T.H_off;
Th_OWID.H_OWID_ma = T.H_off_ma;
Th_OWID.Period_H_Est = 1./T.Param(:,K.H_iPeriod);

Th = synchronize(T_hosp,T_new_adm);
Th = synchronize(Th,Th_OWID);

%%

% 2022.12.15. (december 15, csütörtök), 15:40
if r.Clip
    Th(Th.Date > r.Date_End,:) = [];
end

Idx = find(Th.H_NNK > 0,1,'last');
Idx2 = find(Th.Covid_Miatt > 0,1,'last');
Date_Last_Available_H = Th.Date(max([Idx Idx2]));

weekends = weekday(Th.Date) == 1 | weekday(Th.Date) == 7;

Th.Resp_Inv(weekends) = NaN;
Th.H_NNK(weekends) = NaN;

% Hospitalization
Th.H_NNK_ma = zeros(height(Th),1) + NaN;
Th.H_NNK_ma(1:Idx) = movmean( resolve_missing_weekends(Th.H_NNK(1:Idx),Fill=0), 7);
Th.H_NNK_ma(Th.H_NNK_ma < 0.5) = 0;

% Respiratory 
Th.Resp_Inv_ma = zeros(height(Th),1) + NaN;
Th.Resp_Inv_ma(1:Idx) = movmean( resolve_missing_weekends(Th.Resp_Inv(1:Idx),Fill=0), 7);
Th.Resp_Inv_ma(Th.Resp_Inv_ma < 0.5) = 0;

Th.Covid_Miatt_ma = movmean( resolve_missing_weekends(Th.Covid_Miatt,Fill=0),7 );

% 2023.02.10. (február 10, péntek), 14:09
ldx = Th.Date >= datetime(2022,02,01);
Th.H_NNK_ma(ldx) = Th.Covid_Miatt_ma(ldx) * 1.8;

Idx = find(Th.New_Cases_NNK > 0,1,'last');
Date_Last_Available_New_Cases = Th.Date(Idx);
Th.New_Cases_NNK_ma(1:Idx) = movmean(Th.New_Cases_NNK(1:Idx), 7);

%%

% ldx = Th.Date > datetime(2022,01,31);
% Th.H_NNK_ma(ldx) = Th.Resp_Inv_ma(ldx)*25;

%% 

Date_Start_Est = datetime(2021,01,04);
Date_End_Est = datetime(2021,12,31);

Ldx = Date_Start_Est <= Th.Date & Th.Date <= Date_End_Est;

Date = Th.Date(Ldx);
hosp = Th.H_NNK_ma(Ldx);
new_adm = Th.Hadm_NNK(Ldx);

%use 7 day moving average
new_adm = movmean(new_adm, 7);


Ns = length(hosp);

try
    % the RLS method
    % initial values
    P_LS = 0.1;
    theta_k = 0.3;
    lambda = 0.95;
    theta = theta_k;
    hosp_time=[];
    
    for k=1 : (Ns-1)
        phi_k = hosp(k);
        yk = hosp(k+1) - hosp(k) - new_adm(k);
        P_LS = (lambda * (1/P_LS) + phi_k * phi_k )^(-1);
        theta_k = theta_k + P_LS*phi_k * (yk - phi_k*theta_k);
        theta = [theta; theta_k];
        hosp_time = [hosp_time; -1/theta_k];    
    end % for k= ...
    
    Th.Period_H_Est(Ldx) = min(hosp_time([1 1:end]),12);
    
    K = Epid_Par.GetK;
    k = K.H_iPeriod;

    Parameter_Estimation_OK = true;
catch
    Parameter_Estimation_OK = false;
end

%%

T = synchronize(T,Th,"first");
T.Period_H = T.Period_H_Est;
T.H_iPeriod = 1./T.Period_H;
T.H_off = T.H_NNK;
T.H_off_ma = T.H_NNK_ma;
T.New_Cases_OWID = T.New_Cases_off;
T.New_Cases_off = T.New_Cases_NNK_ma;

if Parameter_Estimation_OK
    T.Param(:,k) = T.H_iPeriod;
    T.Param(:,k) = fillmissing(T.Param(:,k),"previous");
end

T.D_off_OWID = T.D_off;
T.D_off_NNK = cumsum(fillmissing(T.D_NNK,"constant",0));

T.Properties.UserData.Date_Last_Available_H = Date_Last_Available_H;
T.Properties.UserData.Date_Last_Available_New_Cases = Date_Last_Available_New_Cases;

T = removevars(T,["D_NNK","D_off","H_NNK","H_NNK_ma","H_OWID","H_OWID_ma","Period_H_Est"]);

end

function X = resolve_missing_weekends(X,args)
arguments
    X
    args.Fill = 0;
    args.FillEnd = 'repeat';
    args.Method = 'spline';
end
%%
    
    t = 1:numel(X);

    X_isnum = ~isnan(X);

    t_ref = t(X_isnum);
    X_ref = X(X_isnum);
    
    Idx1 = find(X_isnum,1,'first');
    Idx2 = find(X_isnum,1,'last');

    t_ = Idx1:Idx2;
    X(t_) = interp1(t_ref,X_ref,t_,args.Method);

    if strcmp(args.FillEnd,'repeat')
        X(Idx2+1:end) = X(Idx2);
    end

    X(isnan(X)) = args.Fill;

end
