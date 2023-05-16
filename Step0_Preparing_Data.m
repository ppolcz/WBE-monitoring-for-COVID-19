function [T,P,s,T_Ms] = Step0_Preparing_Data(P,r)
arguments
%%
    P
    r.Clip = false;
    r.Date_End = datetime("today");
end
%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2022. July 07. (2022a)

CSV_OWID_TABLE_HUN = "Data/owid-covid-data-HUN-2023-02-10.csv";
XLS_WASTEWATER_DATA = "Data/NNK_Nyers_uj.xlsx";
XLS_HOSPITALIZATION_DATA = "Data/Korhaz_trend.xlsx";
XLS_NEW_HOSP_ADMISSION_2021 = "Data/new_adm_HUN_2021.xlsx";

%% Load parameters

% Read the appropriate parameter table (.xlsx), then, generate the
% parameter trajectories
% [P,Q,~] = Epid_Par.Get(r.CountryCode + r.ParamSetting);


%%

% Load full COVID data from the OWID database (automated download)
T = epid_load_cowid(CSV_OWID_TABLE_HUN, ...
    "Clip",r.Clip,"Date_End",r.Date_End);    

% Merge P and T
T = removevars(T,intersect(P.Properties.VariableNames,T.Properties.VariableNames));
T = synchronize(T,P,'first');
T.Properties.UserData = parsepropval('create',T.Properties.UserData,P.Properties.UserData);

T = renamevars(T,"New_Cases","New_Cases_off");

%% Read szennyviz

% For compatibility reasons
T.Szennyviz_Nyers = zeros(height(T),1) * NaN;
T.Szennyviz = zeros(height(T),1) * NaN;
T.Szennyviz_Nv = zeros(height(T),1) * NaN;

% Optional:
[T,T_Ms] = Step0_Read_Wastewater_data(XLS_WASTEWATER_DATA,T,"Clip",r.Clip,"Date_End",r.Date_End);

% Added on 2022.10.19. (október 19, szerda), 14:03
Date_Last_H = T.Properties.UserData.Date_Last_Available_H;
Date_Last_Szv = T.Properties.UserData.Date_Last_Szennyviz;
T.Properties.UserData.Date_Last_Relevant_New_Cases = max(Date_Last_H - 4,Date_Last_Szv - 2);

%% Read korhaz

% Optional:
T = Step0_Read_Hospitalization_data( ...
    XLS_HOSPITALIZATION_DATA, ...
    XLS_NEW_HOSP_ADMISSION_2021, ...
    T,"Clip",r.Clip,"Date_End",Date_Last_Szv);

% Optional:
T.H_off_ma(T.Date > T.Properties.UserData.Date_Last_Available_H) = NaN;
T.H_off(T.Date > T.Properties.UserData.Date_Last_Available_H) = NaN;
T.D_off(T.Date > T.Properties.UserData.Date_Last_Available_H) = NaN;

Nt = height(T);
Knot_Density = 14;
Nr_Knots = Nt / Knot_Density;
Spline_Order = 5;

H_sp = spap2(Nr_Knots,Spline_Order,T.Day,T.H_off_ma); 
d1_H_sp = fnder(H_sp);
d2_H_sp = fnder(d1_H_sp);
d3_H_sp = fnder(d2_H_sp);
d4_H_sp = fnder(d3_H_sp);

s = v2struct(H_sp,d1_H_sp,d2_H_sp,d3_H_sp,d4_H_sp);

end
