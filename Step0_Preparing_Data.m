function [T,P,s] = Step0_Preparing_Data(P,r)
arguments
%%
    P
    r.Clip = false;
    r.Date_End = datetime(2022,02,06);
end
%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2022. July 07. (2022a)

CSV_OWID_TABLE_HUN = "Data/owid-covid-data-HUN-2023-02-10.csv";
XLS_CG_WW_Daily_Data = "Data/CG_WW_Daily_Data.xls";


%%

% Load full COVID data from the OWID database (automated download)
T = epid_load_cowid(CSV_OWID_TABLE_HUN, ...
    "Clip",r.Clip,"Date_End",r.Date_End);    

% Merge P and T
T = removevars(T,intersect(P.Properties.VariableNames,T.Properties.VariableNames));
T = synchronize(T,P,'first');
T = renamevars(T,"New_Cases","New_Cases_off");

%% Read wastewater data

CG_WW_Daily_Data = readtimetable(XLS_CG_WW_Daily_Data);
CG_WW_Daily_Data = renamevars(CG_WW_Daily_Data,"CG_WW_Daily_Data","WW");

T = synchronize(T,CG_WW_Daily_Data);

Idx = find(~isnan(T.WW),1,"last");
Date_Last_WW = T.Date(Idx);

Date_Last_H = T.Properties.UserData.Date_Last_Available_H;
T.Properties.UserData.Date_Last_WW = Date_Last_WW;
T.Properties.UserData.Date_Last_Relevant_New_Cases = max(Date_Last_H - 4,Date_Last_WW - 2);

%% Spline interpolation for medical data

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
