function [TrMsd,LegEntries,Mtp,Q,Q_orig] = pcz_detect_multiplier(Q,Date,Cmp,Msd,args)
arguments
    Q,Date,Cmp,Msd
    args.Legend = "Official data";
    args.Modify = true;
    args.RemoveLast = 0
end
%%
%  File: trsf_official_new_cases.m
%  Directory: /home/ppolcz/T/_Epid/RecPred_UIO_then_Opt/Ver_2022_05_31_UIO_then_Opt
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 11. (2022a)
%


%% Detect waning basis functions

Q_orig = Q;
for i = height(Q)-1:-1:2
    if Q.Shift(i) == Q.Shift(i-1) && Q.Mtp(i) == Q.Mtp(i-1)
        Q(i,:) = [];
    end
end

% Find index vector
Mat = (string(Q_orig.Properties.RowNames) == string(Q.Properties.RowNames)') .* (1:height(Q));
for i = 2:height(Mat)
    if sum(Mat(i,:)) == 0
        Mat(i,:) = Mat(i-1,:);
    end
end
Idx_Q = sum(Mat,2);

P = epid_get_variant_dominance_pattern(Q,Date([1,end]));

Idx = Date(1) <= P.Date & P.Date <= Date(end);
Pattern = P.Pattern(Idx,:);
Pattern_BIN = Pattern > 0.99 & ~isnan(Msd) & Msd ~= 0;

Lag = Q.Shift;
Mtp = Q.Mtp;

n = numel(Mtp);

Cmp = fillmissing(Cmp,'constant',0);

%%

TrMsd = zeros(numel(Date),n-1);
for i = 1:n-1
    T = timetable(Date,Msd);
    T = fillmissing(lag(T,-Lag(i)),'constant',0);
    
    if args.Modify
        Mtp(i) = floor( ...
            sum(Cmp .* Pattern_BIN(:,i)) / sum(T.Msd .* Pattern_BIN(:,i)) ...
            *10)/10;
    end
    TrMsd(:,i) = T.Msd .* Pattern(:,i) * Mtp(i);
end

% Remove the last data series if applicable
for i = n-1:-1:n-args.RemoveLast
    TrMsd(:,i-1) = TrMsd(:,i-1) + TrMsd(:,i) / Mtp(i) * Mtp(i-1);
    TrMsd(:,i) = [];
end

LegEntries = cell(1,n-1);
for i = 1:n-1
    % if Lag(i) > 0
    %     LegEntries{i} = sprintf('%s $\\times %g$, $%g$-day shifted ($\\leftarrow$)',args.Legend,Mtp(i),Lag(i));
    % elseif Lag(i) < 0
    %     LegEntries{i} = sprintf('%s $\\times %g$, $%g$-day shifted ($\\rightarrow$)',args.Legend,Mtp(i),-Lag(i));
    % else
        LegEntries{i} = sprintf('%s $\\times %g$',args.Legend,Mtp(i));
    % end
end
LegEntries = LegEntries(1:end-args.RemoveLast);
LegEntries = string(strjoin(LegEntries,"; "));

Q.Mtp = Mtp;
Q_orig.Mtp = Mtp(Idx_Q);

% Bar = bar(T.Date,NC_Off(:,1:n-1)',1,'stacked','FaceAlpha',0.5);
% legend(LegEntries,'Interpreter','latex')

end