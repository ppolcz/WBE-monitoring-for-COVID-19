classdef Epid_Par
%%
%  File: Epid_Par.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study18_timedep
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. January 08. (2021b)
%

properties (GetAccess = public, SetAccess = private)
    K,np
end

properties (GetAccess = public, SetAccess = public)
end

methods
    function [o,K,np] = Epid_Par()
        o.K = Epid_Par.GetK;
        K = o.K;
        o.np = numel(fieldnames(K));
        np = o.np;
    end
end

methods (Static)

    function [s,p] = Cas
        import casadi.*

        [K,np] = Epid_Par.GetK;

        s = {};
        s.alpha  = SX.sym('tauL');
        s.zeta   = SX.sym('tauP');
        s.rhoA   = SX.sym('tauA');
        s.rhoI   = SX.sym('tauI');
        s.lambda = SX.sym('tauH');
        s.gamma  = SX.sym('pI');
        s.eta    = SX.sym('pH');
        s.mu     = SX.sym('pD');
        s.delta  = SX.sym('qA');

        p = SX(np,1);
        for fn = string(fieldnames(s))'
            p(K.(fn)) = s.(fn);
        end
    end

    function [s,p] = Cas_New
        import casadi.*

        [K,np] = Epid_Par.GetK;

        s = {};
        s.tauL = SX.sym('tauL');
        s.tauP = SX.sym('tauP');
        s.tauA = SX.sym('tauA');
        s.tauI = SX.sym('tauI');
        s.tauH = SX.sym('tauH');
        s.pI   = SX.sym('pI');
        s.pH   = SX.sym('pH');
        s.pD   = SX.sym('pD');
        s.qA   = SX.sym('qA');

        p = SX(np,1);
        for fn = string(fieldnames(s))'
            p(K.(fn)) = s.(fn);
        end
    end

    function [s,p,S,P] = Symbolic
        [S,P] = Epid_Par.Cas;
        [K,np] = Epid_Par.GetK;

        s = {};
        p = sym(zeros(np,1));
        for fn = string(fieldnames(S))'
            s.(fn) = sym(fn,'real');
            p(K.(fn)) = s.(fn);
        end
    end

    function [s,p,S,P] = Symbolic_New
        [S,P] = Epid_Par.Cas_New;
        [K,np] = Epid_Par.GetK;

        s = {};
        p = sym(zeros(np,1));
        for fn = string(fieldnames(S))'
            s.(fn) = sym(fn,'real');
            p(K.(fn)) = s.(fn);
        end
    end

    function [K,np] = GetK()
        Cnt(0);
        K.L_iPeriod = Cnt;
        K.P_iPeriod = Cnt;
        K.A_iPeriod = Cnt;
        K.I_iPeriod = Cnt;
        K.H_iPeriod = Cnt;
        K.Rel_beta_A = Cnt;
        K.Pr_I = Cnt;
        K.Pr_H = Cnt;
        K.Pr_D = Cnt;
        np = Cnt-1;

        % Masneven:
        K.delta = K.Rel_beta_A;
        ...
        K.alpha = K.L_iPeriod;
        K.zeta = K.P_iPeriod;
        K.rhoA = K.A_iPeriod;
        K.rhoI = K.I_iPeriod;
        K.lambda = K.H_iPeriod;
        ...
        K.gamma = K.Pr_I;
        K.eta = K.Pr_H;
        K.mu = K.Pr_D;

        % Masneven: 2022.07.02. (július  2, szombat), 14:52
        K.qA = K.Rel_beta_A;
        ...
        K.tauL = K.L_iPeriod;
        K.tauP = K.P_iPeriod;
        K.tauA = K.A_iPeriod;
        K.tauI = K.I_iPeriod;
        K.tauH = K.H_iPeriod;
        ...
        K.pI = K.Pr_I;
        K.pH = K.Pr_H;
        K.pD = K.Pr_D;
    end

    % 2022.09.13. (szeptember 13, kedd), 17:29
    function [Q,xls] = Read(CountryCode)
        xls = "Parameters/" + sprintf('Par_%s.xlsx',CountryCode);
        Q = readtable(xls,"ReadRowNames",true,"Sheet","Main");
    end

    function T = Generate_Timetable(Q,xls)

        S = @Epid_Par.Interp_Sigmoid_v3;
        G = @Epid_Par.Interp_Pulse;

        VarMap = {
            'Period_L'        'L_iPeriod'       -1
            'Period_P'        'P_iPeriod'       -1
            'Period_A'        'A_iPeriod'       -1
            'Period_I'        'I_iPeriod'       -1
            'Period_H'        'H_iPeriod'       -1
            ...
            'Rel_beta_A'      'Rel_beta_A'       1
            'Pr_I'            'Pr_I'             1
            'Pr_H'            'Pr_H'             1
            'Pr_D'            'Pr_D'             1
            ...
            'beta0'           'beta0'            1
            'beta_max'        'beta_max'         1
            % 'New_Cases_Mtp'   'New_Cases_Mtp'    1
            % 'New_Cases_Shift' 'New_Cases_Shift'  1
            % 'Szennyviz_Mtp'   'Szennyviz_Mtp'    1
            % 'Szennyviz_Shift' 'Szennyviz_Shift'  1
            };

        % Waning function's table
        W = readtable(xls,"ReadRowNames",true,"Sheet","WaningPulse");

        % Uncertainty of each parameter
        U = readtable(xls,"Sheet","Uncertainty");
        u = U.Variables;

        % Helper cell variable with dates and transition lengths
        DateLen_cell = table2cell(Q(2:end,["Date","TransitionLength"]));
        
        d_Start = Q.Date(1) - Q.TransitionLength(1) - 10;
        d_End = Q.Date(end) + Q.TransitionLength(end) + 10;

        T = table;
        for i = 1:height(VarMap)
            xlsvar = VarMap{i,1};
            var = VarMap{i,2};
            Exp = VarMap{i,3};

            args = [DateLen_cell num2cell(Q.(xlsvar)(2:end))]';

            T.(var) = S(d_Start , Q.(xlsvar)(1) , args{:} , d_End)';

            if Exp == -1
                T.(xlsvar) = T.(var);
                T.(var) = 1./T.(xlsvar);
            end
        end
        T.Param = Epid_Par.Str2Vec(T)';
        T.ParamStd = u .* T.Param / 200;
        T.Date = (d_Start:d_End)';

        T = table2timetable(T,"RowTimes","Date");

        %% Detect waning basis functions

        Qw = Q(Q.wmax > 0,:);

        Bases_Raw = zeros(height(T),height(Qw)-1);
        for i = 1:height(Qw)-1
            Bases_Raw(:,i) = S(d_Start , 0 , Qw.Date(i) , Qw.TransitionLength(i) , 1 , Qw.Date(i+1) , Qw.TransitionLength(i+1) , 0 , d_End)';
        end

        % If wnoms have the same value, the given bases are considered one.
        [S_wnom,S_Unique_Idx,S_Idx] = unique(Qw(1:end-1,["wnom","wmax"]));
        S_n = height(S_wnom);

        Bases = zeros(height(T),S_n + height(W));
        for i = 1:S_n
            Bases(:,i) = sum(Bases_Raw(:,S_Idx == i),2);
        end
        for i = 1:height(W)
            Bases(:,S_n+i) = G(d_Start , 0 , W.PulseDate(i) , W.PulseWidth(i) , 1 , d_End)';
        end

        T.Phases = Bases_Raw;
        T.Waning_Bases = Bases;
        T.Properties.UserData.wnom = [ Qw.wnom(S_Unique_Idx) ; W.wnom ];
        T.Properties.UserData.wmax = [ Qw.wmax(S_Unique_Idx) ; W.wmax ]; 
        T.Properties.UserData.Waning_desc = [ Qw.Properties.RowNames(S_Unique_Idx) ; W.Properties.RowNames ];
        T.Properties.UserData.Phases_desc = Qw.Properties.RowNames(1:end-1);
    
        % 2023.04.21. (április 21, péntek), 11:47
        T = removevars(T,setdiff(unique([ VarMap(:,1) ; VarMap(:,2) ]),VarMap(end-1:end,1)));

    end

    function [T,Q,xls] = Get(CountryCode)
    arguments
        CountryCode = "HUN"
    end
    %%

        xls = "Parameters/" + sprintf('Par_%s.xlsx',CountryCode);
        
        % Parameters' table
        Q = readtable(xls,"ReadRowNames",true,"Sheet","Main");

        T = Epid_Par.Generate_Timetable(Q,xls);
        
    end

    function S = Sigmoid(A,B,N,args)
    arguments
        A,B,N
        args.Support = 5;
    end
        n = numel(A);
        assert(n == numel(B));

        S = zeros(n,N);
        for i = 1:numel(A)
            a = A(i);
            b = B(i);

            t = linspace(-args.Support,args.Support,N);
            s = -1 ./ (1 + exp(-t));
            
            s = s - s(1);
            s = s / s(end);
            S(i,:) = a + s * (b - a);
        end
    end

    function w = Pulse(N,args)
    arguments
        N
        args.Integral (:,1) = 1;
        args.Support = 3;
        args.Offset = 0;
        args.Window (1,:) char {mustBeMember(args.Window,{'hann','gaussian','hamming'})} = 'hann'
    end
        w = zeros(1,N);
        switch args.Window
            case 'hann'
                w = hann(N);
            case 'hamming'
                w = hamming(N);
            case 'gaussian'
                w = gausswin(N,args.Support);
        end
        w = w(:)' - w(1);
        w = w / sum(w) .* args.Integral;
    end

    % function s = Linear(a,b,N)
    %     s = linspace(0,1,N);        
    %     s = a + s * (b - a);
    % end

    function [pp,dd] = Interp_Sigmoid_v1(p,Dates)
        Dates = Dates';
        Dates = Dates(:)';

        assert(all(days(diff(Dates)) > 0),'Dates should be given in an increasing order.')

        N = days(Dates(end) - Dates(1)) + 1;

        t = days(Dates - Dates(1)) + 1;
        dd = Dates(1) + (0:N-1);

        pp = zeros(size(p,1),N);

        for i = 1:numel(Dates)-1
            Range = t(i):t(i+1)-1;
            Ni = t(i+1)-t(i);
            if mod(i,2) == 1
                % Zero order hold
                pp(:,Range) = p(:,(i+1)/2) * ones(1,Ni);
            else
                pp(:,Range) = Epid_Par.Sigmoid(p(:,i/2),p(:,i/2+1),Ni);
            end
        end
        pp(:,end) = p(:,end);
    end

    function [pp,dd] = Interp_Sigmoid_v2(d_Start,v,varargin)
        n = nargin/3-1;

        Dates = [varargin{1:3:3*n}]';
        WinR = [varargin{2:3:3*n}]';
        Vals = [v,varargin{3:3:3*n}];

        Dates = [ [d_Start ; Dates + WinR] , [ Dates - WinR ; varargin{end} ] ];

        [pp,dd] = Epid_Par.Interp_Sigmoid_v1(Vals,Dates);
    end

    % 2022.08.24. (augusztus 24, szerda), 11:22
    function [pp,dd] = Interp_Sigmoid_v3(D0,v0,varargin)
    %%
    %{
        D0 = datetime(2022,06,21);
        v0 = 0.35;
        varargin = { datetime(2022,07,20),15, ...
            0.3 , datetime(2022,08,01),5, ...
            0.2 , datetime(2022,08,31)};
    %}
        % Number of values to interpolate
        n = (numel(varargin)-1) / 3 + 1;
        Dn = varargin{end};

        dd = D0:Dn;
        pp = zeros(size(dd)) + v0;
        N = numel(pp);
        
        % if only a single value is given, i.e., Interp_Sigmoid_v3(D0,v0,D1)
        if n == 1
            return
        end

        d0 = D0;
        r0 = 0;
        for i = 1:n-1
            d1 = varargin{3*i-2};
            r1 = varargin{3*i-1};
            w1 = 2*r1+1;
            v1 = varargin{3*i};

            if d0+r0 <= d1-r1
                rMid = 1;
                dMid = d0+r0 + ceil(days(d1-r1-d0-r0)/2 - 0.5);
            else
                rMid = ceil(days(d0+r0-d1+r1) / 2 - 0.5);
                dMid = d1-r1 + rMid;
            end
    
            if i > 1
                alpha = Epid_Par.Interp_Sigmoid_v3(D0,0,dMid,rMid,1,Dn);
            else
                alpha = ones(size(pp));
            end
    
            Transition = Epid_Par.Sigmoid(v0,v1,w1);
            Idx_center = days(d1-D0)+1;
            Idx = Idx_center-r1:Idx_center+r1;
    
            ldx = 1 <= Idx & Idx <= N;

            pp1 = pp;
            pp1(Idx_center:end) = v1;
            pp1(Idx(ldx)) = Transition(ldx);
                
            pp = (1 - alpha) .* pp + alpha .* pp1;
    
            % plot(dd,pp)
            % grid on

            d0 = d1;
            v0 = v1;
            r0 = r1;
        end

        % plot(dd,pp)
        % grid on
    end

    % Usage: S( Start_Date , value , (Date,trlen , value) , (Date,trlen , value) , End_Date )
    function T = Interp_Sigmoid_Var(VariantName,varargin)
        [pp,Date] = Epid_Par.Interp_Sigmoid_v2(varargin{:});
        T = timetable(pp,'RowTimes',Date,'VariableNames',"Var_" + VariantName);
    end

    function [pp,dd] = Interp_Pulse(d_Start,Offset,varargin)
        n = nargin/3-1;
        d_End = varargin{end};
        Dates = [varargin{1:3:3*n}]';
        WinR = [varargin{2:3:3*n}]';
        WinL = WinR*2+1;
        Vals = [varargin{3:3:3*n}];

        N = days(d_End - d_Start) + 1;
        t = days(Dates - d_Start) + 1;
        dd = d_Start + (0:N-1);
        pp = zeros(1,N) + Offset;

        for i = 1:n
            Range = t(i)-WinR(i):t(i)+WinR(i);
            pp(Range) = pp(Range) + Epid_Par.Pulse(WinL(i),"Integral",Vals(i),"Offset",Offset,"Window","gaussian","Support",5);
        end
    end

    function p = Str2Vec(P)
        [K,np] = Epid_Par.GetK;
        p = zeros(np,numel(P.L_iPeriod));

        fns = fieldnames(P);
        for i = 1:numel(fns)
            if isfield(K,fns{i})
                p(K.(fns{i}),:) = P.(fns{i});
            end
        end
    end

    function P = Vec2Str(p)
        [K,~] = Epid_Par.GetK;

        fns = fieldnames(K);
        for i = 1:numel(fns)
            P.(fns{i}) = p(K.(fns{i}),:);
        end
    end

end

end

%{

alpha = P.L_iPeriod;
pi    = P.P_iPeriod;
rhoA  = P.A_iPeriod;
rhoI  = P.I_iPeriod;
h     = P.H_iPeriod;
gamma = P.Pr_I;
delta = P.Rel_beta_A;
eta   = P.Pr_H;
mu    = P.Pr_D;
vsp   = P.vsp;
orel  = P.orel;
omega = P.omega;

p = SX(np,1);
p(K.L_iPeriod)  = SX.sym('alpha');
p(K.P_iPeriod)  = SX.sym('pi');   
p(K.A_iPeriod)  = SX.sym('rhoA'); 
p(K.I_iPeriod)  = SX.sym('rhoI'); 
p(K.H_iPeriod)  = SX.sym('h');    
p(K.Pr_I)       = SX.sym('gamma');
p(K.Rel_beta_A) = SX.sym('delta');
p(K.Pr_H)       = SX.sym('eta');  
p(K.Pr_D)       = SX.sym('mu');   
p(K.vsp)        = SX.sym('vsp');  
p(K.orel)       = SX.sym('orel'); 
p(K.omega)      = SX.sym('omega');

alpha = SX.sym('alpha');     
pi    = SX.sym('pi');  
rhoA  = SX.sym('rhoA');    
rhoI  = SX.sym('rhoI');    
h     = SX.sym('h'); 
gamma = SX.sym('gamma');     
delta = SX.sym('delta');     
eta   = SX.sym('eta');   
mu    = SX.sym('mu');  
vsp   = SX.sym('vsp');   
orel  = SX.sym('orel');    
omega = SX.sym('omega');     

p(K.L_iPeriod)  = alpha;
p(K.P_iPeriod)  = pi;
p(K.A_iPeriod)  = rhoA;
p(K.I_iPeriod)  = rhoI;
p(K.H_iPeriod)  = h;
p(K.Pr_I)       = gamma;
p(K.Rel_beta_A) = delta;
p(K.Pr_H)       = eta;
p(K.Pr_D)       = mu;
p(K.vsp)        = vsp;
p(K.orel)       = orel;
p(K.omega)      = omega;

%}
