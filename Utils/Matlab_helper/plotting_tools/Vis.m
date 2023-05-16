classdef Vis < handle
%%
%  File: Fig.m
%  Directory: /home/ppolcz/T/_Epid/Utils/Matlab_helper/plotting_tools
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 11. (2022a)
%

properties
    Fig
    tls
    Ax
    lines
    others

    ActNr
end

properties (GetAccess = private, SetAccess = private)
end

properties (GetAccess = public, SetAccess = public)
end

methods (Access = public)

    function o = Vis
        o.Fig = [];
        o.Ax = [];
        o.ActNr = 1001;
    end

    function ax = new_fig(o,fignr,args)
        arguments
            o
            fignr
            args.Name = [];
            args.SaveName = [];
        end

        fig = figure(fignr,"Name",args.Name);
        
        if isempty(o.Fig)
            o.Fig = fig;
        else
            o.Fig = [o.Fig ; fig];
        end

        delete(fig.Children);
        o.Ax(~isvalid(o.Ax)) = [];
        ax = axes('Parent',fig); 
        hold on; 
        grid on; 
        box on;        

    end

end

methods (Static)

    function Ax = clear(Ax,fignr,args)
        arguments
            fignr
            args.Name = [];
            args.SaveName = [];
        end

    end

end

methods (Access = private)

end

end
