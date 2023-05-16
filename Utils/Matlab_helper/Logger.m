classdef Logger
%%
%  Author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2018.02.11. Sunday, 11:17:40
%  Major review on 2020. March 17. (2019b) -- simplified
%  Major review on 2020. May 2. (2019b)
%

%%
properties (GetAccess = public, SetAccess = private)
    dirname;
    file;
    diary_fname;
    latex_fname;
    results_fname;
    fig_fname;
    mat_fname;
    stamp;
    date;
    runID;
end

methods (Access = public)

    function p = Logger(fname, varargin)

        RUN_ID = getenv('RUN_ID');
        p.runID = str2double(RUN_ID);

        p.stamp = datestr(now, 'yyyy-mm-dd_HH:MM');
        p.date = datestr(now, 'yyyy.mm.dd. dddd, HH:MM:SS');

        [~,basename,~] = fileparts(fname);

        p.dirname = [ cd filesep 'results' filesep basename '-' p.stamp '_id' RUN_ID ];
        p.file = @(fname) [p.dirname filesep fname];

        p.diary_fname = p.file('output.txt');
        p.latex_fname = p.file('latex.tex');
        p.results_fname = p.file('results.xlsx');
        p.fig_fname = @(relname) p.file(relname);
        p.mat_fname =  @(relname) p.file([relname '.mat']);

        setenv('RESULTS_DIRNAME',p.dirname);
        setenv('DIARY_FNAME',p.diary_fname);
        setenv('LATEX_FNAME',p.latex_fname);
        setenv('RESULTS_FNAME',p.results_fname);

        if ~exist(p.dirname,'dir')
            mkdir(p.dirname)
        end

        % Start diary (output logging)
        diary(p.diary_fname);
        pcz_info('Output logging (with `diary''): %s', p.diary_fname);

    end

    function stoplog(p)
        diary off
        if exist(p.diary_fname,'file')
            pcz_output2log(p.diary_fname);
            pcz_dispFunction2(' ')
            pcz_info('Logfile `%s'' formatted!', p.diary_fname);
        end
    end

    function store_results_init(obj,sheetnr,varargin)

        Cols = [ {'DateTime','RunID'} varargin ];

        for i = 3:nargin
            if ~ischar(Cols{i})
                assert(~isempty(inputname(i)));
                Cols{i} = inputname(i);
            end
        end

        if ~exist(obj.results_fname, 'file')
            Results = cell2table(cell(0,numel(Cols)), 'VariableNames', Cols);
            writetable(Results,obj.results_fname,'Sheet',sheetnr);
        end

    end

    function store_results(obj,sheetnr,varargin)

        s.DateTime = obj.date;
        s.RunID = str2double(getenv('RUN_ID'));
        for i = 3:nargin
            if isempty(inputname(i))
                if ischar(varargin{i-2})
                    s.(varargin{i-2}) = [ '| ' strtrim(strrep(varargin{i-2},'_',' ')) ':' ];
                else
                    warning('Input %d does not have a name', i);
                end, continue
            end
            s.(inputname(i)) = varargin{i-2};
        end

        if exist(obj.results_fname, 'file')
            Results = readtable(obj.results_fname,'Sheet',1);
        end

        if ~exist('Results', 'var')
            Cols = fieldnames(s)';
            Results = cell2table(cell(0,numel(Cols)), 'VariableNames', Cols);
        else

            % Kiegeszitem az ujonan beszurando sort, ha nem szerepel benne
            % minden oszlop ami mar eddig benne van a tablazatban.
            not_in_s = setdiff(Results.Properties.VariableNames,fieldnames(s));
            for i = 1:numel(not_in_s)
                if isnumeric(Results.(not_in_s{i}))
                    Results.(not_in_s{i}) = 0;
                else
                    Results.(not_in_s{i}) = ' ';
                end
            end
        end

        Results = [ Results ; struct2table(s) ];

        writetable(Results,obj.results_fname,'Sheet',sheetnr);

        pcz_dispFunction('Results stored in `%s'', Sheet %d', obj.results_fname, sheetnr);
        pcz_dispFunction2(evalc('disp(s)'))
    end

end


properties (Constant = true)
    font_axis08 = {'FontSize',8,'FontName','TeX Gyre Schola Math'};
    font_axis10 = {'FontSize',10,'FontName','TeX Gyre Schola Math'};
    font_axis12 = {'FontSize',12,'FontName','TeX Gyre Schola Math'};
    font_axis14 = {'FontSize',14,'FontName','TeX Gyre Schola Math'};
    font_axis18 = {'FontSize',18,'FontName','TeX Gyre Schola Math'};
    font_axis22 = {'FontSize',22,'FontName','TeX Gyre Schola Math'};
    font_axis26 = {'FontSize',26,'FontName','TeX Gyre Schola Math'};
    font_latex18c = {'interpreter','latex','FontSize',18,'HorizontalAlignment','center','FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex8  = {'interpreter','latex','FontSize', 8,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex10 = {'interpreter','latex','FontSize',10,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex12 = {'interpreter','latex','FontSize',12,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex14 = {'interpreter','latex','FontSize',14,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex16 = {'interpreter','latex','FontSize',16,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex18 = {'interpreter','latex','FontSize',18,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex22 = {'interpreter','latex','FontSize',22,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex26 = {'interpreter','latex','FontSize',26,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex = {'interpreter','latex'};
end

methods(Static)
    function ret = latexify_axis(ax, ax_fontsize)
        if isnumeric(ax)
            ax_fontsize = ax;
            ax = gca;
        end

        try
            set(ax.Parent, 'Color', [1 1 1])
        end

        set(ax, ...'FontName','TeX Gyre Schola Math',...
            'GridColor', [0.1 0.1 0.1], 'MinorGridColor', [0.1 0.1 0.1],...
            'XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]);

        set(ax,'TickLabelInterpreter', 'latex');

        xl = get(ax,'XLabel');
        xlFontSize = get(xl,'FontSize');
        xAX = get(ax,'XAxis');
        set(ax,'FontSize', ax_fontsize)
        set(xl, 'FontSize', xlFontSize);

        yl = get(ax,'YLabel');
        ylFontSize = get(yl,'FontSize');
        yAX = get(ax,'YAxis');
        set(yAX,'FontSize', ax_fontsize)
        set(yl, 'FontSize', ylFontSize);

        zl = get(ax,'ZLabel');
        zlFontSize = get(zl,'FontSize');
        zAX = get(ax,'ZAxis');
        set(zAX,'FontSize', ax_fontsize)
        set(zl, 'FontSize', zlFontSize);

        if nargout > 0
            ret = ax;
        end
    end

    function ret = latexify_colorbar(cbar, fontsize)
        set(cbar, 'FontSize', fontsize, 'FontName','TeX Gyre Schola Math',...
            'TickLabelInterpreter', 'latex');
        cbar.Label.Color = [0 0 0];
        cbar.Label.Interpreter = 'latex';

        if nargout > 0
            ret = cbar;
        end
    end

    function ret_ = latexified_labels(ax, fontsize, xl, yl, zl)

        ret = cell(1,nargin - 2);
        if nargin > 2
            ret{1} = xlabel(ax,xl,'interpreter','latex','FontSize', fontsize,'FontName','TeX Gyre Schola Math','Color',[0 0 0]);

        end

        if nargin > 3
            ret{2} = ylabel(ax,yl,'interpreter','latex','FontSize', fontsize,'FontName','TeX Gyre Schola Math','Color',[0 0 0]);
        end

        if nargin > 4
            ret{3} = zlabel(ax,zl,'interpreter','latex','FontSize', fontsize,'FontName','TeX Gyre Schola Math','Color',[0 0 0]);
        end

        if nargout > 0
            ret_ = ret;
        end
    end

end

end
