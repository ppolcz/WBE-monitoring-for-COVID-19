function fig_latexify(fig, position, labels)
%%
%  File: fig_latexify.m
%  Directory: 5_Sztaki20_Main/Models/01_QArm/v4_2021_05_21/Helper
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. May 06. (2021a)
%

% bname = genvarname(fig.Name);
bname = strrep(fig.Name,' ','_');

fname = [ '/home/ppolcz/_/5_Sztaki20_Main/LaTeX/05_Laborszem_2021-05-11/fig/Exported_from_Matlab/' datestr(datetime,29) sprintf('_Fig_%d',fig.Number) '_' bname ];
imgname = [fname '.pdf'];
vecname = [fname '-vec.pdf'];

Font_Size = 14;

if nargin > 1 && ~isempty(position)
    if numel(position) == 2
        position = [ 0 0 position ];
    end
else
    position = fig.Position;
end

Hasznos_Magassag = [ 25 987 ];
hm = Hasznos_Magassag(1);
hM = Hasznos_Magassag(2);

position(1) = 5*fig.Number;
position(2) = min(max(hM - position(4) - 5*fig.Number , hm) , hM - fig.Position(4));

fig.Position = position;

axi = 0;
for i = 1:numel(fig.Children), ax = fig.Children(i); 
if isa(ax,'matlab.graphics.axis.Axes')
    ax.Box = 'on';
    ax.Title.Interpreter = 'latex';

    Logger.latexify_axis(ax,Font_Size)

    axi = axi + 1;
    if nargin > 2
        Logger.latexified_labels(ax,Font_Size,labels{axi,:})
    else
        ax.XLabel.FontSize = Font_Size;
        ax.YLabel.FontSize = Font_Size;
        ax.ZLabel.FontSize = Font_Size;
    end

    ax.Title.FontSize = Font_Size;
end
end

% if ~isfile(vecname)
%     exportgraphics(fig,vecname,ContentType="vector");
% end

% if ~isfile(imgname)
%     print('-dpdf',imgname)
%     pause(0.1);
%     system(sprintf('pdfcrop %s',imgname))    
% end


end