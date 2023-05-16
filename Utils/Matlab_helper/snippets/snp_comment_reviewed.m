function returnval = snp_comment_reviewed(doch,event)
%%
%
%  file:   snp_comment_reviewed.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2016.01.31. Sunday, 12:17:19
%  Reviewed on 2016.02.17. Wednesday, 10:45:50
%  Reviewed on 2016.03.13. Sunday, 12:51:26
%  Reviewed on 2017. August 25. [date be informative]
%  Reviewed on 2017. September 24. [just return (.mlx)]
%  Modified on 2018. March 03. (Modified instead of Reviewed)
%  Reviewed on 2021. September 29. (2021b)
%  Revised on 2022. January 05. (2021b)
%

ret = sprintf('%%  Revised on %s (%s)\n', pcz_fancyDate('informative'), version('-release'));

try
    document = matlab.desktop.editor.getActive;
    selection = document.Selection;
    linenr = selection(1);
    document.insertTextAtPositionInLine(ret,linenr,1);    
catch ex
    getReport(ex)
    returnval = ret;
    
    disp('Error occured, but text copied to clipbord.')

    ret = ret(4:end-1);
    clipboard('copy',ret);
end

end

