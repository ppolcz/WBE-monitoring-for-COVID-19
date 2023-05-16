function [ret] = snp_begin_end_scope(doch,event)
%%
%  File: snp_begin_end_scope.m
%  Directory: 7_ftools/utilities/snippets
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2016.01.31. Sunday, 15:04:58
%  Reviewed on 2021. November 18. (2021b)
%

[a,b] = pcz_generateBeginEndTimer;

try
    document = matlab.desktop.editor.getActive;
    selection = document.Selection;
    linenr = selection(1);

    b = strrep(b, '%%', '');
    b = regexprep(b, 'clear .*', '');
    b = strrep(b, newline, '');
    
    document.insertTextAtPositionInLine(a,linenr,1);
    clipboard('copy',b);
catch ex
    getReport(ex)
end


end
