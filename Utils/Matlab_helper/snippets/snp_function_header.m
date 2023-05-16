function snp_function_header(doch,event)
%%
%  file:   snp_function_header.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2016.01.31. Sunday, 17:38:20
%

try
    active = matlab.desktop.editor.getActive;
    editor = active.JavaEditor;

    f = pcz_resolvePath(active.Filename);
    
    ret = sprintf('function [varargout] = %s(varargin)\n', f.bname);
    
    editor.setCaretPosition(0);
    editor.insertTextAtCaret(ret);
    editor.appendText(sprintf('\nend'));    


    document = matlab.desktop.editor.getActive;
    selection = document.Selection;
    linenr = selection(1);
    document.insertTextAtPositionInLine(ret,linenr,1);    


catch ex
    getReport(ex)
end

end

