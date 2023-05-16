function [ret] = snp_header_comment
%%
%  File: snp_header_comment.m
%  Directory: 7_ftools/utilities/snippets
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. May 19. (2019b)
%

%%

comm = pcz_headerComment;

try
    active = matlab.desktop.editor.getActive;
    editor = active.JavaEditor;
    
    if startsWith(editor.getText,'function')
        caret = editor.lineAndColumnToPosition(2,0);
    else
        caret = editor.lineAndColumnToPosition(1,0);
    end

    editor.setCaretPosition(caret);
    editor.insertTextAtCaret(comm);    
catch ex
    getReport(ex)
    clipboard('copy',comm);
    disp('Text copied to clipboard instead!')
end


end