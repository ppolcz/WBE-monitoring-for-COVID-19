function ret = pcz_print_plot_color
%%
%  File: pcz_print_plot_color.m
%  Directory: 7_ftools/ftools/v12.01/utilities/plotting_tools
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. December 01. (2021b)
%

Colors_Snippet = strjoin({
    'Color_1 = [0 0.4470 0.7410];'
    'Color_2 = [0.8500 0.3250 0.0980];'
    'Color_3 = [0.9290 0.6940 0.1250];'
    'Color_4 = [0.4940 0.1840 0.5560];'
    'Color_5 = [0.4660 0.6740 0.1880];'
    'Color_6 = [0.3010 0.7450 0.9330];'
    'Color_7 = [0.6350 0.0780 0.1840];'
    },newline);

Colors_Snippet_struct = strrep(Colors_Snippet,'Color_','Colors.Color_');

disp(Colors_Snippet)
disp('(Copied to clipboard)')
clipboard("copy",Colors_Snippet);

disp ' ' 
disp(Colors_Snippet_struct)

disp ' ' 
eval(Colors_Snippet_struct)

ret = Colors;

end