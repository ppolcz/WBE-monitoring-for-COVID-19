function [DISTANCE] = pcz_distlfr_report(syslfr1,syslfr2,tol,varargin)
%%
%  File: pcz_distlfr_report.m
%  Directory: 7_ftools/ftools/v12/utilities/symbolical
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. June 08. (2019b)
%

%%

if ischar(tol)
    varargin = [ tol varargin ];
    tol = 0.01;
end

title = '';
if ~isempty(varargin)
    title = sprintf(varargin{:});
end

TMP_ZNEWEagSzRkCGbFsczkg =  pcz_dispFunctionName(title,'',struct('parent',1));

pcz_info('Tolerance for `distlfr'': %g.', tol)

DISTANCE = distlfr(syslfr1,syslfr2);

if nargout == 0
    bool = DISTANCE < tol;
    pcz_dispFunction('LFR distance: %g', DISTANCE)
    pcz_dispFunctionSeparator
    pcz_info(bool, varargin{:}, {'first', 2})
end

pcz_dispFunctionEnd(TMP_ZNEWEagSzRkCGbFsczkg);

end