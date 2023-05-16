function [r_eigs,r_extremum,r_perc] = pcz_posdef_report(A, prec, N, varargin)
%% pcz_posdef_report
%  
%  File: pcz_posdef_report.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. April 26.
%

%%

if nargin < 3 || isempty(N) || ischar(N)
    if nargin >= 3 && ischar(N)
        varargin = [N varargin];
    end
    N = 'not used at the moment';
end

if nargin < 2 || isempty(prec) || ischar(prec)
    if nargin >= 2 && ischar(prec)
        varargin = [prec varargin];
    end
    prec = 10;
    tol = 10^(-prec);
elseif prec < 1
    tol = prec;
    prec = round(-log10(tol));
else
    tol = 10^(-prec);
end    

eigvals = eig(A);

isdef = all([ real(eigvals) == eigvals ; real(eigvals) >= -tol ]);

indices = find(real(eigvals) > -tol);

perc = 1 - numel(indices)/numel(eigvals);

extremum = min(eigvals);


if nargout > 2
    r_perc = perc;
end

if nargout > 1
    r_extremum = extremum;
end

if nargout > 0
    r_eigs = [ eigvals(indices,:) ; eigvals(setdiff(1:numel(eigvals),indices),:) ];
end

if nargout == 0    
    if isdef
        if ~isempty(varargin)
            varargin{1} = [ varargin{1} 'Tolerance: %g,  worst eigenvalue: %g' ];
        else
            varargin{1} = 'Tolerance: %g,  worst eigenvalue: %g';
        end
        
        varargin = [ varargin tol extremum ];
    end
    
    pcz_info(isdef, varargin{:}, {'first',2})
    
    if ~isdef
        pcz_dispFunction('Worst eigenvalue: %g', extremum)
        pcz_dispFunction('Percentage of good eigenvalues: %g%%', perc*100)
        pcz_dispFunction('Tolerance: %g', tol)
        pcz_dispFunction('Nr of bad eigenvalues/all: %d/%d', numel(indices), numel(eigvals))

        % pcz_dispFunctionStackTrace
        % pcz_dispFunctionNeedNewLine
    end
end


end