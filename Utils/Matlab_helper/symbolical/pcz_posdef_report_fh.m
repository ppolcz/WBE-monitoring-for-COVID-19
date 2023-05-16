function pcz_posdef_report_fh(fh, vars, plim, varargin)
%% pcz_posdef_report_fh
%
%  File: pcz_posdef_report_fh.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2019. February 15.
%

if numel(varargin) == 1 && isstruct(varargin{1})

    opts = varargin{1};
    CALLED_FROM_OTHER_REPORT = 1;

else

    opts.title = 'POSDEF check for a function handle';
    opts.postol = 1e-6;
    opts.tolerance = 1e-10;
    opts.RandomPoints = 0;
    opts.VectorArg = 1;
    opts = parsepropval(opts, varargin{:});

    TMP_ZSOAjugbWppryHbRDSsx = pcz_dispFunctionName(opts.title,'',struct('parent',1));
    CALLED_FROM_OTHER_REPORT = 0;

end

%%

np = numel(vars);
[P_v,P_Nr,Mode_str] = P_limvert_to_vert(plim,np);

N = opts.RandomPoints;
if N > 0
    RAND_v = ( rand(np,N) .* (plim * [-1 ; 1] * ones(1,N)) + plim(:,1)*ones(1,N) )';
else
    RAND_v = [];
end

PR_v = [
    P_v
    RAND_v
    ];

if opts.VectorArg
    zero_value = fh(zeros(np,1));
else
    zeros_cell = num2cell(zeros(np,1));
    zero_value = fh(zeros_cell{:});
end
    

pcz_info('Mode: %s, nr. of corners: %d, nr. of random points: %d', Mode_str, P_Nr, N);
pcz_info('Tolerance: %g, positive tolerance: %g.', opts.tolerance, opts.postol)
pcz_info('LMI size: (%dx%d) ', size(zero_value));
% pcz_dispFunction2('X_v = %s', pcz_num2str(X_v));
% pcz_dispFunction ''

feasible = zeros(P_Nr+N,1);
min_eig = zeros(P_Nr+N,1);
max_eig = zeros(P_Nr+N,1);
percentage = zeros(P_Nr+N,1);
iszero = zeros(P_Nr+N,1);

for i = 1:P_Nr+N
    if opts.VectorArg
        matrix = fh(PR_v(i,:)');
    else
        values_cell = num2cell(PR_v(i,:)');
        matrix = fh(values_cell{:});
    end
    
    [eigvals,min_eig(i),percentage(i)] = pcz_posdef_report(matrix, opts.tolerance);
    feasible(i) = min_eig(i) > -opts.tolerance;
    max_eig(i) = max(abs(eigvals));

    if max_eig(i) < opts.postol
        iszero(i) = 1;
    end
end

bool = all(feasible) && ~any(iszero);

for i = 1:P_Nr+N
    if min_eig(i) < 0 && feasible(i)
        pcz_warning(false, 'Least eigv: %10.5d, in %s (%s)', min_eig(i), CorP(i), pcz_num2str(PR_v(i,:)','format', '%g'));
    elseif ~feasible(i)
        pcz_info(false, 'Least eigv: %10.5d, in %s (%s)', min_eig(i), CorP(i), pcz_num2str(PR_v(i,:)','format', '%g'));
    end

    if iszero(i)
        pcz_info(false, 'LMI is almost zero (max eig: %d) in %s (%s)', max_eig(i), CorP(i), pcz_num2str(PR_v(i,:)','format', '%g'));
    end
end

pcz_dispFunctionSeparator

if any(iszero)
    [indices,~,overallperc] = pcz_symzero_report(iszero);
    pcz_info(false, 'LMI is almost zero in %d points out of %d, (%g%%).', numel(indices), P_Nr+N, overallperc);
end

if ~bool
    [indices,~,overallperc] = pcz_symzero_report(~feasible);
    pcz_info(false, 'This LMI is NOT feasible in %d points out of %d, (%g%%).', numel(indices), P_Nr+N, overallperc*100);
else
    pcz_info(true, 'This LMI is feasible along the given tolerance value.');
end

    function ret = CorP(nr)
        if nr <= P_Nr
            ret = 'corner point';
        else
            ret = 'random point';
        end
    end

%%

if ~CALLED_FROM_OTHER_REPORT
    pcz_dispFunctionEnd(TMP_ZSOAjugbWppryHbRDSsx);
end

end

