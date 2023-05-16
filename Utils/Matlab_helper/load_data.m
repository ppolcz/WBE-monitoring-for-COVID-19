function varargout = load_data(fn,dims,tspan)
arguments
    fn
    dims
    tspan = [0 Inf]
end
%%
%  File: load_data_msd.m
%  Directory: /home/ppolcz/_SZTAKI_Doc04/meresek/2022-04-11-10-37
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. April 25. (2022a)
%

K = numel(dims);
varargout = cell(1,K+1);

FileInfo = whos('-file',fn);
s = load(fn,FileInfo(1).name);
d = s.(FileInfo(1).name)';

dims = cumsum([0 dims]);
Msg = sprintf('Size of data = [%d(samples),%d(channels)], requested: %d(channels)!', size(d)-[0 1], dims(end));

if size(d,2) == dims(end)+1
    fprintf('\nLoading `%s`!\n%s\n',fn,Msg);
else
    error('Error loading `%s`!\n%s',fn,Msg);
end

% Load time
t = d(:,1);
d(:,1) = [];

% 2022.04.28. (április 28, csütörtök), 12:49
Ldx = isnan(t);
for i = 1:numel(tspan)/2
    Ldx = Ldx | ( tspan(2*i-1) <= t & t <= tspan(2*i) );
end

% [1] Time
varargout{1} = t(Ldx,:);

% [...] Data
for i = 1:K
    varargout{i+1} = d(Ldx,dims(i)+1:dims(i+1));
end

% [end] Sampling period
Ts = round(mean(diff(t)),10,"significant");
varargout{end+1} = Ts;
fprintf('Sampling time: Ts = %g [sec].\n\n',Ts);

end