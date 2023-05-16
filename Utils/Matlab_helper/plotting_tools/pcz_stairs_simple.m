function [ret] = pcz_stairs_simple(t,x,varargin)
%%
%  File: pcz_stairs.m
%  Directory: /home/ppolcz/_SZTAKI_Doc04/meresek/2022-04-11-10-37
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2022. April 12. (2022a)
%

tt = [t t]';
tt = tt(2:end-1);

xx = [x x]';
xx = xx(1:end-2);

ret = plot(tt,xx,varargin{:});


end