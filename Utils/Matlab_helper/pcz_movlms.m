function [theta] = pcz_movlms(phi,y,wN)
%%
%  File: pcz_movlms.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study31_Szennyviz_analizis
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. April 26. (2022a)
%

N = numel(y);

T = toeplitz(1:N,[1 zeros(1,wN-1)]);

% Kifejtve:
% Phi_Idx = num2cell(T(wN:end,end:-1:1),2);
% Phi_cell = cellfun(@(idx) {phi(idx,:)},Phi_Idx);
% y_cell = cellfun(@(idx) {y(idx,:)},Phi_Idx);
% theta = cellfun(@(phi,y) phi\y, Phi_cell,y_cell);

% Roviden:
theta = cellfun(@(idx) phi(idx,:)\y(idx,:), num2cell(T(wN:end,end:-1:1),2));

end