function [R] = pginput
%% 
%  
%  file:   pginput.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.12.07. Wednesday, 15:10:10
%

i = 1;
while 1
    [X(i),Y(i),button] = ginput(1);
    
    if button == 27, break, end
    
    hold on,
    plot(X(max(i-1,1):i),Y(max(i-1,1):i),'r.--');

    i = i + 1;
end

plot(X([1,i-1]),Y([1,i-1]),'r.--');

R = [X(1:i-1)' Y(1:i-1)'];

end