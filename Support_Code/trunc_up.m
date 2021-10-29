function y = trunc_up(x,d)

% This function will take the value x and truncate it at the d'th decimal place.
% d should be a positive number. It is to show colorbar min bounding tick
% marks cleanly.
%
% Written by Chris Warner 10/11/16

sgn = sign(x);

if(sgn < 0)
    x = abs(x);
end

y = floor(x*10^d)*10^-d + sgn*10^-d;

if(sgn < 0)
    y = -y;
end