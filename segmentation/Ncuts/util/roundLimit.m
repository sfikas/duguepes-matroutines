function x = roundLimit(x,p);
% Timothee Cour, 29-Aug-2006 07:49:15

x = round(x);
x = max(min(x,p),1);
