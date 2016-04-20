function ind = sub2ind2(pq,x,y);
% Timothee Cour, 29-Aug-2006 07:49:15

if nargin<3
    y=x(:,2);
    x=x(:,1);
end
ind = x+pq(1)*(y-1);
