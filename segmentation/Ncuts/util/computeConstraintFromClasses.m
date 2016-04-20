function [C,C12]=computeConstraintFromClasses(classes,indexes1,indexes2,nTot);
% Timothee Cour, 29-Aug-2006 07:49:15

classes=classes(:);
n1=length(classes);
n2=max(classes);
C=sparse(1:n1,classes,1,n1,n2);
sizes=full(sum(C,1));
sizes=1./sizes(:);

C12=sparse(classes,1:n1,sizes(classes),n2,n1);
if nargin <2
    C=[C12,-speye(n2)];    
else    
    C=sparse([classes;(1:n2)'],[indexes1;indexes2],[sizes(classes);-1*ones(n2,1)],n2,nTot);
end

