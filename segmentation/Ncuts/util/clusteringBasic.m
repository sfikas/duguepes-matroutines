function [classes,EigenvectorsRotated,convergence] = clusteringBasic(Eigenvectors,isRandom);
% Timothee Cour, 29-Aug-2006 07:49:15


[n,k]=size(Eigenvectors);

vm = sqrt(sum(Eigenvectors.*Eigenvectors,2));
Eigenvectors = Eigenvectors./(repmat(vm,1,k)+eps);

R=zeros(k);
if isRandom
    R(:,1)=Eigenvectors(1+round(rand(1)*(n-1)),:)';
else
    R(:,1)=Eigenvectors(1,:)';
end
c=zeros(n,1);
for j=2:k
    c=c+abs(Eigenvectors*R(:,j-1));
    [minimum,i]=min(c);
    R(:,j)=Eigenvectors(i,:)';
end

lastObjectiveValue=0;
exitLoop=0;
nbIterationsDiscretisation = 0;
nbIterationsDiscretisationMax = 20;%voir
while exitLoop== 0 
    nbIterationsDiscretisation = nbIterationsDiscretisation + 1 ;   
    
    EigenvectorsRotated = Eigenvectors*R;
    mex_normalizeColumns(EigenvectorsRotated);
    classes = mex_extractMaxima(EigenvectorsRotated);
    temp = mex_XindicatorTimesX(classes,Eigenvectors);
    [U,S,V] = svd(temp,0);

    NcutValue=2*(n-trace(S));
    if abs(NcutValue-lastObjectiveValue) < eps | nbIterationsDiscretisation > nbIterationsDiscretisationMax
        exitLoop=1;
    else
        lastObjectiveValue = NcutValue;
        R=V*U';
    end
end

if nbIterationsDiscretisation > nbIterationsDiscretisationMax
    convergence = 0;
    disp(sprintf('Attention, la discretisation n''a pas converge : NcutValue = %g',NcutValue));
else
    convergence = 1;
end

classes = double(classes);
