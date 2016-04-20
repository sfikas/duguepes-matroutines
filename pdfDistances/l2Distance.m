function y = l2Distance(Rm, Rs2, RPj, m, s2, Pj)
%
%   Calculate distance of two gaussian mixture models.
%   Uses a formula by Kostantinopoulos. That is, assuming p,q are the gmms,
%   and gr() means integral on dx:
%       D = log[ 2*gr(p(x) * q(x)) / gr(p^2(x) + q^2(x)) ]
%
%   s2 einai eite s^2 ean h diastash twn data einai 1, h dxd an h diastash
%   twn data einai d.
%   Pros to paron ypothetei oti to s2 ths eisodou einai scalar opote to
%   metatrepei se pinaka (s2 * I).
%   G.Sfikas 27/7/4
%
%   To s2 einai PINAKAS twra
%   17/10/4
%   See also gaussdistance.
%
%   Update 11/6/2005 : More eps
%   Update 16/11/2005: Twra doulevei gia gmms me diaforetiko arithmo
%   kernels.
%   Update 21/11/2005: Minor changes. Kapoia lathi sta dim/rdim.
%   Update 9/11/2006

Rgaussians = length(RPj);
gaussians = length(Pj);
dim = size(m, 1);
rdim = size(Rm, 1);

sum1 = 0; %Arithmhths
for i=1:gaussians
    for j=1:Rgaussians
        S = s2(:,:,i) + eps*eye(dim);   %(:,i) * eye(dim);CHANGE10/04      %Doulevoume me S,W sto ekshs anti s2,Rs2
        W = Rs2(:,:,j) + eps*eye(rdim);  %(:,j) * eye(dim);
%         V = inv(inv(S) + inv(W));
        V = S - S*inv(S + W)*S;
        mT = (m(:,i)'*inv(S) + Rm(:,j)'*inv(W))*V;
        k = m(:,i)'*inv(S)*(m(:,i) - mT') + Rm(:,j)'*inv(W)*(Rm(:,j) - mT');
        poso = det(V)/ (exp(k) * det(S) * det(W) + eps);
        sum1 = sum1 + Pj(i)*RPj(j)*sqrt(poso);
    end
end
sum2 = 0; %Paronomasths
for i=1:gaussians
    for j=1:gaussians
        S = s2(:,:,i) + eps*eye(dim); %(:,i) * eye(dim);      %Doulevoume me S,W sto ekshs anti s2,Rs2
        W = s2(:,:,j) + eps*eye(dim); %(:,j) * eye(dim);
%         V = inv(inv(S) + inv(W));
        V = S - S*inv(S + W)*S;
        mT = (m(:,i)'*inv(S) + m(:,j)'*inv(W))*V;
        k = m(:,i)'*inv(S)*(m(:,i) - mT') + m(:,j)'*inv(W)*(m(:,j) - mT');
        poso = det(V)/ (exp(k) * det(S) * det(W) + eps);
        sum2 = sum2 + Pj(i)*Pj(j)*sqrt(poso);
    end
end
for i=1:Rgaussians
    for j=1:Rgaussians
        S = Rs2(:,:,i) + eps*eye(rdim); %(:,i) * eye(dim);      %Doulevoume me S,W sto ekshs anti s2,Rs2
        W = Rs2(:,:,j) + eps*eye(rdim); %(:,j) * eye(dim);
%         V = inv(inv(S) + inv(W));
        V = S - S*inv(S + W)*S;
        mT = (Rm(:,i)'*inv(S) + Rm(:,j)'*inv(W))*V;
        k = Rm(:,i)'*inv(S)*(Rm(:,i) - mT') + Rm(:,j)'*inv(W)*(Rm(:,j) - mT');
        poso = det(V)/ (exp(k) * det(S) * det(W) + eps);
        sum2 = sum2 + RPj(i)*RPj(j)*sqrt(poso);
    end
end
y = -(log(2) + log(sum1) - log(sum2))*100;
% y = abs(log((2*sum1)/(sum2+eps)))*100;

return;