function y = kullbackDistance(Rm, Rs2, RPj, m, s2, Pj, RsampleVectors, sampleVectors)
%
%   kullbackdistance(Rm, Rs2, RPj, m, s2, Pj, RsampleVectors,
%   sampleVectors)
%
%   Calculate distance of two gaussian mixture models.
%   Using Symmetric Kullback-Liebler distance:
%       |kl(p||q) + kl(q||p)|
%
%   Check correctness by using:
%
%   kullbackDistance(0, 4, 1, 0, 1, 1, randn(1,5000)*2, randn(1,5000));
%
%   This computes divergence on two gaussians with deviations
%   2 and 1, both centered on 0, and it should result approximately 1.12 .
%
%
% Update 30/10/2005
%    Svhnw ta Rimage/ image orismata..
%   Pleon deigmata tha pernoume me deigmatolhpsia apo to modelo, anti
%   na xrhsimopoieitai plhroforia ths eikonas.
%   Volevei etsi, giati h KL einai h monh apostash pou xreiazetai ta
%   deigmata (logw monte carlo simulation)
% Updated 12/6/2005: Esvhsa 'global inimage Rinimage' (axrhsto?)
%   Vgazw apo ta sxolia ton kwdika undersampling, kai
%   ton allazw gia na doulepsei gia eikones diaf. megethous
%   Exw th symvash: samples = Rsamples = min(MAXsamples, MAXRsamples); 
% Update 20/6/2006
%   No 'RGB' or dimensionality or same kernel number assumption.
%   inimage, Rinimage are considered (:,:,dim) data.
% Update 22/6/2006
%   Ta inimage, Rinimage ginontai sampleVectors kai
%   twra einai pinakas-sthles me dianysmata-samples.
% Update/Fix 28/6/2006
%   Diorthwsh lathos proshmwn stous prostheteous 3, 4.
% Update/Fix 23/10/2006
%   Diorthwsh ston arithmo twn samples. To 2.5% twn samples
%   htan provlhmatiko otan eixame poly liga maximum deigmata.
%
dim = size(Rm, 1);
% Backwards compatibility..
inimage = sampleVectors';
Rinimage = RsampleVectors';
MAXsamples = size(sampleVectors, 2);
MAXRsamples = size(RsampleVectors, 2);
samples = MAXsamples;
Rsamples = MAXRsamples;
% Modify the sample size, for better speed.
% Worse result, though.
if min(MAXsamples, MAXRsamples) < 600
    samples = ceil(min(MAXsamples, MAXRsamples) / 4);
elseif min(MAXsamples, MAXRsamples) < 10000;
    samples = 250;
else
    samples = ceil(min(MAXsamples, MAXRsamples) / 40);
end
Rsamples = samples;

j = 1;
for i = uint32(1:double(MAXsamples/samples):double(MAXsamples))
    temp(j,:) = inimage(i,:);
    j = j+1;
end
inimage = temp;
clear temp;
j = 1;
for i = uint32(1:double(MAXRsamples/Rsamples):double(MAXRsamples))
    temp(j,:) = Rinimage(i,:);
    j = j+1;
end
clear Rinimage;
Rinimage = temp;

y = 0;
M = size(m, 2); % Mixing kernels
N = size(Rm, 2);
clear temp;
% Addend 1 (Entropy)
for j = 1:M
    temp = inimage - ones(length(inimage(:,1)),1) * m(:,j)';
    tempN1 = temp * inv(s2(:,:,j)+eps*eye(dim)) .* temp;
    px_j(:,j) = ((1/(2*pi)).^(dim/2)) * ((det(s2(:,:,j)+eps*eye(dim))).^(-1/2)) * exp(-0.5*sum(tempN1,2));
end
temp3 = px_j*Pj';
sum1 = sum(log(temp3(isfinite(temp3))+eps));
y = y + sum1 / samples;
% Addend 2 (Cross-Entropy)
clear temp;
clear tempN1;
clear px_j;
clear temp3;
for j = 1:N
    temp = inimage - ones(length(inimage(:,1)),1) * Rm(:,j)';
    tempN1 = temp * inv(Rs2(:,:,j)+eps*eye(dim)) .* temp;
    px_j(:,j) = ((1/(2*pi)).^(dim/2)) * ((det(Rs2(:,:,j)+eps*eye(dim))).^(-1/2)) * exp(-0.5*sum(tempN1,2));
end
temp3 = px_j*RPj';
sum1 = sum(log(temp3(isfinite(temp3))+eps));
y = y - sum1 / samples;
% Addend 3 (Cross - Entropy)
clear temp;
clear tempN1;
clear px_j;
clear temp3;
for j = 1:M
    temp = Rinimage - ones(length(Rinimage(:,1)),1) * m(:,j)';
    tempN1 = temp * inv(s2(:,:,j)+eps*eye(dim)) .* temp;
    px_j(:,j) = ((1/(2*pi)).^(dim/2)) * ((det(s2(:,:,j)+eps*eye(dim))).^(-1/2)) * exp(-0.5*sum(tempN1,2));
end
temp3 = px_j*Pj';
sum1 = sum(log(temp3(isfinite(temp3))+eps));
y = y - sum1 / Rsamples;
% Addend 4 (Entropy)
clear temp;
clear tempN1;
clear px_j;
clear temp3;
for j = 1:N
    temp = Rinimage - ones(length(Rinimage(:,1)),1) * Rm(:,j)';
    tempN1 = temp * inv(Rs2(:,:,j)+eps*eye(dim)) .* temp;
    px_j(:,j) = ((1/(2*pi)).^(dim/2)) * ((det(Rs2(:,:,j)+eps*eye(dim))).^(-1/2)) * exp(-0.5*sum(tempN1,2));
end
temp3 = px_j*RPj';
sum1 = sum(log(temp3(isfinite(temp3))+eps));
y = y + sum1 / Rsamples;

if y < 0
    disp('Warning! Sample for KL distance probably too small; computed negative distance.');
end
y = abs(y); 
return;