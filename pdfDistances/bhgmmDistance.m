function y = bhgmmDistance(Rm, Rs2, RPj, m, s2, Pj)
%   Calculate distance of two gaussian mixture models, using
%       Bhatacharryya - based GMM distance.
%
%   Example:
%       bD = bhgmmDistance(m1, s2_1, Pj1, m2, s2_2, Pj2);
%
%   See also l2Distance, kullbackDistance, emdDistance, mahalanobisDistance
%
%   G.Sfikas 22/3/4
%
%   UPDATE HISTORY:
%   ---------------
%   21-10-04 full covariance
%   11-06-05 clear dist
%   30-10-05 Esvisa to kommati me thn KL - htan lathos. 
%                H KL einai ksexwrista twra.
%   13-11-06 Symmazema

clear dist;
Rgaussians = length(RPj);
gaussians = length(Pj);
for i = 1:Rgaussians
    for j = 1:gaussians
    dist((i-1)*gaussians + j) = RPj(i)*Pj(j)* bhattacharyya(Rm(:,i),Rs2(:,:,i),m(:,j),s2(:,:,j));
    end
end
y = sum(dist);
return;

function y = bhattacharyya(m1, s1, m2, s2)
%
% Calculates bhattacharyya distance of two gaussians.
% Prepei na einai idio dim (det(s1+s2)..)
% Updated 11/6/5 (more eps)
% Updated 11/6/5 (to 2 sto klasma tou ln prepei na ginei 2^dim)
% Updated 4/11/5 (allaksa to ln se log.. (giati eixa ftiaksei dikh moy
%                   synarthsh gia ton fys.logarithmo???)
y = 0.125 * (m1 - m2)'*(inv((s1+s2)/2))*(m1 - m2) + 0.5 * log( (det(s1+s2)+eps)/ ((sqrt(det(s1)*det(s2))*(2^size(s1,1))+eps)));
return