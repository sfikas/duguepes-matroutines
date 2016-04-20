function y = emdDistance(Rm, Rs2, RPj, m, s2, Pj)
%
%   Calculate distance of two gaussian mixture models.
%   Use Earth movers / Frechet distance.
%
%   G.Sfikas 11/11/5
%
% For gaussians = 2, Rgaussians = 3, the Aeq matrix 'll be:
%
%      1     1     1     0     0     0     1     0     0     0     0
%      0     0     0     1     1     1     0     1     0     0     0
%      1     0     0     1     0     0     0     0     1     0     0
%      0     1     0     0     1     0     0     0     0     1     0
%      0     0     1     0     0     1     0     0     0     0     1
%      1     1     1     1     1     1     0     0     0     0     0


clear dist;
Rgaussians = length(RPj);
gaussians = length(Pj);
for i = 1:gaussians
    for j = 1:Rgaussians
        dist((i-1)*Rgaussians + j) = frechet(m(:,i),s2(:,:,i), Rm(:,j),Rs2(:,:,j));
        distMX(i,j) = dist((i-1)*Rgaussians + j);
    end
end
sourceEarth = sum(Pj);
targetEarth = sum(RPj);
movableEarth = min(sourceEarth, targetEarth);
%
% Now set optimization constraints
%
for i = 1:gaussians
    Aeq(i, (Rgaussians*(i-1) + 1):Rgaussians*(i-1) + Rgaussians) = ones(1,Rgaussians);
    Aeq((gaussians + 1):(gaussians + Rgaussians), (Rgaussians*(i-1) + 1):(Rgaussians*(i-1) + Rgaussians)) = eye(Rgaussians);
end
Aeq(1:(gaussians+Rgaussians), (gaussians*Rgaussians + 1):(gaussians*Rgaussians + gaussians+Rgaussians)) = eye(gaussians+Rgaussians);
Aeq(gaussians+Rgaussians+1, 1:gaussians*Rgaussians) = ones(1, gaussians*Rgaussians);
%
beq(1:gaussians, 1) = Pj(1, 1:gaussians);
beq(gaussians +1:Rgaussians + gaussians, 1) = RPj(1, 1:Rgaussians);
beq(gaussians + Rgaussians + 1, 1) = movableEarth;
%
% Cost function & lower bound
%
f = zeros(gaussians*Rgaussians + gaussians + Rgaussians, 1);
f(1:gaussians*Rgaussians, 1) = dist(1:gaussians*Rgaussians);
lb = zeros(gaussians*Rgaussians + gaussians + Rgaussians, 1);
%
options = optimset('LargeScale', 'off', 'Simplex', 'on');
[x,fval,exitflag,output,lambda] = linprog(f,[],[],Aeq,beq,lb,[],[],options);
%
% Save the flow matrix.
%
for i = 1:gaussians
    flowMX(i, 1:Rgaussians) = x(1+(i-1)*Rgaussians:i*Rgaussians)';
end
movableEarth;
y = fval / movableEarth;
return

function y = frechet(m1, s1, m2, s2)
%
% Calculates Frechet distance of two gaussians.
% Frechet distance is closed form solution of the EMD in the case
% of two gaussians.
%
y = sqrt(  (m1 - m2)'*(m1 - m2) + trace(s1 + s2 - 2*sqrtm(s1*s2)) );
return
