function [Upzet] = BIDProjection(xi)
%
% Upzet = BIDProjection(xi)
%
% Project vector 'xi' onto space of sum(x) = 1.
%
%
% Note:
%  The original version crashed on this input
%   xi(1) = 2.1196e-027;
%   xi(2) = 1.1342e+016;
%   xi(3) = 0;
% most probably due to machine numerical inacurracy.
%
% By K Blekas plus a minor change
% G.Sfikas 8 Feb 2008
% Fix 20 Mar 2008: sur clumsy fix: il etait 'if ggm > 0 && sum(xi)-ggM == 0'
%                   qui ne marche pas pour non-positives vectors

K= size(xi,2);
% The minor (clumsy) fix:
[ggM ggExcl] = max(xi');
if ggM > 0 && sum(abs(xi)) - ggM == 0
    Upzet = zeros(1, K);
    Upzet(ggExcl) = 1;
    return;
end
%

mes=sum(xi)/K;
            
zet = 1/K + xi - mes;
            
N0 = size(find(zet<0),2) + size(find(zet>1),2);
            
if N0<1
    Upzet = zet;
else
    Eta = ones(1,K);
    notfound = 1;
    VertDist = 1 + sum(xi.^2) - 2 .* xi;
    c=2;
   
    nzi = xi;
    
    while (notfound > 0)
        
            [maxv vtx] = max(VertDist);
            VertDist(vtx) = -1.0;
            
            Eta(vtx) = 0;
            
            tzi = Eta .* nzi;
            
            %for k=1:K
            %    tzi(k) = Eta(k)*nzi(k);
            %end
            
            Tempzet = tzi + (1 - sum(tzi))/(K-c+1);
            
            %for k=1:K
            %    Tempzet(k) = Tempzet(k) * Eta(k);
            %end
            Tempzet = Tempzet .* Eta;
            
            N0 = size(find(Tempzet<0),2) + size(find(Tempzet>1),2);
            
            if N0<1
                Upzet = Tempzet;
                notfound=0;
            end
            c = c + 1;
    end
    
end