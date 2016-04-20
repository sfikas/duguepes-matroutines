function [polarity, l1, l2] = computePolarity(inimage, scale)
% Computes polarity for given MxN matrix 'inimage' which
% normally should be an L* component. Also computes l1,l2
% eigenvalue matrices (see paper).
%
% [polarity, l1, l2] = computePolarity(inimage, scale)
%
% 'polarity'    Polarity matrix
% 'l1','l2'     Matrix of principal & secondary eigenvalues
%                on gradient covariance (see paper)
%
% 'inimage'     Input. Must be L* component.
% 'scale'       _standard deviation_ of the used gaussian kernel.
%               It can be scalar (constant scale) or a matrix (variant
%               scale)
%
% G.Sfikas 11/4/2006
% update: May 3 2006
% comments: May 8 2006
%

%Is scale constant or not?
if isscalar(scale) == 0
    polarity = zeros(size(inimage,1), size(inimage,2));
    l1 = zeros(size(inimage,1), size(inimage,2));
    l2 = zeros(size(inimage,1), size(inimage,2));
    for u = 1:8
        vScale = (u-1)/2;
        [vPolarity, vL1, vL2] = computePolarity(inimage, vScale);
        vMap = (scale == vScale);
        polarity = polarity + vPolarity.*vMap;
        l1 = l1 + vL1.*vMap;
        l2 = l2 + vL2.*vMap;        
    end
    return;
end
%Compute gradients
[Lx Ly] = gradient(double(inimage));
%Compute gradients autocorellation
acXX = Lx.^2;
acXY = Lx.*Ly;
acYY = Ly.^2;
%Smooth 'em
acXX = convolution2D(acXX, scale);
acXY = convolution2D(acXY, scale);
acYY = convolution2D(acYY, scale);
%Compute eigenvalues, principal eigenvector.
% Ta ypologizw 'xeirokinita' dhl kanw prakseis apo
% |A-lambda*I| = 0 klp. (pio grhgoro)
tr = acXX + acYY;
diakrinousa(:,:) = real( 0.5 * sqrt(tr.^2-4*(acXX.*acYY-acXY.^2)));
l1 = 0.5*tr + diakrinousa;
l2 = 0.5*tr - diakrinousa;
phi2(:,:,1) = -1./(acXY+eps);
phi2(:,:,2) = 1./(acYY-l1-2*eps);
%
% Project onto principal eigenvectors.
% NOTE: Trick!!
%
% Epeidh oi provoles pou prepei na ginoun einai ypervolika polles
% (peripou cx*cy*dx*dy opou cx,cy diastaseis parathirou synelikshs
%  kai dx,dy diastaseis eikonas) kai to xeirotero oti anagkastika prepei 
% na ginoun epanalhptika (me for) 
%
clear i;
discrete_angles = 32;
minimum_angle = angle(exp(i*(pi/discrete_angles)));
evectorFakeComplex = phi2(:,:,1) + i*phi2(:,:,2);
evectorFakeAngle = mod(round(angle(evectorFakeComplex) ./ minimum_angle), discrete_angles);
% Remedy 'zero' angle.
evectorFakeAngle = evectorFakeAngle + discrete_angles*(evectorFakeAngle==0);
finalEplus = zeros(size(evectorFakeAngle,1), size(evectorFakeAngle,2));
finalEminus = zeros(size(evectorFakeAngle,1), size(evectorFakeAngle,2));
for k = 1:discrete_angles
    currAngle = exp(i*k*(pi/discrete_angles));
    proj(:,:,k) = Lx.*real(currAngle) + Ly.*imag(currAngle);
    Eplus(:,:,k) = proj(:,:,k).*(proj(:,:,k)>=0);
    Eminus(:,:,k) = Eplus(:,:,k) - proj(:,:,k);
    Eplus(:,:,k) = convolution2D(Eplus(:,:,k), scale);
    Eminus(:,:,k) = convolution2D(Eminus(:,:,k), scale);
    %
    finalEplus = finalEplus + (evectorFakeAngle==k).*Eplus(:,:,k);
    finalEminus = finalEminus + (evectorFakeAngle==k).*Eminus(:,:,k);
end     
polarity = abs(finalEplus - finalEminus)./(finalEplus+finalEminus+eps);

return;