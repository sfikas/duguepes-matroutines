function [classes,X,lambda,Xr,W,C,timing] = ncut_multiscale(image,nsegs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                  %
%   Multiscale Normalized Cuts Segmentation Code   %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% inputs: 
% image: image to segment (size pxq or pxqx3)
% nsegs: number of segments requested
% outputs:
% classes: image regions (size pxq)
% X: eigenvectors (size pxqxnsegs)
% lamda: eigenvalues
% Xr: rotated eigenvectors (computed during discretisation)
% W: multiscale affinity matrix
% C: multiscale constraint matrix
% timing: timing information
% 
% source code available at http://www.seas.upenn.edu/~timothee
% Authors: Timothee Cour, Florence Benezit, Jianbo Shi
% Related publication:
% Timothee Cour, Florence Benezit, Jianbo Shi. Spectral Segmentation with
% Multiscale Graph Decomposition. IEEE International Conference on Computer
% Vision and Pattern Recognition (CVPR), 2005.
% 
% Please cite the paper and source code if you are using it in your work.


image=im2double(image);

%% compute multiscale affinity matrix W and multiscale constraint matrix C
t= cputime;
[W,C]=compute_W_C_multiscale(image);
disp(['compute W,C time: ',num2str(cputime-t)]);

%% compute constrained normalized cuts eigenvectors
disp('starting multiscale ncut...');
t = cputime;
if ~isempty(C)
    [X,lambda,timing] = computeNcutConstraint_projection(W,C,nsegs);
else
    [X,lambda,timing] = computeKFirstEigenvectors(W,nsegs);
end
disp(['ncut time: ',num2str(cputime-t)]);
disp(['(including time for W*X operations: ',num2str(timing.A_times_X),' )']);

%% compute discretisation
[p,q,r]=size(image);
indPixels = (1:p*q)';
X = reshape(X(indPixels,:),p,q,nsegs);
t =cputime;
[classes,Xr] =discretisation(X);
disp(['discretize time: ',num2str(cputime-t)]);

