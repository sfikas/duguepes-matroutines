function [options,nbEigenvectors] = getDefaultOptionsEigs(n,nbEigenvectors);
% Timothee Cour, 29-Aug-2006 07:49:15

options.issym = 1;
options.maxit = 100;

options.tol = 1e-3;% 1e-3; 1e-4; 1e-2

options.v0 = ones(n,1);
%options.v0 = rand(n,1);
options.p = min(3*nbEigenvectors+3,n);%TODO: more ? (since arpack is bottleneck)
% options.p = min(2*nbEigenvectors+3,n);%TODO: more ? (since arpack is bottleneck)
% options.p = min(2*nbEigenvectors,n);
% voir : 2.5*nbEigenvectors ??
%options.p = min(8,n);

options.fastMode = 1;
options.computeX = 1;
options.warningConvergence = 1;
options.sigma = 'LA';
options.n = n;

nbEigenvectors = min(nbEigenvectors , options.p-1);
