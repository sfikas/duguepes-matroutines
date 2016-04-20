function [C,C12]=computeMultiscaleConstraint_1scale(p1,q1,p2,q2,indexes1,indexes2,nTot);
% input: parameters for 2 consecutive layers
% output:
% C: rows of the global multiscale constraint matrix corresponding to those 2 layers; 
% indexes1 and indexes2 are used for indexing into the global multiscale
% constraint matrix
% C12: interpolation matrix between the 2 layers
% Timothee Cour, 29-Aug-2006 07:49:15


classes=mex_constraint_classes(p1,q1,p2,q2);
[C,C12]=computeConstraintFromClasses(classes,indexes1,indexes2,nTot);
