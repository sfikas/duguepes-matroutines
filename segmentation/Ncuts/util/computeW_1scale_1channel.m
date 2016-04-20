function W=computeW_1scale_1channel(image,layer,dataW);
% compute 1 scale of multiscale image affinity matrix W for gray scale
% image
% input:
% image:pxq image
% layer:struct containing all layer-specific information such as :
% layer.p,layer.q,layer.radius,layer.scales,layer.weight,layer.location
% layer.mode2:'F' | 'IC' | 'mixed' |'hist'[not used in this function but in another]
% dataW:struct containing other parameters such as:
% dataW.sigmaI,dataW.edgeVariance
% Timothee Cour, 29-Aug-2006 07:49:15


[wi,wj] = cimgnbmap_lower([layer.p,layer.q],layer.radius,1);

if isfield(layer,'mode2') && ~isempty(layer.mode2)
    mode2=layer.mode2;
else
    mode2='mixed';    
end

if ismember(mode2,{'F','mixed'})
    sigmaI=(std(image(:)) + 1e-10 )* dataW.sigmaI;
end
if ismember(mode2,{'IC','mixed'})
    [emag,ephase]=computeEdges(image,layer.scales);
    edgeVariance=max(emag(:)) * dataW.edgeVariance/sqrt(0.5);
end

switch mode2
    case 'F'
        W=computeW_multiscale_F(layer,image,sigmaI,wi,wj);
    case 'IC'
        W=computeW_multiscale_IC(layer,emag,edgeVariance,ephase,wi,wj);
    case 'mixed'
        W=computeW_multiscale_mixed(layer,image,sigmaI,emag,edgeVariance,ephase,wi,wj);
    otherwise
        error('?');
end
W = W*layer.weight;


function [emag,ephase]=computeEdges(image,scale);
filter_par = [4,scale,30,3];  %[9,30,4]
threshold=1e-14;
[x,y,gx,gy,par,threshold,emag,ephase,g,FIe,FIo,mago] = quadedgep_optimized(image,filter_par,threshold);


function W=computeW_multiscale_F(layer,image,sigmaI,wi,wj);
options.mode='multiscale_option';
options.mode2='F';
options.F=image/sigmaI;
options.location=layer.location;
W = mex_affinity_option(wi,wj,options);

function W=computeW_multiscale_IC(layer,emag,edgeVariance,ephase,wi,wj);
options.mode='multiscale_option';
options.mode2='IC';
options.emag=emag/edgeVariance;
options.ephase=ephase;
options.isPhase=1;
options.location=layer.location;
W = mex_affinity_option(wi,wj,options);

function W=computeW_multiscale_mixed(layer,image,sigmaI,emag,edgeVariance,ephase,wi,wj);
options.mode='multiscale_option';
options.mode2='mixed';
options.F=image/sigmaI;
options.emag=emag/edgeVariance;
options.ephase=ephase;
options.isPhase=1;
options.location=layer.location;
W = mex_affinity_option(wi,wj,options);



