function W=computeW_1scale(image,layer,dataW);
% compute 1 scale of multiscale image affinity matrix W
% input:
% image: pxq or pxqxk (for example, RGB image)
% layer: parameters for current layer
% dataW: parameters for affinity matrix W
% Timothee Cour, 29-Aug-2006 07:49:15


%{
if isfield(layer,'mode2') && ~isempty(layer.mode2)
    if strcmp(layer.mode2,'hist')
        [wi,wj] = cimgnbmap_lower([layer.p,layer.q],layer.radius,1);
        W=computeW_multiscale_hist(layer,image,wi,wj);
        W = W*layer.weight;

        return;
    end
end
%}

% for each layer in image, compute corresponding partial affinity matrix
[p,q,r]=size(image);
W=computeW_1scale_1channel(image(:,:,1),layer,dataW);
for j=2:r,
    Wj = computeW_1scale_1channel (image(:,:,j),layer,dataW);
    W = min(W,Wj);
end


%{
function W=computeW_multiscale_hist(layer,image,wi,wj);
nbins=200;
sigmaHist=2;

[p,q,r]=size(image);
n=p*q;
[H,map,Fp]=computeFeatureHistogram(image,nbins,2);
Fp=reshape2(Fp,n);
mex_normalizeColumns(Fp);
Fp=Fp.*repmat(1./std(Fp),n,1);
Fp=Fp*sigmaHist;

options.mode='multiscale_hist';
options.F=Fp;
options.hist=-H;
options.map=map;
options.ephase=[];
options.isPhase=0;
options.location=layer.location;
W = mex_affinity_option(wi,wj,options);
%}
