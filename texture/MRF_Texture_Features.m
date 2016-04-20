function [X] =  MRF_Texture_Features(Im,wnd,Gauss_wnd,Gauss_sigma,pca_dim)
% Im is the original RGB or grayscale image of dimensions [DIMR DIMC DIMZ]
% Gauss_wnd is the size of the Gaussian blurring window
% Gauss_sigma is the std of the Gaussian blur
% pca_dim is the dimension of the PCA decomposition of the MRF features
% X has dimensions [DIMR DIMC DIMZ]

[DIMR DIMC DIMZ] = size(Im);

h = fspecial('gaussian',Gauss_wnd,Gauss_sigma);
Im2 = imfilter(Im,h,'replicate','full');


[DIMR DIMC DIMZ] = size(Im2);

X = zeros (DIMR,DIMC,DIMZ*wnd*wnd);
DD=zeros(DIMR,DIMC,DIMZ*wnd*wnd);
wnd2 = floor(wnd/2);

k=0;
for i=(wnd2+1):(DIMR-wnd2)
    for j=(wnd2+1):(DIMC-wnd2)
        xx = Im2((i-wnd2):(i+wnd2),(j-wnd2):(j+wnd2),:);
        k=k+1;
        DD(i,j,:) = xx(:);
    end
end

D=zeros(DIMR*DIMC,DIMZ*wnd*wnd);
for k=1:DIMZ*wnd*wnd
    yy=reshape(DD(:,:,k),DIMR,DIMC);
    D(:,k) = yy(:);
end
clear DD;

%C=cov(D);
%[Evec,Eval] = pcacov(C);
[Evec,NewCoords,Eval,tsquare] = princomp(D);
Evec2=Evec(:,1:pca_dim);


%% Visualize the Evectors %%%%
% for k=1:wnd*wnd
%     yy=Evec(:,k);
%     yy2=reshape(yy,wnd,wnd);
%     figure;imagesc(yy2);colormap gray;
% end
    

% for k=1:DIMR*DIMC
%     NewCoords(k,:) = Evec2'*( D(k,:) - mean(D))';
% end



Y = zeros (DIMR,DIMC,pca_dim);
for k=1:pca_dim
    Y(:,:,k)=reshape(NewCoords(:,k),DIMR,DIMC);
end

[DIMRnew DIMCnew DIMZ] = size(Im);


X = zeros (DIMRnew,DIMCnew,pca_dim);
for k=1:pca_dim
    X(:,:,k) = Y((wnd2+1):(DIMR-wnd2),(wnd2+1):(DIMC-wnd2),k);
end


clear D;
clear NewCoords;
clear Evec;
clear Eval;
clear xx;
clear yy;

