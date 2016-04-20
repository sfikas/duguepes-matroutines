function Ncut_Segmentation(FullPathName,ImageFileName,NB_CLUSTERS,bTextureFeatures,ImResizeFactor)
%%% filename is the file containing the image
%% NB_CLUSTERS is a vector containing the various number of segments
%% bTextureFeatures indicates if MRF texture features or RGB features are to be considered

image = imread(FullPathName);


image = imresize(image,ImResizeFactor);
figure;imshow(image,[]);

if bTextureFeatures
    disp('Using MRF textures');
    %%% MRF texture features %%%
    wnd = 7;
    Gauss_wnd = 7;
    Gauss_sigma = 2.0;
    pca_dim =8;
    ImF = MRF_Texture_Features(image,wnd,Gauss_wnd,Gauss_sigma,pca_dim);
else
    disp('Using normal features');
    %%%%%%%%   RGB or Lab features   %%%%%
     %ImF = RGB2Lab(double(image));
     ImF = double(image);
end

%image=rgb2gray(image);
[DIMR,DIMC,DIMZ]=size(ImF);


%%%% Create matrix D of dim (DIMRxDIMC) x DIMZ 
%%%% whose rows correspond to pixel features 

D=[];
for k=1:DIMZ
    It=ImF(:,:,k);
    D=[D;It(:)'];
end
D=D';


for k=1:length(NB_CLUSTERS)

    close all;
    
    NbClusters = NB_CLUSTERS(k);
    
    disp(['image size : ',mat2str([DIMR,DIMC,DIMZ]) ]);

    disp('Ncut Segmentation...');


   
    %[IDX,Centers]=DCM_GMM(D,ImF,NbClusters,bSemanticSegmentation);
    [IDX,X,lambda,Xr,W,C,timing] = ncut_multiscale(ImF,NbClusters);
    %[IDX,Centers] = kmeans(D,NbClusters,'maxiter',30,'emptyaction','singleton','start','cluster','replicates',5);


    %%% Show the Segmented Image
    SegMap = reshape(IDX,DIMR,DIMC);
    figure;imshow(SegMap,[]);

    fnImg = strtok(ImageFileName,'.');
    
    fn1 = 'Results\Ncut_MRF\Neighborhood_21x21\Segments';
    
    
    fn2=sprintf('%s%d\',fn1,NbClusters);
    
    if ~isdir(fn2)
        mkdir(fn2);
    end
    
    fnImgOut = sprintf('%s\\%s.%s',fn2,fnImg,'jpg');
    imwrite(mat2gray(SegMap),fnImgOut);

    fnMatOut = sprintf('%s\\%s.mat',fn2,fnImg);
    save(fnMatOut,'SegMap');

    fprintf('Segmentation of image %s into %d Clusters by Ncut finished\n\n\n',ImageFileName,NbClusters);



end %% for NbClusters%%


            
       





