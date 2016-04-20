%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contents in G.Sfikas library ('sfikasLibrary')  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Update I: 24 Feb 2009
% Last update : 20 Jul 2009
%
% Not included in this list:
% /edgemap                  (Martin)
% /edgemap/superpixels      (Yi Ma)
% /lightspeed               (Minka)
% /AnalyzeToolbox           (Medical imaging toolbox)
% /mixtureLearning/nombre   (Ipse, u.c.)
%
%
% =========================================================================
% /                             General functions
% =========================================================================
%
% buildRetrievalIndex           Creates a list of all files of a given
%                               extension contained on a given folder 
%                               (including its subfolders).
% deterministicKmeans           A set of centroids is returned for the
%                               given dataset, found using k-means. The
%                               initialization depends on the data, hence
%                               the name 'deterministic'.
% multistartKmeans              Use several random initializations for
%                                k-means.
% gaussianValue                 Computes the value of a normal distribution
%                               for a given datum or set of data.
% studentValue                  Computes the value of a Student-t
%                               distribution for a given datum of set of
%                               data.
% logGaussianValue              Compute the log of "gaussianValue" (more
%                                stable)
% logStudentValue               Compute the log of "studentValue" (more
%                                stable)
% mahalanobis                   Compute the mahalanobis norm for a set of
%                               vectors. i.e X'*inv(A)*X (note the inv on
%                               the A)
% squaredist                    Compute X'A*X, like 'mahanalobis'.
% makeMovie                     Show the contents of a 3d image (eg, an
%                               MRI) as an AVI movie file.
% model2image & seg2image       Change the extension and path of a given file.
% lab2rgb & rgb2lab             Convert (X,Y,3) matrices from lab to rgb
%                                   and from rgb to lab.
% xrgb2lab                      An old version of rgb2lab. May be required
%                                   by some old code.
% randGmm                       Samples from a GMM.
% imnoiseSNR                    Add noise to given signal.
%                                   Noise strengh is entered either
%                                   in decibels or noise variance.
% medoid                    	Compute medoid out of set of vectors.
% imRAG                         Compute adjacency graph for a K-class
%                                   image.
%
% =========================================================================
% /matrixManipulation
% =========================================================================
%
% convertJxN                    Converts a (X,Y,J)-sized to a (J,X*Y)
%                               matrix. Useful if you want to pass data to
%                               some training algorithm.
% convolution2D                 Convolutes a (X,Y)-sized matrix with a 2d
%                               kernel.
% maxVote                       Use a maximum-vote 3x3 filter on input
%                                   segmentation.
% smoothUsingVariantScale       Smooth an image using gaussian kernels of
%                               spatially variant scale.
% translation                   Translate (move) a 2d matrix by a given
%                                   offset. 
%
% =========================================================================
% /mixtureLearning
% =========================================================================
%
% gaussianMixEmFit              Learn a Gaussian MM.
% studentMixEmFit               Learn a Student-t MM.
% gaussianMixGreedyEmFit        Learn a Gaussian MM using Greedy EM.
% studentMixGreedyEmFit         Learn a Student-t MM using Greedy EM.
% VARIATIONAL/                  Variation methodology applications
%   gaussianMixBayesian         Learn a Gaussian MM with priors.
%   studentMixBayesian          Learn a Student-t MM with priors on all
%                               parameters except for the degrees of 
%                               freedom.
%   studentMixBayesianXP        Learn a Student-t MM with priors on all
%                               parameters except for the degrees of
%                               freedom _and_ the weights.
% MARKOV/
%   gaussianMixBayesianContinuousLp
%                               Learn a model with continuous line process.
%                                (CVPR08 proposal)
%   gaussianMixBayesianLp       Learn a model with discrete line process.
%                                (MICCAI08 proposal)
%   gaussianMixDCASV            Learn a model with class- and directional- 
%                                adaptive priors (Nikou07 TIP paper)
%
% NOMBRE/
%   gaussianMixNombre
%                               Learn a model which uses a spatial MRF
%                               and can find the number of classes
%                               automatically
%   gaussianMixNombre2          Same as "gaussianMixNombre" - but without
%                               automatic number of classes selection
%                               (ie almost like "gaussianMixContinuousLp")
%   gaussianMixNombre3          Same as "gaussianMixNombre" - but without
%                               an MRF incorporated 
%                               (ie almost like "gaussianMixBayesian")
% =========================================================================
% /pdfDistances
% =========================================================================
%
% bhgmmDistance                 Bhatacharryya-based distance for GMMs.
% emdDistance                   Earth movers distance.
% kullbackDistance              Symmetric kullback-liebler distance.
% l2Distance                    L2 distance for GMMs. [Sfikas04]
% mahalanobisDistance           Quadratic distance for blobworld region
%                               descriptors.
%
% =========================================================================
% /segmentation
% =========================================================================
%
% buildSegmentation             Builds a segmentation for given 2D image,
%                                using a variety of methods.
% BoundaryDetectionError
% GlobalConsistencyError
% probabilisticRandIndex
% RandIndex                     Fast version of 'probabilisticRandIndex'.
% VariationOfInformation
%
% =========================================================================
% /texture
% =========================================================================
%
% MRF_texture_features              Compute MRF texture feature vectors,
%                                       by default 8-variate.
% computeBlobworldFeatureVectors    Compute Blobworld feature vectors,
%                                       ie smooth Lab, Pol-Ani-Con
%                                       and x-y. (8-variate).
%
% =========================================================================
% /mexRoutines
% =========================================================================
%
% BIDProjection                 Project input vector 'x' onto space 
%                                sum(x) = 1 && x > 0.
%                               Currently uses method presented in
%                                   Sfikas et al [MLSP workshop 09]
%                               Formerly the one in
%                                   Blekas et al [TIP 05]
% xConjugateProjection          Deprecated. A failed attempt for
%                                   a better projection, for the same
%                                   problem treated by BIDProjection.
%
%