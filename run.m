% Before running, make sure the data is available in the data folder
% The DME volumes are expected to be in /data/DME folder.
% The normal volumes are expected to be in /data/normal folder.

clc; clear all; close all;

%% Initalize - some settings and variables to run the algorithm

do_plot = 0;
% should we do preprocessing to arrive at X matrix. 
% if not, we use the saved X - X_normal.mat and X_patient.mat
do_preprocess = 0;
% should we use the saved X values - X_normal.mat and X_patient.mat?
saved_results = 1;
% Varies based on the dataset
num_vols = 16;
num_frames = num_vols*128;
% K for PCA dimensions
pca_k = 900;
% K for the GMM model
gmm_k = 10;
% Threshold for the mahalnobis distance to separate the outliers.
% the more we set, the lesser outliers (diseased) frames.
mahal_thresh = 1100;
% Number of frames to train with. Used for cross validation
train_frames = 12*128; 


%% Read data and pre-process them for flattening

if (do_preprocess)
    X = preprocess('data/DME/patient',2);
end


%% Do PCA to lower dimension

if (saved_results)
    X = load('X_normal');
    X = X.X;
end
X = X';
[X_res, x_mapping] = do_pca(X,pca_k);
disp('Done with loading X - normal cases');

% Choose K based on this ratio. Above 0.9 is good
%sum(x_mapping.lambda(1:K))/sum(x_mapping.lambda);
%plot3(X_res(:,1),X_res(:,2),X_res(:,3),'rx');

if (saved_results)
    Y = load('X_patient');
    Y = Y.X;
end
Y = Y';
[Y_res, y_mapping] = do_pca(Y,pca_k);
disp('done with loading Y - DME cases');


%% Now fit a GMM with the data

options = statset('Display','final','MaxIter',1000);

obj = gmdistribution.fit(X_res(1:train_frames,:),gmm_k, ...
                        'Options',options, ...
                        'Regularize',0.001, ...
                        'Start', 'randSample', ...
                        'CovType', 'diagonal');


%% Cross validation using volumens 13-16

clear x_outliers min_mahal;
n = 1;
for i = train_frames+1:num_frames
   %frame_dists(n) = norm(X_res(i,:) - gmm_means);
   [x_outliers(n) min_mahal(n)] = check_outlier(...
                                  obj, X_res(i,:), gmm_k, mahal_thresh);
   n = n+1;
end

fprintf('num outliers for X val data (vols 13-16) is: %d out of %d \n', ...
        sum(x_outliers), length(x_outliers));


%% Classification of diseased volumes
% To run this step, run the GMM fit step for all volumes.

clear y_outliers min_mahal;
n = 1;
for i = 1:num_frames
    [y_outliers(n),min_mahal(n)] = check_outlier(...
                                   obj, Y_res(i,:),gmm_k, mahal_thresh);
    n = n+1;
end

fprintf('Diseased frames are: %d out of %d \n', sum(y_outliers), length(y_outliers));






