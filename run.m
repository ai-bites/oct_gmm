% Before running, make sure the data is available in the data folder
% The DME volumes are expected to be in /data/DME folder.
% The normal volumes are expected to be in /data/normal folder.

clc; clear all; close all;
dataset = 'SERI';

%% Initalize - some settings and variables to run the algorithm

do_plot = 0;
% should we do preprocessing to arrive at X matrix. 
% if not, we use the saved X - X_normal.mat and X_patient.mat
do_preprocess = 1;
% should we use the saved X values - X_normal.mat and X_patient.mat?
saved_results = 0;
% Varies based on the dataset
num_vols = 16;
frames_per_vol = 128; % for SERI
% frames_per_vol = 97;    % for duke dataset
total_frames = num_vols*frames_per_vol;
% K for PCA dimensions
pca_k = 900; % for SERI
% pca_k = 500;   % for duke dataset
% K for the GMM model
gmm_k = 5;
% Threshold for the mahalnobis distance to separate the outliers.
% the more we set, the lesser outliers (diseased) frames.
mahal_thresh = 1200; % for our dataset
% mahal_thresh = 900;    % for duke dataset
% Number of frames to train with. Used for cross validation
train_frames = 12*frames_per_vol; 
% train_frames = 1019;


%% Read data and pre-process them for flattening

if (do_preprocess)
    [X,X_fvols] = preprocess('data/normal/normal',num_vols, dataset);
    [Y,Y_fvols] = preprocess('data/DME/patient',num_vols, dataset);
end

% save('X_patient_1.mat','Y','-v7.3');
% save('X_normal_1.mat','X','-v7.3');

%% Do PCA to lower dimension

if (saved_results)
    X = load('X_normal');
    X = X.X;
end
X = X';
[X_res, x_mapping] = do_pca(X,pca_k);
disp('Done with loading X - normal cases');

% Choose K based on this ratio. Above 0.9 is good
%sum(x_mapping.lambda(1:pca_k))/sum(x_mapping.lambda);
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

obj = gmdistribution.fit(X_res,gmm_k, ...
                        'Options',options, ...
                        'Regularize',0.001, ...
                        'Start', 'randSample', ...
                        'CovType', 'diagonal');


%% Cross validation using volumens 13-16


outliers = {};
x_outliers = [];
n = 1; % count number of outliers
frame_no = 1; % frame number in the volume
vol = 1;
for i = train_frames+1:total_frames
   %frame_dists(n) = norm(X_res(i,:) - gmm_means);
   [x_outliers(n) min_mahal(n)] = check_outlier(...
                                  obj, X_res(i,:), gmm_k, mahal_thresh);
   n = n+1;
   frame_no = frame_no+1;
    
    % reset for next volume
    if (frame_no > frames_per_vol) 
        outliers{vol} = x_outliers;
        frame_no = 1;
        vol = vol+1;
        n = 1;
        x_outliers = [];
    end
end

% Analyse the outliers obtained
diseased_frames = 0;
diseased_vols = 0;
for i = 1:length(outliers)
    if (sum(outliers{i} > 0))
       diseased_vols = diseased_vols+1;
       diseased_frames = diseased_frames + sum(outliers{i});
       fprintf('diseased vol number is: %d \n', i);
    end
end
fprintf('Diseased volumes are: %d, total diseased frames are: %d \n',...
        diseased_vols, diseased_frames);

%% Classification of diseased volumes
% To run this step, run the GMM fit step for all volumes.

outliers = {};
y_outliers = [];
n = 1; % count number of outliers
frame_no = 1; % frame number in the volume
vol = 1;
for i = 1:total_frames
    [y_outliers(n)] = check_outlier(...
                                   obj, Y_res(i,:),gmm_k, mahal_thresh);
    n = n+1;
    frame_no = frame_no+1;
    % reset for next volume
    if (frame_no > frames_per_vol) 
        outliers{vol} = y_outliers;
        frame_no = 1;
        vol = vol+1;
        n = 1;
        y_outliers = [];
    end
end

% Analyse the outliers obtained
diseased_frames = 0;
diseased_vols = 0;
for i = 1:length(outliers)
    if (sum(outliers{i} > 0))
       diseased_vols = diseased_vols+1;
       diseased_frames = diseased_frames + sum(outliers{i});
       fprintf('diseased vol number is: %d \n', i);
    end
end
fprintf('Diseased volumes are: %d, total diseased frames are: %d \n',...
        diseased_vols, diseased_frames);

