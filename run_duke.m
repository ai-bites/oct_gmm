
clc; clear all; close all;
dataset = 'DUKE';

%%

do_plot = 0;
% should we do preprocessing to arrive at X matrix. 
% if not, we use the saved X - X_normal.mat and X_patient.mat
do_preprocess = 1;
% should we use the saved X values - X_normal.mat and X_patient.mat?
saved_results = 0;
% Varies based on the dataset
num_vols = 15;
% frames_per_vol = 97;
% total_frames = num_vols*frames_per_vol;
% K for PCA dimensions
pca_k = 500;
% K for the GMM model
gmm_k = 10;
% Threshold for the mahalnobis distance to separate the outliers.
% the more we set, the lesser outliers (diseased) frames.
mahal_thresh = 1100;
% Number of frames to train with. Used for cross validation
train_frames = 1019;

%% Read data and pre-process them for flattening

if (do_preprocess)
    [X,X_fvols] = preprocess('duke_data/normal/normal',num_vols, dataset);
    [Y,Y_fvols] = preprocess('duke_data/DME/patient',7, dataset);
end


%% Do PCA to lower dimension

X = X';
[X_res, x_mapping] = do_pca(X,pca_k);
disp('Done with PCA on X - normal cases');

% Choose K based on this ratio. Above 0.9 is good
%sum(x_mapping.lambda(1:pca_k))/sum(x_mapping.lambda);

Y = Y';
[Y_res, y_mapping] = do_pca(Y,pca_k);
disp('done with PCA on Y - DME cases');

%% Now fit a GMM with the data

options = statset('Display','final','MaxIter',1000);

% either all of X or first 12 volumes
%train_data = X_res(1:train_frames,:);
train_data = Y_res;
obj = gmdistribution.fit(train_data,gmm_k, ...
                        'Options',options, ...
                        'Regularize',0.001, ...
                        'Start', 'randSample', ...
                        'CovType', 'diagonal');


%% Cross validation using volumens 13-15

mahal_thresh = [600 800 900 1000 1100 1200 1300];

for ii = 1:length(mahal_thresh)
    outliers = {};
    x_outliers = [];
    n = 1; % count number of outliers
    frame_no = 1; % frame number in the volume
    vol = 1;
    for i = train_frames+1:1407
        [x_outliers(n) min_mahal(n)] = check_outlier(...
                                      obj, X_res(i,:), gmm_k, mahal_thresh(ii));
        n = n+1;
        frame_no = frame_no+1;
        [a b frames_in_vol] = size(X_fvols{vol});
        % reset for next volume
        if (frame_no > frames_in_vol)
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
    for j = 1:length(outliers)
        if (sum(outliers{j} > 0))
           diseased_vols = diseased_vols+1;
           diseased_frames = diseased_frames + sum(outliers{j});
           fprintf('diseased vol number is: %d \n', j);
        end
    end
    fprintf('Diseased volumes are: %d, total diseased frames are: %d \n',...
            diseased_vols, diseased_frames);

end % try for different mahal distances


%% Classification of diseased volumes

mahal_thresh = [600 625 650 675 700 725 750 775 800];

for ii = 1:length(mahal_thresh)
    outliers = {};
    y_outliers = [];
    n = 1; % count number of outliers
    frame_no = 1; % frame number in the volume
    vol = 1;

    % for 7 volumes, 679
    for i = 1:679
        [y_outliers(n)] = check_outlier(...
                                       obj, Y_res(i,:),gmm_k, mahal_thresh(ii));
        n = n+1;
        frame_no = frame_no+1;
        [a b frames_in_vol] = size(Y_fvols{vol});
        % reset for next volume
        if (frame_no > frames_in_vol) 
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
           % fprintf('diseased vol number is: %d \n', i);
        end
    end
    fprintf('Diseased volumes are: %d, total diseased frames are: %d \n',...
            diseased_vols, diseased_frames);

end % try for different mahal distances

