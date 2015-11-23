% Before running, make sure the data is available in the data folder
% The DME volumes are expected to be in /data/DME folder.
% The normal volumes are expected to be in /data/normal folder.

clc; clear all; close all;
dataset = 'SERI';

%% Initalize - some settings and variables to run the algorithm

do_plot = 0;
% should we do preprocessing to arrive at X matrix. 
% if not, we use the saved X - X_normal.mat and X_patient.mat
do_preprocess = 0;
% should we use the saved X values - X_normal.mat and X_patient.mat?
saved_results = 1;
% Varies based on the dataset
num_vols = 1;
frames_per_vol = 128; 
total_frames = num_vols*frames_per_vol;
% K for PCA dimensions
pca_k = 300;
% Cross validation using volumens 13-16
% try different K values for GMM to decide which is best
GMM_Ks = [2 3 5 8 10 12 15];
% Threshold for the mahalnobis distance to separate the outliers.
% the more we set, the lesser outliers (diseased) frames.
mahal_thresh = chi2inv(0.95, pca_k); % for our dataset
% Number of frames to train with. Used for cross validation
train_frames = 12*frames_per_vol; 
% Number of iterations to run for each GMM K value to 
% decide which K is best. The higher, the longer wait!
num_iter = 1;


%% Read data and pre-process them for flattening

if (do_preprocess)
    [X_vols] = preprocess('data/normal/normal',num_vols, dataset, 'intensity');
    %[Y_vols] = preprocess('data/DME/patient',num_vols, dataset, 'intensity');
else
    % load already pre-processed and stored data
    X_vols = load('X_seri_lbp');
    X_vols = X_vols.X_vols;
    Y_vols = load('Y_seri_lbp');
    Y_vols = Y_vols.Y_vols;
end

%% train using the preprocessed volumes
% to find the best K value for GMM

clc;
[avg_diseased, used_vols, x_test_vols] = train_gmm(X_vols, pca_k, GMM_Ks, 'intensity', num_iter, dataset);

% Using the best average values, choose the best K value for final GMM
fprintf('Best K for GMM is: %d \n', GMM_Ks(find(avg_diseased == min(avg_diseased)))); 
best_gmm_k = GMM_Ks(find(avg_diseased == min(avg_diseased)));


%% Fit the final model with all used volumes

X = get_features(used_vols, 'intensity', dataset);
[X_train, X_mapping] = do_pca(X',pca_k);

options = statset('Display','final','MaxIter',1000);
% final model with the best K value for GMM
% train with all training data
obj = gmdistribution.fit(X_train,best_gmm_k, ...
                            'Options',options, ...
                            'Regularize',0.001, ...
                            'Start', 'randSample', ...
                            'CovType', 'diagonal', ...
                            'Replicates', 5, ...
                            'Start', 'randSample');

diseased_vol_thresh = get_disease_thresh(obj, X_train, used_vols, best_gmm_k, mahal_thresh);

% test on the normal and diseased test data
[dis_outliers, norm_outliers] = test_gmm(obj, x_test_vols, Y_vols, ...
            dataset, pca_k, best_gmm_k, mahal_thresh, diseased_vol_thresh, 'intensity');


% Method 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do PCA to lower dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if (saved_results)
%     X = load('X_normal_1');
%     X = X.X;
%     disp('Done with loading X - normal cases');
% end
% [X_res, x_mapping] = do_pca(X',pca_k);
% % Choose K based on this ratio. Above 0.9 is good
% fprintf('Variance retained after PCA on X is: %f \n',sum(x_mapping.lambda(1:pca_k))/sum(x_mapping.lambda));
% 
% if (saved_results)
%     Y = load('X_patient_1');
%     Y = Y.Y; 
%     disp('done with loading Y - DME cases');
% end
% [Y_res, y_mapping] = do_pca(Y',pca_k);
% fprintf('Variance retained after PCA on Y is: %f \n',sum(y_mapping.lambda(1:pca_k))/sum(y_mapping.lambda));
% 
% 
% %% Now fit a GMM with the data
% % Cross validate for different values of GMM K values
% % and choose which is best 
% 
% clc;
% options = statset('Display','off','MaxIter',1000);
% 
% for ii = 1:length(GMM_Ks)
%     gmm_k = GMM_Ks(ii);
%     
%     obj = gmdistribution.fit(X_res(1:train_frames,:),gmm_k, ...
%                             'Options',options, ...
%                             'Regularize',0.001, ...
%                             'Start', 'randSample', ...
%                             'CovType', 'diagonal', ...
%                             'Replicates', 5, ...
%                             'Start', 'randSample');
%     outliers = {};
%     x_outliers = [];
%     n = 1; % count number of outliers
%     frame_no = 1; % frame number in the volume
%     vol = 1;
%     for i = train_frames+1:total_frames
%        %frame_dists(n) = norm(X_res(i,:) - gmm_means);
%        [x_outliers(n)] = check_outlier(obj, X_res(i,:), gmm_k, mahal_thresh);
%        n = n+1;
%        frame_no = frame_no+1;
% 
%         % reset for next volume
%         if (frame_no > frames_per_vol) 
%             outliers{vol} = x_outliers;
%             frame_no = 1;
%             vol = vol+1;
%             n = 1;
%             x_outliers = [];
%         end
%     end
% 
%     % Analyse the outliers obtained
%     diseased_frames = 0;
%     diseased_vols = 0;
%     for i = 1:length(outliers)
%         fprintf('Outliers in Vol %d are: %d \n', i, sum(outliers{i}));
%         if (sum(outliers{i}) > diseased_vol_thresh)
%            diseased_vols = diseased_vols+1;
%            diseased_frames = diseased_frames + sum(outliers{i});
%            %fprintf('diseased vol number is: %d \n', i);
%         end
%     end
%     
%     fprintf('For GMM k value: %d \n', gmm_k);
%     fprintf('Total diseased volumes are: %d, total diseased frames are: %d \n',...
%             diseased_vols, diseased_frames);
%     
%     total_diseased(ii) = diseased_frames;
% end
% 
% % Choose the best K value for GMM. 
% % We choose the K that returns the least number of diseased frames
% % If more than 2 Ks are best, chose the  max.
% best_gmm_k = max(GMM_Ks(find(total_diseased == min(total_diseased))));
% 
% 
% %% Classification of diseased volumes
% % To run this step, run the GMM fit step for all volumes.
% 
% clc;
% best_gmm_k = 8;
% outliers = {};
% y_outliers = [];
% n = 1; % count number of outliers
% frame_no = 1; % frame number in the volume
% vol = 1;
% options = statset('Display','final','MaxIter',1000);
% 
% % final model with the best K value for GMM
% % train with all training data
% obj = gmdistribution.fit(X_res,best_gmm_k, ...
%                             'Options',options, ...
%                             'Regularize',0.001, ...
%                             'Start', 'randSample', ...
%                             'CovType', 'diagonal', ...
%                             'Replicates', 5, ...
%                             'Start', 'randSample');
% 
% for i = 1:total_frames
%     [y_outliers(n)] = check_outlier(obj, Y_res(i,:),gmm_k, mahal_thresh);
%     n = n+1;
%     frame_no = frame_no+1;
%     % reset for next volume
%     if (frame_no > frames_per_vol)
%         outliers{vol} = y_outliers;
%         frame_no = 1;
%         vol = vol+1;
%         n = 1;
%         y_outliers = [];
%     end
% end
% 
% % Analyse the outliers obtained
% diseased_frames = 0;
% diseased_vols = 0;
% for i = 1:length(outliers)
%     if (sum(outliers{i}) > diseased_vol_thresh)
%        diseased_vols = diseased_vols+1;
%        diseased_frames = diseased_frames + sum(outliers{i});
%        fprintf('diseased vol number is: %d \n', i);
%     end
% end
% fprintf('Diseased volumes are: %d, total diseased frames are: %d \n',...
%         diseased_vols, diseased_frames);
% 
