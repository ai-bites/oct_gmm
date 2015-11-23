
clc; clear all; close all;
dataset = 'DUKE';

%%

do_plot = 0;
% should we do preprocessing to arrive at X matrix. 
% if not, we use the saved X - X_normal.mat and X_patient.mat
do_preprocess = 0;
% should we use the saved X values - X_normal.mat and X_patient.mat?
saved_results = 0;
% Varies based on the dataset
num_vols = 15;
% frames_per_vol = 97;
% total_frames = num_vols*frames_per_vol;
% K for PCA dimensions
pca_k = 300;
% Cross validation using volumens 13-15
% try different K values for GMM to decide which is best
GMM_Ks = [2 3 5 8 10 12 15];
% Threshold for the mahalnobis distance to separate the outliers.
% the more we set, the lesser outliers (diseased) frames.
mahal_thresh = chi2inv(0.95, pca_k);
% Number of frames to train with. Used for cross validation
train_frames = 1116;
% If there are diseased frames are more than this threshold
% then the volume is diseased
% diseased_vol_thresh = 10;

%% Read data and pre-process them for flattening

if (do_preprocess)
    disp('doing pre-processing');
    X_vols = preprocess('duke_data/normal/normal',num_vols, dataset, 'intensity');
    Y_vols = preprocess('duke_data/DME/patient',num_vols, dataset, 'intensity');
else
    % load already pre-processed and stored data
    X_vols = load('X_duke_lbp');
    X_vols = X_vols.X_vols;
    Y_vols = load('Y_duke_lbp');
    Y_vols = Y_vols.Y_vols;
end


%% Choose the best model parameter
% In case of GMM, it is to find the best K value for GMM

clc;
num_iter = 1;
[avg_diseased, used_vols, x_test_vols] = train_gmm(X_vols, pca_k, GMM_Ks, 'intensity',num_iter, dataset);

% Using the best average values, choose the best K value for final GMM
fprintf('Best K for GMM is: %d \n', GMM_Ks(find(avg_diseased == min(avg_diseased)))); 
best_gmm_k = max(GMM_Ks(find(avg_diseased == min(avg_diseased))));


%% Fit the final model with all the used volumes

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

%% Now test with the obtained model 

test_gmm(obj, x_test_vols, Y_vols, ...
            dataset, pca_k, best_gmm_k, mahal_thresh, diseased_vol_thresh, 'intensity');


%%
% 
% %% Do PCA to lower dimension
% 
% X = X';
% [X_res, x_mapping] = do_pca(X,pca_k);
% % Choose K based on this ratio. Above 0.9 is good
% fprintf('Variance retained in X after PCA is: %f \n', ...
%     sum(x_mapping.lambda(1:pca_k))/sum(x_mapping.lambda));
% 
% Y = Y';
% [Y_res, y_mapping] = do_pca(Y,pca_k);
% disp('done with PCA on Y - DME cases');
% fprintf('Variance retained in Y after PCA is: %f \n', ...
%     sum(y_mapping.lambda(1:pca_k))/sum(y_mapping.lambda));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Now fit a GMM with the data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clc; clear total_diseased;
% % either all of X or first 12 volumes
% train_data = X_res(1:train_frames,:);
% options = statset('Display','off','MaxIter',1000);
% 
% for ii = 1:length(GMM_Ks)
%     gmm_k = GMM_Ks(ii);
%     obj = gmdistribution.fit(train_data,gmm_k, ...
%                         'Options',options, ...
%                         'Regularize',0.001, ...
%                         'Start', 'randSample', ...
%                         'CovType', 'diagonal', ...
%                         'Replicates', 5, ...
%                         'Start', 'randSample');
%     
%     outliers = {};
%     x_outliers = [];
%     n = 1; % count number of outliers
%     frame_no = 1; % frame number in the volume
%     vol = 1;
%     for i = train_frames+1:1407
%         [x_outliers(n)] = check_outlier(...
%                                  obj, X_res(i,:), gmm_k, mahal_thresh);
%         n = n+1;
%         frame_no = frame_no+1;
%         [a b frames_in_vol] = size(X_fvols{vol});
%         % reset for next volume
%         if (frame_no > frames_in_vol)
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
%     for j = 1:length(outliers)
%         if (sum(outliers{j}) > diseased_vol_thresh)
%            diseased_vols = diseased_vols+1;
%            diseased_frames = diseased_frames + sum(outliers{j});
%            fprintf('diseased vol number is: %d \n', j);
%         end
%     end
%     
%     fprintf('For GMM k value: %d \n', gmm_k);
%     fprintf('Diseased volumes are: %d, total diseased frames are: %d \n',...
%             diseased_vols, diseased_frames);
% 
%     total_diseased(ii) = diseased_frames;
% 
% end
% 
% % Choose the best K value for GMM. 
% % We choose the K that returns the least number of diseased frames
% % If more than 2 Ks are best, chose the  min.
% best_gmm_k = min(GMM_Ks(find(total_diseased == min(total_diseased))));
% fprintf('BEST K value for GMM is: %d \n', best_gmm_k);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Classification of diseased volumes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clc;
% % final model with the best K value for GMM
% % train with all training data
% options = statset('Display','final','MaxIter',1000);
% obj = gmdistribution.fit(X_res,best_gmm_k, ...
%                             'Options',options, ...
%                             'Regularize',0.001, ...
%                             'Start', 'randSample', ...
%                             'CovType', 'diagonal', ...
%                             'Replicates', 5, ...
%                             'Start', 'randSample');
% outliers = {};
% y_outliers = [];
% n = 1; % count number of outliers
% frame_no = 1; % frame number in the volume
% vol = 1;
% 
% % for 7 volumes, 679
% for i = 1:679
%     [y_outliers(n)] = check_outlier(...
%                           obj, Y_res(i,:),gmm_k, mahal_thresh);
%     n = n+1;
%     frame_no = frame_no+1;
%     [a b frames_in_vol] = size(Y_fvols{vol});
%     % reset for next volume
%     if (frame_no > frames_in_vol) 
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
%        % fprintf('diseased vol number is: %d \n', i);
%     end
% end
% fprintf('Diseased volumes are: %d, total diseased frames are: %d \n',...
%         diseased_vols, diseased_frames);
