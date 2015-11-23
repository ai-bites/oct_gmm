function [avg_disesaed, used_vols, test_vols] = train_gmm(volumes, pca_k, GMM_Ks, method, num_iterations, dataset)
% Trains different GMM models and finds the best GMM model suited for the
% provided pre-processed volumes.

    num_vols = length(volumes);
    num_used = 11;
    mahal_thresh = chi2inv(0.95, pca_k);
    options = statset('Display','off','MaxIter',1000);
    
    % try different K values for GMM to arrive at the best one
    for ii = 1:length(GMM_Ks)
        % choose 11 volumes at random from all volumes
        n = randperm(num_vols);
        for i = 1:num_used
            used_vols{i} = volumes{n(i)};
        end
        a = 1;
        for uu = num_used+1:num_vols
            test_vols{a} = volumes{n(uu)};
            a = a+1;
        end
        
        for iter = 1:num_iterations
            % Choose 8 vols randomly for GMM fitting
            nn = randperm(num_used);
            for j = 1:8
                train_vols{j} = used_vols{nn(j)};
            end
            for jj = 1:3
                xval_vols{jj} = used_vols{nn(jj)};
            end
            [d1 d2 frames_per_vol] = size(xval_vols{1});

            % Get the features as vectors stacked into a matrix
            if (strcmp(method,'intensity'))
                X = get_features(train_vols, 'intensity', dataset);
                xval_X = get_features(xval_vols,'intensity', dataset);
            end
            if (strcmp(method,'lbp'))
                X = get_features(train_vols, 'lbp', dataset);
                xval_X = get_features(xval_vols,'lbp', dataset);
            end
            
            [m n] = size(X);
            [dummy xval_frames] = size(xval_X);
            % do PCA to reduce dimensionality
            D = [X xval_X];
            [D_red, D_mapping] = do_pca(D',pca_k);
            fprintf('Variance retained after PCA on X is: %f \n', ...
                     sum(D_mapping.lambda(1:pca_k))/sum(D_mapping.lambda));
            fprintf('size of reduced data is: %d * %d \n', size(D_red));
            train_data = D_red(1:n,:);
            xval_data = D_red(n+1:end,:);

            % fit GMM to the training data
            gmm_k = GMM_Ks(ii);
            obj = gmdistribution.fit(train_data,gmm_k, ...
                                    'Options',options, ...
                                    'Regularize',0.001, ...
                                    'Start', 'randSample', ...
                                    'CovType', 'diagonal', ...
                                    'Replicates', 1, ...
                                    'Start', 'randSample');

            % Cross validation using 3 volumes
            outliers = {};
            x_outliers = [];
            n = 1; % count number of outliers
            frame_no = 1; % frame number in the volume
            vol = 1;
            for i = 1:xval_frames
                [x_outliers(n)] = check_outlier(obj, xval_data(i,:), gmm_k, mahal_thresh);
                n = n+1;
                frame_no = frame_no+1;
                [dum1 dum2 frames_in_vol] = size(xval_vols{vol});
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
            diseased_frames(iter) = 0;
            for a = 1:length(outliers)
                   diseased_frames(iter) = diseased_frames(iter) + sum(outliers{a});
            end
            fprintf('Total diseased frames are: %d \n',diseased_frames(iter));
            
        end % iterations for each GMM K value
        avg_disesaed(ii) = sum(diseased_frames)/num_iterations;
    end % for each GMM N value
end