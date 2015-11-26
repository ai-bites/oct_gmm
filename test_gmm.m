function [actual_idx, predicted_idx] = test_gmm(obj, x_test_vols, Y_cropped_vols, ...
                    dataset, pca_k, best_gmm_k, mahal_thresh, diseased_vol_thresh, method)


    % some pre-requisites - count number of frames
    % in test data
    num_x = length(x_test_vols);
    num_y = length(Y_cropped_vols);
    num_xframes = 0;
    num_yframes = 0;
    % result of testing expected. 1 for diseased, 0 for normal
    actual_idx = [ones(num_y,1) ; zeros(num_x, 1)];
    % we don't really get the index corresponding to the volumes
    predicted_idx = zeros(num_y+num_x,1);
    % size of the normal volumes tested
    for a = 1:num_x
        [d1 d2 x_d3] = size(x_test_vols{a});
        num_xframes = num_xframes+x_d3;
    end
    % size of the diseased volumes tested
    for b = 1:num_y
        [d1 d2 y_d3] = size(Y_cropped_vols{b});
        num_yframes = num_yframes+y_d3;
    end

    Y = get_features(Y_cropped_vols, method, dataset);
    X = get_features(x_test_vols, method, dataset);
    [m,n] = size(Y);
    
%     data = [Y X];
%     [data_red] = do_pca(data', pca_k);
%     Y_test = data_red(1:n,:);
%     X_test = data_red(n+1:end,:);
    Y_test = do_pca(Y', pca_k);
    X_test = do_pca(X', pca_k);

    %----------------------------------------------------------------------
    % All DME volumes
    %----------------------------------------------------------------------
    dis_outliers = {};
    y_outliers = [];
    n = 1; % count number of outliers
    frame_no = 1; % frame number in the volume
    vol = 1;
    
    for i = 1:num_yframes
        [y_outliers(n)] = check_outlier(...
                              obj, Y_test(i,:),best_gmm_k, mahal_thresh);
        n = n+1;
        frame_no = frame_no+1;
        [dum1 dum2 frames_in_vol] = size(Y_cropped_vols{vol});
        % reset for next volume
        if (frame_no > frames_in_vol) 
            dis_outliers{vol} = y_outliers;
            frame_no = 1;
            vol = vol+1;
            n = 1;
            y_outliers = [];
        end
    end
    
    % Analyse the outliers obtained
    diseased_frames = 0;
    diseased_vols = 0;
    for i = 1:length(dis_outliers)
        if (sum(dis_outliers{i}) > diseased_vol_thresh)
           diseased_vols = diseased_vols+1;
           diseased_frames = diseased_frames + sum(dis_outliers{i});
           % fprintf('diseased vol number is: %d \n', i);
        end
    end
    fprintf('TEST: Diseased volumes are: %d, total diseased frames are: %d \n',...
            diseased_vols, diseased_frames);
    predicted_idx(1:diseased_vols,1) = 1; 
    
    %----------------------------------------------------------------------
    % All Normal test volumes
    %----------------------------------------------------------------------
    norm_outliers = {};
    x_outliers = [];
    n = 1; % count number of outliers
    frame_no = 1; % frame number in the volume
    vol = 1;
    
    % for volumes other than X val normal volumes
    for i = 1:num_xframes
        [x_outliers(n)] = check_outlier(...
                              obj, X_test(i,:),best_gmm_k, mahal_thresh);
        n = n+1;
        frame_no = frame_no+1;
        [a b frames_in_vol] = size(x_test_vols{vol});
        % reset for next volume
        if (frame_no > frames_in_vol) 
            norm_outliers{vol} = x_outliers;
            frame_no = 1;
            vol = vol+1;
            n = 1;
            x_outliers = [];
        end
    end
    % Analyse the outliers obtained in normal cases
    diseased_frames = 0;
    diseased_vols = 0;
    for i = 1:length(norm_outliers)
        if (sum(norm_outliers{i}) > diseased_vol_thresh)
           diseased_vols = diseased_vols+1;
           diseased_frames = diseased_frames + sum(norm_outliers{i});
           % fprintf('diseased vol number is: %d \n', i);
        end
    end
    fprintf('TEST: Normal volumes as diseased are: %d, total diseased frames are: %d \n',...
            diseased_vols, diseased_frames);
    predicted_idx((num_y+1):(num_y+1+diseased_vols-1), 1) = 1;
    
end