function [threshold] = get_disease_thresh(obj, X_train, used_vols, best_gmm_k, mahal_thresh)


    outliers = {};
    x_outliers = [];
    n = 1; % count number of outliers
    frame_no = 1; % frame number in the volume
    vol = 1;
    [m n] = size(X_train);
    for i = 1:m
        [x_outliers(n)] = check_outlier(obj, X_train(i,:), best_gmm_k, mahal_thresh);
        n = n+1;
        frame_no = frame_no+1;
        [dum1 dum2 frames_in_vol] = size(used_vols{vol});
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
    for a = 1:length(outliers)
        diseased_frames = diseased_frames + sum(outliers{a});
    end
    threshold = round(diseased_frames / length(used_vols));
    fprintf('Threshold for diseased frames is: %d \n', threshold);
    
end