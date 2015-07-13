function [ isOut,min_mahal] = check_outlier( obj, frame, gmm_k, mahal_thresh)
% CHECK_OUTLIER Checks if a given frame is an outlier or is within the GMM.
%
% obj   - the objective obtained after fitting the GMM model
% frame - the frame to be checked to see if it falls within the GMM model
% num_gauss - number of gaussians in the model
% isIn  - 0|1. 1 if it is a outlier else 0


    % Assume it is out
    isOut = 1;

    [rx ry] = size(frame);
    if (rx ~= 1 && ry ~= 1) 
        error('Frame is not vectorised.');
    end

    % Distance to each gaussian in the fit.
    if (rx == 1)
        mahal_dist = mahal(obj,frame(1,:));
    else
        mahal_dist = mahal(obj,frame(:,1));
    end
    min_mahal = min(mahal_dist);
    % mahal_dist = sqrt((xn - m)' * inv(S) * (xn - m));
    
    for j = 1:gmm_k
        curr_mean = obj.mu(j,:);
        curr_sigma = obj.Sigma(:,:,j);
        
        % dm(j) = ((frame - curr_mean)*pinv(curr_sigma)*(frame - curr_mean)')
        % fprintf('Mahal dist computed: %f \n' , mahal_dist(j));
        if (mahal_dist(j) < mahal_thresh)
            %disp('not an outlier');
            isOut = 0;
            return;
        else
            %disp('it is oulier');
            isOut = 1;
        end
    end

end 

