function [ isOut ] = check_outlier( obj, frame, gmm_k, mahal_thresh)
% CHECK_OUTLIER Checks if a given frame is an outlier or is within the GMM.
%
% obj   - the objective obtained after fitting the GMM model
% frame - the frame to be checked to see if it falls within the GMM model
% gmm_k - number of gaussians in the model
% mahal_thresh - threshold distance for Mahalnobis distance
% isOut  - 0|1. 1 if it is a outlier else 0


    [rx ry] = size(frame);
    if (rx ~= 1 && ry ~= 1) 
        error('Frame is not vectorised.');
    end

    % Distance to each gaussian in the fit.
    % mahal_dist = sqrt((xn - m)' * inv(S) * (xn - m));
    if (rx == 1)
        mahal_dist = mahal(obj,frame(1,:));
    else
        mahal_dist = mahal(obj,frame(:,1));
    end
    
    % Check to see if it is an outlier
    if (min(mahal_dist) > mahal_thresh)
        isOut = 1;
    else
        isOut = 0;
    end
    
%     for j = 1:gmm_k
%         if (mahal_dist(j) < mahal_thresh)
%             %disp('not an outlier');
%             isOut = 0;
%             return;
%         else
%             %disp('it is oulier');
%             isOut = 1;
%         end
%     end

end 

