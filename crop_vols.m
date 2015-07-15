function [cropped_vols] = crop_vols(volumes)
% Given a struct of volumes, crops and removes the zero rows 
% from the frames in the volumes and returns the cropped volumes in a
% struct.
%
% volumes      - A struct with several volumes of data
% cropped_vols - A Struct of resulting cropped volumes

    num_vols = length(volumes);
    
    % Decide on the upper and lower bound for the ROI
    k = 1;
    for i = 1:num_vols
        if (i == 11) continue; end
        [d1 d2 d3] = size(volumes{i});
        for j = 1:d3
            img = volumes{i}(:,:,j);
            % sum along the columns and find locations
            % where intensities are not equal to 0.
            % We crop these locations
            dummy = find(sum(img')' ~= 0);
            min_idxs(k) = min(dummy);
            max_idxs(k) = max(dummy);
            k = k+1;
        end
    end
    
    min_crop = min(min_idxs);
    max_crop = max(max_idxs);
    
    % we know the index where to crop. So now crop all volumes
    for i = 1:num_vols
        cropped_vols{i} = volumes{i}(min_crop:max_crop,:,:);
    end

end