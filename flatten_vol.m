function [res] = flatten_vol(vol, do_plot)
% Flattens all the frames in volumetric data and returns the 
% flattened volume.
%
% do_plot - pick first frame from the volume and plot it b4 and after
%           flattening
% vol     - Actual volume to flatten 
% res     - volume returned after flattening

    [d1 d2 d3] = size(vol);
    res = zeros(d1,d2,d3);
    
    if (do_plot) 
        figure;
        subplot(1,2,1); 
        imshow(vol(:,:,1),[]); 
    end
    
    % One frame from the volume at a time
    for j = 1:d3
        curr_img = vol(:,:,j);

        % Threshold 
        thresh_at = 50;
        mask = curr_img > thresh_at;
        %imshow(mask,[]);

        % Apply median filter
        img_med = medfilt2(mask,[5 5]);
        % imshow(img_med,[]);

        % close and then open 
        sec = strel('disk', 20);
        seo = strel('disk', 5);
        img_morph = imopen(imclose(img_med, sec), seo);

        % Fit a 2nd order polynomial
        [m n] = size(img_morph);
        [I J] = find(img_morph > 0);
        p = polyfit(J,I,2);
        x = 1:1:n;
        f = polyval(p,x);

        % Do some warping to flatten
        f_round = round(f);
        f_max = max(f_round); 
        f_min = min(f_round);
        new_mask = zeros(m,n);
        new_img = zeros(m,n);
        % Move each column of the image
        % down by a distance based on the polynomial fit
        for i = 1:n
            curr_col = img_morph(:,i);
            [idx val] = find(curr_col == 1);
            new_mask(idx+(f_max - f_round(i)),i) = val;
            img_val = curr_img(idx,i);
            new_img(idx+(f_max - f_round(i)),i) = img_val;
        end
        res(:,:,j) = new_img;
    end
    
    if (do_plot) 
        subplot(1,2,2); 
        imshow(res(:,:,1), []); 
    end
    
end