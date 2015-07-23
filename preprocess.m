function [X,flat_vols] = preprocess(location,num_vols, dataset, method)
% Does the preprocessing on the mat files in the location.
%
% location - folder+filename. Eg, 'data/DME/patient' where 
%            patient is the name of mat files indexed patient1.mat,
%            patient2.mat, etc.
% num_vols - number of volumes to preprocess and form X


    % Read all the needed volumes
    for i = 1:num_vols
        fn = strcat(location,num2str(i));
        vol = load(fn);
        volumes{i} = vol.vol;
        fprintf('loaded vol %i \n', i);
    end
    
    % resize images in the volume
    for i = 1:num_vols
        [d1 d2 d3] = size(volumes{i});
        for j = 1:d3
            reduced_vols{i}(:,:,j) = imresize(volumes{i}(:,:,j), 0.5);
        end
    end
    
    % Flatten each volume
    for i = 1:num_vols
       flat_vols{i} = flatten_vol(reduced_vols{i},0); 
       fprintf('flattened vol %i \n', i);
    end
    disp('Done with flattening the volumes');
    
    % Filter using non-local Means filter
    if (strcmp(lower(method),'lbp'))
        options.kernelratio = 6;
        options.windowratio = 6;
        for i = 1:num_vols
            [d1 d2 d3] = size(flat_vols{i});
            for ii = 1:d3
                I = (flat_vols{i}(:,:,ii));
                max_I = max(max(I));
                I = I / max_I;
                J = NLMF(I,options);
                flat_vols{i}(:,:,ii) = J*max_I;
            end
        end
        disp('Done with NL means filtering images');
    end
    
    % Crop the volumes
    [cropped_vols] = crop_vols(flat_vols);
    disp('Done with cropping the volumes');

    switch lower(method)
        % We are dealing with intensity
        case 'intensity'
            % Now vectorise frames and form X matrix
            n = 1;
            for i = 1:num_vols
               [d1 d2 d3] = size(cropped_vols{i});
               for j = 1:d3
                   %img = imresize(cropped_vols{i}(:,:,j), 0.5);
                   img = cropped_vols{i}(:,:,j);
                   X(:,n) = img(:);
                   n = n+1;
               end
            end
        % We are dealing with the texture
        case 'lbp'
            n = 1;
            for i = 1:length(volumes)
               [d1 d2 d3] = size(cropped_vols{i});
               for j = 1:d3
                  img = cropped_vols{i}(:,:,j);
                  J = lbp(img, 2, 16, 0, 'h');
                  X(:,n) = J(:);
                  n = n+1;
               end
            end
        % None of the above
        otherwise
            disp('not the right method for pre-processing');
            X = 0;
            return;
    end
    
    disp('Done with forming the X matrix');

