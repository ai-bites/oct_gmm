function [X,flat_vols] = preprocess(location,num_vols, dataset)
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

    % Flatten each volume
    for i = 1:num_vols
       flat_vols{i} = flatten_vol(volumes{i},0); 
       fprintf('flattened vol %i \n', i);
    end
    disp('Done with flattening the volumes');

    % Crop the volumes
    [cropped_vols] = crop_vols(flat_vols);
    disp('Done with cropping the volumes');

    % Now vectorise frames and form X matrix
    n = 1;
    for i = 1:num_vols
       [d1 d2 d3] = size(cropped_vols{i});
       for j = 1:d3
           img = imresize(cropped_vols{i}(:,:,j), 0.5);
           X(:,n) = img(:);
           n = n+1;
       end
    end
    disp('Done with forming the X matrix');

end