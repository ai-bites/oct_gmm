function [flat_vols] = preprocess(location, num_vols, dataset, method)
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
            if (ii == 1)
               figure; imshow(I,[]);
               title('Flat b4 filter');
               figure; imshow(J,[]);
               title('Flat after filter');
            end
        end
    end
    disp('Done with NL means filtering images');
    
    % Crop the volumes
    % [cropped_vols] = crop_vols(flat_vols);
    % disp('Done with cropping the volumes');
    
    
    
    
    
end

