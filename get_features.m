function [X] = get_features(volumes, method, dataset)


    switch lower(method)
        % We are dealing with intensity
        case 'intensity'
            % Now vectorise frames and form X matrix
            n = 1;
            for i = 1:length(volumes)
               [d1 d2 d3] = size(volumes{i});
               for j = 1:d3
                   %img = imresize(cropped_vols{i}(:,:,j), 0.5);
                   if (strcmp(dataset,'DUKE') && strcmp(method,'intensity'))
                       img = volumes{i}(:,1:256,j);
                   else
                       img = volumes{i}(:,:,j);
                   end
                   X(:,n) = img(:);
                   n = n+1;
               end
            end
        % We are dealing with the texture
        case 'lbp'
            n = 1;
            for i = 1:length(volumes)
               [d1 d2 d3] = size(volumes{i});
               for j = 1:d3
                  img = volumes{i}(:,:,j);
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

end