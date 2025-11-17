%
%   calculate density based on segmented nucleus image
%

function [den, x, y, n_cell, imsize, valid_idx] = calculate_density(fname, den_reg)
    label_nucleus = imread(fname);
    
    info = imfinfo(fname);
    
    imsize = [info.Height info.Width];
    
    % calculate centroid and nucleus size
    centroids = regionprops(label_nucleus, 'Centroid', 'EquivDiameter');
    centroids = struct2table(centroids);
    
    % remove nucleus that are too big or too small
    nuc_lwr = mean(centroids.EquivDiameter) - 2*std(centroids.EquivDiameter);
    nuc_upr = mean(centroids.EquivDiameter) + 2*std(centroids.EquivDiameter);
    valid_idx = find(centroids.EquivDiameter > nuc_lwr & centroids.EquivDiameter < nuc_upr);
    centroids = centroids(valid_idx, :);
    
    % calculate distance between centroids
    D = squareform(pdist(centroids.Centroid));
    
    % calculate the number of cells in its neighboring region
    n_cell = sum(D < den_reg, 1);
    
    x = centroids.Centroid(:,1);
    y = centroids.Centroid(:,2);
    
    % interpolate the densitys by a grid
    [Xq, Yq] = ndgrid(1:1:imsize(1), 1:1:imsize(2));
    vq = griddata(x, y, n_cell, Xq, Yq, 'natural');
    
    den = vq';
end