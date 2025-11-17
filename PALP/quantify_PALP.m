%
%   quantify_PALP_assay.m
%
clearvars
clc
close all

ROOTDIR = "E:/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath D:\Spatiotemporal_analysis\code
addpath D:\Spatiotemporal_analysis\code\util

warning('off', 'MATLAB:MKDIR:DirectoryExists');
warning('off', 'images:imfindcircles:warnForSmallRadius');

%% CONFIGURATIONS
prefix = "PALP-control-0202";
path = ROOTDIR + "PALP/" + prefix + "/";

ftype = "HD_align-";
fstart = 1;
fend = 8;

roi_sz = 4;

%% quantify PALP signals
close all

summary = table;

for ii = [fstart:fend]
    rep = sprintf('%02d', ii);

    % obtain the coordinate of induction sites and create mask on the sites
    t = 1;
    fname = path + "/" + ftype + rep + "/" + ftype + rep + "_t"+t+"c2";
    I = load_image_file(t, fname + ".tif");

    [centers, radii] = imfindcircles(I(:,:,1), [4 9], 'Sensitivity', 0.95);

    imsize = size(I(:, :, 1));
    center_x = centers(:, 1);
    center_y = centers(:, 2);
    [X, Y] = meshgrid(1:imsize(2), 1:imsize(1));
    mask = false(imsize);
    radii(:) = roi_sz;
    for i = 1:numel(radii)
        mask = mask | hypot(X - center_x(i), Y - center_y(i)) <= radii(i);
    end

    figure
    B = labeloverlay(histeq(I), mask);
    imshow(B, 'Border', 'tight')

    % obtain the labels for the induction sites
    t = 1;
    I = load_image_file(t, path + "/" + ftype + rep + "/" + ftype + rep + "_t"+t+"c1_ORG.tif");

    stats = regionprops("table", mask, I, ["MeanIntensity","Centroid"]);
    [~, index] = sort(stats.MeanIntensity, 'descend');
    data = [repelem(ii, length(index))', (1:length(index))', stats.Centroid(index, 1), stats.Centroid(index, 2)];
    
    arrayfun(@(i)text(data(i,3), data(i,4), sprintf("%d",data(i,2)), 'color', 'w', 'FontSize', 20), 1:size(data,1));

    induction_site = array2table(data, "VariableNames", ["rep", "site_index" "x" "y"]);

    % save the labeled image for reference
    labeled_im = fname + "_labeled.tif";
    exportgraphics(gcf, labeled_im);

    % load time 2 after induction
    t = 2;
    oxi_im = load_image_file(t, path + "/" + ftype + rep + "/" + ftype + rep + "_t2c1_ORG.tif");
    red_im = load_image_file(t, path + "/" + ftype + rep + "/" + ftype + rep + "_t2c2_ORG.tif");

    norm_oxi_im = oxi_im;
    norm_red_im = red_im;

    for i = 1:height(induction_site)
    % for i = 1
        center_x = induction_site.x(i);
        center_y = induction_site.y(i);
        [X, Y] = meshgrid(1:imsize(2), 1:imsize(1));
        mask = false(imsize);
        mask = mask | hypot(X - center_x, Y - center_y) <= roi_sz;

        stats = regionprops(mask, norm_oxi_im, 'PixelValues');
        induction_site.o(i) = mean(stats.PixelValues);

        stats = regionprops(mask, norm_red_im, 'PixelValues');
        induction_site.r(i) = mean(stats.PixelValues);

        induction_site.or(i) = induction_site.o(i)/(induction_site.o(i)+induction_site.r(i));
    end

    summary = [summary; induction_site];
end
% close all;

%% quantify over time

close all

tstart = 1;
tend = 4;

summary = table;

% for ii = [fstart:fend]
for ii = [1:4 6:8]
    rep = sprintf('%02d', ii);
    % rep = "";


    % obtain the coordinate of induction sites and create mask on the sites
    t = 1;
    fname = path + "/" + ftype + rep + "/" + ftype + rep + "_t"+t+"c2";
    I = load_image_file(t, fname + ".tif");

    % [centers, radii] = imfindcircles(I(:,:,1), [4 9], 'Sensitivity', 0.95);
    [centers, radii] = imfindcircles(I(:,:,1), [2 8], 'Sensitivity', 0.95); % for control-0202

    imsize = size(I(:, :, 1));
    center_x = centers(:, 1);
    center_y = centers(:, 2);
    [X, Y] = meshgrid(1:imsize(2), 1:imsize(1));
    mask = false(imsize);
    radii(:) = roi_sz;
    for i = 1:numel(radii)
        mask = mask | hypot(X - center_x(i), Y - center_y(i)) <= radii(i);
    end

    figure
    B = labeloverlay(histeq(I), mask);
    imshow(B, 'Border', 'tight')

    % obtain the labels for the induction sites
    t = 1;
    I = load_image_file(t, path + "/" + ftype + rep + "/" + ftype + rep + "_t"+t+"c1_ORG.tif");

    stats = regionprops("table", mask, I, ["MeanIntensity","Centroid"]);
    [~, index] = sort(stats.MeanIntensity, 'descend');
    data = [repelem(ii, length(index))', (1:length(index))', stats.Centroid(index, 1), stats.Centroid(index, 2)];
    
    arrayfun(@(i)text(data(i,3), data(i,4), sprintf("%d",data(i,2)), 'color', 'w', 'FontSize', 20), 1:size(data,1));

    induction_site = array2table(data, "VariableNames", ["rep", "site_index" "x" "y"]);

    % save the labeled image for reference
    labeled_im = fname + "_labeled.tif";
    exportgraphics(gcf, labeled_im);

    % load time 2 after induction
    for t = tstart:tend
        oxi_im = load_image_file(t, path + "/" + ftype + rep + "/" + ftype + rep + "_t" + sprintf("%02d", t) + "c1_ORG.tif");
        red_im = load_image_file(t, path + "/" + ftype + rep + "/" + ftype + rep + "_t" + sprintf("%02d", t) + "c2_ORG.tif");
    
        norm_oxi_im = oxi_im;
        norm_red_im = red_im;
    
        for i = 1:height(induction_site)
        % for i = 1
            center_x = induction_site.x(i);
            center_y = induction_site.y(i);
            [X, Y] = meshgrid(1:imsize(2), 1:imsize(1));
            mask = false(imsize);
            mask = mask | hypot(X - center_x, Y - center_y) <= roi_sz;
    
            stats = regionprops(mask, norm_oxi_im, 'PixelValues');
            induction_site.o(i) = mean(stats.PixelValues);
    
            stats = regionprops(mask, norm_red_im, 'PixelValues');
            induction_site.r(i) = mean(stats.PixelValues);
    
            induction_site.or(i) = induction_site.o(i)/(induction_site.o(i)+induction_site.r(i));

            induction_site.t(i) = t;
        end
    
        summary = [summary; induction_site];
    end
end

%% FUNCTIONS
function I = load_image_file(t, file_path)
    try
        I = imread(file_path);
    catch
        file_path = strrep(file_path, "t"+t, "t0"+t);                                                                                                       
        I = imread(file_path);
    end
end

function I = min_max_normalization(I, Imin, Imax, a, b)
    I = (I - Imin) / (Imax - Imin) * (b-a) + a;
end

function quantile_normalization(images)
    flattened_im = cellfun(@im)
end

function norm_I = percentile_normalization(I, lower_p, upper_p)
    p_low = prctile(I(:), lower_p);
    p_high = prctile(I(:), upper_p);

    % I_clip = max(min(I(:), p_high), p_low);
    norm_I = (I - p_low) ./ (p_high - p_low+1e-20);
    norm_I = reshape(norm_I, size(I));
end