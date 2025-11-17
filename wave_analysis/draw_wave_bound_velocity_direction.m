%
%   set_wave_boundary.m
%
clearvars
clc
close all

addpath("D:/Spatiotemporal_analysis/code/util/");

ROOTDIR = "E:/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% load image for 36 h

cy5_file = "D:\Spatiotemporal_analysis\wave_mCh300FBS10-1\Well10\c2_cy5_adj\mCh300FBS10-1_s10t66c2_ORG_ff_pff_bc_crop_final_bw.tif";
bw_file  = "D:\Spatiotemporal_analysis\wave_mCh300FBS10-1\Well10\c2_cy5_adj\mCh300FBS10-1_s10t66c2_ORG_ff_pff_bc_crop_final_bw.tif";
out_prefix = "D:\Spatiotemporal_analysis\wave_mCh300FBS10-1\Well10\c2_cy5_adj\mCh300FBS10-1_s10t66c2_ORG_ff_pff_bc_crop_final_bw";
im_range = [1 1]; %[14 21];
xstart = 0;
ystart = 0;
imsize = [2200 3700];

cy5 = cell(im_range(2), 1);
bw = cell(im_range(2), 1);

prev_bw = zeros(imsize(1), imsize(2));
for i = im_range(1):im_range(2)
    cy5{i} = imread(cy5_file, i); %, 'Info', info);
    bw{i}  = imread(bw_file, i); %, 'Info', info);
end

%%
bound_type = "gauss";

if strcmp(bound_type, "nuc")
    % perform boundary calculation by dilating and connecting the nucleus
    [bdy, Bo] = wave_boundary_nuc(bw, im_range);
elseif strcmp(bound_type, "gauss")
    % perform boundary calculation by using a gaussian filter
    sigma = 15;% 15; %15; %20; %' 15 20;
    binarize_threshold = 0.75; %0.75; % 0.75; %2.5; %0.75; %7.5; %5; %7.5;
    medfilt_sz = 12; %25;
    
    [bdy, Bo] = wave_boundary_gaussfilt(bw, im_range, sigma, binarize_threshold, medfilt_sz);
    
end


%% save files
bo_outfile = out_prefix + "_bo.mat";
bdy_outfile = out_prefix + "_bdy.mat";
save(bo_outfile, 'Bo');
save(bdy_outfile, 'bdy');

%%
outfile = out_prefix + "_bo.pdf";

figure
imshow(zeros(imsize)+255, 'Border', 'tight')

hold on

for i = im_range(1):1:im_range(2)
    bo = Bo{i};
    for j = 1:length(bo)
        b = bo{j};
        if ~isempty(b)
            plot(b(:, 2), b(:, 1), 'k-', 'LineWidth',2)
        end
    end

end
hold off

exportgraphics(gcf, outfile)

%% construct contours

% foi = find(~cellfun(@(x) isempty(bwboundaries(full(x))), bdy)); %
foi = find(~cellfun(@(x) isempty(x), Bo)); %
C = bdy(foi);
% filteredImage = zeros(imsize, imsize, length(C));
filteredImage = zeros(imsize(1), imsize(2), length(C));
for i = 2:2:length(C)
    filteredImage(:,:,i) = full(C{i});  % assign each non-empty array
end

Z = length(C)+1-sum(filteredImage,3);   % construct contours (Z dim)

% rf = 2; % resize_factor
% Z = imcrop(Z, [sz/2-sz/(rf*2) sz/2-sz/(rf*2) sz/rf sz/rf]);
% Z = imcrop(Z, [sz/2-sz/(rf*2) sz/5*2-sz/(rf*2) sz/rf sz/rf]);
% Z = imcrop(Z, [50 400 sz/rf sz/rf]);

%% calculate vector
% vector_per_side = 13; %200;
vector_per_side = 20;
sigma = 1.5; %3.5;

[X, Y] = meshgrid(1:1:size(Z, 2), 1:1:size(Z, 1));

ws = size(Z, 2)/vector_per_side;
[Xq, Yq] = meshgrid(1:ws:size(Z, 2), 1:ws:size(Z, 1));
Vq = interp2(X, Y, Z, Xq, Yq);

[Gmag, Gdir] = imgradient(Vq, 'central');
Gmag(Gmag == 0) = nan;
Gmag = 1./Gmag;
U = Gmag .* cosd(Gdir);
V = -1 .* Gmag .* sind(Gdir);

U(isnan(U)) = 0;
V(isnan(V)) = 0;

% smoothen component
UG = imgaussfilt(U, sigma);
VG = imgaussfilt(V, sigma);


magnitude = sqrt(UG.^2 + VG.^2);
normUG = UG ./ magnitude; % Normalize U
normVG = VG ./ magnitude; % Normalize V

normUG(isnan(normUG)) = 0;
normVG(isnan(normVG)) = 0;

%% plot the contours
figure
imshow(Z, [min(Z(:))-1 max(Z(:))+1], 'Border', 'tight')

hold on
quiver(Xq, Yq, UG, VG, 0.5, 'color', 'k', 'LineWidth', 2)
hold off


%% plot the speed
foi = find(~cellfun(@(x) isempty(x), Bo));
C = bdy(foi);
filteredImage = zeros(imsize(1), imsize(2), length(C));
for i = 1:1:length(C)
    filteredImage(:,:,i) = full(C{i});  % assign each non-empty array
end

Z = length(C)+1-sum(filteredImage,3);   % construct contours (Z dim)

% calculate vector
vector_per_side = 30;
sigma = 1.5; %3.5;

[X, Y] = meshgrid(1:1:size(Z, 2), 1:1:size(Z, 1));

ws = size(Z, 2)/vector_per_side;
[Xq, Yq] = meshgrid(1:ws:size(Z, 2), 1:ws:size(Z, 1));
Vq = interp2(X, Y, Z, Xq, Yq);

[Gmag, Gdir] = imgradient(Vq, 'central');
Gmag(Gmag == 0) = nan;
Gmag = 1./Gmag;
U = Gmag .* cosd(Gdir);
V = -1 .* Gmag .* sind(Gdir);

U(isnan(U)) = 0;
V(isnan(V)) = 0;

% smoothen component
UG = imgaussfilt(U, sigma);
VG = imgaussfilt(V, sigma);

magnitude = sqrt(UG.^2 + VG.^2);
mag_smooth = imgaussfilt(magnitude, 3.5);
%% plot the contours
figure
% imshow(Z, [0 length(C)+1], 'Border', 'tight')
surf(Xq, Yq, zeros(size(mag_smooth)), mag_smooth, 'EdgeColor', 'none')
view(2)
axis ij square off
clim([0 1.2])

%% 
figure
imshow(zeros(size(Z))+255, 'Border', 'tight')
hold on
quiver(Xq, Yq, UG, VG, 0.5, 'color', 'k', 'LineWidth', 2)
hold off

%% Compute gradient direction using atan2
Gdir = atan2(VG, UG); % atan2(Y, X), note the negative VG to correct orientation

Gdir_deg = rad2deg(Gdir);

% Convert from [-180, 180] range to [0, 360]
Gdir_corrected = mod(Gdir_deg, 360);

% Flatten and remove NaNs
Gdir_rad = deg2rad(Gdir_corrected(:));
% Gdir_rad = Gdir_rad(~isnan(Gdir_rad));

% Plot polar histogram
figure;
polarhistogram(Gdir_rad, 18);  % 36

%% plot the wave outlines

cy5_fr = im_range(2)-1;

all_sigma = 5:5:20; % 7.5:2.5:15;
all_bin_thresh = [0.5 0.75 1 2.5 5]; %0.5:0.1:1 %[0.5 0.75 1 2.5 5];

figure
imshow(zeros(size(cy5{cy5_fr}))+255, 'Border', 'tight')

hold on
for i = im_range(1):2:im_range(2)
    bo = Bo{i};
    for j = 1:length(bo)
        if ~isnan(bo{j})
            % plot(bo(:, 2), bo(:, 1), '-', "LineWidth", 4, 'Color', 'black');
            plot(bo{j}(:, 2), bo{j}(:, 1), '-', "LineWidth", 3, 'Color', 'black');
        end
    end
end
hold off


%% FUNCTION
function [bdy, Bo] = wave_boundary_gaussfilt(bw, im_range, sigma, binarize_threshold, medfilt_sz)
    sm_area_thres = 3000;
    sm_hole_thres = 10000; %50000;

    bdy = cell(im_range(2), 1);
    Bo = cell(im_range(2), 1);

    run_medfilt = 1;
    if ~exist('medfilt_sz', 'var')
        run_medfilt = 0;
    end

    for i = im_range(1):im_range(2)
        B = imgaussfilt(double(bw{i}), sigma);

        BW = imbinarize(B, binarize_threshold);    
        BW = bwareaopen(BW, sm_area_thres);         % Remove small objects and speckles
        BW = ~bwareaopen(~BW, sm_hole_thres); 

        % if run_medfilt
        %     % BW = imclose(BW, strel('disk', 10));
        %     BW = medfilt2(BW, [medfilt_sz medfilt_sz]);
        % end


        bdy{i} = sparse(BW);
       
        bo = bwboundaries(BW);
        for j = 1:length(bo)
            if length(bo{j}) < 300
                bo{j} = [];
            end
        end
        Bo{i} = bo;
    end
end

function [bdy, Bo] = wave_boundary_nuc(cy5, im_range)
    sm_nucl_thres = 5; %15 20
    lg_nucl_thres = 500;
    sm_area_thres = 30000;
    sm_hole_thres = 10000; %50000;
    
    dilate_sz = 15; %20;
    close_sz = 20; %15;
    % erode_sz = dilate_sz - 15;
    % medfilt_sz = 30; %25;
    
    % bdy keeps the mask for dead cell population
    bdy = cell(im_range(2), 1);
    Bo = cell(im_range(2), 1);
    
    % for each image, outline the wave
    for i = im_range(1):im_range(2)
        % if i <= im_range(1)
        if i == 14
            sm_nucl_thres = 1;
            sm_area_thres = 3000;
        else
            sm_nucl_thres = 5;
            sm_area_thres = 30000;
        end
    
        I = cy5{i};
    
        % remove any labels with area > lg_nucl_thres, which is assumed to be debris
        L = bwlabel(I);
        obj = regionprops('table', L, 'Area'); % calculate the area for each component
        if any(obj.Area > lg_nucl_thres)
            I(ismember(L,find(obj.Area > lg_nucl_thres))) = 0; % rm by setting to 0
        end
    
        % remove any labels with area < sm_nucl_thres, which is assumed to be debris
        L = bwlabel(I);
        obj = regionprops('table', L, 'Area'); % calculate the area for each component
        if any(obj.Area < sm_nucl_thres)
            I(ismember(L,find(obj.Area < sm_nucl_thres))) = 0; % rm by setting to 0
        end
    
        % Dilate, Close, Rm sm obj, Rm holes, Erode, Rm sm obj
        BW = imdilate(I, strel('disk', dilate_sz));   % Dilate to connect the dead cells
        BW = imclose(BW, strel('disk', close_sz));  % Close the gaps in the dead cells
        BW = bwareaopen(BW, sm_area_thres);         % Remove small objects and speckles
        BW = ~bwareaopen(~BW, sm_hole_thres);       % Remove holes
        BW = medfilt2(BW, [medfilt_sz medfilt_sz]);
        
        % store the outline in bdy
        bdy{i} = sparse(BW);
        Bo{i} = bwboundaries(BW);
    end
end
