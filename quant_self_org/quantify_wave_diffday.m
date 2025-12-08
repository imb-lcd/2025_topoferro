%
%   quantify_wave_selforg.m
%
clearvars
clc
close all

addpath("D:/Spatiotemporal_analysis/code/util/");

ROOTDIR = "E:/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% load image

% % for 12 h
cy5_file = "E:/movies_diff_day/seed05_day10_s24_diffday-1023-all_s24tAllc2_ORG_rgb.tif";
bw_file  = "E:/movies_diff_day/seed05_day10_s24_diffday-1023-all_s24tAllc2_ORG_bw.tif";
im_range = [15 25];
outpath = "E:/movies_diff_day/Figure1/24h_diffday-1023-1d_s24/";
outfile_prefix = "diffday-1023-1d_s02t" + im_range(1) + "-" + im_range(2) + "c2_ORG";
xstart = 2400; % 1; 
ystart = 3200; % 1;

% % for 24 h
% cy5_file = "E:/movies_diff_day/seed10_day15_s02_diffday-1023-1d_s02tAllc2_ORG_rgb.tif";
% bw_file  = "E:/movies_diff_day/seed10_day15_s02_diffday-1023-1d_s02tAllc2_ORG_bw.tif";
% im_range = [25 35];
% outpath = "E:/movies_diff_day/Figure1/36h_diffday-1023-1d_s02/";
% outfile_prefix = "diffday-1023-1d_s02t" + im_range(1) + "-" + im_range(2) + "c2_ORG";
% xstart = 650; %510; % 1; 
% ystart = 3070; % 1;

% % % for 24 h second choice
% cy5_file = "E:/wave_diffday-1118-1d/Well14/c2_cy5_adj/diffday-1118-1d_s14tAllc2_ORG_adj_rgb.tif";
% % bw_file  = "E:/wave_diffday-1118-1d/Well14/c2_cy5_adj/diffday-1118-1d_s14tAllc2_ORG_adj_bw2.tif";
% bw_file  = "E:/wave_diffday-1118-1d/Well14/c2_cy5_adj/diffday-1118-1d_s14tAllc2_ORG_adj_bw_wv1only.tif";
% im_range = [20 30];
% outpath = "E:/movies_diff_day/Figure1/36h_diffday-1118-1d_s14/";
% outfile_prefix = "diffday-1023-1d_s02t" + im_range(1) + "-" + im_range(2) + "c2_ORG";
% xstart = 3670; % 3670+530;
% ystart = 3360; % 3360+270;

% for 48 h
% cy5_file = "D:/Spatiotemporal_analysis/wave_mCh300FBS10-1/Well10/c2_cy5_adj/mCh300FBS10-1_s10t30-66c2_ORG_ff_pff_bc_rgb.tif";
% bw_file = "D:/Spatiotemporal_analysis/wave_mCh300FBS10-1/Well10/c2_cy5_adj/mCh300FBS10-1_s10t30-66c2_ORG_bw.tif";
% im_range = [22 32]; %[14 21];
% outpath = "E:/movies_diff_day/Figure1/48h_mCh300FBS10-1_s10/";
% outfile_prefix = "mCh300FBS10-1_s10t" + im_range(1) + "-" + im_range(2) + "c2_ORG";
% xstart = 610; %510; 
% ystart = 320; 

% for all
sz = 1400; %877;%  1400;
imsize = [sz sz]; %4977; %5076;

cy5 = cell(im_range(2), 1);
bw = cell(im_range(2), 1);

for i = im_range(1):im_range(2)
    cy5{i} = imread(cy5_file, i); %, 'Info', info);
    bw{i}  = imread(bw_file, i); %, 'Info', info);

    cy5{i} = imcrop(cy5{i}, [xstart, ystart, imsize(2)-1, imsize(1)-1]);
    bw{i}  = imcrop(bw{i},  [xstart, ystart, imsize(2)-1, imsize(1)-1]);
end

%% draw wave boundaries
bound_type = "gauss";

if strcmp(bound_type, "nuc")
    % perform boundary calculation by dilating and connecting the nucleus
    [bdy, Bo] = wave_boundary_nuc(bw, im_range);
elseif strcmp(bound_type, "gauss")
    % perform boundary calculation by using a gaussian filter
    sigma = 15; %15; %20; %' 15 20;
    binarize_threshold = 0.75; % 0.75; %2.5; %0.75; %7.5; %5; %7.5;
    medfilt_sz = 12; %12; %25;
    
    [bdy, Bo] = wave_boundary_gaussfilt(bw, im_range, sigma, binarize_threshold, medfilt_sz);
end

%% construct contours

foi = find(~cellfun(@(x) isempty(x), Bo)); %
C = bdy(foi);

filteredImage = zeros(imsize(1), imsize(2), length(C));
for i = 1:1:length(C) % 1:1:length(C)
    filteredImage(:,:,i) = full(C{i});  % assign each non-empty array
end

Z = length(C)+1-sum(filteredImage,3);   % construct contours (Z dim)

figure
imshow(Z, [min(Z(:))-1 max(Z(:))], 'Border', 'tight')

% %% calculate nematic order vector field

vector_per_side = 10; %10; % 15
sigma = 1; %1.5; %1.; %1.5; % 2;

[X, Y] = meshgrid(1:1:size(Z, 2), 1:1:size(Z, 1));

ws = size(Z, 2)/vector_per_side;
[Xq, Yq] = meshgrid(ws/2:ws:size(Z, 2), ws/2:ws:size(Z, 1));
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

% % recalculate to make the magnitude the same
% mag = sqrt(UG.^2 + VG.^2);
% normUG = UG ./ mag; % Normalize U
% normVG = VG ./ mag; % Normalize V
% 
% normUG(isnan(normUG)) = 0;
% normVG(isnan(normVG)) = 0;

%% plot the quivers
savefigure = 1;

figure 
imshow(zeros(size(Z))+255, 'Border', 'tight')
axis ij
hold on
quiver(Xq, Yq, UG, VG, 1, 'color', 'k', 'LineWidth', 2)
% q = quiver(Xq, Yq, UG, VG, 1, 'color', 'k', 'LineWidth', 2, 'ShowArrowHead', 'off');
hold off

if savefigure
    % outfile = outpath + outfile_prefix + "_quiver_noarrow_new.pdf";
    outfile = outpath + outfile_prefix + "_quiver_arrow_new.pdf";
    exportgraphics(gcf, outfile);
end

%% plot the contours and nematic order
savefigure = 0;

figure
imshow(Z, [min(Z(:))-1 max(Z(:))], 'Border', 'tight')

hold on
% % quiver(Xq, Yq, normUG, normVG, 0.5, 'color', 'k', 'LineWidth', 2)
q = quiver(Xq, Yq, UG, VG, 0.5, 'color', 'k', 'LineWidth', 2);
hold off

if savefigure 
    outfile = outpath + outfile_prefix + "_contour.pdf";
    exportgraphics(gcf, outfile);
end

% plot the polarhistogram while quiver is still opened
% Extract the vector components from the quiver handle
Uq = q.UData;
Vq = q.VData;

% Compute each vector’s angle (in radians), wrapping to [0,2π)
theta = atan2(-Vq(:), Uq(:));           
theta = mod(theta, 2*pi);              

% Remove any zero-length vectors
mag = sqrt(UG.^2 + VG.^2);
theta = theta(mag>0);
%%
savefigure = 0;
% Choose number of bins and plot
nbins = 24;   % 10 degree bins
figure;
polarhistogram(theta, nbins, 'Normalization', 'pdf');
rlim([0 0.8])
ax = gca;
ax.ThetaTick = 0:45:315;     
if savefigure
    outfile = outpath + outfile_prefix + "_polarhistogram.pdf";
    exportgraphics(gcf, outfile);
end


%% calculate cirular dispersion
% from the Circular Statistics Toolbox
R = circ_r(theta); % calculate mean resultant length    
circ_var = 1 - R;        

[p, ~] = circ_rtest(theta);

disp(circ_var)
disp(p)

%% plot boundary on cy5
savefigure = 1;

for i = im_range(1):2:length(Bo)
    figure
    imshow(cy5{i}, 'Border', 'tight')
    hold on
    bo = Bo{i};
    for j = 1:length(bo)
        if ~isnan(bo{j})
            plot(bo{j}(:, 2), bo{j}(:, 1), '-', 'Color', 'white', 'LineWidth', 2)
        end
    end
    hold off
    outname = "diffday-1023-1d_s02t" + sprintf("%02d", i) + "c2_ORG";
    if savefigure 
        outfile = outpath + outname + "cy5_boundary.pdf";
        exportgraphics(gcf, outfile);
    end
end



%% plot boundary only
savefigure = 1;

figure
imshow(zeros(size(Z))+255, 'Border', 'tight')
hold on
for i = im_range(1):2:length(Bo)
    bo = Bo{i};
    for j = 1:length(bo)
        if ~isnan(bo{j})
            plot(bo{j}(:, 2), bo{j}(:, 1), '-', 'Color', 'k', 'LineWidth', 2)
        end
    end
end
hold off

if savefigure 
    outfile = outpath + outfile_prefix + "_boundary.pdf";
    exportgraphics(gcf, outfile);
end

%% plot speed
savefigure = 1;

vector_per_side = 30;
sigma = 1.5; %3.5;    

filteredImage = zeros(imsize(1), imsize(2), length(C));
for i = 1:2:length(C)
    filteredImage(:,:,i) = full(C{i});  % assign each non-empty array
end

Z = length(C)+1-sum(filteredImage,3);   % construct contours (Z dim)

[X, Y] = meshgrid(1:1:size(Z, 2), 1:1:size(Z, 1));

ws = size(Z, 2)/vector_per_side;
[Xq, Yq] = meshgrid(ws/2:ws:size(Z, 2), ws/2:ws:size(Z, 1));
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

mag = sqrt(UG.^2 + VG.^2);
mag_smooth = imgaussfilt(mag, 3.5);

figure
surf(Xq, Yq, zeros(size(mag_smooth)), mag_smooth, 'EdgeColor', 'none')
view(2)
axis ij square off
clim([0 1.4])

if savefigure
    outfile = outpath + outfile_prefix + "_velocity.pdf";
    exportgraphics(gcf, outfile);
end

%% plot vector field only
savefigure=1

vector_per_side = 15;
sigma = 1.5; %3.5;    

filteredImage = zeros(imsize(1), imsize(2), length(C));
for i = 1:2:length(C)
    filteredImage(:,:,i) = full(C{i});  % assign each non-empty array
end

Z = length(C)+1-sum(filteredImage,3);   % construct contours (Z dim)

[X, Y] = meshgrid(1:1:size(Z, 2), 1:1:size(Z, 1));

ws = size(Z, 2)/vector_per_side;
% [Xq, Yq] = meshgrid(ws/2:ws:size(Z, 2), ws/2:ws:size(Z, 1));
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


mag = sqrt(UG.^2 + VG.^2);
normUG = UG ./ mag; % Normalize U
normVG = VG ./ mag; % Normalize V

normUG(isnan(normUG)) = 0;
normVG(isnan(normVG)) = 0;

figure
imshow(zeros(size(Z))+255, 'Border', 'tight')

hold on
quiver(Xq, Yq, normUG, normVG, 0.5, 'color', 'k', 'LineWidth', 2)
% q = quiver(Xq, Yq, UG, VG, 0.5, 'color', 'k', 'LineWidth', 2);
hold off

if savefigure
    % outfile = outpath + outfile_prefix + "_vector_wMagnitude.pdf";
    outfile = outpath + outfile_prefix + "_vector_sameMagnitude.pdf";
    exportgraphics(gcf, outfile);
end


%% plot velocity using violin plot
[~, p1] = vartest2(h12_mag(:), h24_mag(:));
[~, p2] = vartest2(h24_mag(:), h48_mag(:));

disp("p12-24: " + p1 + "; p24-48" + p2);

figure
hold on
% swarmchart(repelem(1, height(h12_mag(:))), h12_mag(:), 100,'filled', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', [0.8 0.8 0.8])
% swarmchart(repelem(2, height(h24_mag(:))), h24_mag(:), 100,'filled', 'MarkerFaceColor', [0.55 0.55 0.55], 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', [0.55 0.55 0.55])
% swarmchart(repelem(3, height(h48_mag(:))), h48_mag(:), 100,'filled', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', [0.3 0.3 0.3])

violinplot(repelem(1, height(h12_mag(:))), h12_mag(:))
violinplot(repelem(2, height(h24_mag(:))), h24_mag(:))
violinplot(repelem(3, height(h48_mag(:))), h48_mag(:))


boxchart(repelem(1, height(h12_mag(:))), h12_mag(:), 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 2);
boxchart(repelem(2, height(h24_mag(:))), h24_mag(:), 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 2);
boxchart(repelem(3, height(h48_mag(:))), h48_mag(:), 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 2);
hold off
xticks([1 2 3])
xticklabels({"12h", "24h", "48h"})
ylabel('Velocity (a.u)')
ylim([0 1.2])

%% plot velocity using barplot
[~, p1] = vartest2(h12_mag(:), h24_mag(:));
[~, p2] = vartest2(h24_mag(:), h48_mag(:));

disp("p12-24: " + p1 + "; p24-48" + p2);

sd12 = std(h12_mag(:))^2;
sd24 = std(h24_mag(:))^2;
sd48 = std(h48_mag(:))^2;

figure

bar([sd12 sd24 sd48])

b.FaceColor = 'flat';
b.CData = [0.8 0.8 0.8; 0.55 0.55 0.55; 0.3 0.3 0.3];



xticks([1 2 3])
xticklabels({"12h", "24h", "48h"})
ylabel('Velocity variance')


%% plot barplot for circular variance

% cirvar = [0.921 0.726 0.250];
cirvar = [0.904 0.921 0.250];

figure
hold on
b = bar(cirvar)
b.FaceColor = 'flat';
b.CData = [0.8 0.8 0.8; 0.55 0.55 0.55; 0.3 0.3 0.3];
hold off
ylim([0 1])

xticks([1 2 3])
xticklabels({"12h", "24h", "48h"})
ylabel('Circular variance')



%% FUNCTION
function [bdy, Bo] = wave_boundary_gaussfilt(bw, im_range, sigma, binarize_threshold, medfilt_sz)
    sm_area_thres = 7000;
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
        if run_medfilt
            % BW = imclose(BW, strel('disk', 10));
            BW = medfilt2(BW, [medfilt_sz medfilt_sz]);
        end


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