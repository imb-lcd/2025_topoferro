%
%   quantify cell elongation over self-organization
%

clearvars 
clc
close all

ROOTDIR = "E:/";

addpath D:/Spatiotemporal_analysis/code/
addpath D:/Spatiotemporal_analysis/code/util/

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% calculate elongation of the cells
rng(2);

% prefix = ["iNAP 200K_s1t" "c2_ORG"];
% path = "D:/Spatiotemporal_analysis/iNAP_cellpose/iNAP_mhalf_cellpose/";

prefix = ["iNAP_200K_s1t" "c2_ORG"];
path = "D:/Spatiotemporal_analysis/iNAP/iNAP_200k_iNAP1/";

framelist = [1 16 32 48];
nframe = length(framelist);

cell_shape = table;

max_cells = 10000; % the highest number of cells being segmented 

for fr = 1:nframe
    mask_file = path + "/" + prefix(1) + sprintf("%03d", framelist(fr)) + prefix(2) + "_cp_masks.png";
    mask = imread(mask_file);


    stats = regionprops(mask, {'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Area', 'Perimeter'});
    stats = struct2table(stats);

    stats.AxisRatio = stats.MajorAxisLength ./ stats.MinorAxisLength;
    stats.AspectRatio = stats.Area ./ stats.Perimeter;
    stats.Frame = repelem(framelist(fr), height(stats))';
    
    stats = stats(randperm(height(stats), max_cells), :);

    cell_shape = [cell_shape; stats];
end

%% statistical difference
[~, p12] = ttest2(cell_shape(cell_shape.Frame == 1,:).AxisRatio,  cell_shape(cell_shape.Frame == 16,:).AxisRatio, 'Tail', 'left', 'Vartype', 'unequal');
[~, p23] = ttest2(cell_shape(cell_shape.Frame == 16,:).AxisRatio, cell_shape(cell_shape.Frame == 32,:).AxisRatio, 'Tail', 'left', 'Vartype', 'unequal');
[~, p34] = ttest2(cell_shape(cell_shape.Frame == 32,:).AxisRatio, cell_shape(cell_shape.Frame == 48,:).AxisRatio, 'Tail', 'left', 'Vartype', 'unequal');

disp(p12 + " " + p23 + " " + p34)

%% plot the figures
figure
hold on
% swarmchart(categorical(cell_shape.Frame), cell_shape.AxisRatio, 30, 'MarkerFaceColor', '#808080', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', '#808080')
violinplot(categorical(cell_shape.Frame), cell_shape.AxisRatio)
boxchart(categorical(cell_shape.Frame), cell_shape.AxisRatio, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k');
xticklabels({"12", "24", "36", "48"})
xlabel('Time (h)')
ylabel('Aspect ratio')
ylim([0 11])
% title('AxisRatio')
hold off


%% calculate defect centers and aligned region
prefix = ["iNAP_200k_s1t" "c2_ORG"];

fr = 1;
framelist = [48];

% load mask file
path = "D:/Spatiotemporal_analysis/iNAP/iNAP_200k_iNAP1/";

mask_file = path + "/" + prefix(1) + sprintf("%03d", framelist(fr)) + prefix(2) + "_cp_masks.png";
mask = imread(mask_file);

% load entropy file
prefix = ["iNAP_200k_s1t" "c3_ORG"];

path = "D:/Spatiotemporal_analysis/iNAP/iNAP_200k_dic/";

ent_file = path + "/" + prefix(1) + sprintf("%03d", framelist(fr)) + prefix(2) + "_ent.mat";
ent = struct2array(load(ent_file));
ent = ent(2:end-1, 2:end-1);

ent_cutoff = mean(ent(:), 'omitnan') - std(ent(:), 'omitnan');
% -2.8012; -2.5617;

%% obtain statistics
mask_cropped = imcrop(mask, [26 26 5948-2 5948-2]);

stats = regionprops(mask_cropped, ent, {'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'MeanIntensity'});
stats = struct2table(stats);

stats.AxisRatio = stats.MajorAxisLength ./ stats.MinorAxisLength;
stats.Frame = repelem(framelist(fr), height(stats))';

%%
% load defect file
prefix = ["iNAP_200k_s1t" "c3_ORG"];

path = "D:/Spatiotemporal_analysis/iNAP/iNAP_200k_dic/";

phalf_file = path + "/" + prefix(1) + sprintf("%03d", framelist(fr)) + prefix(2) + "_phalf.txt";
phalf = readmatrix(phalf_file);

mhalf_file = path + "/" + prefix(1) + sprintf("%03d", framelist(fr)) + prefix(2) + "_mhalf.txt";
mhalf = readmatrix(mhalf_file);

figure;
imshow(mask, [0 1])
hold on
% scatter(mhalf(:, 1), mhalf(:, 2), 'blue', 'filled')
scatter(phalf(:, 1), phalf(:, 2), 'red', 'filled')
hold off

%% get mask based on phalf and mhalf

ws = 300;

mhalf_mask = zeros(size(mask));

for i = 1:height(mhalf)
    x = mhalf(i, 1);
    y = mhalf(i, 2);

    x_min = max(1, round(x - ws/2));
    x_max = min(size(mask, 2), round(x + ws/2));
    y_min = max(1, round(y - ws/2));
    y_max = min(size(mask, 1), round(y + ws/2));

    mhalf_mask(y_min:y_max, x_min:x_max) = 1;
end

phalf_mask = zeros(size(mask));
for i = 1:height(phalf)
    x = phalf(i, 1);
    y = phalf(i, 2);

    x_min = max(1, round(x - ws/2));
    x_max = min(size(mask, 2), round(x + ws/2));
    y_min = max(1, round(y - ws/2));
    y_max = min(size(mask, 1), round(y + ws/2));

    phalf_mask(y_min:y_max, x_min:x_max) = 1;
end

figure
nexttile
imshow(mhalf_mask)
nexttile
imshow(phalf_mask)

%% label cells as within defects and if they are aligned.
stats.mhalf_labels = false(height(stats), 1);
stats.phalf_labels = false(height(stats), 1);
stats.labels = repmat("others", height(stats), 1);
for i = 1:height(stats)
    x = round(stats.Centroid(i, 1));
    y = round(stats.Centroid(i, 2));

    if y < size(mask_cropped, 2) && x < size(mask_cropped, 1)
        if mhalf_mask(y, x) == 1
            stats.mhalf_labels(i) = true;
            stats.labels(i) = "mhalf";
        end
        if phalf_mask(y, x) == 1
            stats.phalf_labels(i) = true;
            stats.labels(i) = "phalf";
        end
        if stats.MeanIntensity(i) < ent_cutoff
            stats.labels(i) = "aligned";
        end
    end
end

stats.labels = categorical(stats.labels);

%% plot and calculate ttests
[~, p12] = ttest2(stats(stats.labels == "mhalf",:).AxisRatio,  stats(stats.labels == "aligned",:).AxisRatio, 'Tail', 'left', 'Vartype', 'unequal');
[~, p23] = ttest2(stats(stats.labels == "phalf",:).AxisRatio,  stats(stats.labels == "aligned",:).AxisRatio, 'Tail', 'left', 'Vartype', 'unequal');
[~, p34] = ttest2(stats(stats.labels == "others",:).AxisRatio, stats(stats.labels == "aligned",:).AxisRatio, 'Tail', 'left', 'Vartype', 'unequal');

disp(p12 + " " + p23 + " " + p34)


stats.labels = categorical(stats.labels, {'mhalf', 'phalf', 'aligned', 'others'}, 'Ordinal', true);

figure
hold on
% swarmchart(stats.labels, stats.AxisRatio, 30, 'MarkerFaceColor', '#808080', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', '#808080')
violinplot(stats.labels, stats.AxisRatio)
boxchart(stats.labels, stats.AxisRatio, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k');
ylabel("Aspect ratio")