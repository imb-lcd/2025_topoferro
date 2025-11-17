%
%   nuc_to_cyt_ratio_analysis_defects.m
%   Analysis and plots for nuc to cyt ratio analyses
%

clearvars
clc
close all

ROOTDIR = "E:/IF/";

addpath D:/Spatiotemporal_analysis/code/
addpath D:/Spatiotemporal_analysis/code/util/

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%%   load nc ratio files

if_type = categorical("NRF2_KI696"); % or "NRF2"

% for YAP-1208-100k 550/2

if if_type == "YAP"
    all_prefix = ["YAP-1018-day1-tile" "YAP-1018-day2-tile" "YAP-1018-day3-tile" "YAP-1208-100k"];
    
    load('E:\IF\IF_results\YAP-1018-dayAll_nc.mat');     % load nc ratio on the nucleus positions
    load('E:\IF\IF_results\YAP-1018-dayAll_nc_nuc.mat'); % load nc ratio as a list in the order of segmented nucleus ids

    mycmap = create_bwr_cmap();

    def_pos = {[], [3317 359; 1479 240], [1229.5 3221.5; 851.5 1757.5; 3725.5 4589.5; 2783.5 3467.5; 761.5 599.5; 4499.5 497.5], ...
        [1356 1992; 1943.5 2663.5; 1481.5 3371.5; 3005.5 2291.5; 2999.5 2945.5; 1457.5 179.5; 4313.5 2501.5]};
    phalf_pos = {[], [], [], [1356 1992; 3005.5 2291.5; 2999.5 2945.5; 4313.5 2501.5]};
    mhalf_pos = {[], [], [], [1943.5 2663.5; 1481.5 3371.5; 1457.5 179.5]};

    ws = 550/2; 

    framelist = [10 6 1 1];

    y_lim = [0 3.5];
    c_lim = [0.4 1.6];

elseif if_type == "NRF2"
    all_prefix = ["NRF2-1125-100k-day1-tile" "NRF2-1125-100k-day2-tile" "NRF2-1125-100k-day3-tile"];
    
    load('E:\IF\IF_results\NRF2-1125-100k-dayAll_nc.mat');     % load nc ratio on the nucleus positions
    load('E:\IF\IF_results\NRF2-1125-100k-dayAll_nc_nuc.mat'); % load nc ratio as a list in the order of segmented nucleus ids

    mycmap = [0 0 0; parula(256)];

    def_pos = {[], [1659 4610; 3191 4009; 3590 1485; 193 1085; 2246 4610; 2637 3610; 4521,3235], [551 331; 963 4440; 3005 1104; 4589 8]};
    phalf_pos = {[], [1659 4610; 3191 4009; 3590 1485; 4521 3235], []};
    mhalf_pos = {[], [193 1085; 2246 4610; 2637 3610], []};

    ws = 585/2; 
    
    framelist = [7 8 7]; % [6 9 7];

    y_lim = [0 22];
    c_lim = [3.5 8];
elseif if_type == "YAP_TRULI"
    all_prefix = ["YAP-0625-TRULI"]; 

    load('E:\IF\IF_results\YAP-0625-TRULI_nc.mat');     % load nc ratio on the nucleus positions
    load('E:\IF\IF_results\YAP-0625-TRULI_nc_nuc.mat'); % load nc ratio as a list in the order of segmented nucleus ids
    
    mycmap = create_bwr_cmap();

    def_pos = {[361 3266; 3970,3039; 97,4974; 1552,904; 2116,1979; 4775,345]};
    phalf_pos = {[361 3266; 3970,3039;]};
    mhalf_pos = {[97,4974; 1552,904; 2116,1979; 4775,3457]};

    ws = 585/2; 
    
    framelist = [6];
    c_lim = [0.4 1.6];
    y_lim = [0 3.5];
elseif if_type == "NRF2_KI696"
    all_prefix = ["NRF2-0625-KI696"]; 

    load('E:\IF\IF_results\NRF2-0625-KI696_nc.mat');     % load nc ratio on the nucleus positions
    load('E:\IF\IF_results\NRF2-0625-KI696_nc_nuc.mat'); % load nc ratio as a list in the order of segmented nucleus ids
    
    mycmap = [0 0 0; parula(256)];

    def_pos = {[2678,192; 1408,1554; 3933,2193]};
    phalf_pos = {[2678 192; 3804 2288]};
    mhalf_pos = {[ 1408 1554 ;  3933 2193]};

    ws = 585/2; 
    
    framelist = [8];
    y_lim = [0 22];
    c_lim = [3.5 8];
elseif if_type == "YAP-0613-TRULI"
    all_prefix = ["YAP-0613-TRULI"]; 

    load('E:\IF\IF_results\YAP-0613-TRULI_nc.mat');     % load nc ratio on the nucleus positions
    load('E:\IF\IF_results\YAP-0613-TRULI_nc_nuc.mat'); % load nc ratio as a list in the order of segmented nucleus ids
    
    mycmap = create_bwr_cmap();

    def_pos = {[266 2607;415,4543; 4899,1534]};
    phalf_pos = {[266 2607; 266 2607]};
    mhalf_pos = {[415,4543; 4899,1534]};

    ws = 585/2; 
    
    framelist = [4];
    c_lim = [0.4 1.6];
    y_lim = [0 3.5];
elseif if_type == "NRF2-0613-KI696-200k"
    all_prefix = ["NRF2-0613-KI696-200k"]; 

    load('E:\IF\IF_results\NRF2-0613-KI696_nc.mat');     % load nc ratio on the nucleus positions
    load('E:\IF\IF_results\NRF2-0613-KI696_nc_nuc.mat'); % load nc ratio as a list in the order of segmented nucleus ids
    
    mycmap = [0 0 0; parula(256)];

    def_pos = {[1304.5,2528.5; 3800.5,2499; 3296,202]};
    phalf_pos = {[1304.5,2528.5; 3800.5,2499]};
    mhalf_pos = {[3296,202; 3296,202]};

    ws = 585/2; 
    
    framelist = [6];
    y_lim = [0 22];
    c_lim = [3.5 8];
end

%% draw defects on nc ratio nuc plot

for p = 1:length(all_prefix)
    prefix = all_prefix{p};

    nc = all_nc{p};

    def = def_pos{p};
    phalf = phalf_pos{p};
    mhalf = mhalf_pos{p};

    figure
    imshow(nc, c_lim, 'Border', 'tight')
    colormap(mycmap)
end

%% load density

all_den = cell(length(all_prefix), 1);

for p = 1:length(all_prefix)
    prefix = all_prefix(p);
    den_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/" + prefix + "_zMaxc3_ORG_den220.mat";
    load(den_file);

    all_den{p} = den;
end

%% load entropy
all_ent = cell(length(all_prefix), 1);

ent_cutoff = nan(length(all_prefix), 1);

for p = 1:length(all_prefix)
    prefix = all_prefix(p);
    
    ent_file = ROOTDIR + "/" + prefix + "/c2_DIC_nematics/" + prefix + "_z" + sprintf("%02d", framelist(p)) + "c2_ORG_ent.mat";
    load(ent_file);

    all_ent{p} = H;
end

%% load defects
d = [];

if if_type == "YAP"
    pp = 4;
elseif if_type == "NRF2"
    pp = 2;
else 
    pp = 1;
end

for p = pp
    prefix = all_prefix{p};
    nc_nuc = all_nc_nuc{p};
    den = all_den{p};
    ent = imresize(all_ent{p}, [5000 5000]);
    
    seg_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/segment/" + prefix + "_zMaxc3_ORG_stardist.tif";
    seg = imread(seg_file);
    
    num_nuc = max(seg(:));
    
    % Compute centroids and round to nearest pixel
    % centroids = round(cell2mat(struct2cell(regionprops(seg, 'Centroid')))');
    props = regionprops(seg, 'Centroid');
    centroids = cat(1, props.Centroid);
    centroids = round(centroids);
    
    % Extract x, y from matrix m
    x = centroids(:,1);
    y = centroids(:,2);
    
    if strcmp(if_type, "NRF2")
        curr_def_pos = def_pos{p}';
    else
        curr_def_pos = def_pos{p};
        curr_phalf_pos = phalf_pos{p};
        curr_mhalf_pos = mhalf_pos{p};
    end
    
    % Initialize logical mask for all ROIs
    in_roi = false(size(nc_nuc,1), 1);
    in_phalf_roi = false(size(nc_nuc,1), 1);
    in_mhalf_roi = false(size(nc_nuc,1), 1);
    
    for i = 1:size(curr_def_pos, 1)
        x_min = curr_def_pos(i,1) - ws;
        x_max = curr_def_pos(i,1) + ws;
        y_min = curr_def_pos(i,2) - ws;
        y_max = curr_def_pos(i,2) + ws;
    
        % Check if each point is within the current ROI
        in_roi = in_roi | ((x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max));
    end
    
    for i = 1:size(curr_phalf_pos, 1)
        x_min = curr_phalf_pos(i,1) - ws;
        x_max = curr_phalf_pos(i,1) + ws;
        y_min = curr_phalf_pos(i,2) - ws;
        y_max = curr_phalf_pos(i,2) + ws;
    
        % Check if each point is within the current ROI
        in_phalf_roi = in_phalf_roi | ((x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max));
    end

    for i = 1:size(curr_mhalf_pos, 1)
        x_min = curr_mhalf_pos(i,1) - ws;
        x_max = curr_mhalf_pos(i,1) + ws;
        y_min = curr_mhalf_pos(i,2) - ws;
        y_max = curr_mhalf_pos(i,2) + ws;
    
        % Check if each point is within the current ROI
        in_mhalf_roi = in_mhalf_roi | ((x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max));
    end
    
    curr_den = den(sub2ind(size(den), y, x));
    curr_ent = ent(sub2ind(size(ent), y, x));

    d = [d; x y curr_den curr_ent nc_nuc double(in_roi) double(in_phalf_roi) double(in_mhalf_roi)];
end

d = array2table(d, 'VariableNames', {'X', 'Y', 'den', 'ent' 'nc', 'defects', 'phalf', 'mhalf'});

%% plot defects vs non-defects

[~, p] = ttest2(d(d.defects == 0,:).nc, d(d.phalf == 1,:).nc, 'Tail', 'left', 'Vartype',  'unequal')
[~, p] = ttest2(d(d.defects == 0,:).nc, d(d.mhalf == 1,:).nc, 'Tail', 'left', 'Vartype',  'unequal')

figure
nexttile
hold on
swarmchart(categorical(d.defects), d.nc, 15, 'Color', '#D95319');
boxchart(categorical(d.defects), d.nc, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k');
hold off
if if_type == "YAP"
    yline(1, '--')
end

xticklabels({"non-defects", "defects"})
ylim(y_lim)

nexttile
hold on
swarmchart(categorical(d.phalf), d.nc, 15, 'Color', '#D95319');
boxchart(categorical(d.phalf), d.nc, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k');
hold off
if if_type == "YAP"
    yline(1, '--')
end
xticklabels({"non-phalf", "phalf"})
ylim(y_lim)

nexttile
hold on
swarmchart(categorical(d.mhalf), d.nc, 15, 'Color', '#D95319');
boxchart(categorical(d.mhalf), d.nc, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k');
hold off
if if_type == "YAP"
    yline(1, '--')
end
xticklabels({"non-mhalf", "mhalf"})
ylim(y_lim)

%% set cutoffs for den and ent
p2 = 2; 
p3 = 3;

if if_type == "YAP"
    den_lwr = 22.5408;
    den_upr = 35.8129;
    ent_lwr = prctile(all_ent{p2}(:), 25);
    ent_upr = prctile(all_ent{p2}(:), 75);
elseif if_type == "NRF2"
    den_lwr = 23.0555;
    den_upr = 31.5634;
    ent_lwr = prctile(all_ent{p2}(:), 25);
    ent_upr = prctile(all_ent{p2}(:), 75);
elseif if_type == "YAP_TRULI"
    den_lwr = prctile(all_den{1}(:), 25);
    den_upr = prctile(all_den{1}(:), 75);
    ent_lwr = prctile(all_ent{1}(:), 25);
    ent_upr = prctile(all_ent{1}(:), 75);
elseif if_type == "NRF2_KI696"
    den_lwr = prctile(all_den{1}(:), 25);
    den_upr = prctile(all_den{1}(:), 75);
    ent_lwr = prctile(all_ent{1}(:), 25);
    ent_upr = prctile(all_ent{1}(:), 75);
elseif if_type == "YAP-0613-TRULI"
    den_lwr = prctile(all_den{1}(:), 25);
    den_upr = prctile(all_den{1}(:), 75);
    ent_lwr = prctile(all_ent{1}(:), 25);
    ent_upr = prctile(all_ent{1}(:), 75);
elseif if_type == "NRF2-0613-KI696-200k"
    den_lwr = prctile(all_den{1}(:), 25);
    den_upr = prctile(all_den{1}(:), 75);
    ent_lwr = prctile(all_ent{1}(:), 25);
    ent_upr = prctile(all_ent{1}(:), 75);
end

disp("den: " +den_lwr + " "+ den_upr + "; ent: " + ent_lwr + " " + ent_upr);

% categorize the cells by spatial properties labels
label = repmat("NA", height(d), 1);
label(d.defects == 1) = "defects";
label(d.phalf == 1) = "phalf";
label(d.mhalf == 1) = "mhalf";
label(d.defects == 0 & d.den > den_upr & d.ent > ent_upr) = "Hden-MAln";
label(d.defects == 0 & d.den > den_upr & d.ent < ent_lwr) = "Hden-Aln";
label(d.defects == 0 & d.den < den_lwr & d.ent > ent_upr) = "Lden-MAln";
label(d.defects == 0 & d.den < den_lwr & d.ent < ent_lwr) = "Lden-Aln";

d.label = categorical(label);

%% assign randomly select rows for background
rng(3, "twister");

rand_label = zeros(height(d), 1); 

ridx = randperm(height(d), sum((d.label=='phalf') | (d.label=='mhalf')));
rand_label(ridx) = 1;

d.rand_label = rand_label;

%% set new label for spatial properties in defect's range

def_den_upr = mean(d(d.defects == 1, :).den, 'omitnan') + 2*std(d(d.defects == 1, :).den, 'omitnan');
def_den_lwr = mean(d(d.defects == 1, :).den, 'omitnan') - 2*std(d(d.defects == 1, :).den, 'omitnan');

def_ent_upr = mean(d(d.defects == 1, :).ent, 'omitnan') + 2*std(d(d.defects == 1, :).ent, 'omitnan');
def_ent_lwr = mean(d(d.defects == 1, :).ent, 'omitnan') - 2*std(d(d.defects == 1, :).ent, 'omitnan');

d.def_sp_label = strings(height(d),1);

def_sp_idx = d.den > def_den_lwr & d.den < def_den_upr & d.ent > def_ent_lwr & d.ent < def_ent_upr & d.label ~= "defect";

d.def_sp_label(def_sp_idx) = "def_sp_range";

d.def_sp_label = categorical(d.def_sp_label);

%% calculate t tests
% compare with random label
% [~, p_nc1p] = ttest2(d(d.rand_label == 1,:).nc, d(d.label == "phalf",:).nc,    'Tail', 'left', 'Vartype',  'unequal');
% [~, p_nc1m] = ttest2(d(d.rand_label == 1,:).nc, d(d.label == "mhalf",:).nc,    'Tail', 'left', 'Vartype',  'unequal');
[~, p_nc1] = ttest2(d(d.rand_label == 1,:).nc, d(d.defects == 1,:).nc,         'Tail', 'left', 'Vartype',  'unequal');
% [~, p_nc2] = ttest2(d(d.rand_label == 1,:).nc, d(d.label == "Lden-MAln",:).nc, 'Tail', 'left', 'Vartype', 'unequal');
% [~, p_nc3] = ttest2(d(d.rand_label == 1,:).nc, d(d.label == "Lden-Aln",:).nc,  'Tail', 'left', 'Vartype',  'unequal');
% 
% [~, p_nc4] = ttest2(d(d.rand_label == 1,:).nc, d(d.label == "Hden-MAln",:).nc, 'Tail', 'left', 'Vartype', 'unequal');
% [~, p_nc5] = ttest2(d(d.rand_label == 1,:).nc, d(d.label == "Hden-Aln",:).nc,  'Tail', 'left', 'Vartype', 'unequal');
% % 
% disp("t test with random: " +p_nc1+" "+p_nc2+" "+p_nc3+" "+p_nc4+" "+p_nc5)

[~, p_nc2] = ttest2(d(d.rand_label == 1,:).nc, d(d.label == "Lden-MAln" | d.label == "Lden-Aln",:).nc, 'Tail', 'right', 'Vartype', 'unequal');
[~, p_nc3] = ttest2(d(d.rand_label == 1,:).nc, d(d.label == "Hden-MAln" | d.label == "Hden-Aln",:).nc, 'Tail', 'right', 'Vartype', 'unequal');
disp("t test with random: " +p_nc1+" "+p_nc2+" "+p_nc3)

% compare lden-maln vs lden-aln
[~, p_nc6] = ttest2(d(d.label == "Lden-MAln",:).nc, d(d.label == "Lden-Aln",:).nc, 'Tail', 'right', 'Vartype',  'unequal');
[~, p_nc7] = ttest2(d(d.label == "Hden-MAln",:).nc, d(d.label == "Hden-Aln",:).nc, 'Tail', 'right', 'Vartype',  'unequal');

disp("t test maln vs aln; lden: " + p_nc6 + "; hden: " + p_nc7);

% compare df_sp_range with defect and random
[~, p_nc8] = ttest2(d(d.def_sp_label == "def_sp_range",:).nc, d(d.defects == 1,:).nc,  'Tail', 'left', 'Vartype',  'unequal');
[~, p_nc9] = ttest2(d(d.def_sp_label == "def_sp_range",:).nc, d(d.rand_label == 1,:).nc,    'Tail', 'right', 'Vartype',  'unequal');

disp("t test def_sp_range vs defect and rand: " + p_nc8 + " " + p_nc9);

% compare LD with HD, defect
[~, p_nc10] = ttest2(d(d.label == "Lden-MAln" | d.label == "Lden-Aln",:).nc, d(d.label == "Hden-MAln" | d.label == "Hden-Aln",:).nc,    'Tail', 'right', 'Vartype',  'unequal');
[~, p_nc11] = ttest2(d(d.label == "Lden-MAln" | d.label == "Lden-Aln",:).nc, d(d.defects == 1,:).nc,    'Tail', 'right', 'Vartype',  'unequal');
disp("t test LD: HD: " + p_nc10 + "; defect: " + p_nc11);

%% compare with single data
rand_median = 1.3167; %1.7544; 
rand_median = 5.7166; % NRF2

[~, p1] = ttest(d(d.defects == 1,:).nc, rand_median, 'Tail', 'left');
[~, p2] = ttest(d(d.label == "Lden-MAln" | d.label == "Lden-Aln",:).nc, rand_median, 'Tail', 'left');
[~, p3] = ttest(d(d.label == "Hden-MAln" | d.label == "Hden-Aln",:).nc, rand_median, 'Tail', 'left');

disp("t test to a single data value: " + p1 + " " + p2 + " " + p3)

%% calculate cohen's D
cohen1 = meanEffectSize(d(d.label == "Lden-MAln" | d.label == "Lden-Aln",:).nc, d(d.label == "Hden-MAln" | d.label == "Hden-Aln",:).nc, "Effect", "cohen");
cohen2 = meanEffectSize(d(d.label == "Lden-MAln" | d.label == "Lden-Aln",:).nc, d(d.defects == 1,:).nc, "Effect", "cohen");
cohen3 = meanEffectSize(d(d.label == "Lden-MAln" | d.label == "Lden-Aln",:).nc, d(d.rand_label == 1,:).nc, "Effect", "cohen");

disp("cohen's with LD: "+cohen1.Effect+" "+cohen2.Effect+" "+cohen3.Effect)

cohen4 = meanEffectSize(d(d.label == "Lden-MAln",:).nc, d(d.label == "Lden-Aln",:).nc, 'Effect', 'cohen');
cohen5 = meanEffectSize(d(d.label == "Hden-MAln",:).nc, d(d.label == "Hden-Aln",:).nc, 'Effect', 'cohen');

disp("cohen's aln-maln within density: "+cohen4.Effect+" "+cohen5.Effect)

cohen6 = meanEffectSize(d(d.label == "phalf",:).nc, d(d.label == "mhalf",:).nc, 'Effect', 'cohen');

disp("cohen's defect: "+cohen6.Effect)


%%
group_phalf = d(d.label == "phalf",:);
group_mhalf = d(d.label == "mhalf",:);

group_defrange = d(d.def_sp_label == "def_sp_range",:);

group_ldenmaln = d(d.label == "Lden-MAln",:);
group_ldenaln  = d(d.label == "Lden-Aln",:);

group_hdenmaln = d(d.label == "Hden-MAln",:);
group_hdenaln  = d(d.label == "Hden-Aln",:);

group_rand = d(d.rand_label == 1, :);

%%
figure
clf
hold on
swarmchart(ones(height(group_defrange), 1)*-1, group_defrange.nc, 'filled', 'MarkerFaceColor', '#808080')
swarmchart(ones(height(group_phalf), 1)*0,     group_phalf.nc,   'filled', 'MarkerFaceColor', '#808080')
swarmchart(ones(height(group_mhalf), 1)*1,     group_mhalf.nc,   'filled', 'MarkerFaceColor', '#808080')
swarmchart(ones(height(group_ldenmaln), 1)*2, group_ldenmaln.nc, 'filled', 'MarkerFaceColor', '#808080')
swarmchart(ones(height(group_ldenaln), 1)*3,  group_ldenaln.nc,  'filled', 'MarkerFaceColor', '#808080')
swarmchart(ones(height(group_hdenmaln), 1)*4, group_hdenmaln.nc, 'filled', 'MarkerFaceColor', '#808080')
swarmchart(ones(height(group_hdenaln), 1)*5,  group_hdenaln.nc,  'filled', 'MarkerFaceColor', '#808080')
swarmchart(ones(height(group_rand), 1)*6,     group_rand.nc,     'filled', 'MarkerFaceColor', '#808080')

boxchart(ones(height(group_defrange), 1)*-1, group_defrange.nc, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 1.5);
boxchart(ones(height(group_phalf), 1)*0,     group_phalf.nc,   'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 1.5);
boxchart(ones(height(group_mhalf), 1)*1,     group_mhalf.nc,   'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 1.5);
boxchart(ones(height(group_ldenmaln), 1)*2,  group_ldenmaln.nc, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 1.5);
boxchart(ones(height(group_ldenaln), 1)*3,   group_ldenaln.nc,  'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 1.5);
boxchart(ones(height(group_hdenmaln), 1)*4,  group_hdenmaln.nc, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 1.5);
boxchart(ones(height(group_hdenaln), 1)*5,   group_hdenaln.nc,  'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 1.5);
boxchart(ones(height(group_rand), 1)*6,      group_rand.nc,     'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'LineWidth', 1.5);
hold off
if if_type == "YAP" || if_type == "YAP_TRULI"
    yline(1, '--', 'LineWidth', 1.5)
else
    yline(median(group_rand.nc, 'omitnan'), '--', 'LineWidth', 1.5)

    y_upr = prctile(d.nc, 75) + 1.5*(prctile(d.nc, 75) - prctile(d.nc, 25));
    y_lwr = max(0, prctile(d.nc, 25) - 1.5*(prctile(d.nc, 75) - prctile(d.nc, 25)));

    ylim([y_lwr y_upr])
end
xticks([-1:6])
xticklabels({'Def range', 'phalf', 'mhalf', 'LDenMAln', 'LDenAln', 'HDenMAln', 'HDenAln', 'Rand'});
xtickangle(0)

%% show split violin plot
figure
hold on

% defects
violinplot(ones(height(group_phalf),1)*1,  group_phalf.nc, DensityDirection="negative")
violinplot(ones(height(group_mhalf),1)*1,  group_mhalf.nc,  DensityDirection="positive")

boxchart(ones(height(group_phalf),1)*0.95, group_phalf.nc, 'BoxFaceColor', "#0072BD", 'BoxFaceAlpha', 0.7, 'BoxWidth', 0.04, 'MarkerStyle', 'none');
boxchart(ones(height(group_mhalf),1)*1.05, group_mhalf.nc, 'BoxFaceColor', '#D95319', 'BoxFaceAlpha', 0.7, 'BoxWidth', 0.04, 'MarkerStyle', 'none');

% low density
violinplot(ones(height(group_ldenmaln),1)*2,  group_ldenmaln.nc, DensityDirection="negative")
violinplot(ones(height(group_ldenaln),1)*2,   group_ldenaln.nc,  DensityDirection="positive")

boxchart(ones(height(group_ldenmaln),1)*1.95, group_ldenmaln.nc, 'BoxFaceColor', "#0072BD", 'BoxFaceAlpha', 0.7, 'BoxWidth', 0.04, 'MarkerStyle', 'none');
boxchart(ones(height(group_ldenaln),1)*2.05,  group_ldenaln.nc,  'BoxFaceColor', '#D95319', 'BoxFaceAlpha', 0.7, 'BoxWidth', 0.04, 'MarkerStyle', 'none');

% high density
violinplot(ones(height(group_hdenmaln),1)*3,  group_hdenmaln.nc, DensityDirection="negative")
violinplot(ones(height(group_hdenaln),1)*3,   group_hdenaln.nc,  DensityDirection="positive")

boxchart(ones(height(group_hdenmaln),1)*2.95, group_hdenmaln.nc, 'BoxFaceColor', "#0072BD", 'BoxFaceAlpha', 0.7, 'BoxWidth', 0.04, 'MarkerStyle', 'none');
boxchart(ones(height(group_hdenaln),1)*3.05,  group_hdenaln.nc,  'BoxFaceColor', '#D95319', 'BoxFaceAlpha', 0.7, 'BoxWidth', 0.04, 'MarkerStyle', 'none');

% random
violinplot(ones(height(group_rand), 1)*4,  group_rand.nc)
boxchart(ones(height(group_rand), 1)*4,    group_rand.nc,  'BoxFaceColor', '#D95319', 'BoxFaceAlpha', 0.7, 'BoxWidth', 0.1, 'MarkerStyle', 'none');

if if_type == "YAP" || if_type == "YAP_TRULI" || if_type == "YAP-0613-TRULI"
    yline(1, '--', 'LineWidth', 1.5)
else
    yline(median(group_rand.nc, 'omitnan'), '--', 'LineWidth', 1.5)
    % ylim([2.5 9])
    y_upr = prctile(d.nc, 75) + 1.5*(prctile(d.nc, 75) - prctile(d.nc, 25));
    y_upr = round(y_upr/10)*10;
    y_lwr = max(0, prctile(d.nc, 25) - 1.5*(prctile(d.nc, 75) - prctile(d.nc, 25)));
    y_lwr = floor(y_lwr/10)*10;

    ylim([y_lwr 25])
end
xlim([0 5])
xticks([1:4])
xticklabels({'defects', 'low density', 'high density', 'rand'});
hold off

%% FUNCTIONS
function [lwr, upr] = iqr_outlier(d)
    q1 = prctile(d, 25);
    q3 = prctile(d, 75);

    r = iqr(d);

    lwr = q1 - 1.5 * r;
    upr = q3 + 1.5 * r;
end

function bwr = create_bwr_cmap()
    % Create bwr colormap
    blue = [0, 0, 1];
    white = [1, 1, 1];
    red = [1, 0, 0];
    
    nSteps = 25;
    
    bwr = [linspace(blue(1), white(1), nSteps)' linspace(blue(2), white(2), nSteps)' linspace(blue(3), white(3), nSteps)'; ...
            linspace(white(1), red(1), nSteps)' linspace(white(2), red(2), nSteps)' linspace(white(3), red(3), nSteps)'];
    bwr = [0 0 0; bwr];
end