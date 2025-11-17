%
%   nuc_to_cyt_ratio_analysis.m
%   Analysis and plots for nuc to cyt ratio analyses
%

clearvars
clc
close all

ROOTDIR = "E:/IF/";

addpath D:/Spatiotemporal_analysis/code/
addpath D:/Spatiotemporal_analysis/code/util/

addpath D:/Matlab_FileExchange/scatter_kde

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% load ratio summary 
if_type = categorical("NRF2"); % or "YAP

if if_type == "YAP"
    % all_prefix = ["YAP-1018-day1-tile" "YAP-1018-day2-tile" "YAP-1018-day3-tile"]; % "YAP-1208-100k"];
    all_prefix = ["YAP-0625-Control" "YAP-0625-TRULI"]; % "YAP-1208-100k"];
    ylab = "Percent of YAP Nucleus Localization";
    mycmap = generate_bwr();
    c_lim = [0.4 1.6];
elseif if_type == "NRF2"
    % all_prefix = ["NRF2-1125-100k-day1-tile" "NRF2-1125-100k-day2-tile" "NRF2-1125-100k-day3-tile"];
    all_prefix = ["NRF2-0625-Control" "NRF2-0625-KI696"]; % "YAP-1208-100k"];
    ylab = "Percent of NRF2 Nucleus Intensity";
    mycmap = [0 0 0; parula(256)];
    c_lim = [3.5 8];
end

mycol = [0.8 0.8 0.8];
mycol = [mycol; mycol-0.4; mycol-0.8];

den_bin = 5;
coh_bin = 0.01;

days = cell(length(all_prefix), 1);

for p = 1:length(all_prefix)
    prefix = all_prefix(p);

    ratio_file = ROOTDIR + "/" + prefix + "/c1_sig_adj/" + prefix + "_ratio_summary_ilastik_filter_sp.txt";
    d = readtable(ratio_file, "Delimiter", "\t");
    d.den_bin = round(d.den/den_bin)*den_bin;
    d.coh_bin = round(d.coh/coh_bin)*coh_bin;
    d(d.filter == 1 | isnan(d.ratio_top50), :) = [];

    days{p} = d;
end

%% defie ranges for kernel density and boxplots
x_range = [0 60]; 
y_range_ratio = [0 3]; 

if strcmp(if_type, "YAP")
    y_range_intensity = [0 25]; 
elseif strcmp(if_type, "NRF2")
    y_range_intensity = [0 12.5]; 
end

%% plot kernel density plots
figure
tiledlayout(length(days), 3, 'TileSpacing', 'compact', 'Padding', 'compact')
sgtitle(if_type)
for i = 1:length(days)
    d = days{i};
    sp = d.den;

    ax1 = nexttile;
    scatter_kde(sp, d.ratio);
    yline(1, "--")
    ylabel("ratio")
    ylim(y_range_ratio)
    title("ratio")

    ax2 = nexttile;
    scatter_kde(sp, d.nuc_intensity);
    ylabel("nuc intensity")
    ylim(y_range_intensity)
    title("nuc intensity")
    
    ax3 = nexttile;
    scatter_kde(sp, d.cyto_intensity);
    ylabel("cyto intensity")
    ylim(y_range_intensity)
    title("cyto intensity")

    linkaxes([ax1, ax2, ax3], 'x');
    xlim(x_range)
end

%% plot boxplots

bin_range = x_range(1):den_bin:x_range(2);

figure
tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact')
sgtitle(if_type)
for i = 1:length(days)
    d = days{i};

    missing_bins = setdiff(bin_range, unique(d.den_bin));
    y_dummy = NaN(size(missing_bins));

    x_all = [d.den_bin; missing_bins'];

    ax1 = nexttile;
    y_all = [d.ratio; y_dummy'];
    hold on
    swarmchart(categorical(x_all), y_all, 7, 'Color', '#D95319');
    boxchart(categorical(x_all), y_all, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k');
    yline(1,"--")
    hold off
    ylabel("ratio")
    ylim(y_range_ratio)
    title("ratio")
    
    ax2 = nexttile;
    y_all = [d.nuc_intensity; y_dummy'];
    hold on
    swarmchart(categorical(x_all), y_all, 7, 'Color', '#D95319');
    boxchart(categorical(x_all), y_all, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k');
    hold off
    ylabel("nuc intensity")
    ylim(y_range_intensity)
    title("nuc intensity")

    ax3 = nexttile;
    y_all = [d.cyto_intensity; y_dummy'];
    hold on
    swarmchart(categorical(x_all), y_all, 7, 'Color', '#D95319');
    boxchart(categorical(x_all), y_all, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k');
    hold off
    ylabel("cyto intensity")
    ylim(y_range_intensity)
    title("cyto intensity")
end

%% plot nc ratio between control and perturbation

figure

xticklabel = strings(2, 1);
hold on
for p = 1:length(all_prefix)
    d = days{p};
    if if_type == "YAP"
        violinplot(ones(height(d), 1)*p, d.ratio, 'DensityScale', 'count')
        boxchart(ones(height(d), 1)*p, d.ratio, 'MarkerStyle', 'none', 'BoxWidth', 0.3);
    elseif if_type == "NRF2"
        violinplot(ones(height(d), 1)*p, d.nuc_intensity, 'DensityScale', 'count')
        boxchart(ones(height(d), 1)*p, d.nuc_intensity, 'MarkerStyle', 'none', 'BoxWidth', 0.3);
        ylim([0 25])
    end
    xticklabel(p) = all_prefix(p);
end
xticks([1 2])
xticklabels(xticklabel)
hold off

ylim

%% calculate cohen's d

d1 = days{1};
d2 = days{2};

cohen = meanEffectSize(d1.ratio, d2.ratio, 'Effect', 'cohen');


%% load density

all_den = cell(length(all_prefix), 1);

for p = 1:length(all_prefix)
    prefix = all_prefix(p);
    den_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/" + prefix + "_zMaxc3_ORG_den220.mat";
    load(den_file);

    all_den{p} = den;
end

%% determine density cutoffs for percentage of barplots

all_d = all_den{2}(:);

d_lwr = prctile(all_d(:), 25); % mean(all_d, 'omitnan') - 1*std(all_d, 'omitnan');
d_upr = prctile(all_d(:), 75); % mean(all_d, 'omitnan') + 1*std(all_d, 'omitnan');

disp("den lwr: " + d_lwr + "; upr: " + d_upr)
% YAP:  den lwr: 25.9245; upr: 34.7454 
% NRF2: den lwr: 23.0555; upr: 31.5634

%% generate percentages barplot for YAP and NRF2

r_lwr = 0.9; % 0.95;
r_upr = 1.1; % 1.05;

all_nc = cell(length(all_prefix), 1);

figure
tiledlayout(length(all_prefix), 1)
sgtitle(if_type)
for i = 1:length(days)
    d = days{i};
    all_nc{i} = nan(3,2);

    r_cut = median(d.nuc_intensity);
    
    % calculate percentage for low density
    d_L = d(d.den < d_lwr, :);
    if strcmp(if_type, "YAP")
        [cyt, nuc] = calculate_percent_localization(if_type, d_L.ratio, [r_lwr r_upr]);
        all_nc{i}(1, :) = [cyt nuc];
    elseif strcmp(if_type, "NRF2")
        [low, high] = calculate_percent_localization(if_type, d_L.nuc_intensity, r_cut);
        all_nc{i}(1, :) = [low high];
    end

    % calculate percentage for mid density
    d_M = d(d.den >= d_lwr & d.den < d_upr, :);
    if strcmp(if_type, "YAP")
        [cyt, nuc] = calculate_percent_localization(if_type, d_M.ratio, [r_lwr r_upr]);
        all_nc{i}(2, :) = [cyt nuc];
    elseif strcmp(if_type, "NRF2")
        [low, high] = calculate_percent_localization(if_type, d_M.nuc_intensity, r_cut);
        all_nc{i}(2, :) = [low high];
    end

    % calculate percentage for high density
    d_H = d(d.den >= d_upr, :);
    if strcmp(if_type, "YAP")
        [cyt, nuc] = calculate_percent_localization(if_type, d_H.ratio, [r_lwr r_upr]);
        all_nc{i}(3, :) = [cyt nuc];
    elseif strcmp(if_type, "NRF2")
        [low, high] = calculate_percent_localization(if_type, d_H.nuc_intensity, r_cut);
        all_nc{i}(3, :) = [low high];
    end
    
    % record number of cells counted
    counts = [height(d_L) height(d_M) height(d_H)];


    nexttile
    bar(all_nc{i}, 'FaceColor', mycol(i, :))
    ylim([0 1])

    disp("day " + i + "; " + counts(1) + " " + counts(2) + " " + counts(3) + " " )
end
all_nc{1}
all_nc{2}

%% plot changes of percentages over time
nc_plot = nan(2,3);

for i = [1 3]
    nc = all_nc{i};

    % for low density
    nc_plot(1, i) = nc(1, 1);
    
    % for high density
    if i == 1
        nc_plot(2, i) = nc(2, 1);
    elseif i > 1
        % nc_plot(2, i) = mean([nc(2, 1) nc(3, 1)]);
        nc_plot(2, i) = nc(3, 1);
    end
end

figure
hold on
plot([1 3], nc_plot(1, [1 3]), '.-', 'Color', mycol(2, :), 'LineWidth', 2)
plot([1 3], nc_plot(2, [1 3]), '.-', 'Color', mycol(3, :), 'LineWidth', 2)
hold off
xticks([1:3])
xlim([0.5 3.5])
ylim([0 1])
xlabel("Day")
ylabel(ylab)

%% calculate chisq test of independence

x = [282 208; 564 622; 175 376];

row_sum = sum(x, 2);
col_sum = sum(x, 1);
grand_total = sum(x, 'all');

% Expected counts under independence
expected = (row_sum * col_sum) / grand_total;

% Chi-squared statistic
chi2_stat = sum((x - expected).^2 ./ expected, 'all');

% Degrees of freedom = (rows-1)*(cols-1) = (2-1)*(3-1) = 2
df = (size(x,1)-1)*(size(x,2)-1);

% p-value
p = 1 - chi2cdf(chi2_stat, df)

%% plot nc ratio on nucleus

mycmap = generate_ylbu();
c_lim = [3.5 8];

for p = [1 2]
    prefix = all_prefix{p};

    nc = all_nc{p};
    
    figure
    imshow(nc, c_lim, 'Border', 'tight')
    colormap(mycmap)
    hold on
end

%% FUNCTIONS
function bwr = generate_bwr()
    blue = [0, 0, 1];
    white = [1, 1, 1];
    red = [1, 0, 0];
    
    nSteps = 25;
    
    bwr = [linspace(blue(1), white(1), nSteps)' linspace(blue(2), white(2), nSteps)' linspace(blue(3), white(3), nSteps)'; ...
            linspace(white(1), red(1), nSteps)' linspace(white(2), red(2), nSteps)' linspace(white(3), red(3), nSteps)'];
    bwr = [0 0 0; bwr];
end

function ylbu = generate_ylbu()
    blue = [28, 117, 255];
    yellow = [249, 237, 50];
    
    nSteps = 25;
    
    ylbu = [linspace(blue(1), yellow(1), nSteps)' linspace(blue(2), yellow(2), nSteps)' linspace(blue(3), yellow(3), nSteps)'];

    ylbu = ylbu / 255;
    ylbu = [0 0 0; ylbu];
end

function [low, high] = calculate_percent_localization(if_type, r, threshold)
    if strcmp(if_type, "YAP")
        lwr = threshold(1);
        upr = threshold(2);

        % calculate percent in cytoplasm
        low = numel( r(r>upr)) / ( numel(r)- numel(r(r>=lwr & r<=upr)) );
        
        % calculate percent in nucleus
        high = numel( r(r<lwr)) / ( numel(r)- numel(r(r>=lwr & r<=upr)) );
        disp(numel(r)- numel(r(r>=lwr & r<=upr)) + " " + numel( r(r>upr)) + " " + numel( r(r<lwr)));
        
    elseif strcmp(if_type, "NRF2")
        r_cut = threshold(1);

        % calculate low intensity
        low = numel( r(r<r_cut)) / numel(r);

        % calculate high intensity
        high = numel( r(r>=r_cut)) / numel(r);
    end
end