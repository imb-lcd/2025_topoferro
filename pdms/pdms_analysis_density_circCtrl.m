%
%   pdms_analysis_density_circCtrl.m
%
clearvars
clc
close all

addpath("D:/Spatiotemporal_analysis/code/util/");

ROOTDIR = "D:/Spatiotemporal_analysis/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATION
all_prefix = ["pdms-0404", "pdms-0404", "pdms-0418", "pdms-0418"];

welllist = [10 12 20 23];
fr = 1;
imsize = 5848;

n_prefix = length(all_prefix);
n_wells  = length(welllist);

%% load pdms positions
pdms_xy = [1580 1755;
           1950 1940
           1750 1830
           1910 2110];
pdms_r  = 2075/2;

%% load orientation and density
den_all  = cell(n_prefix, n_wells);

for p = 1:n_prefix
    prefix = all_prefix(p);
    well = sprintf("%02d", welllist(p));
    disp(prefix + " " + well)

    PATH = ROOTDIR + "/wave_" + prefix + "/Well" + well + "/";
    den_path = PATH + "/c1_mCherry_density/density_reg120/";
        
    den_file = den_path + prefix+"_s"+well+"t"+sprintf("%02d",fr)+"c1_ORG_den120.mat";
    load(den_file);
    den = imresize(den, [imsize imsize]);

    den_all{p} = den;
end

%% load initiation of arms

arm_file = "D:\Spatiotemporal_analysis\wave_0_analyses\wave_pdms\" + ...
    "pdms_arm_ini_summary.txt";

arm = readtable(arm_file, "Delimiter", "\t");
arm.experiment = categorical(arm.experiment);

%% calculate density
showfigure = false;

sigma = 15;
binarize_threshold = 0.5; 

num_images = 33;

wv_den_mean = nan(n_prefix, 7, num_images);
wv_den_std  = nan(n_prefix, 7, num_images);

ctrl_den_mean = nan(n_prefix, 7, num_images);
ctrl_den_std  = nan(n_prefix, 7, num_images);

% precompute PDMS grid
[Xgrid, Ygrid] = meshgrid(1:imsize, 1:imsize);

for p = 1:n_prefix
    prefix = all_prefix(p);

    well = sprintf("%02d", welllist(p));
    disp(prefix + " " + well)
    
    % define pdms mask
    pdms_center = pdms_xy(p, :) + pdms_r;   % [x_center, y_center]
    pdms_mask = (Xgrid - pdms_center(1)).^2 + (Ygrid - pdms_center(2)).^2 <= pdms_r^2;


    fname = ROOTDIR + "/wave_" + prefix + "/Well" + well + "/c2_cy5_adj/" + ...
        prefix + "_s"+well+"tAllc2_ORG_bw_adj.tif";

    info = imfinfo(fname);

    % keep two versions of density
    den_raw = den_all{p};       % keeps PDMS region
    den = den_raw;              % used for arm/control analysis
    den(pdms_mask) = nan;       % exclude PDMS after frame 10

    curr_arm = arm(arm.experiment == prefix & arm.well == welllist(p), :);

    for i = 1:height(curr_arm)
        disp(i)
        a_ini = [curr_arm.arm_x(i) curr_arm.arm_y(i)];
        a_idx = curr_arm.arm(i);
        
        wv_arm_file = ROOTDIR + "/wave_" + prefix + "/Well" + well + "/c2_cy5_adj/" + ...
            prefix + "_s"+well+"tAllc2_ORG_bw_adj_arm" + a_idx + ".tif";

        for k = 3:num_images
            %%% define wave region
            bw = imread(fname, k);
            
            [~, Bo] = wave_boundary_gaussfilt(bw, 1, sigma, binarize_threshold);
            bo = Bo{1};
            wv_mask = poly2mask(bo(:, 2), bo(:, 1), imsize, imsize);


            if k <= 10
                wv_region = wv_mask & ~isnan(den_raw);

                wv_vals = den_raw(wv_region);

                wv_den_mean(p, i, k) = mean(wv_vals, 'omitnan');
                wv_den_std(p, i, k)  = std(wv_vals,  'omitnan');

                % no arm-based control region before frame 10
                ctrl_den_mean(p, i, k) = nan;
                ctrl_den_std(p, i, k)  = nan;

                if showfigure
                    den_pdms = den_raw;
                    den_pdms(~pdms_region) = nan;

                    figure
                    imshow(den_pdms, [1 160], 'Border', 'tight')
                    colormap(gca, jet)
                    title("Density inside PDMS region")
                end

                continue
            end

            % get the wave arm mask
            bw = imread(wv_arm_file, k-10);
            [~, Bo] = wave_boundary_gaussfilt(bw, 1, sigma, binarize_threshold);
            bo = Bo{1};

            arm_mask = poly2mask(bo(:, 2), bo(:, 1), imsize, imsize);

            %%% density inside the wave
            wv_region = arm_mask & ~isnan(den);

            wv_vals = den(wv_region);

            wv_den_mean(p, i, k) = mean(wv_vals, 'omitnan');
            wv_den_std(p, i, k)  = std(wv_vals,  'omitnan');

            %%% control region outside the wave. 
            [arm_y, arm_x] = find(arm_mask);

            if isempty(arm_x)
                continue
            end

            % distance from initiation center to all pixels in wave mask
            dist_to_ini = sqrt((arm_x - a_ini(1)).^2 + ...
                               (arm_y - a_ini(2)).^2);

            % farthest distance defines control-circle radius
            ctrl_r = max(dist_to_ini);

                        % circle centered at initiation site
            ctrl_circle = (Xgrid - a_ini(1)).^2 + ...
                          (Ygrid - a_ini(2)).^2 <= ctrl_r^2;

            % control region = inside circle, not wave, not PDMS/NaN
            ctrl_region = ctrl_circle & ~wv_mask & ~isnan(den);

            ctrl_vals = den(ctrl_region);

            ctrl_den_mean(p, i, k) = mean(ctrl_vals, 'omitnan');
            ctrl_den_std(p, i, k)  = std(ctrl_vals,  'omitnan');

            if showfigure
                den_wv = den;
                den_wv(~wv_region) = nan;
            
                den_ctrl = den;
                den_ctrl(~ctrl_region) = nan;
            
                figure
                tiledlayout(1, 3)
            
                nexttile
                imshow(den, [1 160], 'Border', 'tight')
                colormap(gca, jet)
                title("Density")
            
                nexttile
                imshow(den_wv, [1 160], 'Border', 'tight')
                colormap(gca, jet)
                title("Density in selected arm")
            
                nexttile
                imshow(den_ctrl, [1 160], 'Border', 'tight')
                colormap(gca, jet)
                title("Density in control region")
            end
        end
    end
end


%% plot densiyt over time

savefile = true;

% x-axis frames
frames = 1:num_images;

sm_factor_wv = 2;
sm_factor_ctrl = 2;

% mean of means across prefix and arms
wv_mean_overall   = squeeze(mean(wv_den_mean, [1 2], 'omitnan'));
ctrl_mean_overall = squeeze(mean(ctrl_den_mean, [1 2], 'omitnan'));


frames = 1:num_images;

wv_mean_raw = wv_mean_overall;      % keep original for reference
wv_mean_plot = wv_mean_overall;     % this will be used for plotting

% region to replace around kink
idx_replace = 10:11;

% anchor points before and after kink
idx_anchor = [8 9 12 13];

% only use anchors that are valid
valid_anchor = idx_anchor(~isnan(wv_mean_overall(idx_anchor)));

if numel(valid_anchor) >= 2
    wv_mean_plot(idx_replace) = interp1( ...
        frames(valid_anchor), ...
        wv_mean_overall(valid_anchor), ...
        frames(idx_replace), ...
        'pchip');
end


wv_mean_overall   = smoothdata(wv_mean_overall, 'movmean', sm_factor_wv);
ctrl_mean_overall = smoothdata(ctrl_mean_overall, 'movmean', sm_factor_ctrl);

% mean of standard deviations across prefix and arms
wv_shade = squeeze(std(wv_den_mean, 0, [1 2], 'omitnan'));
ctrl_shade = squeeze(std(ctrl_den_mean, 0, [1 2], 'omitnan'));


wv_shade_plot = wv_shade;     % this will be used for plotting

% region to replace around kink
idx_replace = 11:14;

% anchor points before and after kink
idx_anchor = [9 10 15 16];

% only use anchors that are valid
valid_anchor = idx_anchor(~isnan(wv_shade(idx_anchor)));

if numel(valid_anchor) >= 2
    wv_shade_plot(idx_replace) = interp1( ...
        frames(valid_anchor), ...
        wv_shade(valid_anchor), ...
        frames(idx_replace), ...
        'pchip');
end

wv_shade   = smoothdata(wv_shade_plot, 'movmean', sm_factor_wv+10);
ctrl_shade = smoothdata(ctrl_shade, 'movmean', sm_factor_ctrl);


% valid points
wv_valid   = ~isnan(wv_mean_overall)   & ~isnan(wv_shade);
ctrl_valid = ~isnan(ctrl_mean_overall) & ~isnan(ctrl_shade);

% plot
figure
hold on

% wave shaded region
x_wv = frames(wv_valid);
y_wv = wv_mean_overall(wv_valid);
s_wv = wv_shade(wv_valid);

fill([x_wv fliplr(x_wv)], ...
     [(y_wv - s_wv)' fliplr((y_wv + s_wv)')], ...
     [0.75 0.75 0.75], ...
     'FaceAlpha', 0.35, ...
     'EdgeColor', 'none');

% control shaded region
x_ctrl = frames(ctrl_valid);
y_ctrl = ctrl_mean_overall(ctrl_valid);
s_ctrl = ctrl_shade(ctrl_valid);

fill([x_ctrl fliplr(x_ctrl)], ...
     [(y_ctrl - s_ctrl)' fliplr((y_ctrl + s_ctrl)')], ...
     [0 0 0.65], ...
     'FaceAlpha', 0.25, ...
     'EdgeColor', 'none');

% mean lines
plot(frames, wv_mean_plot, 'k-', 'LineWidth', 2);
plot(frames, ctrl_mean_overall, 'b-', 'LineWidth', 2);

xlabel('Frame')
ylabel('Density')

ylim([60 110])


if savefile
    outfile = "D:\Spatiotemporal_analysis\wave_0_analyses\wave_pdms\" + ...
        "pdms_density_time_ctrl";
    exportgraphics(gcf, outfile + ".jpg");
    exportgraphics(gcf, outfile + ".pdf");
    savefig(outfile + ".fig");
end

%% FUNCTION
function [bdy, Bo] = wave_boundary_gaussfilt(bw, nframes, sigma, binarize_threshold, medfilt_sz)
    sm_area_thres = 3000;
    sm_hole_thres = 10000; %50000;

    bdy = cell(nframes, 1);
    Bo = cell(nframes, 1);

    run_medfilt = 1;
    if ~exist('medfilt_sz', 'var')
        run_medfilt = 0;
    end

    for i = 1:nframes
        B = imgaussfilt(double(bw), sigma);

        BW = imbinarize(B, binarize_threshold);    
        BW = bwareaopen(BW, sm_area_thres);         % Remove small objects and speckles
        BW = ~bwareaopen(~BW, sm_hole_thres); 
        BW = imclose(BW, strel('disk', 10)); 

        
        if run_medfilt
            
            BW = medfilt2(BW, [medfilt_sz medfilt_sz]);
        end

        bdy{i} = sparse(BW);
       
        bo = bwboundaries(BW);
        areas = cellfun(@(b) polyarea(b(:,2), b(:,1)), bo);
        
        % Find index of largest area
        [~, idx] = max(areas);
        
        % Keep only the largest boundary
        bo = bo{idx};

        bo(:, 2) = smoothdata(bo(:, 2), 'sgolay', 75)';
        bo(:, 1) = smoothdata(bo(:, 1), 'sgolay', 75)';

        Bo{i} = bo;
    end
end

function myparula = generate_black_parula(steps_black)
    if ~exist('steps_black')
        steps_black = 25;
    end
    % steps_black = 200;

    % Number of steps for each segment
    n_black_to_blue = steps_black;
    n_parula = 412;
    
    % Get the parula colormap
    parula_map = parula(n_parula);
    
    % Identify the dark blue of parula (first color)
    parula_start = parula_map(1, :);  % usually a dark blue
    
    % Interpolate from black to parula's blue
    black = [0/256 0/256 0/256];
    ramp = [linspace(black(1), parula_start(1), n_black_to_blue)', ...
            linspace(black(2), parula_start(2), n_black_to_blue)', ...
            linspace(black(3), parula_start(3), n_black_to_blue)'];
    
    % Concatenate the colormaps
    myparula = [ramp; parula_map];
end