%
%   quantify_deadcell_percent.m
%   quantify percentages of dead cells
%
clearvars 
clc
close all

ROOTDIR = "E:/";

addpath D:/Spatiotemporal_analysis/code/
addpath D:/Spatiotemporal_analysis/code/util/

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% load segment and cy5 images

all_prefix = ["seed05_day10_s24_diffday-1023-all", "seed10_day15_s02_diffday-1023-all" "mCh300FBS10-1"];
well = [24 2 10 14];
frames = [34 40 66 41];

nframes = max(frames);

th = [20000, 30000, 0]; % third is processed separately

num_death = nan(nframes, 4);
pct_death = nan(nframes, 4);
num_cell  = nan(nframes, 4);

%%

for p = 1:length(all_prefix)
% for p = 4
    prefix = all_prefix(p);

    seg_file = ROOTDIR + "/movies_diff_day/" + ...
        prefix + "_s" + sprintf("%02d", well(p)) + "tAllc2_ORG_stardist.tif";

    cy5_file = ROOTDIR + "/movies_diff_day/" + ...
        prefix + "_s" + sprintf("%02d", well(p)) + "tAllc2_ORG.tif";
        
    % for i = 1:frames(p)
    for i = 1:frames(p)
        seg = imread(seg_file, i);
        cy5 = imread(cy5_file, i);

        idx = setdiff(unique(seg), 0); % get all the ids in the label

        stats = regionprops(seg, cy5, {'MeanIntensity', 'Centroid'});
        stats = stats(idx);

        mean_int = [stats(:).MeanIntensity];
        xy = round(reshape([stats(:).Centroid], 2, []))';

        is_dead = mean_int > th(p);
        is_cell = true(size(is_dead));

        if i == 1
            num_death_1 = sum(is_dead & is_cell);
            num_death(i, p) = 0;

            num_cell(i, p) = max(sum(is_cell) - num_death_1, 0);
        else
            curr_death = max(sum(is_dead & is_cell) - num_death_1, 0);
            num_death(i, p) = max(curr_death, num_death(i-1, p));

            curr_cell = max(sum(is_cell) - num_death_1, 0);
            num_cell(i, p) = max(curr_cell, num_cell(i-1, p));
        end
 
        % pct_death(i, p) = num_death(i) / num_cell(i);
    end
end

%% specific process for third wave mCh300FBS10-1

t_start = 28;
for p = 3
    seg_file = "D:/Spatiotemporal_analysis/wave_mCh300FBS10-1/Well10/c2_cy5_adj/mCh300FBS10-1_s10tAllc2_ORG_stardist.tif";

    max_cell = 65535;

    seg = imread(seg_file, t_start); 
    num_death_1 = max(seg(:));

    for t = 1:t_start-1
        num_death(t, p) = 0;
        num_cell(t, p) = max_cell;

        pct_death(t, p) = num_death(t) / num_cell(t) ;
    end

    for t = t_start+1:frames(p)
        seg = imread(seg_file, t);

        num_death(t, p) = max(seg(:)) - num_death_1;
        num_cell(t, p) = max_cell + max_cell/3;

        pct_death(t, p) = num_death(i)/ num_cell(t);
        
    end 
end

%% test show figure
max_cell = num_cell(:, p);

figure
nexttile
plot(1:nframes, num_death(:, p))

nexttile
plot(1:nframes, num_cell(:, p))

nexttile
plot(1:nframes, num_death(:, p) ./ max_cell)
ylim([0 1])

nexttile
plot(1:nframes, smoothdata(num_death(:, p) ./ max_cell, 'movmean', 5))
ylim([0 1])

%% show final figure
% figure
hold on
for p = 4
    if p < 3
        plot(1:nframes, num_death(:, p) ./ num_cell(:, p))
    else
        plot(1:nframes, smoothdata(num_death(:, p) ./ num_cell(:, p), 'movmean', 3))
    end
    
end
hold off
ylim([0 1])