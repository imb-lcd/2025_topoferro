%
% Characterize the general properties of waves
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath D:\Spatiotemporal_analysis\code
addpath D:\Spatiotemporal_analysis\code\util

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATION

PATH = ROOTDIR + "/wave_0_analyses/";

den_reg = 120;
MCHCH = "c1";
DENPATH = MCHCH + "_mCherry_density/";

CY5CH = "c2";
CY5ADJPATH = CY5CH + "_cy5_adj/";

DICCH = "c3";
NEMPATH = DICCH + "_DIC_nematics/";
ORIENTPATH = NEMPATH + "orientation_files/";

IMSIZE = 5000;

%% Load wave information table

fname = PATH + "all_wave_info.txt";

wvinfo = readtable(fname);

% exclude waves that will not be analyzed
wvinfo(wvinfo.include == 0, :) = [];

%% Load stream info
stream = cell(height(wvinfo),1);

bo = cell(height(wvinfo), 1);
nframe = cell(height(wvinfo), 1);

wv_xlims = nan(height(wvinfo),2);
wv_ylims = nan(height(wvinfo),2);

all_spatial_wv_num = [];

for ii = 1:height(wvinfo)
    fname_prefix = ROOTDIR + "wave_" + wvinfo.Set(ii) + "/" + ...
        "Well" + sprintf('%02d', wvinfo.well(ii)) + "/wave" + wvinfo.wave(ii) + "_analyses/" + ...
        set_filename(wvinfo.Set(ii), wvinfo.well(ii), wvinfo.start_frame(ii)+"-"+wvinfo.end_frame(ii), "");

    % load streamline info
    fname = fname_prefix + "_streamline_table.mat";
    stream{ii} = struct2array(load(fname));

    % load wave boundary info
    fname = fname_prefix + "_bo.mat";
    bo{ii} = struct2array(load(fname));
    nframe{ii} = length(bo{ii});

    % load wv lims
    [wv_xlims(ii,:), wv_ylims(ii,:)] = get_wave_pos_limits(bo{ii}, 'side');
    if ismember(ii, [3, 44])
        wv_xlims(ii, :) = wv_xlims(ii, :) - 250;
    end
end

% set min_frame
min_frame = min(cell2mat(nframe));

%% Boundary for the example wave 6

bo_file = "D:\Spatiotemporal_analysis\propagation_examples\laln_hden_mCh300FBS10-1_s05t28-41c2_figures\mCh300FBS10-1_s05t28-41_ORG_bo.mat";
full_Bo = struct2array(load(bo_file));

[wv_xlims, wv_ylims] = get_wave_pos_limits(bo{6}, 'side');

imsize = [5000 5000];

figure
imshow(zeros(imsize)+255);

hold on
for i = 1:3:length(full_Bo)
    Bo = full_Bo{i};
    plot(Bo(:,2), Bo(:,1), 'Color', 'k', 'LineWidth', 2);
end
xlim(wv_xlims)
ylim(wv_ylims)
hold off

%% calculate area of contour

iso_area = nan(length(full_Bo), 1);

figure
imshow(zeros(5000, 5000)+255)
hold on

centerx = (max(full_Bo{1}(:,1)) - min(full_Bo{1}(:,1)))/2 + min(full_Bo{1}(:,1));
centery = (max(full_Bo{1}(:,2)) - min(full_Bo{1}(:,2)))/2 + min(full_Bo{1}(:,2));

for i = 1:length(full_Bo)
    iso_area(i) = polyarea(full_Bo{i}(:,1), full_Bo{i}(:,2));

    % plot circle
    r = sqrt(iso_area(i)/pi);
    th = 0:pi/50:2*pi;
    xunit = r*cos(th) + centery;
    yunit = r*sin(th) + centerx;
    plot(xunit, yunit, '-.k');

    % plot boundary
    plot(full_Bo{i}(:,2), full_Bo{i}(:,1), 'Color', 'k', 'LineWidth', 2);
end
hold off

xlim(wv_xlims);
ylim(wv_ylims-200);

%% load density, entropy and angles
max_den = 174.8893;
min_den = 1.0038;

all_den = cell(height(wvinfo), 1);
% all_ent = cell(height(wvinfo), 1);

for ii = 1:height(wvinfo)
% for ii = 6
    % load den
    den_path = ROOTDIR + "wave_" + wvinfo.Set(ii) + "/Well" + sprintf('%02d', wvinfo.well(ii)) + "/" + DENPATH + "/density_reg" + den_reg + "/";
    den_name = den_path + set_filename(wvinfo.Set(ii), wvinfo.well(ii), wvinfo.spatial_frame(ii), MCHCH) + "_den" + den_reg + "_" + IMSIZE + ".mat";

    all_den{ii} = struct2array(load(den_name));

    % load ent
    % ent_path = ROOTDIR + "wave_" + wvinfo.Set(ii) + "/Well" + sprintf('%02d', wvinfo.well(ii)) + "/" + ORIENTPATH + "/";
    % ent_name = ent_path + set_filename(wvinfo.Set(ii), wvinfo.well(ii), wvinfo.spatial_frame(ii), DICCH) + "_ent.mat";
    % 
    % all_ent{ii} = struct2array(load(ent_name));
end

%% compare wae entropy with isotropic entropy

all_wv_ent = cell(height(wvinfo), 1);
all_iso_ent = cell(height(wvinfo), 1);

for ii = 1:height(wvinfo)
    single_Bo = bo{ii};

    % figure
    % wave entropy area
    wv_mask = poly2mask(single_Bo{end}(:,2), single_Bo{end}(:,1), 4949, 4949);
    wv_ent = double(all_ent{ii}) .* wv_mask;
    wv_ent(wv_ent == 0) = nan;

    all_wv_ent{ii} = wv_ent;
    % nexttile
    % imshow(wv_ent,  [min(ent(:)) max(ent(:))], colormap=jet)
    
    % isotropic entropy
    iso_area = polyarea(single_Bo{end}(:,1), single_Bo{end}(:,2));
    centerx = (max(single_Bo{1}(:,1)) - min(single_Bo{1}(:,1)))/2 + min(single_Bo{1}(:,1));
    centery = (max(single_Bo{1}(:,2)) - min(single_Bo{1}(:,2)))/2 + min(single_Bo{1}(:,2));

    r = sqrt(iso_area/pi);
    th = 0:pi/50:2*pi;
    xunit = r*cos(th) + centery;
    yunit = r*sin(th) + centerx;
    
    iso_mask = poly2mask(xunit, yunit, 4949, 4949);
    
    iso_ent = double(all_ent{ii}) .* iso_mask;
    iso_ent(iso_ent == 0) = nan;
    % nexttile
    % imshow(iso_ent, [min(ent(:)) max(ent(:))], colormap=jet)

    all_iso_ent{ii} = iso_ent;
end

samp_all_ent = datasample(all_ent{6}(:), length(wv_ent(:)));
% 
% % statistics
% [~, p] = ttest2(wv_ent(:), iso_ent(:));
% 
figure
boxplot([wv_ent(:), iso_ent(:), samp_all_ent(:)])


%% stream
for ii = 1:height(wvinfo)
    if ii == 2 || ii == 46 || ii == 47
        continue
    end
    ent = all_ent{ii};
    ent = imresize(ent, [5000 5000]);
    for i = 1:length(stream{ii})
    % for i = 1
        s = stream{ii}{i};

        xpos = fix(stream{ii}{i}.XData);
        ypos = fix(stream{ii}{i}.YData);
        
        stream{ii}{i}.Entropy = ent(sub2ind(size(ent), ypos, xpos));
    end
end

%%
stream{2} = [];
stream{46} = [];
stream{47} = [];
%%
all_stream = vertcat(stream{:});
all_stream = vertcat(all_stream{:});
all_stream(isnan(all_stream.Frame),:) = [];

start_stream = all_stream(all_stream.Frame == 1, :);
prop_stream = all_stream((all_stream.Frame ~= 1) & (all_stream.SmoothSpeed >= 45), :);
stop_stream = all_stream(all_stream.SmoothSpeed < 45, :);

%%
sample_size = 100000;

samp_all_ent = datasample(all_stream.Entropy, sample_size);
samp_start_ent = datasample(start_stream.Entropy, sample_size);
samp_prop_ent = datasample(prop_stream.Entropy, sample_size);
samp_stop_ent = datasample(stop_stream.Entropy, sample_size);

%%
% data = [all_stream.Entropy' rmini_stream.Entropy' samp_ent'];
% grp = [zeros(1, length(all_stream.Entropy)) ones(1, length(rmini_stream.Entropy)) (ones(1, length(samp_ent))*2) ];

% data = [all_stream.Entropy' start_stream.Entropy' prop_stream.Entropy' stop_stream.Entropy' samp_ent'];
data = [all_stream.Entropy' start_stream.Entropy' prop_stream.Entropy' stop_stream.Entropy'];
grp = [zeros(1, length(all_stream.Entropy)) ones(1, length(start_stream.Entropy)) ones(1, length(prop_stream.Entropy))*2 ones(1, length(stop_stream.Entropy))*3];

% data = [samp_all_ent' samp_start_ent' samp_prop_ent' samp_stop_ent'];
% grp = [zeros(1, sample_size) ones(1, sample_size) 2*ones(1, sample_size) 3*ones(1, sample_size)];

figure
swarmchart(grp, data)

%% direction along alignment vs each point from initiation
all_stream = vertcat(stream{:}); % vertcat(stream{6});
% all_stream = vertcat(stream{6});
% all_stream = stream;

for i = 1:length(all_stream)
    s = all_stream{i};
    all_stream{i}.Iso_angdiff = atan(abs(s.YData-s.YData(end))./abs(s.XData-s.XData(end)));
end

all_stream = vertcat(all_stream{:});

figure
bin_sz = 7; %7; % 10;
% rlimit =  2.8*10^4; %12*10^5; %2.8*10^4; % 12*10^5;

nexttile
polarhistogram(all_stream.AngDiff, bin_sz, 'Normalization', 'pdf')
rlim([0 2])

nexttile
polarhistogram(all_stream.Iso_angdiff, bin_sz, 'Normalization', 'pdf')
rlim([0 2])

% continuation - calculate polarhistogram for 360
bin_sz = 12;

figure(1)
clf
nexttile
angles = mod(all_stream.AngDiff360 + pi/2, pi) - pi/2;
polarhistogram(angles, bin_sz, 'Normalization', 'pdf')
rlim([0 1])
rticks([0 0.25 0.5 0.75 1])

nexttile
angles = mod(all_stream.Iso_angdiff + pi/2, pi) - pi/2;
polarhistogram(angles, bin_sz, 'Normalization', 'pdf')
rlim([0 1])
rticks([0 0.25 0.5 0.75 1])
%% compare wave density with isotropic density

sample_size = 10000;

all_wv_den = cell(height(wvinfo), 1);
all_iso_den = cell(height(wvinfo), 1);

for ii = 1:height(wvinfo)
% for ii = 1
    single_Bo = bo{ii};

    % figure
    % wave density area
    wv_mask = poly2mask(single_Bo{end}(:,2), single_Bo{end}(:,1), 5000, 5000);
    wv_den = all_den{ii} .* wv_mask;
    wv_den(wv_den == 0) = nan;
    wv_den = wv_den(:);
    wv_den(isnan(wv_den)) = [];

    all_wv_den{ii} = datasample(wv_den, sample_size);

    % nexttile
    % imshow(wv_den,  [min(den(:)) max(den(:))], colormap=jet)
 
    % isotropic density
    iso_area = polyarea(single_Bo{end}(:,1), single_Bo{end}(:,2));
    centerx = (max(single_Bo{1}(:,1)) - min(single_Bo{1}(:,1)))/2 + min(single_Bo{1}(:,1));
    centery = (max(single_Bo{1}(:,2)) - min(single_Bo{1}(:,2)))/2 + min(single_Bo{1}(:,2));

    r = sqrt(iso_area/pi);
    th = 0:pi/50:2*pi;
    xunit = r*cos(th) + centery;
    yunit = r*sin(th) + centerx;
    
    iso_mask = poly2mask(xunit, yunit, 5000, 5000);
    
    iso_den = all_den{ii} .* iso_mask;
    iso_den(iso_den == 0) = nan;
    iso_den = iso_den(:);
    iso_den(isnan(iso_den)) = [];

    all_iso_den{ii} = datasample(iso_den, sample_size);

    % nexttile
    % imshow(iso_den, [min(den(:)) max(den(:))], colormap=jet)
end
%%
flat_wv_den = vertcat(all_wv_den{:});
flat_iso_den = vertcat(all_iso_den{:});

% statistics
[~, p] = ttest2(flat_wv_den(:), flat_iso_den(:));

% figure
% boxplot([flat_wv_den(:), flat_iso_den(:)])

figure
hold on
% violinplot(repelem(1, length(flat_wv_den(:)))', flat_wv_den(:))
boxchart(repelem(1, length(flat_wv_den(:)))', flat_wv_den(:), 'BoxWidth', 0.3, 'MarkerStyle', 'none');
% violinplot(repelem(2, length(flat_iso_den(:)))', flat_iso_den(:))
boxchart(repelem(2, length(flat_iso_den(:)))', flat_iso_den(:), 'BoxWidth', 0.3, 'MarkerStyle', 'none');
xticks([1 2])
xlim([0.5 2.5])
yticks([0 28.27433388, 56.54866776, 84.82300165, 113.0973355, 141.3716694])
yticklabels([0 1 2 3 4 5])
ylim([0 141.3716694])
hold off


%%
all_wv_den = cell(height(wvinfo), 1);
all_iso_den = cell(height(wvinfo), 1);
% figure
for ii = 1:height(wvinfo)
% for ii = 6
    curr_bo = bo{ii};

    wv_den = cell(length(bo), 1);
    iso_den = cell(length(bo), 1);

    prev_mask = false(5000, 5000);
    prev_iso_mask = false(5000, 5000);

    den = all_den{ii};
    for b = 1:1:length(curr_bo)
    % for b = 1:3:length(curr_bo)
        wv_mask = poly2mask(curr_bo{b}(:,2), curr_bo{b}(:,1), 5000, 5000);
        wv_mask = wv_mask .* ~prev_mask;

        curr_wv_den = wv_mask .* den;
        wv_den{b} = process_den(curr_wv_den);

        % nexttile
        % imshow(wv_mask)
        prev_mask = prev_mask + wv_mask;

        % isotropic 
        iso_area = polyarea(curr_bo{b}(:,1), curr_bo{b}(:,2));
        centerx = (max(curr_bo{1}(:,1)) - min(curr_bo{1}(:,1)))/2 + min(curr_bo{1}(:,1));
        centery = (max(curr_bo{1}(:,2)) - min(curr_bo{1}(:,2)))/2 + min(curr_bo{1}(:,2));
    
        r = sqrt(iso_area/pi);
        th = 0:pi/50:2*pi;
        xunit = r*cos(th) + centery;
        yunit = r*sin(th) + centerx;
        
        iso_mask = poly2mask(xunit, yunit, 5000, 5000);
        iso_mask = iso_mask .* ~prev_iso_mask;

        % nexttile
        % imshow(iso_mask)

        curr_iso_den = iso_mask .* den;
        iso_den{b} = process_den(curr_iso_den);

        prev_iso_mask = prev_iso_mask + iso_mask; 
    end

    all_wv_den{ii} = wv_den;
    all_iso_den{ii} = iso_den;
end

%%
nsample = 5000;

figure
hold on
for b = 1:1:length(bo)
    wv_sample = [];
    iso_sample = [];

    for ii = 1:height(wvinfo)
        curr_wv_den  = all_wv_den{ii}{b};
        curr_iso_den = all_iso_den{ii}{b};

        if length(curr_wv_den) < 1  || length(curr_iso_den) < 1
            continue;
        end

        den_sample = datasample(curr_wv_den, nsample);
    
        wv_sample = [wv_sample; den_sample];

        den_sample = datasample(curr_iso_den, nsample);

        iso_sample = [iso_sample; den_sample];
    end

    violinplot(repelem(b, length(wv_sample)), wv_sample)
    boxchart(repelem(b, length(wv_sample)), wv_sample, 'MarkerStyle', 'none');

    violinplot(repelem(b+1, length(iso_sample)), iso_sample)
    boxchart(repelem(b+1, length(iso_sample)), iso_sample, 'MarkerStyle', 'none');
end
hold off
    
%%
wv_den_median = cell(10, 1);
iso_den_median = cell(10, 1);
figure
hold on
for b = 1:1:10
    curr_wv_den_median = nan(height(wvinfo), 1);
    curr_iso_den_median = nan(height(wvinfo), 1);
    for ii = 1:height(wvinfo)
        curr_wv_den  = all_wv_den{ii}{b};
        curr_iso_den = all_iso_den{ii}{b};


        curr_wv_den_median(ii) = median(curr_wv_den);
        curr_iso_den_median(ii) = median(curr_iso_den);
    end

    violinplot(repelem(b, height(wvinfo)), curr_wv_den_median)
    boxchart(repelem(b, height(wvinfo)), curr_wv_den_median, 'MarkerStyle', 'none');

    violinplot(repelem(b+1, height(wvinfo)), curr_iso_den_median)
    boxchart(repelem(b+1, height(wvinfo)), curr_iso_den_median, 'MarkerStyle', 'none');
    b
end
hold off

%% scatter plot
figure
hold on
for b = 1:1:10
    curr_wv_den_median = nan(height(wvinfo), 1);
    curr_iso_den_median = nan(height(wvinfo), 1);
    for ii = 1:height(wvinfo)
        curr_wv_den  = all_wv_den{ii}{b};
        curr_iso_den = all_iso_den{ii}{b};


        curr_wv_den_median(ii) = median(curr_wv_den);
        curr_iso_den_median(ii) = median(curr_iso_den);
    end
    swarmchart(repelem(b, height(wvinfo)), curr_wv_den_median,  'r', 'XJitter', 'rand')
    swarmchart(repelem(b, height(wvinfo)), curr_iso_den_median, 'b', 'XJitter', 'rand')
end
hold off

%% FUNCITON
function den = process_den(den)
   den(den == 0) = nan;
   den = den(:);
   den(isnan(den)) = [];
end