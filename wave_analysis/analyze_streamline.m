%
% data analysis for streamlines
%

clearvars
clc
close all


ROOTDIR = "D:/Spatiotemporal_analysis/";
    MATLAB_FX = "D:/Matlab_FileExchange/";

addpath(ROOTDIR + "/code");
addpath(ROOTDIR + "/code/util");

addpath(MATLAB_FX + "scatter_kde");

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
    fname = fname_prefix + "_streamline_table_intpl.mat";
    stream{ii} = struct2array(load(fname));

    % load wave boundary info
    fname = fname_prefix + "_bo.mat";
    bo{ii} = struct2array(load(fname));
    nframe{ii} = length(bo{ii});

    % % load wv lims
    [wv_xlims(ii,:), wv_ylims(ii,:)] = get_wave_pos_limits(bo{ii}, 'side');
    % if ismember(ii, [3, 44])
    %     wv_xlims(ii, :) = wv_xlims(ii, :) - 250;
    % end

    % load wave numbers
    % all_spatial_wv_num = [all_spatial_wv_num; repelem(ii, size(stream_intpl{ii},1))'];
end

% set min_frame
min_frame = min(cell2mat(nframe));

%% load interpolated stream info
stream_intpl = cell(length(stream), 1);

for ii = 1:height(wvinfo)
    fname_prefix = ROOTDIR + "wave_" + wvinfo.Set(ii) + "/" + ...
        "Well" + sprintf('%02d', wvinfo.well(ii)) + "/wave" + wvinfo.wave(ii) + "_analyses/" + ...
        set_filename(wvinfo.Set(ii), wvinfo.well(ii), wvinfo.start_frame(ii)+"-"+wvinfo.end_frame(ii), "");

    % load streamline info
    fname = fname_prefix + "_streamline_table_intpl.mat";
    stream_intpl{ii} = struct2array(load(fname));
end

%% interpolate the streamlines for the same units
% 
% Interpolate to be ws size
% ws = median of wave speed
% compress or expand the streamlines to 10 windows per frame
%

ws = 25;

stream_intpl = cell(length(stream), 1);

for ii = 1:height(wvinfo)

    stream_intpl{ii} = cell(length(stream{ii}), 1);

    for i = 1:length(stream{ii})

        stream_intpl{ii}{i} = nan(nframe{ii}*ws, 5);

        s = stream{ii}{i};
        last_record_frame = nan;
        for j = 1:nframe{ii}
            if any(s.Frame == j)
                sj = s(s.Frame == j, :);
                last_record_frame = j;

                if height(sj) == 1
                    % if there is only one row, no interpolation can be
                    % done. Data will be repeated for the window size
                    stream_intpl{ii}{i}((j-1)*ws+1:j*ws, :) = [(1:ws)+ws*(j-1); repelem(sj.AngDiff, ws); repelem(sj.Coh, ws); repelem(sj.Density, ws); repelem(sj.SmoothSpeed, ws)]';

                else
                    % interpolate the data to fit the window size
                    x = 1:height(sj);
                    xi = linspace(1, height(sj), ws);

                    interp_dir = interp1(x, sj.AngDiff, xi);
                    interp_coh = interp1(x, sj.Coh, xi);
                    interp_den = interp1(x, sj.Density, xi);
                    interp_spd = interp1(x, sj.SmoothSpeed, xi);

                    stream_intpl{ii}{i}((j-1)*ws+1:j*ws, :) = [(1:ws)+ws*(j-1); fliplr(interp_dir); fliplr(interp_coh); fliplr(interp_den); fliplr(interp_spd)]';

                end
            else
                % if there is no frame information, assume speed is zero
                % also assume spatial prop are the same as the last frame
                sj = s(s.Frame == last_record_frame, :);

                if height(sj) == 1
                    stream_intpl{ii}{i}((j-1)*ws+1:j*ws, :) = [(1:ws)+ws*(j-1); repelem(sj.AngDiff, ws); repelem(sj.Coh, ws); repelem(sj.Density, ws); repelem(sj.SmoothSpeed, ws)]';

                else
                    x = 1:height(sj);
                    xi = linspace(1, height(sj), ws);

                    interp_dir = interp1(x, sj.AngDiff, xi);
                    interp_coh = interp1(x, sj.Coh, xi);
                    interp_den = interp1(x, sj.Density, xi);

                    stream_intpl{ii}{i}((j-1)*ws+1:j*ws, :) = [(1:ws)+ws*(j-1); fliplr(interp_dir); fliplr(interp_coh); fliplr(interp_den); repelem(0, ws)]';

                end
            end
        end

        stream_intpl{ii}{i} = array2table(stream_intpl{ii}{i});
        stream_intpl{ii}{i}.Properties.VariableNames = [{'Index'}, {'AngDiff'}, {'Coh'}, {'Density'}, {'SmoothSpeed'}];

    end
end


%% re-interpolate the streamlines onto plane
disp("Interpolate streamlines onto speed plane")

plotfigure = 1;

spd_smooth_alpha = 5;

% for ii = 1:height(wvinfo)
for ii = 1
    disp(ii)
    wave_bound = false(IMSIZE, IMSIZE);

    for j = 1:numel(bo{ii})
        boundary = bo{ii}{j};
        wave_bound = wave_bound | poly2mask(boundary(:,2), boundary(:,1), IMSIZE, IMSIZE);
    end

    % interpolate the wave and plot
    s = vertcat(stream{ii}{:});

    x = s.XData;
    y = s.YData;
    v = s.SmoothSpeed;

    [xq, yq] = meshgrid(1:1:IMSIZE);

    vq = griddata(x, y, v, xq, yq, 'nearest');

    vg = imgaussfilt(vq, spd_smooth_alpha);

    speed = block_nondead_cells(wave_bound, IMSIZE, IMSIZE, vg);

    speed_fname = ROOTDIR + "wave_" + wvinfo.Set(ii) + "/" + ...
        "Well" + sprintf('%02d', wvinfo.well(ii)) + "/wave" + wvinfo.wave(ii) + "_analyses/" + ...
        set_filename(wvinfo.Set(ii), wvinfo.well(ii), wvinfo.start_frame(ii)+"-"+wvinfo.end_frame(ii), "") + "_smoothspeed_" + IMSIZE + "_2.mat";
    % save(speed_fname, 'speed');

    % plot the interpolated speed
    if plotfigure
        figure(1)
        s = pcolor(xq, yq, speed);
        set(s, 'FaceColor', 'interp', 'EdgeColor', 'none');
        set(gca, 'xtick', [], 'ytick', [])
        axis square tight ij
        colormap parula
    end

end

%%
all_stream_intpl = vertcat(stream_intpl{:});

%%
a = all_stream_intpl.AngDiff;
c = all_stream_intpl.Coh;
d = all_stream_intpl.Density;
spd = all_stream_intpl.SmoothSpeed;

nan_idx = isnan(a) | isnan(c) | isnan(d) | isnan(spd);

a = a(~nan_idx);
c = c(~nan_idx);
d = d(~nan_idx);
spd = spd(~nan_idx);

an = normalize(a);
cn = normalize(c);
dn = normalize(d);

%%
corrcoef([a c d spd])
partialcorr([a c d spd])

%%
mdl = fitlm([an dn an.*dn], normalize(spd))

%%
se = mdl.Coefficients.SE(2:end)

figure
x = 1:length(mdl.Coefficients.Estimate(2:end));
bar(x, abs([mdl.Coefficients.Estimate(2:end)]));
hold on
er = errorbar(x, abs([mdl.Coefficients.Estimate(2:end)]), se, se);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
% xticklabels({"along alignment", "coh", "dens", "a x c", "a x d", "c x d"})
xticklabels({"along alignment", "dens", "a x d"})
ylabel("Absolute Coefficient")
ylim([0 0.5])

%%
ra = an(randperm(length(an)));
rd = dn(randperm(length(dn)));
rad = ra .* rd;

rmdl = fitlm([ra rd rad], spd)

se = rmdl.Coefficients.SE(2:end)

figure
x = 1:length(rmdl.Coefficients.Estimate(2:end));
bar(x, abs([rmdl.Coefficients.Estimate(2:end)]));
hold on
er = errorbar(x, abs([rmdl.Coefficients.Estimate(2:end)]), se, se);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
% xticklabels({"along alignment", "coh", "dens", "a x c", "a x d", "c x d"})
xticklabels({"along alignment", "dens", "a x d"})
ylabel("Absolute Coefficient")

%%
%
% Generate overall heatmap and streamline plots
%
%% calculate total rows, total mean and total standard deivation
total_rows = sum(cellfun(@(x) size(x,1), stream_intpl));

all_stream_intpl = vertcat(stream_intpl{:});
all_stream_intpl = vertcat(all_stream_intpl{:});

aln_mean = mean(all_stream_intpl(:,2), 'omitnan'); % 0.7107
aln_std = std(table2array(all_stream_intpl(:,2)), 'omitnan');
aln_max = max(all_stream_intpl(:,2));
aln_min = 0;

den_mean = mean(all_stream_intpl(:,4), 'omitnan'); % 69.7602
den_std = std(table2array(all_stream_intpl(:,4)), 'omitnan');
den_max = max(all_stream_intpl(:,4));
den_min = min(all_stream_intpl(:,4));

%% Select the window with the largest nematics change
% select frames with the largest change per window size

all_maxdiff_aln_pos = cell(height(wvinfo), 1);
all_maxdiff_aln_neg = cell(height(wvinfo), 1);

all_maxdiff_den_pos = cell(height(wvinfo), 1);
all_maxdiff_den_neg = cell(height(wvinfo), 1);

all_maxdiff_both_pos = cell(height(wvinfo), 1);
all_maxdiff_both_neg = cell(height(wvinfo), 1);

all_maxdiff_both_flat = cell(height(wvinfo), 1);

for ii = 1:height(wvinfo)
% for ii = 22

    maxdiff_aln_pos = nan(height(stream_intpl{ii}), 2);
    maxdiff_aln_neg = nan(height(stream_intpl{ii}), 2);

    maxdiff_den_pos = nan(height(stream_intpl{ii}), 2);
    maxdiff_den_neg = nan(height(stream_intpl{ii}), 2);

    maxdiff_both_pos = nan(height(stream_intpl{ii}), 2);
    maxdiff_both_neg = nan(height(stream_intpl{ii}), 2);

    maxdiff_both_flat = nan(height(stream_intpl{ii}), 2);

    for i = 1:height(stream_intpl{ii})
    % for i = 101
        s = stream_intpl{ii}{i};

        % calculate max diff for along alignment and flattest density
        diffs = table2array(diff(s(:, 2))); % Differences between each consecutive row
        movdiffs2 = filter(ones(1, ws), 1, diffs); % sum the differences
        [maxdiff_aln_pos(i, 2), maxdiff_aln_pos(i, 1)] = max(movdiffs2);
        [maxdiff_aln_neg(i, 2), maxdiff_aln_neg(i, 1)] = min(movdiffs2);

        % diffs = diff(s(:, 4)); % Differences between each consecutive row
        % movdiffs4 = filter(ones(1, ws), 1, diffs); % sum the differences
        % 
        % [ sorted2, sorted_idx2 ] = sort(movdiffs2, 'descend');
        % [~, ranks2] = sort(sorted_idx2);
        % [ sorted4, sorted_idx4 ] = sort(abs(movdiffs4), 'descend');
        % [~, ranks4] = sort(sorted_idx4);
        % 
        % [~, maxdiff_aln_pos(i,1)] = min(ranks2+ranks4);
        % maxdiff_aln_pos(i, 2) = movdiffs2(maxdiff_aln_pos(i,1));
 
        % calculate max diff for density
        diffs = table2array(diff(s(:, 4))); % Differences between each consecutive row
        movdiffs4 = filter(ones(1, ws), 1, diffs); % sum the differences
        [maxdiff_den_pos(i, 2), maxdiff_den_pos(i, 1)] = max(movdiffs4);
        [maxdiff_den_neg(i, 2), maxdiff_den_neg(i, 1)] = min(movdiffs4);
        % diffs = diff(s(:, 2)); % Differences between each consecutive row
        % movdiffs2 = filter(ones(1, ws), 1, diffs); % sum the differences
        % 
        % [ sorted4, sorted_idx4 ] = sort(movdiffs4, 'descend');
        % [~, ranks4] = sort(sorted_idx4);
        % [ sorted2, sorted_idx2 ] = sort(abs(movdiffs2), 'descend');
        % [~, ranks2] = sort(sorted_idx2);
        % 
        % [~, maxdiff_den_pos(i,1)] = min(ranks2+ranks4);
        % maxdiff_den_pos(i, 2) = movdiffs4(maxdiff_den_pos(i,1));

        % calculate max diff for along alignment AND density
        s2 = s(:, 2);
        s2 = (s2 - aln_mean) ./ aln_std;

        s4 = s(:, 4);
        s4 = (s4 - den_mean) ./ den_std;

        diffs1 = table2array(diff(s2)); % Differences between each consecutive row
        diffs2 = table2array(diff(s4)); % Differences between each consecutive row
        movdiffs1 = filter(ones(1, ws), 1, diffs1); % sum the differences
        movdiffs2 = filter(ones(1, ws), 1, diffs2); % sum the differences
        movdiffs = movdiffs1 + movdiffs2;
        [maxdiff_both_pos(i, 2), maxdiff_both_pos(i, 1)] = max(movdiffs);
        [maxdiff_both_neg(i, 2), maxdiff_both_neg(i, 1)] = min(movdiffs);

        % calculate min diff for along alignment AND density
        [maxdiff_both_flat(i, 2), maxdiff_both_flat(i, 1)] = min(abs(movdiffs));
        
    end
    
    % remove any positions that are not negative or positive
    maxdiff_aln_pos(maxdiff_aln_pos(:,2)<=0, :) = nan;
    maxdiff_aln_neg(maxdiff_aln_neg(:,2)>=0, :) = nan;

    maxdiff_den_pos(maxdiff_den_pos(:,2)<=0, :) = nan;
    maxdiff_den_neg(maxdiff_den_neg(:,2)>=0, :) = nan;

    maxdiff_both_pos(maxdiff_both_pos(:,2)<=0, :) = nan;
    maxdiff_both_neg(maxdiff_both_neg(:,2)>=0, :) = nan;

    % store all waves
    all_maxdiff_aln_pos{ii} = maxdiff_aln_pos;
    all_maxdiff_aln_neg{ii} = maxdiff_aln_neg;

    all_maxdiff_den_pos{ii} = maxdiff_den_pos;
    all_maxdiff_den_neg{ii} = maxdiff_den_neg;

    all_maxdiff_both_pos{ii} = maxdiff_both_pos;
    all_maxdiff_both_neg{ii} = maxdiff_both_neg;

    all_maxdiff_both_flat{ii} = maxdiff_both_flat;
end

%% select spatial frames based on the max differences

present_frame = 6*ws;

all_aln_pos_ws = cell(height(wvinfo), 1);
all_aln_neg_ws = cell(height(wvinfo), 1);

all_den_pos_ws = cell(height(wvinfo), 1);
all_den_neg_ws = cell(height(wvinfo), 1);

all_both_pos_ws = cell(height(wvinfo), 1);
all_both_neg_ws = cell(height(wvinfo), 1);

all_both_flat_ws = cell(height(wvinfo), 1);

all_spd_aln_pos_ws = cell(height(wvinfo), 1);
all_spd_aln_neg_ws = cell(height(wvinfo), 1);
all_spd_den_pos_ws = cell(height(wvinfo), 1);
all_spd_den_neg_ws = cell(height(wvinfo), 1);

all_spd_both_pos_ws = cell(height(wvinfo), 1);
all_spd_both_neg_ws = cell(height(wvinfo), 1);
all_spd_both_flat_ws = cell(height(wvinfo), 1);

for ii = 1:height(wvinfo) 
    aln_pos_ws = nan(length(stream_intpl{ii}), present_frame*2);
    aln_neg_ws = nan(length(stream_intpl{ii}), present_frame*2);
    
    den_pos_ws = nan(length(stream_intpl{ii}), present_frame*2);
    den_neg_ws = nan(length(stream_intpl{ii}), present_frame*2);

    both_pos_ws = nan(length(stream_intpl{ii}), present_frame*2);
    both_neg_ws = nan(length(stream_intpl{ii}), present_frame*2);

    both_flat_ws = nan(length(stream_intpl{ii}), present_frame*2);

    spd_aln_pos_ws = nan(length(stream_intpl{ii}), present_frame);
    spd_aln_neg_ws = nan(length(stream_intpl{ii}), present_frame);
    
    spd_den_pos_ws = nan(length(stream_intpl{ii}), present_frame);
    spd_den_neg_ws = nan(length(stream_intpl{ii}), present_frame);

    spd_both_pos_ws = nan(length(stream_intpl{ii}), present_frame);
    spd_both_neg_ws = nan(length(stream_intpl{ii}), present_frame);

    spd_both_flat_ws = nan(length(stream_intpl{ii}), present_frame);

    maxdiff_aln_pos = all_maxdiff_aln_pos{ii};
    maxdiff_aln_neg = all_maxdiff_aln_neg{ii};

    maxdiff_den_pos = all_maxdiff_den_pos{ii};
    maxdiff_den_neg = all_maxdiff_den_neg{ii};

    maxdiff_both_pos = all_maxdiff_both_pos{ii};
    maxdiff_both_neg = all_maxdiff_both_neg{ii};

    maxdiff_both_flat = all_maxdiff_both_flat{ii};

    for i = 1:length(stream_intpl{ii})
    % for i = 101
        s = stream_intpl{ii}{i};

        % construct data structure based on max positive aln diff
        aln_pos_ws_aln = find_maxdiff_idx(s, maxdiff_aln_pos(i,1), present_frame, height(s), 2);
        aln_pos_ws_aln = (aln_pos_ws_aln - aln_mean) ./ aln_std; 
        aln_pos_ws_den = find_maxdiff_idx(s, maxdiff_aln_pos(i,1), present_frame, height(s), 4);
        aln_pos_ws_den = (aln_pos_ws_den - den_mean) ./ den_std; 
        aln_pos_ws(i, :) = [aln_pos_ws_aln aln_pos_ws_den];

        spd_aln_pos_ws(i, :) = find_maxdiff_idx(s, maxdiff_aln_pos(i,1), present_frame, height(s), 5);

        % construct data structure based on max negative aln diff
        aln_neg_ws_aln = find_maxdiff_idx(s, maxdiff_aln_neg(i,1), present_frame, height(s), 2);
        aln_neg_ws_aln = (aln_neg_ws_aln - aln_mean) ./ aln_std; 
        aln_neg_ws_den = find_maxdiff_idx(s, maxdiff_aln_neg(i,1), present_frame, height(s), 4);
        aln_neg_ws_den = (aln_neg_ws_den - den_mean) ./ den_std; 
        aln_neg_ws(i, :) = [aln_neg_ws_aln aln_neg_ws_den];

        spd_aln_neg_ws(i, :) = find_maxdiff_idx(s, maxdiff_aln_neg(i,1), present_frame, height(s), 5);
    
        % construct data structure based on max positive den diff
        den_pos_ws_den = find_maxdiff_idx(s, maxdiff_den_pos(i,1), present_frame, height(s), 4);
        den_pos_ws_den = (den_pos_ws_den - den_mean) ./ den_std; 
        den_pos_ws_aln = find_maxdiff_idx(s, maxdiff_den_pos(i,1), present_frame, height(s), 2);
        den_pos_ws_aln = (den_pos_ws_aln - aln_mean) ./ aln_std; 
        den_pos_ws(i, :) = [den_pos_ws_aln den_pos_ws_den];

        spd_den_pos_ws(i, :) = find_maxdiff_idx(s, maxdiff_den_pos(i,1), present_frame, height(s), 5);

        % construct data structure based on max negative den diff
        den_neg_ws_den = find_maxdiff_idx(s, maxdiff_den_neg(i,1), present_frame, height(s), 4);
        den_neg_ws_den = (den_neg_ws_den - den_mean) ./ den_std; 
        den_neg_ws_aln = find_maxdiff_idx(s, maxdiff_den_neg(i,1), present_frame, height(s), 2);
        den_neg_ws_aln = (den_neg_ws_aln - aln_mean) ./ aln_std; 
        den_neg_ws(i, :) = [den_neg_ws_aln den_neg_ws_den];

        spd_den_neg_ws(i, :) = find_maxdiff_idx(s, maxdiff_den_neg(i,1), present_frame, height(s), 5);

        % construct data structure based on max positive both diff
        both_pos_ws_aln = find_maxdiff_idx(s, maxdiff_both_pos(i,1), present_frame, height(s), 2);
        both_pos_ws_aln = (both_pos_ws_aln - aln_mean) ./ aln_std; 
        both_pos_ws_den = find_maxdiff_idx(s, maxdiff_both_pos(i,1), present_frame, length(s), 4);
        both_pos_ws_den = (both_pos_ws_den - den_mean) ./ den_std; 
        both_pos_ws(i, :) = [both_pos_ws_aln both_pos_ws_den];

        spd_both_pos_ws(i, :) = find_maxdiff_idx(s, maxdiff_both_pos(i,1), present_frame, height(s), 5);

        % construct data structure based on max negative den diff
        both_neg_ws_aln = find_maxdiff_idx(s, maxdiff_both_neg(i,1), present_frame, height(s), 2);
        both_neg_ws_aln = (both_neg_ws_aln - aln_mean) ./ aln_std; 
        both_neg_ws_den = find_maxdiff_idx(s, maxdiff_both_neg(i,1), present_frame, height(s), 4);
        both_neg_ws_den = (both_neg_ws_den - den_mean) ./ den_std; 
        both_neg_ws(i, :) = [both_neg_ws_aln both_neg_ws_den];

        spd_both_neg_ws(i, :) = find_maxdiff_idx(s, maxdiff_both_neg(i,1), present_frame, height(s), 5);

        % construct data structure based on smallest aln and den change
        both_neg_ws_aln = find_maxdiff_idx(s, maxdiff_both_neg(i,1), present_frame, height(s), 2);
        both_neg_ws_aln = (both_neg_ws_aln - aln_mean) ./ aln_std; 
        both_neg_ws_den = find_maxdiff_idx(s, maxdiff_both_neg(i,1), present_frame, height(s), 4);
        both_neg_ws_den = (both_neg_ws_den - den_mean) ./ den_std; 
        both_neg_ws(i, :) = [both_neg_ws_aln both_neg_ws_den];

        spd_both_neg_ws(i, :) = find_maxdiff_idx(s, maxdiff_both_neg(i,1), present_frame, height(s), 5);


        % construct data structure based on smallest to zero for aln and den change
        both_flat_ws_aln = find_maxdiff_idx(s, maxdiff_both_flat(i,1), present_frame, height(s), 2);
        both_flat_ws_aln = (both_flat_ws_aln - aln_mean) ./ aln_std; 
        both_flat_ws_den = find_maxdiff_idx(s, maxdiff_both_flat(i,1), present_frame, height(s), 4);
        both_flat_ws_den = (both_flat_ws_den - den_mean) ./ den_std; 
        both_flat_ws(i, :) = [both_flat_ws_aln both_flat_ws_den];

        spd_both_flat_ws(i, :) = find_maxdiff_idx(s, maxdiff_both_flat(i,1), present_frame, height(s), 5);


        % aln_pos_ws(i,:) = find_maxdiff_idx(s, maxdiff_aln_pos(i,1), present_frame, length(s), 2);
        % aln_neg_ws(i,:) = find_maxdiff_idx(s, maxdiff_aln_neg(i,1), present_frame, length(s), 2);
        % den_pos_ws(i,:) = find_maxdiff_idx(s, maxdiff_den_pos(i,1), present_frame, length(s), 4);
        % den_neg_ws(i,:) = find_maxdiff_idx(s, maxdiff_den_neg(i,1), present_frame, length(s), 4);
    end

    all_aln_pos_ws{ii} = aln_pos_ws;
    all_aln_neg_ws{ii} = aln_neg_ws;
    
    all_den_pos_ws{ii} = den_pos_ws;
    all_den_neg_ws{ii} = den_neg_ws;

    all_both_pos_ws{ii} = both_pos_ws;
    all_both_neg_ws{ii} = both_neg_ws;

    all_both_flat_ws{ii} = both_flat_ws;

    all_spd_aln_pos_ws{ii} = spd_aln_pos_ws;
    all_spd_aln_neg_ws{ii} = spd_aln_neg_ws;

    all_spd_den_pos_ws{ii} = spd_den_pos_ws;
    all_spd_den_neg_ws{ii} = spd_den_neg_ws;

    all_spd_both_pos_ws{ii} = spd_both_pos_ws;
    all_spd_both_neg_ws{ii} = spd_both_neg_ws;

    all_spd_both_flat_ws{ii} = spd_both_flat_ws;
end

all_aln_pos_ws = vertcat(all_aln_pos_ws{:});
all_aln_neg_ws = vertcat(all_aln_neg_ws{:});

all_den_pos_ws = vertcat(all_den_pos_ws{:});
all_den_neg_ws = vertcat(all_den_neg_ws{:});

all_both_pos_ws = vertcat(all_both_pos_ws{:});
all_both_neg_ws = vertcat(all_both_neg_ws{:});

all_both_flat_ws = vertcat(all_both_flat_ws{:});

all_spd_aln_pos_ws = vertcat(all_spd_aln_pos_ws{:});
all_spd_aln_neg_ws = vertcat(all_spd_aln_neg_ws{:});

all_spd_den_pos_ws = vertcat(all_spd_den_pos_ws{:});
all_spd_den_neg_ws = vertcat(all_spd_den_neg_ws{:});

all_spd_both_pos_ws = vertcat(all_spd_both_pos_ws{:});
all_spd_both_neg_ws = vertcat(all_spd_both_neg_ws{:});

all_spd_both_flat_ws = vertcat(all_spd_both_flat_ws{:});


%% number of windows to contain NA
nan_ws = ws*2; % select for streamlines windows that are not with large window

rmnan_aln_pos_ws = all_aln_pos_ws; 
% rmnan_aln_neg_ws = all_aln_neg_ws; 
rmnan_den_pos_ws = all_den_pos_ws;
% rmnan_den_neg_ws = all_den_neg_ws;
rmnan_both_pos_ws = all_both_pos_ws;
% rmnan_both_neg_ws = all_both_neg_ws;
rmnan_both_flat_ws = all_both_flat_ws;

rmnan_spd_aln_pos_ws = all_spd_aln_pos_ws;
% rmnan_spd_aln_neg_ws = all_spd_aln_neg_ws;
rmnan_spd_den_pos_ws = all_spd_den_pos_ws;
% rmnan_spd_den_neg_ws = all_spd_den_neg_ws;

rmnan_spd_both_pos_ws = all_spd_both_pos_ws;
% rmnan_spd_both_neg_ws = all_spd_both_neg_ws;
rmnan_spd_both_flat_ws = all_spd_both_flat_ws; 

rmnan_aln_pos_ws(sum(isnan(all_aln_pos_ws), 2) >= nan_ws, :) = [];
% rmnan_aln_neg_ws(sum(isnan(all_aln_neg_ws), 2) >= nan_ws, :) = [];

rmnan_den_pos_ws(sum(isnan(all_den_pos_ws), 2) >= nan_ws, :) = [];
% rmnan_den_neg_ws(sum(isnan(all_den_neg_ws), 2) >= nan_ws, :) = [];

rmnan_both_pos_ws(sum(isnan(all_both_pos_ws), 2) >= nan_ws, :) = [];
% rmnan_both_neg_ws(sum(isnan(all_both_neg_ws), 2) >= nan_ws, :) = [];

rmnan_both_flat_ws(sum(isnan(all_both_flat_ws), 2) >= nan_ws, :) = [];

rmnan_spd_aln_pos_ws(sum(isnan(all_aln_pos_ws), 2) >= nan_ws, :) = [];
% rmnan_spd_aln_neg_ws(sum(isnan(all_aln_neg_ws), 2) >= nan_ws, :) = [];
rmnan_spd_den_pos_ws(sum(isnan(all_den_pos_ws), 2) >= nan_ws, :) = [];
% rmnan_spd_den_neg_ws(sum(isnan(all_den_neg_ws), 2) >= nan_ws, :) = [];

rmnan_spd_both_pos_ws(sum(isnan(all_both_pos_ws), 2) >= nan_ws, :) = []; 
% rmnan_spd_both_neg_ws(sum(isnan(all_both_neg_ws), 2) >= nan_ws, :) = [];
rmnan_spd_both_flat_ws(sum(isnan(all_both_flat_ws), 2) >= nan_ws, :) = [];

%% cluster for aln pos
% select rows with changes larger than a specific value
threshold = 2;
flat_thres = 1;

% for aln_pos_ws
% rowdiff = max(rmnan_aln_pos_ws(:, 1:present_frame), [], 2) - min(rmnan_aln_pos_ws(:, 1:present_frame), [], 2);
% rowflat = max(rmnan_aln_pos_ws(:, present_frame+1:end), [], 2);
% [M, I] = max(rmnan_aln_pos_ws(:, 1:present_frame), [], 2);
% 
% sp_ws = rmnan_aln_pos_ws((rowdiff>threshold) & (rowflat<flat_thres) & (I>present_frame/2), :);
% spd_sp_ws = rmnan_spd_aln_pos_ws((rowdiff>threshold) & (rowflat<flat_thres) & (I>present_frame/2), :);

% aln1 = mean(rmnan_aln_pos_ws(:, 1:present_frame/2), 2);
% aln2 = mean(rmnan_aln_pos_ws(:, present_frame/2+1:present_frame), 2);
% den1 = mean(rmnan_aln_pos_ws(:, present_frame+1:present_frame*1.5), 2);
% den2 = mean(rmnan_aln_pos_ws(:, present_frame*1.5+1:end), 2);
% 
% sp_ws = rmnan_aln_pos_ws((aln1 < 0) & (aln2 > 0) & (den1 < 0.25) & (den2 < 0.25), :);
% spd_sp_ws = rmnan_spd_aln_pos_ws((aln1 < 0) & (aln2 > 0) & (den1 < 0.25) & (den2 < 0.25), :);
% 
% wv_nums = all_spatial_wv_num;
% wv_nums(sum(isnan(all_aln_pos_ws), 2) >= nan_ws, :) = [];
% wv_nums = wv_nums((aln1 < 0) & (aln2 > 0) & (den1 < 0.25) & (den2 < 0.25), :);

% for den_pos_ws
% rowdiff = max(rmnan_den_pos_ws(:, present_frame+1:end), [], 2) - min(rmnan_den_pos_ws(:, present_frame+1:end), [], 2);
% rowflat = max(rmnan_den_pos_ws(:, 1:present_frame), [], 2);
% [M, I] = max(rmnan_den_pos_ws(:, present_frame+1:end), [], 2);
% 
% sp_ws = rmnan_den_pos_ws((rowdiff>threshold) & (rowflat<flat_thres) & (I>present_frame/2), :);
% spd_sp_ws = rmnan_spd_den_pos_ws((rowdiff>threshold) & (rowflat<flat_thres) & (I>present_frame/2), :);
% 
% wv_nums = all_spatial_wv_num;
% wv_nums(sum(isnan(all_den_pos_ws), 2) >= nan_ws, :) = [];
% wv_nums = wv_nums((rowdiff>threshold) & (rowflat<flat_thres) & (I>present_frame/2), :);

aln1 = mean(rmnan_den_pos_ws(:, 1:present_frame/2), 2);
aln2 = mean(rmnan_den_pos_ws(:, present_frame/2+1:present_frame), 2);
den1 = mean(rmnan_den_pos_ws(:, present_frame+1:present_frame*1.5), 2);
den2 = mean(rmnan_den_pos_ws(:, present_frame*1.5+1:end), 2);
rowdiff = max(rmnan_den_pos_ws(:, present_frame+1:end), [], 2) - min(rmnan_den_pos_ws(:, present_frame+1:end), [], 2);
[aln_M, I] = max(rmnan_den_pos_ws(:, 1:present_frame), [], 2);

sp_ws = rmnan_den_pos_ws((aln1 < 0.25) & (aln2 < 0.25) & (den1 < 0) & (den2 > 0), :);
spd_sp_ws = rmnan_spd_den_pos_ws((aln1 < 0.25) & (aln2 < 0.25) & (den1 < 0) & (den2 > 0), :); % & (max(den2,[],2) > 0.5), :);

% wv_nums = all_spatial_wv_num;
% wv_nums(sum(isnan(all_den_pos_ws), 2) >= nan_ws, :) = [];
% wv_nums = wv_nums((aln1 < 0.25) & (aln2 < 0.25) & (den1 < 0) & (den2 > 0), :); % & (max(den2,[],2) > 0.5), :);


% % for both_pos_ws
% aln_high = mean(rmnan_both_pos_ws(:, present_frame/2:present_frame), 2);
% den_high = mean(rmnan_both_pos_ws(:, present_frame*1.5+1:end), 2);
% aln_low = mean(rmnan_both_pos_ws(:, 1:present_frame/2), 2);
% den_low = mean(rmnan_both_pos_ws(:, present_frame+1:present_frame*1.5), 2);
% 
% sp_ws = rmnan_both_pos_ws((aln_high > 0) & (den_high > 0) & (aln_low < 0) & (den_low < 0), :);
% spd_sp_ws = rmnan_spd_both_pos_ws((aln_high > 0) & (den_high > 0) & (aln_low < 0) & (den_low < 0), :);
% 
% wv_nums = all_spatial_wv_num;
% wv_nums(sum(isnan(all_both_pos_ws), 2) >= nan_ws, :) = [];
% wv_nums = wv_nums((aln_high > 0) & (den_high > 0) & (aln_low < 0) & (den_low < 0), :);

% % both flat : two halves are < 0
% aln_flat1 = mean(rmnan_both_flat_ws(:, 1:present_frame/2), 2);
% den_flat1 = mean(rmnan_both_flat_ws(:, present_frame+1:present_frame*1.5), 2);
% aln_flat2 = mean(rmnan_both_flat_ws(:, present_frame/2+1:present_frame), 2);
% den_flat2 = mean(rmnan_both_flat_ws(:, present_frame*1.5+1:end), 2);
% aln_min = min(rmnan_both_flat_ws(:, 1:present_frame), [], 2);
% den_min = min(rmnan_both_flat_ws(:, present_frame+1:end), [], 2);
% 
% sp_ws = rmnan_both_flat_ws((aln_flat1 < -0) & (den_flat1 < -0) & (aln_flat2 < -0) & (den_flat2 < -0) & (aln_min < -1)  & (den_min < -1), :);
% spd_sp_ws = rmnan_spd_both_flat_ws((aln_flat1 < -0) & (den_flat1 < -0) & (aln_flat2 < -0) & (den_flat2 < -0) & (aln_min < -1)  & (den_min < -1), :);
% 
% wv_nums = all_spatial_wv_num;
% wv_nums(sum(isnan(all_both_flat_ws), 2) >= nan_ws, :) = [];
% wv_nums = wv_nums((aln_flat1 < -0) & (den_flat1 < -0) & (aln_flat2 < -0) & (den_flat2 < -0) & (aln_min < -1)  & (den_min < -1), :);
% 

%% plot clustergram
cg_input = sp_ws;
cg_input(cg_input<-2) = -2;
cg_input(cg_input>2) = 2;

cg = clustergram(cg_input, 'Colormap', redbluecmap, 'Cluster', 'Column', 'ImputeFun', @knnimpute, 'linkage', 'ward');

%% save cluster figure
outpath = "D:\Spatiotemporal_analysis\propagation_all\";
f = cg.plot;
% outname = wvinfo.Set(ii) + "_s" + sprintf('%02d', wvinfo.well(ii)) + "_wave" + wvinfo.wave(ii) + "_3_" + outtype + "_heatmap.jpg";
% outname = "waveAll" + "_" + outtype + "_rmnan_heatmap.jpg";
outname = "haln_lden_allwave_heatmap.jpg";
exportgraphics(gcf, outpath + outname);

%% show speed based on cluster
spd_ordered = spd_sp_ws;
spd_ordered = spd_ordered(str2num(cell2mat(flipud(cg.RowLabels))),:);

figure
% imagesc(ordered_spatial_spd, [0 mean(spatial_spd(:),'omitnan')+2*std(spatial_spd(:),[],'omitnan')])
imagesc(spd_ordered, [0 200])
set(gca,'xtick',[], 'ytick',[],'xticklabel',[],'yticklabel',[])

outname = "haln_lden_allwave_speed_map.jpg";
exportgraphics(gcf, outpath + outname);


%% plot streamline plot for each cluster
savefigure = 1;
colorOrder = ["#0072BD", "#EDB120", '#A5A5A5', "#7E2F8E"];
% outpath = PATH + "/streamline_plot/hclust_streamline_frame6_wave27/" + outtype + "/";
% outpath = PATH + "/hclust_allframe_wave27/" + outtype + "/";
% mkdir(outpath)

plot_range = 6*ws;

for i = 1:length(clustgroup)

    cidx = str2num(cell2mat(cgroup{i}.RowLabels));

    % den_clust = aln_pos_ws(cidx, present_frame+1:end);
    den_clust = sp_ws(cidx, present_frame+1:end);
    
    den_mean = mean(den_clust, 'omitnan');
    den_std1 = den_mean - std(den_clust, 'omitnan');
    den_std2 = den_mean + std(den_clust, 'omitnan');
    den_mean = smoothdata(den_mean(~isnan(den_mean)), 2);
    den_std1 = smoothdata(den_std1(~isnan(den_std1)), 2);
    den_std2 = smoothdata(den_std2(~isnan(den_std2)), 2);

    % aln_clust = aln_pos_ws(cidx, 1:present_frame);
    aln_clust = sp_ws(cidx, 1:present_frame);

    aln_mean = mean(aln_clust, 'omitnan');
    aln_std1 = aln_mean - std(aln_clust, 'omitnan');
    aln_std2 = aln_mean + std(aln_clust, 'omitnan');
    aln_mean = smoothdata(aln_mean(~isnan(aln_mean)), 2);
    aln_std1 = smoothdata(aln_std1(~isnan(aln_std1)), 2);
    aln_std2 = smoothdata(aln_std2(~isnan(aln_std2)), 2);

    % speed_clust = spd_aln_pos_ws(cidx, :);
    speed_clust = spd_sp_ws(cidx, :);

    speed_mean = mean(speed_clust, 'omitnan');
    speed_std1 = speed_mean - std(speed_clust, 'omitnan');
    speed_std2 = speed_mean + std(speed_clust, 'omitnan');
    speed_mean = smoothdata(speed_mean(~isnan(speed_mean)), 2);
    speed_std1 = smoothdata(speed_std1(~isnan(speed_std1)), 2);
    speed_std2 = smoothdata(speed_std2(~isnan(speed_std2)), 2);

    % speed_mean = speed_mean(:, 1:110);
    % speed_std1 = speed_std1(:, 1:110);
    % speed_std2 = speed_std2(:, 1:110);

    % plot_range = max(size(speed_mean));

    figure
    % nexttile
    hold on
    yyaxis left
    plot(1:plot_range, speed_mean, '-', 'Color', colorOrder(1), 'LineWidth', 2);
    patch([1:plot_range fliplr(1:plot_range)], [speed_std1 fliplr(speed_std2)], 1, 'FaceColor', colorOrder(1), 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    xlim([1 plot_range])
    xticks([1:ws:plot_range])
    xticklabels({0 ((1:ws:plot_range)-1)/25+1})
    
    ylim([0 250])
    xlabel('Frame')
    ylabel('Speed')

    yyaxis right
    plot(1:plot_range, den_mean, '--', 'Color', colorOrder(2), 'LineWidth', 2);
    patch([1:plot_range fliplr(1:plot_range)], [den_std1 fliplr(den_std2)], 1, 'FaceColor', colorOrder(2), 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    ylim([-2 2]);
    ylabel('Density')

    plot(1:plot_range, aln_mean, '--', 'Color', colorOrder(3), 'LineWidth', 2);
    patch([1:plot_range fliplr(1:plot_range)], [aln_std1 fliplr(aln_std2)], 1, 'FaceColor', colorOrder(3), 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    %ylim([-0.2 pi/2])

    hold off

    % if savefigure
    %     % outname = outpath +  wvinfo.Set(ii) + "_s" + sprintf('%02d', wvinfo.well(ii)) + "_wave" + wvinfo.wave(ii) + "_streamline_" + outtype + "_clust" + i + ".jpg";
    %     outname = outpath +  "waveAll_streamline_" + outtype + "_rmnan_clust" + i + ".jpg";
    %     exportgraphics(gcf, outname)
    % end

end
% close all

%% plot streamline plot from selected streamlines based on densit/aln changes
colorOrder = ["#0072BD", "#EDB120", '#A5A5A5', "#7E2F8E"];
plot_range = 6*ws;

den_clust = sp_ws(:, present_frame+1:end);

den_mean = mean(den_clust, 'omitnan');
den_std1 = den_mean - std(den_clust, 'omitnan');
den_std2 = den_mean + std(den_clust, 'omitnan');
den_mean = smoothdata(den_mean(~isnan(den_mean)), 2);
den_std1 = smoothdata(den_std1(~isnan(den_std1)), 2);
den_std2 = smoothdata(den_std2(~isnan(den_std2)), 2);

% aln_clust = aln_pos_ws(cidx, 1:present_frame);
aln_clust = sp_ws(:, 1:present_frame);

aln_mean = mean(aln_clust, 'omitnan');
aln_std1 = aln_mean - std(aln_clust, 'omitnan');
aln_std2 = aln_mean + std(aln_clust, 'omitnan');
aln_mean = smoothdata(aln_mean(~isnan(aln_mean)), 2);
aln_std1 = smoothdata(aln_std1(~isnan(aln_std1)), 2);
aln_std2 = smoothdata(aln_std2(~isnan(aln_std2)), 2);

speed_clust = spd_sp_ws(:, :);

speed_mean = mean(speed_clust, 'omitnan');
speed_std1 = speed_mean - std(speed_clust, 'omitnan');
speed_std2 = speed_mean + std(speed_clust, 'omitnan');
% speed_mean = smoothdata(speed_mean(~isnan(speed_mean)), 2, 'movmean', 50); 
% speed_std1 = smoothdata(speed_std1(~isnan(speed_std1)), 2, 'movmean', 50);
% speed_std2 = smoothdata(speed_std2(~isnan(speed_std2)), 2, 'movmean', 50);
speed_mean = smoothdata(speed_mean(~isnan(speed_mean)), 2);
speed_std1 = smoothdata(speed_std1(~isnan(speed_std1)), 2);
speed_std2 = smoothdata(speed_std2(~isnan(speed_std2)), 2);

figure
% nexttile
hold on
yyaxis left
plot(1:plot_range, speed_mean, '-', 'Color', colorOrder(1), 'LineWidth', 2);
patch([1:plot_range fliplr(1:plot_range)], [speed_std1 fliplr(speed_std2)], 1, 'FaceColor', colorOrder(1), 'FaceAlpha', 0.2, 'EdgeColor', 'none')
xlim([1 plot_range])
xticks([1:ws:plot_range])
xticklabels({0 ((1:ws:plot_range)-1)/25+1})

ylim([0 250])
xlabel('Frame')
ylabel('Speed')

yyaxis right
plot(1:plot_range, den_mean, '--', 'Color', colorOrder(2), 'LineWidth', 2);
patch([1:plot_range fliplr(1:plot_range)], [den_std1 fliplr(den_std2)], 1, 'FaceColor', colorOrder(2), 'FaceAlpha', 0.2, 'EdgeColor', 'none')
ylim([-1.5 2.5]);
ylabel('Density')

plot(1:plot_range, aln_mean, '--', 'Color', colorOrder(3), 'LineWidth', 2);
patch([1:plot_range fliplr(1:plot_range)], [aln_std1 fliplr(aln_std2)], 1, 'FaceColor', colorOrder(3), 'FaceAlpha', 0.2, 'EdgeColor', 'none')
%ylim([-0.2 pi/2])

hold off

outname = "laln_hden_allwave_streamline_plot.jpg";
exportgraphics(gcf, outpath + outname);
outname = "laln_hden_allwave_streamline_plot.svg";
saveas(gcf, outpath + outname, 'svg');


%% FUNCTION
function [varargout] = block_nondead_cells(wave_bound, vector_per_side, IMSIZE, varargin)
    varargout = cell((nargin)-3, 1);
    for i = 1:(nargin-3)
        lg = imresize(varargin{i}, [IMSIZE IMSIZE], 'bilinear');
        lg(~wave_bound) = NaN;
        varargout{i} = imresize(lg, [vector_per_side vector_per_side], "bilinear");
    end
end

function maxdiff_s = find_maxdiff_idx(s, maxdiff, present_frame, len_s, col)
    s_idx = maxdiff-present_frame/2;
    e_idx = maxdiff+present_frame/2-1;

    if isnan(maxdiff)
        maxdiff_s = nan(1, present_frame);
    else
        s_nan = max(0, 1-s_idx);
        s_idx = max(1, s_idx);
    
        e_nan = max(0, e_idx-len_s);
        e_idx = min(len_s, e_idx);

        maxdiff_s = [nan(1, s_nan) s(s_idx:e_idx, col)' nan(1, e_nan)];
    end
end

% function maxdiff_s = find_maxdiff_idx2(s, maxdiff, present_frame, len_s, col1, col2)
%     s_idx = maxdiff-present_frame/2;
%     e_idx = maxdiff+present_frame/2-1;
% 
%     if isnan(maxdiff)
%         maxdiff_s = nan(1, present_frame);
%     else
%         s_nan = max(0, 1-s_idx);
%         s_idx = max(1, s_idx);
% 
%         e_nan = max(0, e_idx-len_s);
%         e_idx = min(len_s, e_idx);
% 
%         maxdiff_s = [nan(1, s_nan) s(s_idx:e_idx, col)' nan(1, e_nan)];
%     end
% end
