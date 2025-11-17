%
% analyze intrpl wave
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
%ROOTDIR = "/home/N417/Jen-Hao/Spatiotemporal/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath(ROOTDIR + "/code");
addpath(ROOTDIR + "/code/util");

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


%% load wave information table
perturb = "AA2a"; % "AA2a" "AY";

wvinfo_fname = "perturb_wave_info_all.txt";

wvinfo = readtable(PATH + wvinfo_fname);


% exclude waves that will not be analyzed
wvinfo = wvinfo(contains(wvinfo.perturb, perturb)| contains(wvinfo.perturb, "AA2b"), :);
wvinfo(wvinfo.include == 0, :) = [];

%% load the smooth interpolated wave
disp("load the interpolated data")
wv_data = cell(height(wvinfo), 1);
all_data = [];
for ii = 1:height(wvinfo)
    prefix = wvinfo.Set(ii);
    w = wvinfo.well(ii);
    wv = wvinfo.wave(ii);
    fr = wvinfo.spatial_frame(ii);

    frame_range = sprintf("%02d", wvinfo.start_frame(ii)) + "-" + ...
            sprintf("%02d", wvinfo.end_frame(ii));

    path = ROOTDIR + "wave_" + prefix + "/" + "Well" + sprintf('%02d', w) + "/";

    % ws = IMSIZE/reduced_sz;

    wv_name = path + "/wave" + wv + "_analyses/" + ...
            set_filename(prefix, w, frame_range, '') + "_wvdata_intpl100.mat"; % _old.mat";

    wv_data{ii} = struct2array(load(wv_name));

    all_data = [all_data; wv_data{ii}];
end

wave = rmmissing(all_data);

x = wave.angdiff;
y = wave.density;
v = wave.speed;

idx = isnan(x) | isnan(y) | isnan(v);
x = x(~idx);
y = y(~idx);
v = v(~idx);

v_max = max(v);
v_min = min(v);

v = (v - v_min) ./ (v_max - v_min);

myparula = generate_black_parula(); 
%% plot correlation

figure
nexttile
scatter_kde(x, v); %, 20, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('along alignment');
ylabel('velocity')
ylim([30 250])
axis square

nexttile
scatter_kde(y, v); %, 20, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('density');
ylabel('velocity');
ylim([30 250])
axis square

%% plot scatter plot
figure
scatter(x, y, 20, v, 'filled', 'MarkerFaceAlpha', 0.5)

xlim([0 0.85]);
ylim([30 130]);
xlabel('along alignment');
ylabel('density')
axis square

%% calculate interpolated plot of directional discordance, density, and colored by speed
% myparula = generate_black_parula_red_old();
figure(2)
clf

for nbin = 23
    for upr = 103
        lwr = 30; % 50; %30;

        sigma = 3;

        [xq, yq] = meshgrid(linspace(0, 1.6, nbin), linspace(lwr, upr, nbin));
        [~, ~, vq] = griddata(x, y, v, xq, yq, 'nearest');

        vg = imgaussfilt(vq, sigma);

        %% plot the interpolated wave
        % figure
        nexttile
        pcolor(xq, yq, vg, 'FaceColor', 'interp', 'EdgeColor', 'none');
        % imagesc(xq(1,:), yq(:,1), flipud(vg));
        xticks(0:pi/8:pi/2)
        xticklabels({'0', '\pi/8','\pi/4','3\pi/8','\pi/2'})

        axis square tight
        colormap(myparula);
        colorbar

        title(nbin + "+" + upr)
    end
end
% title(nbin)
% end

%% plot for each wave
nbin = 25;
lwr = 30; %30;
upr = 90; %90;
sigma = 2;

[xq, yq] = meshgrid(linspace(0, 1.6, nbin), linspace(lwr, upr, nbin));

figure
tiledlayout(3, height(wvinfo))
for ii = 1:height(wvinfo)
    wv = wv_data{ii};

    x = wv.angdiff;
    y = wv.density-30;
    v = wv.speed;
    

    idx = isnan(x) | isnan(y) | isnan(v);
    x = x(~idx);
    y = y(~idx);
    v = v(~idx);
    % v(v>200) = 200;
v = (v - v_min) ./ (v_max - v_min);

    nexttile(ii)
    scatter(x, y, 20, v, 'filled', 'MarkerFaceAlpha', 0.5)
    % xlim([0 0.85]);
    ylim([lwr upr]); %ylim([30 130]);
    axis square

    [~, ~, vq] = griddata(x, y, v, xq, yq, 'natural');
    
    nexttile(height(wvinfo)+ii)
    pcolor(xq, yq, vq, 'FaceColor', 'interp', 'EdgeColor', 'none');
    xticks(0:pi/8:pi/2)
    xticklabels({'0', '\pi/8','\pi/4','3\pi/8','\pi/2'})
    axis square
    % clim([v_min v_max])
    clim([0.2 0.5])

    [~, ~, vq] = griddata(x, y, v, xq, yq, 'nearest');
    vg = imgaussfilt(vq, sigma);

    nexttile(height(wvinfo)*2+ii)
    pcolor(xq, yq, vg, 'FaceColor', 'interp', 'EdgeColor', 'none');
    xticks(0:pi/8:pi/2)
    xticklabels({'0', '\pi/8','\pi/4','3\pi/8','\pi/2'})
    axis square off
    % clim([v_min v_max])
    clim([0.2 0.5])
    % colorbar
    title( wvinfo.index(ii))

end

%% calculate coefficients
x = wave.angdiff;
y = wave.density;
v = wave.speed*2;

idx = isnan(x) | isnan(y) | isnan(v);
x = x(~idx);
y = y(~idx);
v = v(~idx);

a = normalize(x);
d = normalize(y);
ad = a .* d;

mdl = fitlm([a d ad], normalize(v))

se = mdl.Coefficients.SE(2:end);
%%
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
ylim([0 0.6])

%% plot all coeffs together

coeff_file = "D:\Spatiotemporal_analysis\wave_0_analyses\wave_relative_contribution.txt";
coeffs = readtable(coeff_file, 'Delimiter', '\t');

d = table2array(coeffs(:, 2:5));
se = table2array(coeffs(:, 6:end));

figure
b = bar(d, 'grouped')
hold on
x = nan(3, 4);

for i = 1:4
    x(:,i) = b(i).XEndPoints;

    er = errorbar(x(:, i), d(:, i), se(:, i), 'k', 'LineStyle', 'none');
end
hold off

xticklabels({"along alignment", "dens", "a x d"})
ylabel("Absolute Coefficient")
ylim([0 0.6])


%% FUNCTION
function myparula = generate_black_parula(steps_black)
    if ~exist('steps_black')
        steps_black = 25;
    end

    % Number of steps for each segment
    n_black_to_blue = steps_black;
    n_parula = 412;
    
    % Get the parula colormap
    parula_map = parula(n_parula);
    
    % Identify the dark blue of parula (first color)
    parula_start = parula_map(1, :);  % usually a dark blue
    
    % Interpolate from black to parula's blue
    black = [50/256 50/256 50/256];
    ramp = [linspace(black(1), parula_start(1), n_black_to_blue)', ...
            linspace(black(2), parula_start(2), n_black_to_blue)', ...
            linspace(black(3), parula_start(3), n_black_to_blue)'];
    
    % Concatenate the colormaps
    myparula = [ramp; parula_map];
end

function myredparula = generate_black_parula_red()
    % Parameters
    data_min = 50;
    data_split = 135;
    data_max = 200;
    n_colors = 256;
    
    % Define warm sub-ranges (in data units)
    r_flat_yellow  = 20;  % 150–160
    r_y2o = 10;           % 160–170
    r_o2r = 15;           % 170–185
    r_r2m = 15;           % 185–200
    r_warm_total = r_flat_yellow + r_y2o + r_o2r + r_r2m;
    
    % Compute total proportions
    r_cold = data_split - data_min;     % 50–140 → 90
    r_total = r_cold + r_warm_total;
    ratio_cold = r_cold / r_total;
    ratio_warm = r_warm_total / r_total;
    
    n1 = round(n_colors * ratio_cold); % cold: black → parula
    n2 = n_colors - n1;                % warm: yellow → magenta
    
    % === Section 1: Black to Parula ===
    steps_black = 25;
    myparula_full = generate_black_parula(steps_black);
    cmap1 = interp1(linspace(0,1,size(myparula_full,1)), myparula_full, linspace(0,1,n1));

    % Split into 4 transitions (yellow→orange, orange→red, red→magenta, magenta→white)
n_sections = 4;
segment_lengths = round(linspace(0, n2, n_sections + 1)); % 5 edges

% Define color anchors
c1 = [1, 1, 0];  % Yellow
c2 = [1, 0.5, 0]; % Orange
c3 = [1, 0, 0];   % Red
c4 = [1, 0, 1];   % Magenta
c5 = [1, 1, 1];   % White

% Interpolate each segment
seg1 = [linspace(c1(1), c2(1), segment_lengths(2)-segment_lengths(1))', ...
        linspace(c1(2), c2(2), segment_lengths(2)-segment_lengths(1))', ...
        linspace(c1(3), c2(3), segment_lengths(2)-segment_lengths(1))'];

seg2 = [linspace(c2(1), c3(1), segment_lengths(3)-segment_lengths(2))', ...
        linspace(c2(2), c3(2), segment_lengths(3)-segment_lengths(2))', ...
        linspace(c2(3), c3(3), segment_lengths(3)-segment_lengths(2))'];

seg3 = [linspace(c3(1), c4(1), segment_lengths(4)-segment_lengths(3))', ...
        linspace(c3(2), c4(2), segment_lengths(4)-segment_lengths(3))', ...
        linspace(c3(3), c4(3), segment_lengths(4)-segment_lengths(3))'];

seg4 = [linspace(c4(1), c5(1), segment_lengths(5)-segment_lengths(4))', ...
        linspace(c4(2), c5(2), segment_lengths(5)-segment_lengths(4))', ...
        linspace(c4(3), c5(3), segment_lengths(5)-segment_lengths(4))'];

% Combine the full extended segment
cmap2 = [seg1; seg2; seg3; seg4];
 myredparula = [cmap1; cmap2];
end


function myredparula = generate_black_parula_red_new()
    % Parameters
    data_min = 50;
    data_split = 160;
    data_max = 200;
    n_colors = 256;
    
    % Define warm sub-ranges (in data units)
    r_flat_yellow  = 20;  % 150–160
    r_y2o = 10;           % 160–170
    r_o2r = 15;           % 170–185
    r_r2m = 15;           % 185–200
    r_warm_total = r_flat_yellow + r_y2o + r_o2r + r_r2m;
    
    % Compute total proportions
    r_cold = data_split - data_min;     % 50–140 → 90
    r_total = r_cold + r_warm_total;
    ratio_cold = r_cold / r_total;
    ratio_warm = r_warm_total / r_total;
    
    n1 = round(n_colors * ratio_cold); % cold: black → parula
    n2 = n_colors - n1;                % warm: yellow → magenta
    
    % === Section 1: Black to Parula ===
    steps_black = 25;
    myparula_full = generate_black_parula(steps_black);
    cmap1 = interp1(linspace(0,1,size(myparula_full,1)), myparula_full, linspace(0,1,n1));
    
    % === Section 2: Warm Colors ===
    % Allocate warm steps proportionally
    n_flat_yellow = round(n2 * r_flat_yellow / r_warm_total);
    n_y2o = round(n2 * r_y2o / r_warm_total);
    n_o2r = round(n2 * r_o2r / r_warm_total);
    n_r2m = n2 - n_flat_yellow - n_y2o - n_o2r;
    
    % Define RGB points
    yellow  = [1, 1, 0];
    orange  = [1, 0.5, 0];
    red     = [1, 0, 0];
    magenta = [1, 0, 1];
    
    % Create segments
    flat_yellow = repmat(yellow, n_flat_yellow, 1);
    
    yellow_to_orange = [linspace(yellow(1), orange(1), n_y2o)', ...
                        linspace(yellow(2), orange(2), n_y2o)', ...
                        linspace(yellow(3), orange(3), n_y2o)'];
    
    orange_to_red = [linspace(orange(1), red(1), n_o2r)', ...
                     linspace(orange(2), red(2), n_o2r)', ...
                     linspace(orange(3), red(3), n_o2r)'];
    
    red_to_magenta = [linspace(red(1), magenta(1), n_r2m)', ...
                      linspace(red(2), magenta(2), n_r2m)', ...
                      linspace(red(3), magenta(3), n_r2m)'];
    
    % Combine warm section
    cmap2 = [flat_yellow; yellow_to_orange; orange_to_red; red_to_magenta];
        % Combine both
    myredparula = [cmap1; cmap2];
end

function myredparula = generate_black_parula_red_old()
    % Set full data range
    data_min = 50;
    data_split = 135;
    data_max = 200;
    
    % Total number of colors (choose based on your resolution)
    n_colors = 256;
    
    % Proportional lengths based on data ranges
    range1 = data_split - data_min; % = 90
    range2 = data_max - data_split; % = 60
    ratio1 = range1 / (range1 + range2); % ~0.6
    ratio2 = range2 / (range1 + range2); % ~0.4
    
    n1 = round(n_colors * ratio1); % number of colors for parula section
    n2 = n_colors - n1;            % number of colors for yellow→red section
    
    % Colormap section 1: parula
    % cmap1 = parula(n1);
    % cmap1 = generate_black_parula();
    % steps_black = 25; % or however many you want
    myparula_full = generate_black_parula(); % your function
    % Resample to n1 total colors
    cmap1 = interp1(linspace(0,1,size(myparula_full,1)), myparula_full, linspace(0,1,n1));

    
    % Colormap section 2: yellow → orange → red
    yellow = [1, 1, 0];
    orange = [1, 0.5, 0];
    red    = [1, 0, 0];
    magenta = [1, 0, 1];
    
    % n_half = floor(n2 / 2);
    % yellow_to_orange = [linspace(yellow(1), orange(1), n_half)', ...
    %                     linspace(yellow(2), orange(2), n_half)', ...
    %                     linspace(yellow(3), orange(3), n_half)'];
    % 
    % orange_to_red = [linspace(orange(1), red(1), n2 - n_half)', ...
    %                  linspace(orange(2), red(2), n2 - n_half)', ...
    %                  linspace(orange(3), red(3), n2 - n_half)'];
    % 
    % cmap2 = [yellow_to_orange; orange_to_red];
    % Divide n2 into three segments
n_y2o = round(n2 / 3);
n_o2r = round(n2 / 3);
n_r2m = n2 - n_y2o - n_o2r;

yellow_to_orange = [linspace(yellow(1), orange(1), n_y2o)', ...
                    linspace(yellow(2), orange(2), n_y2o)', ...
                    linspace(yellow(3), orange(3), n_y2o)'];

orange_to_red = [linspace(orange(1), red(1), n_o2r)', ...
                 linspace(orange(2), red(2), n_o2r)', ...
                 linspace(orange(3), red(3), n_o2r)'];

red_to_magenta = [linspace(red(1), magenta(1), n_r2m)', ...
                  linspace(red(2), magenta(2), n_r2m)', ...
                  linspace(red(3), magenta(3), n_r2m)'];

cmap2 = [yellow_to_orange; orange_to_red; red_to_magenta];
    
    % Combine both
    myredparula = [cmap1; cmap2];
end


%%
% %% load the streamline and boundary information
% rng(4);
% 
% stream = cell(height(wvinfo),1);
% bo = cell(height(wvinfo), 1);
% nframe = cell(height(wvinfo), 1);
% % samp_n = 100;
% for ii = 1:height(wvinfo)
% % for ii = [2 3 4 5]
% 
%     start_frame = sprintf('%02d', wvinfo.start_frame(ii));
%     end_frame   = sprintf('%02d', wvinfo.end_frame(ii));
%     % load streamline info
%     fname = ROOTDIR + "wave_" + wvinfo.Set(ii) + "/" + ...
%         "Well" + sprintf('%02d', wvinfo.well(ii)) + "/wave" + wvinfo.wave(ii) + "_analyses/" + ...
%         set_filename(wvinfo.Set(ii), wvinfo.well(ii), start_frame+"-"+end_frame, "") + "_streamline_table.mat";    
% 
%     stream{ii} = struct2array(load(fname));
% 
%     % curr_s = struct2array(load(fname));
%     % % if (height(curr_s) > samp_n)
%     % if ii == 1
%     %     p = randperm(height(curr_s), samp_n);
%     % 
%     %     stream{ii} = curr_s(p);
%     % else
%     %     stream{ii} = curr_s;
%     % end
% 
% 
%     % load wave boundary info
%     fname = ROOTDIR + "wave_" + wvinfo.Set(ii) + "/" + ...
%         "Well" + sprintf('%02d', wvinfo.well(ii)) + "/wave" + wvinfo.wave(ii) + "_analyses/" + ...
%         set_filename(wvinfo.Set(ii), wvinfo.well(ii), start_frame+"-"+end_frame, "") + "_bo.mat";
% 
%     bo{ii} = struct2array(load(fname));
% 
%     nframe{ii} = length(bo{ii});
% end
% 
% 
% %% plot scatter plot of directional discordance, density, and colored by speed
% nbin = 40; 
% 
% wave = vertcat(stream{:});
% wave = vertcat(wave{:});
% 
% x = wave.AngDiff;
% y = wave.Density;
% v = wave.SmoothSpeed;
% 
% idx = isnan(x) | isnan(y) | isnan(v);
% x = x(~idx);
% y = y(~idx);
% v = v(~idx);
% 
% v(v>300) = 200;
% 
% % v = normalize(v);
% % %%
% figure
% scatter(x, y, 100, v, 'filled', 'MarkerFaceAlpha', 0.01)
% ylim([30 130])
% axis square
% 
% % %% interpolate wave speed on density and along cell orientation
% 
% 
% 
% % y(y>130) = 130;
% y(y>90) = 90;
% y(y<30) = 30;
% 
% 
% 
% [xq, yq] = meshgrid(linspace(0, 1.6, nbin), linspace(30, 130, nbin));
% 
% [Xq,Yq,vq] = griddata(x, y, v, xq, yq, 'nearest');
% 
% vg = imgaussfilt(vq, 5);
% %
% figure
% s = pcolor(xq, yq, vg);
% set(s, 'FaceColor', 'interp', 'EdgeColor', 'none');
% % xticks(0:pi/8:pi/2)
% % xticklabels({'0', '\pi/
% % 8','\pi/4','3\pi/8','\pi/2'})
% % ylim([20 160])
% xlabel('along alignment');
% ylabel('density');
% axis square tight
% myparula = generate_black_parula(25);
% colormap(myparula);
% colorbar
% % clim([50 150])
% % title(rr)
% 
% 
% 
% %% interpolate for each wave
% nbin = 40;
% 
% figure
% tiledlayout(1, 6)
% for i = 1:height(wvinfo)
% 
%     wave = vertcat(stream{i});
%     wave = vertcat(wave{:});
% 
%     x = wave.AngDiff;
%     y = wave.Density;
%     v = wave.SmoothSpeed;
% 
%     idx = isnan(x) | isnan(y) | isnan(v);
%     x = x(~idx);
%     y = y(~idx);
%     v = v(~idx);
% 
% nexttile
% scatter(x, y, 50, v, 'filled', 'MarkerFaceAlpha', 0.05)
% ylim([30 130])
% axis square
% continue
% 
% 
%     v(v>300) = 300;
%     y(y>130) = 130;
%     y(y<30) = 30;
% 
% 
% 
%     [xq, yq] = meshgrid(linspace(0, 1.6, nbin), linspace(30, 130, nbin));
% 
%     [Xq,Yq,vq] = griddata(x, y, v, xq, yq, 'nearest');
% 
%     vg = imgaussfilt(vq, 5);
%     %
%     nexttile
%     s = pcolor(xq, yq, vg);
%     set(s, 'FaceColor', 'interp', 'EdgeColor', 'none');
%     % xticks(0:pi/8:pi/2)
%     % xticklabels({'0', '\pi/8','\pi/4','3\pi/8','\pi/2'})
%     % ylim([20 160])
%     xlabel('along alignment');
%     ylabel('density');
%     axis square tight
%     myparula = generate_black_parula(25);
%     colormap(myparula);
%     colorbar
%     % clim([50 150])
% 
% end
% 
% %%
% idx = find(y > 120);
% tempy = y(idx);
% tempv = v(idx);
% 
% rng(1)
% tempy1 = tempy + (rand(1,length(tempy)).*20-10)';
% tempv1 = tempv + (rand(1,length(tempv)).*20-10)';
% tempy2 = tempy + (rand(1,length(tempy)).*20-10)';
% tempv2 = tempv + (rand(1,length(tempv)).*20-10)';
% tempy3 = tempy + (rand(1,length(tempy)).*20-10)';
% tempv3 = tempv + (rand(1,length(tempv)).*20-10)';
% tempy4 = tempy + (rand(1,length(tempy)).*20-10)';
% tempv4 = tempv + (rand(1,length(tempv)).*20-10)';
% 
% newy = [y;tempy1;tempy2;tempy3;tempy4;];
% newv = [v;tempv1;tempv2;tempv3;tempv4;];
% 
% figure
% nexttile
% scatter_kde(x, v); %, 20, 'filled', 'MarkerFaceAlpha', 0.5);
% xlabel('along alignment');
% ylabel('velocity')
% ylim([30 240])
% axis square
% 
% nexttile
% scatter_kde(newy, newv); %, 20, 'filled', 'MarkerFaceAlpha', 0.5);
% xlabel('density');
% ylabel('velocity');
% ylim([30 240])
% axis square
