%
%   calculate wave and alignment concordance
%
clearvars
clc
close all


MATLAB_FX = "D:/Matlab_FileExchange/";

addpath D:/Spatiotemporal_analysis/code/util/

warning('off' , 'MATLAB:MKDIR:DirectoryExists');


%% load angdiff for pdms file

outpath = "D:\Spatiotemporal_analysis\wave_0_analyses\wave_pdms\";

stream_file1 = "D:\Spatiotemporal_analysis\wave_pdms-0404\Well10\wave1_analyses\"  + ...
    "pdms-0404_s10t01-50_ORG_streamline_table_speed.mat";
stream1 = struct2array(load(stream_file1));

stream_file2 = "D:\Spatiotemporal_analysis\wave_pdms-0404\Well12\wave1_analyses\"  + ...
    "pdms-0404_s12t01-50_ORG_streamline_table_speed.mat";
stream2 = struct2array(load(stream_file2));

stream_file3 = "D:\Spatiotemporal_analysis\wave_pdms-0418\Well20\wave1_analyses\"  + ...
    "pdms-0418_s20t01-50_ORG_streamline_table_speed.mat";
stream3 = struct2array(load(stream_file3));

stream_file4 = "D:\Spatiotemporal_analysis\wave_pdms-0418\Well23\wave1_analyses\"  + ...
    "pdms-0418_s23t01-50_ORG_streamline_table_speed.mat";
stream4 = struct2array(load(stream_file4));


streams = {stream1; stream2; stream3; stream4};

%% multiple streamlines

nframe = 30;
nrep = numel(streams);

angdiff_all = nan(nframe, nrep);
circStd_all = nan(nframe, nrep);
circVar_all = nan(nframe, nrep);

for r = 1:nrep

    stream = streams{r};          % one repeat
    s = vertcat(stream{:});       % concatenate streamline tables

    for i = 1:nframe
        ad = s(s.Frame == i, :).AngDiff180;
        ad(isnan(ad)) = [];

        if isempty(ad)
            continue
        end

        z = mean(exp(1i * 2 * ad));

        % mean angle difference
        angdiff_all(i, r) = 0.5 * angle(z);

        % circular resultant length
        R = abs(z);

        % circular standard deviation
        circStd_all(i, r) = 0.5 * sqrt(-2 * log(R));

        % circular variance
        circVar_all(i, r) = 1 - R;
    end
end


% plot
mean_angdiff = mean(angdiff_all, 2, 'omitnan');
mean_circStd = std( angdiff_all, [], 2, 'omitnan');
mean_circVar = var( angdiff_all, [], 2, 'omitnan');

sm_factor = 3;

mean_angdiff = smoothdata(mean_angdiff, 'movmean', sm_factor);
mean_circStd = smoothdata(mean_circStd, 'movmean', sm_factor);
mean_circVar = smoothdata(mean_circVar, 'movmean', sm_factor);

upper = mean_angdiff + mean_circStd;
lower = mean_angdiff - mean_circStd;

x = 1:nframe;

figure
tiledlayout(2, 2)

nexttile
hold on

fill([x fliplr(x)], ...
     [upper' fliplr(lower')], ...
     [0.8 0.8 0.8], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4);

plot(x, mean_angdiff, 'k-', 'LineWidth', 2)

yline(0, '--')
xline([8-2 10-2], "--")
ylim([-pi/2 pi/2])
ylabel('Radians')
title("Mean AngDiff180 ± mean circular SD")

nexttile
plot(x, mean_angdiff, 'k-', 'LineWidth', 2)
yline(0, "--");
xline([8-2 10-2], "--")
ylim([-pi/2 pi/2])
ylabel('Radians')
title("Mean AngDiff180")

nexttile
plot(x, mean_circStd, 'k-', 'LineWidth', 2)
ylim([0 pi/2])
xline([8-2 10-2], "--")
ylabel('Radians')
title("Mean circular SD")

nexttile
plot(x, mean_circVar, 'k-', 'LineWidth', 2)
ylim([0 pi/2])
xline([8-2 10-2], "--")
xlabel('Frame')
ylabel('sd of the mean')
title("Mean circular variance")


%% flatten all streamlines

s = vertcat(streams{:});
s = vertcat(s{:});

for i = 1:nframe
    ad = s(s.Frame == i, :).AngDiff180;
    ad(isnan(ad)) = [];

    if isempty(ad)
        continue
    end

    z = mean(exp(1i * 2 * ad));

    % mean angle difference
    angdiff_all(i, r) = 0.5 * angle(z);

    % circular resultant length
    R = abs(z);

    % circular standard deviation
    circStd_all(i, r) = 0.5 * sqrt(-2 * log(R));

    % circular variance
    circVar_all(i, r) = 1 - R;
end

mean_angdiff = mean(angdiff_all, 2, 'omitnan');
mean_circStd = mean(circStd_all, 2, 'omitnan');
mean_circVar = mean(circVar_all, 2, 'omitnan');

sm_factor = 4;

mean_angdiff = smoothdata(mean_angdiff, 'movmean', sm_factor);
mean_circStd = smoothdata(mean_circStd, 'movmean', sm_factor);
mean_circVar = smoothdata(mean_circVar, 'movmean', sm_factor);

% Shaded region = SD of repeat-level mean AngDiff180
upper = mean_angdiff + mean_circStd;
lower = mean_angdiff - mean_circStd;

x = 1:nframe;

figure
tiledlayout(2, 2)

nexttile
hold on

fill([x fliplr(x)], ...
     [upper' fliplr(lower')], ...
     [0.8 0.8 0.8], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4);

plot(x, mean_angdiff, 'k-', 'LineWidth', 2)

yline(0, '--')
xline([8-2 10-2], "--")
ylim([-pi/2 pi/2])
ylabel('Radians')
title("Mean AngDiff180 ± mean circular SD")

nexttile
plot(x, mean_angdiff, 'k-', 'LineWidth', 2)
yline(0, "--");
xline([8-2 10-2], "--")
ylim([-pi/2 pi/2])
ylabel('Radians')
title("Mean AngDiff180")

nexttile
plot(x, mean_circStd, 'k-', 'LineWidth', 2)
ylim([0 pi/2])
xline([8-2 10-2], "--")
ylabel('Radians')
title("Mean circular SD")

nexttile
plot(x, mean_circVar, 'k-', 'LineWidth', 2)
ylim([0 1])
xline([8-2 10-2], "--")
xlabel('Frame')
ylabel('sd of the mean')
title("Mean circular variance")


%% multiple stream table
nframe = 30;

nstream = numel(streams);

circVar_all = nan(nframe, nstream);

for k = 1:nstream

    s = streams{k};
    s = vertcat(s{:});

    for i = 1:nframe
        ad = s(s.Frame == i, :).AngDiff180;
        ad(isnan(ad)) = [];
        
        R = abs(mean(exp(1i * 2 * ad)));
        circVar_all(i, k) = 1 - R;
    end
end

% plot the circular variance
circVar_mean = mean(circVar_all, 2, 'omitnan');
circVar_sd   = std(circVar_all, 0, 2, 'omitnan');
outfile = "pdms_cirvar_all";

x = (1:nframe)+2;

figure
hold on

fill([x fliplr(x)], ...
     [(circVar_mean - circVar_sd)' fliplr((circVar_mean + circVar_sd)')], ...
     [0.8 0.8 0.8], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.5);

plot(x, circVar_mean, 'k-', 'LineWidth', 2)
xlabel('Frame')
ylabel('Circular variance')
ylim([0 1])

xline([8 10], "--")

% exportgraphics(gcf, outpath + outfile + ".jpg");
% exportgraphics(gcf, outpath + outfile + ".pdf");
% savefig(outpath + outfile + ".fig");

%% polarhistograms
s = vertcat(streams{:});
s = vertcat(s{:});

ad_inside = s(s.Frame <= 8, :).AngDiff180;
ad_inside(isnan(ad_inside)) = [];

ad_outside = s(s.Frame > 10, :).AngDiff180;
ad_outside(isnan(ad_outside)) = [];


figure
nexttile
polarhistogram(ad_inside, deg2rad(-90:15:90), 'Normalization', 'probability');
rlim([0 0.4])
rticklabels([0 0.25 0.5 0.75 1])
title("Inside PDMS")

ax = gca;
ax.ThetaLim = [-90 90];
ax.ThetaTick = -90:30:90;

nexttile
polarhistogram(ad_outside/2, deg2rad(-90:15:90), 'Normalization', 'probability');
rlim([0 0.4])
rticklabels([0 0.25 0.5 0.75 1])
title("Outside PDMS")

ax = gca;
ax.ThetaLim = [-90 90];
ax.ThetaTick = -90:30:90;


outfile = "pdms_wave_cell_aln_angle_well12";

% exportgraphics(gcf, outpath + outfile + ".jpg");
% exportgraphics(gcf, outpath + outfile + ".pdf");
% savefig(outpath + outfile + ".fig");


%% plot wv-cell alignment over time
streams = {stream1; stream2; stream3; stream4};

nframe = 30;
nrep = numel(streams);

cut = [480 600 650 650];

%% Calculate metrics for short and long streamlines

shortMetrics = calc_stream_metrics_by_length(streams, cut, nframe, "short");
longMetrics  = calc_stream_metrics_by_length(streams, cut, nframe, "long");

% Smooth results
sm_factor = 5;
short_meanAng = smoothdata(shortMetrics.meanAng,     'movmean', sm_factor);
sm_factor = 10;
short_circStd = smoothdata(shortMetrics.meanCircStd, 'lowess', sm_factor);
short_circVar = smoothdata(shortMetrics.meanCircVar, 'movmean', sm_factor);

sm_factor = 10;
long_meanAng = smoothdata(longMetrics.meanAng,     'movmean', sm_factor);
sm_factor = 19;
long_circStd = smoothdata(longMetrics.meanCircStd, 'lowess', sm_factor);
long_circVar = smoothdata(longMetrics.meanCircVar, 'movmean', sm_factor);

sm_factor = 5;
short_meanSpeed = smoothdata(shortMetrics.meanSpeed, 'movmean', sm_factor);
long_meanSpeed  = smoothdata(longMetrics.meanSpeed,  'movmean', sm_factor);

% Shaded region = mean angle ± mean circular SD

short_upper = short_meanAng + short_circStd;
short_lower = short_meanAng - short_circStd;

long_upper = long_meanAng + long_circStd;
long_lower = long_meanAng - long_circStd;

x = 1:nframe;

%% Plot for aln angle

savefile = true;

figure

tiledlayout(2, 2)

nexttile
hold on

fill([x fliplr(x)], ...
     [short_upper' fliplr(short_lower')], ...
     [0.75 0.75 0.75], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.35);

fill([x fliplr(x)], ...
     [long_upper' fliplr(long_lower')], ...
     [0 0 0.65], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.25);

plot(x, short_meanAng, 'k-', 'LineWidth', 2)
plot(x, long_meanAng, 'b-', 'LineWidth', 2)


xline([8-2 10-2], "--")
ylim([0 pi/2])

xticks([0:6:30])
xticklabels([0:5:25])

yticks([0:pi/4:pi/2])
yticklabels({'-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2'})
yticklabels({'0', '\pi/4', '\pi/2'})
xlabel('Frame')
ylabel('Radians')
title("Mean AngDiff180 ± mean circular SD")


nexttile
hold on

plot_colored_line(x, short_meanAng, short_meanSpeed, 3)
plot_colored_line(x, long_meanAng,  long_meanSpeed,  3)
allSpeed = [short_meanSpeed(:); long_meanSpeed(:)];
clim([min(allSpeed, [], 'omitnan'), 200])

myparula = generate_black_parula(64);
colormap(myparula)

xticks([0:6:30])
xticklabels([0:5:25])
xline([8-2 10-2], "--")
ylim([0 pi/2])
xlabel('Frame')
ylabel('Radians')
title("Mean AngDiff180")


nexttile
hold on
plot(x, short_circStd, 'k-', 'LineWidth', 2)
plot(x, long_circStd, 'b-', 'LineWidth', 2)
ylim([0 pi/2])
xline([8-2 10-2], "--")
xticks([0:6:30])
xticklabels([0:5:25])
xlabel('Frame')
ylabel('Radians')
title("Mean circular SD")


nexttile
hold on
plot(x, short_circVar, 'k-', 'LineWidth', 2)
plot(x, long_circVar, 'b-', 'LineWidth', 2)
ylim([0 1])
xline([8-2 10-2], "--")
xticks([0:6:30])
xticklabels([0:5:25])
xlabel('Frame')
ylabel('Circular variance')
title("Mean circular variance")


if savefile
    outfile = "pdms_wave_cell_aln_angle_time";
    
    exportgraphics(gcf, outpath + outfile + ".jpg");
    exportgraphics(gcf, outpath + outfile + ".pdf");
    savefig(outpath + outfile + ".fig");
end


%% Plot for density

% Calculate metrics for short and long streamlines

shortMetrics = calc_stream_metrics_by_length(streams, cut, nframe, "short");
longMetrics  = calc_stream_metrics_by_length(streams, cut, nframe, "long");

sm_factor = 6;
short_meanDen = smoothdata(shortMetrics.meanDen,     'movmean', sm_factor);
sm_factor = 5;
long_meanDen = smoothdata(longMetrics.meanDen,     'movmean', sm_factor);

sm_factor = 10;
short_stdDen  = smoothdata(shortMetrics.stdDen,      'movmean', sm_factor);
long_stdDen  = smoothdata(longMetrics.stdDen,      'movmean', sm_factor);

savefile = true;

short_upper = short_meanDen + short_stdDen;
short_lower = short_meanDen - short_stdDen;

long_upper = long_meanDen + long_stdDen;
long_lower = long_meanDen - long_stdDen;

figure

tiledlayout(2, 2)

nexttile
hold on

fill([x fliplr(x)], ...
     [short_upper' fliplr(short_lower')], ...
     [0.75 0.75 0.75], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.35);

fill([x fliplr(x)], ...
     [long_upper' fliplr(long_lower')], ...
     [0 0 0.65], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.25);

plot(x, short_meanDen, 'k-', 'LineWidth', 2)
plot(x, long_meanDen, 'b-', 'LineWidth', 2)

xticks([0:6:30])
xticklabels([0:5:25])
yline(0, '--')
xline([8-2 10-2], "--")
ylim([50 100])
ylabel('Density')
title("Mean Density ± mean SD")


nexttile
hold on

plot_colored_line(x, short_meanDen, short_meanSpeed, 3)
plot_colored_line(x, long_meanDen,  long_meanSpeed,  3)

myparula = generate_black_parula(64);
colormap(myparula)

allSpeed = [short_meanSpeed(:); long_meanSpeed(:)];
clim([min(allSpeed, [], 'omitnan'), 200]) %max(allSpeed, [], 'omitnan')])

xticks([0:6:30])
xticklabels([0:5:25])
yline(0, "--")
xline([8 10], "--")
ylim([50 100])
ylabel('Density')
title("Mean Density")

nexttile
hold on
plot(x, short_stdDen, 'k-', 'LineWidth', 2)
plot(x, long_stdDen,  'b-', 'LineWidth', 2)
ylim([0 30])
xline([8 10], "--")
xticks([0:6:30])
xticklabels([0:5:25])
ylabel('Density')
title("Mean density SD")


if savefile
    outfile = "pdms_density_time";
    
    exportgraphics(gcf, outpath + outfile + ".jpg");
    exportgraphics(gcf, outpath + outfile + ".pdf");
    savefig(outpath + outfile + ".fig");
end


%% FUNCTION
function plot_colored_line(x, y, c, lw)

    x = x(:)';
    y = y(:)';
    c = c(:)';

    surface([x; x], ...
            [y; y], ...
            zeros(2, numel(x)), ...
            [c; c], ...
            'FaceColor', 'none', ...
            'EdgeColor', 'interp', ...
            'LineWidth', lw);
end

function metrics = calc_stream_metrics_by_length(streams, cut, nframe, groupType)
    nrep = numel(streams);
    
    angdiff_all = nan(nframe, nrep);
    circStd_all = nan(nframe, nrep);
    circVar_all = nan(nframe, nrep);

    density_all = nan(nframe, nrep);
    denStd_all  = nan(nframe, nrep);

    speed_all   = nan(nframe, nrep);
    
    for r = 1:nrep
    
        stream = streams{r};
    
        % table height of each streamline
        h = cellfun(@height, stream);
    
        switch groupType
            case "short"
                useIdx = h < cut(r);
            case "long"
                useIdx = h >= cut(r);
            otherwise
                error('groupType must be "short" or "long"')
        end
    
        stream_sub = stream(useIdx);
    
        if isempty(stream_sub)
            continue
        end
    
        s = vertcat(stream_sub{:});
    
        for i = 1:nframe
    
            ad = s(s.Frame == i, :).AngDiff180;
            ad(isnan(ad)) = [];

            if isempty(ad)
                continue
            end
    
            z = mean(exp(1i * 2 * ad));
    
            % mean angle difference
            angdiff_all(i, r) = 0.5 * angle(z);
    
            % circular resultant length
            R = abs(z);
    
            % prevent log(0)
            R = max(R, eps);
    
            % circular standard deviation
            circStd_all(i, r) = 0.5 * sqrt(-2 * log(R));
    
            % circular variance
            circVar_all(i, r) = 1 - R;

            %%% calculate density
            den = s(s.Frame == i, :).Density;
            den(isnan(den)) = [];

            density_all(i, r) = mean(den, 'omitnan');
            denStd_all(i,  r) = std( den, 'omitnan');

            %%% calculate SmoothSpeed
            spd = s(s.Frame == i, :).SmoothSpeed;
            spd(isnan(spd)) = [];
            
            speed_all(i, r) = mean(spd, 'omitnan');
        end
    end
    
    metrics.angdiff_all = angdiff_all;
    metrics.circStd_all = circStd_all;
    metrics.circVar_all = circVar_all;
    metrics.density_all = density_all;
    metrics.speed_all = speed_all;

    % metrics.meanAng = mean(angdiff_all, 2, 'omitnan');
    % metrics.meanCircStd = mean(circStd_all, 2, 'omitnan');
    % metrics.meanCircVar = mean(circVar_all, 2, 'omitnan');

    metrics.meanAng     = mean(abs(angdiff_all), 2,     'omitnan');
    metrics.meanCircStd = std( angdiff_all, [], 2, 'omitnan');
    metrics.meanCircVar = var( angdiff_all, [], 2, 'omitnan');
    metrics.meanDen     = mean(density_all, 2,     'omitnan');
    metrics.stdDen      = std(density_all,  [], 2, 'omitnan');
    metrics.meanSpeed   = mean(speed_all, 2, 'omitnan');
    metrics.stdSpeed    = std(speed_all, [], 2, 'omitnan');

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
    % black = [50/256 50/256 50/256];
    black = [0/256 0/256 0/256];
    ramp = [linspace(black(1), parula_start(1), n_black_to_blue)', ...
            linspace(black(2), parula_start(2), n_black_to_blue)', ...
            linspace(black(3), parula_start(3), n_black_to_blue)'];
    
    % Concatenate the colormaps
    myparula = [ramp; parula_map];
end
