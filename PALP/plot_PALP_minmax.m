%
%    Analyze PALP signal differences
%
clearvars
clc
close all

MATLAB_FX = "D:/Matlab_FileExchange/";

addpath D:/Spatiotemporal_analysis/code/util/

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%%
align_type = "Align"; % "HD_align";
unaln_type = "Unalign"; %["pos_def" "neg_def"];

atype = [align_type unaln_type];

prefix = ["PALP-control-0202" "PALP-control-0708" "PALP-control-0709" "PALP-control2-0801" ...
          "PALP-AA-0801" "PALP-AY-1019"];

path = "E:/PALP/";

bwr = mybwr();

%%
align_data = cell(length(prefix), 1);
unaln_data = cell(length(prefix), 1);

all_minmax = nan(length(prefix), 4);

for p = 1:length(prefix)
    disp(prefix(p))
    fname = "E:/PALP/" + prefix(p) + "_summary.xlsx"; 

    curr_minmax = [Inf 0 Inf 0]; % o min o max; r min max

    sheets = sheetnames(fname);

    align_data{p} = table;
    for t = align_type
        disp(t)
        if ismember(t, sheets)
            d = readtable(fname, 'Sheet', t, 'VariableNamingRule', 'preserve');
            d.Polar = d.or1;
            d.Periphery = d.or2;
            d.note = string(d.note);
            d(d.note ~= "NaN" & d.note ~= "" & ~ismissing(d.note), :) = [];
            
    
            align_data{p} = [align_data{p}; d];
    
            curr_minmax(1) = min([curr_minmax(1) d.o1' d.o2']); % for o min
            curr_minmax(2) = max([curr_minmax(2) d.o1' d.o2']); % for o max
            curr_minmax(3) = min([curr_minmax(3) d.r1' d.r2']); % for r min
            curr_minmax(4) = max([curr_minmax(4) d.r1' d.r2']); % for r max
        end
    end

    unaln_data{p} = table;
    for t = unaln_type
        disp(t)
        if ismember(t, sheets)
            d = readtable(fname, 'Sheet', t, 'VariableNamingRule', 'preserve');
            d.Polar = d.or1;
            d.Periphery = d.or2;
            d.note = string(d.note);
            d(d.note ~= "NaN" & d.note ~= "" & ~ismissing(d.note), :) = [];

            unaln_data{p} = [unaln_data{p}; d];
            curr_minmax(1) = min([curr_minmax(1) d.o1' d.o2']); % for o min
            curr_minmax(2) = max([curr_minmax(2) d.o1' d.o2']); % for o max
            curr_minmax(3) = min([curr_minmax(3) d.r1' d.r2']); % for r min
            curr_minmax(4) = max([curr_minmax(4) d.r1' d.r2']); % for r max
        end
    end
    
    all_minmax(p, :) = curr_minmax;
end


%% plot distirbution of all o and r
figure
tiledlayout(6, 4, "TileSpacing", "compact", "Padding", "compact");
for p = 1:length(prefix)
    d = align_data{p};

    nexttile
    histogram(d.o1, 5)
    xlim([0 6]*10^4)
    ylim([0 25])
    nexttile
    histogram(d.o2, 5)
    xlim([0 6]*10^4)
    ylim([0 25])

    nexttile
    histogram(d.r1, 5)
    xlim([0 4000])
    ylim([0 22])
    nexttile
    histogram(d.r2, 5)
    xlim([0 4000])
    ylim([0 22])
end

%% run min-max normalization for all data
for p = 1:length(prefix)
    align_data{p} = run_minmax(align_data{p}, all_minmax);
    unaln_data{p} = run_minmax(unaln_data{p}, all_minmax);
end

% adjust for bias due to different microsope
median_align6 = median(align_data{6}.diff);
median_unaln6 = median(unaln_data{6}.diff);

align_data{6}.diff = rescale(align_data{6}.diff, min(align_data{5}.diff), max(align_data{5}.diff));% - median_align6;
unaln_data{6}.diff = rescale(unaln_data{6}.diff, min(unaln_data{5}.diff), max(unaln_data{5}.diff));% - median_unaln6;

align_data{6}.diff = align_data{6}.diff + (median_align6 - median(align_data{6}.diff));
unaln_data{6}.diff = unaln_data{6}.diff + (median_unaln6 - median(unaln_data{6}.diff));


%% HD align
% Create a custom colormap bwr
nSteps = 25;

figure
for p = 1:length(prefix)
    % d = align_data{p};
    d = unaln_data{p};
    
    max_diff = max(abs(d.diff));

    d.norm_diff = round(interp1([-max_diff, 0, max_diff], [1, nSteps, 2*nSteps], d.diff, 'linear', 'extrap'));

    % Number of rows in the table
    numRows = size(d, 1);

    mycolor = [linspace(0, 1, numRows)', linspace(0, 1, numRows)', ones(numRows, 1)];
    
    nexttile
    hold on
    for i = 1:height(d)
        mycolor = bwr(d.norm_diff(i), :);
        edge_color = 'w';
        if d.diff(i) > 0
            edge_color = 'b';
            
        elseif d.diff(i) < 0
            edge_color = 'r';
        end

        
        plot([1, 2], [d.Polar(i) d.Periphery(i)], 'o-', "Color", mycolor, 'MarkerFaceColor', mycolor, 'MarkerEdgeColor', edge_color, 'LineWidth', 1);
    end
    hold off
    [~, pval] = ttest(d.Polar, d.Periphery, 'Tail', 'right');
    title({prefix(p), pval})
    xlim([0.5 2.5]);
    ylim([0 1])
end

%% regroup data

d = table( ...
    [repelem("Ctrl Aln",   height(align_data{1}.diff))'; ... 
     repelem("Ctrl Aln",   height(align_data{2}.diff))'; ... 
     repelem("Ctrl Aln",   height(align_data{3}.diff))'; ... 
     repelem("Ctrl Aln",   height(align_data{4}.diff))'; ... 
     repelem("Ctrl Unaln", height(unaln_data{1}.diff))'; ... 
     repelem("Ctrl Unaln", height(unaln_data{2}.diff))'; ... 
     repelem("Ctrl Unaln", height(unaln_data{3}.diff))'; ... 
     repelem("Ctrl Unaln", height(unaln_data{4}.diff))'; ... 
     repelem("AA Aln",     height(align_data{5}.diff))'; ... 
     repelem("AA Unaln",   height(unaln_data{5}.diff))'; ... 
     repelem("AY Aln",     height(align_data{6}.diff))'; ... 
     repelem("AY Unaln",   height(unaln_data{6}.diff))'], ...
    [align_data{1}.diff; ...
     align_data{2}.diff; ...
     align_data{3}.diff; ...
     align_data{4}.diff; ...
     unaln_data{1}.diff; ...
     unaln_data{2}.diff; ...
     unaln_data{3}.diff; ...
     unaln_data{4}.diff; ...
     align_data{5}.diff; ...
     unaln_data{5}.diff; ...
     align_data{6}.diff; ...
     unaln_data{6}.diff], ...
     'VariableNames', {'Exp', 'Diff'});

d.Exp = categorical(d.Exp, ["Ctrl Aln", "Ctrl Unaln", "AY Aln", "AY Unaln", "AA Aln", "AA Unaln"], 'Ordinal', true);

%% plot the data

figure
hold on
% swarmchart(categorical(d.Exp), d.Diff)
violinplot(categorical(d.Exp), d.Diff)
boxchart(categorical(d.Exp), d.Diff, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'BoxWidth', 0.15);
yline(0, '--')
% yticks([-0.75 -0.5 -0.25 0 0.25 0.5 0.75])
ylim([-0.75 0.85])
hold off

%% calculate p-vals

[~, pval0] = ttest2(d(d.Exp == "Ctrl Aln", :).Diff, d(d.Exp == "Ctrl Unaln", :).Diff, 'Tail', 'Right', 'Vartype', 'unequal');

[~, pval1] = ttest(d(d.Exp == "Ctrl Aln", :).Diff,   0, 'Tail', 'right');
[~, pval2] = ttest(d(d.Exp == "Ctrl Unaln", :).Diff, 0, 'Tail', 'right');
[~, pval3] = ttest(d(d.Exp == "AA Aln", :).Diff, 0, 'Tail', 'right');
[~, pval4] = ttest(d(d.Exp == "AA Unaln", :).Diff, 0, 'Tail', 'right');
[~, pval5] = ttest(d(d.Exp == "AY Aln", :).Diff, 0, 'Tail', 'right');
[~, pval6] = ttest(d(d.Exp == "AY Unaln", :).Diff, 0, 'Tail', 'right');

disp(pval0)
disp(pval1 + " " + pval2 + " " + pval3 + " " + pval4 + " " + pval5 + " " + pval6)


%% plot over time for control 0202
close all 

p = 1;

fname = "E:/PALP/PALP-control-0202_time_summary.xlsx";
ad = readtable(fname, 'Sheet', "Align", 'VariableNamingRule', 'preserve');

time_points = unique(ad.t);

a_mean_or1 = zeros(size(time_points));
a_std_or1  = zeros(size(time_points));
a_mean_or2 = zeros(size(time_points));
a_std_or2  = zeros(size(time_points));

for t = 1:length(time_points)
    or1 = ad.o1(ad.t == t) / max(ad.o1(ad.t == 2) + ad.r1(ad.t == 2));
    a_mean_or1(t) = mean(or1);
    a_std_or1(t)  = std(or1);

    or2 = ad.o2(ad.t == t) / max(ad.o1(ad.t == 2) + ad.r1(ad.t == 2));
    a_mean_or2(t) = mean(or2);
    a_std_or2(t)  = std(or2);
end

ud = readtable(fname, 'Sheet', "Unalign", 'VariableNamingRule', 'preserve');

u_mean_or1 = zeros(size(time_points));
u_std_or1  = zeros(size(time_points));
u_mean_or2 = zeros(size(time_points));
u_std_or2  = zeros(size(time_points));

for t = 1:length(time_points)
    or1 = ud.o1(ud.t == t) / max(ad.o1(ad.t == 2) + ad.r1(ad.t == 2));
    u_mean_or1(t) = mean(or1);
    u_std_or1(t)  = std(or1);

    or2 = ud.o2(ud.t == t) / max(ad.o1(ad.t == 2) + ad.r1(ad.t == 2));
    u_mean_or2(t) = mean(or2);
    u_std_or2(t)  = std(or2);
end

polar_color = [0 0.447 0.741];
perip_color = [0.929 0.694 0.125];

figure
nexttile
hold on
patch([time_points; flipud(time_points)], [a_mean_or2-a_std_or2; flipud(a_mean_or2+a_std_or2)], polar_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
patch([time_points; flipud(time_points)], [a_mean_or1-a_std_or1; flipud(a_mean_or1+a_std_or1)], perip_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(time_points, a_mean_or2, '-', 'Color', polar_color, 'LineWidth', 2);
plot(time_points, a_mean_or1, '-', 'Color', perip_color, 'LineWidth', 2);
ylim([0 1])
xticks([1:4])
xticklabels({"0", "20", "40", "60"})
xlabel("Time (s)")
axis square
hold off

nexttile
hold on
patch([time_points; flipud(time_points)], [u_mean_or2-u_std_or2; flipud(u_mean_or2+u_std_or2)], perip_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
patch([time_points; flipud(time_points)], [u_mean_or1-u_std_or1; flipud(u_mean_or1+u_std_or1)], polar_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(time_points, u_mean_or2, '-', 'Color', perip_color, 'LineWidth', 2);
plot(time_points, u_mean_or1, '-', 'Color', polar_color, 'LineWidth', 2);
ylim([0 1])
xticks([1:4])
xticklabels({"0", "20", "40", "60"})
xlabel("Time (s)")
axis square
hold off


%% FUNCTION
function d = run_minmax(d, all_minmax)

    all_min_o = min(all_minmax(:, 1));
    all_max_o = max(all_minmax(:, 2));
    all_min_r = min(all_minmax(:, 3));
    all_max_r = max(all_minmax(:, 4));

    all = [d.o1; d.r1; d.o2; d.r2];
    o = [d.o1; d.o2];
    r = [d.r1; d.r2];
    or = [d.or1; d.or2];

    d.o1_norm = min_max_normalization(d.o1, all_min_o, all_max_o, '', ''); % min(o), max(o),  % min(all), max(all), 
    d.r1_norm = min_max_normalization(d.r1, all_min_r, all_max_r, '', ''); % min(r), max(r),  % min(all), max(all), 
    d.o2_norm = min_max_normalization(d.o2, all_min_o, all_max_o, '', ''); % min(o), max(o),  % min(all), max(all), 
    d.r2_norm = min_max_normalization(d.r2, all_min_r, all_max_r, '', ''); % min(r), max(r),  % min(all), max(all), 
    
    d.or1_norm = d.o1_norm ./ (d.o1_norm + d.r1_norm + 1e-8);
    d.or2_norm = d.o2_norm ./ (d.o2_norm + d.r2_norm + 1e-8);

    d.Polar = d.or1_norm;
    d.Periphery = d.or2_norm;
    d.diff = d.or1_norm - d.or2_norm;

    % d.Polar = d.or1;
    % d.Periphery = d.or2;
    % d.diff = d.or1 - d.or2;    
end

function I = min_max_normalization(I, Imin, Imax, a, b)
    I = (I - Imin) / (Imax - Imin); % * (b-a) + a;
end

function bwr = mybwr()
    blue = [0, 0, 1];
    white = [1, 1, 1];
    red = [1, 0, 0];
    
    nSteps = 25;
    
    bwr = [linspace(blue(1), white(1), nSteps)' linspace(blue(2), white(2), nSteps)' linspace(blue(3), white(3), nSteps)'; ...
            linspace(white(1), red(1), nSteps)' linspace(white(2), red(2), nSteps)' linspace(white(3), red(3), nSteps)'];
    bwr = flipud(bwr);
end
