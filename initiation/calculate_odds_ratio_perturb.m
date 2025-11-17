%
% Odds of initiation occuring at cellular patterns
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath D:\Spatiotemporal_analysis\code
addpath D:\Spatiotemporal_analysis\code\util

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATIONS
ini_size = 150;

pattern = ["phalf" "mhalf" "splay" "bend" "aligned"];

bin_sz = 25;

%% photoinduction
% 
% load initiation from photoinduction
% 

prefix = ["photoinduct-1" "photoinduct-2" "photoinduct-3" "photoinduct-4"];
all_prefix = "photoinduct-all";
imsize = {[5000 5000], [5000 5000], [5000 5000], [5000 5000]};

bg_rep = 10; %1;
res_rep = 11; %1;

outpath = ROOTDIR + "/initiation_statistics/photo/photos_pattern_odds_bg"+bg_rep+"_res"+res_rep+"/";
% outpath = ROOTDIR + "/initiation_statistics/perturb-photo/photos_pattern_odds_bg"+bg_rep+"_res"+res_rep+"/";
mkdir(outpath);

% load initiations
ini_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix + "_initiation_details_pattern_per5_cov10_sp150.txt";
ini = readtable(ini_file, 'Delimiter', '\t');

% load random selected non-initiation background file
bg_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix + "_randbg" + bg_rep + "_details_pattern_per5_cov10_sp" + ini_size + ".txt";
bg = readtable(bg_file, 'Delimiter', '\t');

% load resistant files
res_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix + "_resist" + res_rep + "_details_pattern_per5_cov10_sp" + ini_size + ".txt";
res = readtable(res_file, 'Delimiter', '\t');

[ctrl_ini, ctrl_ini_n, ctrl_bg, ctrl_bg_n, ctrl_res, ctrl_res_n] = load_by_perturb('', ini, bg, res, prefix, imsize, ini_size);

% %% calculate odds ratio
% % Control
% ctrl_ini_or = calculate_odds_ratio(pattern, ctrl_ini, ctrl_ini_n, ctrl_bg, ctrl_bg_n);
% ctrl_res_or = calculate_odds_ratio(pattern, ctrl_res, ctrl_res_n, ctrl_bg, ctrl_bg_n);

%% perturbation
% 
%   load all photoinduction files for TRULI, KI696, MY and GNE
%
prefix1 = ["photo-truli-KI-1" "photo-truli-KI-4" "photo-truli-ki-0730-1" "photo-truli-ki-0730-4"];
all_prefix1 = "photo-truli-KI-all";
imsize1 = {[4800 4800], [5120 5120], [4800 4800], [5120 5120]};

prefix2 = ["yapi-0816" "yapi-0819" "yapi-0827z4"];
all_prefix2 = "yapi-all";
imsize2 = {[5120 5120], [5120 5120], [5120 5120]};

bg_rep = 1;
res_rep = 1;

outpath = ROOTDIR + "/initiation_statistics/perturb-photo/photos_pattern_odds_bg"+bg_rep+"_res"+res_rep+"/";
mkdir(outpath);

% load initiatios
ini1_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix1 + "_initiation_details_pattern_per5_cov10_sp150.txt";
ini1 = readtable(ini1_file, 'Delimiter', '\t');

% load random selected non-initiation background file
bg1_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix1 + "_randbg" + bg_rep + "_details_pattern_per5_cov10_sp" + ini_size + ".txt";
bg1 = readtable(bg1_file, 'Delimiter', '\t');

% load resistant files
res1_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix1 + "_resist" + res_rep + "_details_pattern_per5_cov10_sp" + ini_size + ".txt";
res1 = readtable(res1_file, 'Delimiter', '\t');

% load initiatios
ini2_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix2 + "_initiation_details_pattern_per5_cov10_sp150.txt";
ini2 = readtable(ini2_file, 'Delimiter', '\t');

% load random selected non-initiation background file
bg2_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix2 + "_randbg" + bg_rep + "_details_pattern_per5_cov10_sp" + ini_size + ".txt";
bg2 = readtable(bg2_file, 'Delimiter', '\t');

% load resistant files
res2_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix2 + "_resist" + res_rep + "_details_pattern_per5_cov10_sp" + ini_size + ".txt";
res2 = readtable(res2_file, 'Delimiter', '\t');

%% load each perturbation

% Control 1
[ctrl1_ini, ctrl1_ini_n, ctrl1_bg, ctrl1_bg_n, ctrl1_res, ctrl1_res_n] = load_by_perturb('Control', ini1, bg1, res1, prefix1, imsize1, ini_size);

% TRULI
[truli_ini, truli_ini_n, truli_bg, truli_bg_n, truli_res, truli_res_n] = load_by_perturb('TRULI', ini1, bg1, res1, prefix1, imsize1, ini_size);

% KI696
[ki696_ini, ki696_ini_n, ki696_bg, ki696_bg_n, ki696_res, ki696_res_n] = load_by_perturb('KI696', ini1, bg1, res1, prefix1, imsize1, ini_size);

% Control 2
[ctrl2_ini, ctrl2_ini_n, ctrl2_bg, ctrl2_bg_n, ctrl2_res, ctrl2_res_n] = load_by_perturb('Control', ini2, bg2, res2, prefix2, imsize2, ini_size);

% MY
[my_ini, my_ini_n, my_bg, my_bg_n, my_res, my_res_n] = load_by_perturb('MY', ini2, bg2, res2, prefix2, imsize2, ini_size);

% GNE
[gne_ini, gne_ini_n, gne_bg, gne_bg_n, gne_res, gne_res_n] = load_by_perturb('GNE', ini2, bg2, res2, prefix2, imsize2, ini_size);

% % Control
% ctrl_ini = [ctrl1_ini; ctrl2_ini];
% ctrl_bg  = [ctrl1_bg;  ctrl2_bg];
% ctrl_res = [ctrl1_res; ctrl2_res];
% 
% ctrl_ini_n = ctrl1_ini_n + ctrl2_ini_n;
% ctrl_bg_n  = ctrl1_bg_n  + ctrl2_bg_n;
% ctrl_res_n = ctrl1_res_n + ctrl2_res_n;

% % All Control
% % Combine photoinduction-1-4 and ctrls of perturbation experiments
% ctrl_ini = [ctrl_ini; ctrl1_ini; ctrl2_ini];
% ctrl_bg  = [ctrl_bg;  ctrl1_bg;  ctrl2_bg];
% ctrl_res = [ctrl_res; ctrl1_res; ctrl2_res];
% 
% ctrl_ini_n = ctrl_ini_n + ctrl1_ini_n + ctrl2_ini_n;
% ctrl_bg_n  = ctrl_bg_n +  ctrl1_bg_n  + ctrl2_bg_n;
% ctrl_res_n = ctrl_res_n + ctrl1_res_n + ctrl2_res_n;

%% calculate odds ratio

% Control
ctrl_ini_or = calculate_odds_ratio(pattern, ctrl_ini, ctrl_ini_n, ctrl_bg, ctrl_bg_n);
ctrl_res_or = calculate_odds_ratio(pattern, ctrl_res, ctrl_res_n, ctrl_bg, ctrl_bg_n);

% Control 1
ctrl1_ini_or = calculate_odds_ratio(pattern, ctrl1_ini, ctrl1_ini_n, ctrl1_bg, ctrl1_bg_n);
ctrl1_res_or = calculate_odds_ratio(pattern, ctrl1_res, ctrl1_res_n, ctrl1_bg, ctrl1_bg_n);

% Control 2
ctrl2_ini_or = calculate_odds_ratio(pattern, ctrl2_ini, ctrl2_ini_n, ctrl2_bg, ctrl2_bg_n);
ctrl2_res_or = calculate_odds_ratio(pattern, ctrl2_res, ctrl2_res_n, ctrl2_bg, ctrl2_bg_n);

truli_ini_or = calculate_odds_ratio_by_sp(truli_ini, truli_ini_n, truli_bg, truli_bg_n);
truli_res_or = calculate_odds_ratio_by_sp(truli_res, truli_res_n, truli_bg, truli_bg_n);

gne_ini_or   = calculate_odds_ratio_by_sp(gne_ini, gne_ini_n, gne_bg, gne_bg_n);
gne_res_or   = calculate_odds_ratio_by_sp(gne_res, gne_res_n, gne_bg, gne_bg_n);

ki696_ini_or = calculate_odds_ratio_by_sp(ki696_ini, ki696_ini_n, ki696_bg, ki696_bg_n);
ki696_res_or = calculate_odds_ratio_by_sp(ki696_res, ki696_res_n, ki696_bg, ki696_bg_n);

%%
%{'LL', 'LH', 'HL', 'HH'};
or = truli_ini_or;

figure
hold on
for sp = 2:height(or)
    errorbar(sp, or(sp,1), or(sp,2), or(sp,3), 'bo')
    plot(sp, or(sp,1), '.', 'MarkerSize', 30)
end
hold off
yline(1, '--', 'Color', '#808080')
% ylim([0 3.5])
% xlim([1 10])
xticks(2:5)


ylabel("Odds ratio of occurrence")

%% calculate odds ratio

% Control
ctrl_ini_or = calculate_odds_ratio(pattern, ctrl_ini, ctrl_ini_n, ctrl_bg, ctrl_bg_n);
ctrl_res_or = calculate_odds_ratio(pattern, ctrl_res, ctrl_res_n, ctrl_bg, ctrl_bg_n);

% Control 1
ctrl1_ini_or = calculate_odds_ratio(pattern, ctrl1_ini, ctrl1_ini_n, ctrl1_bg, ctrl1_bg_n);
ctrl1_res_or = calculate_odds_ratio(pattern, ctrl1_res, ctrl1_res_n, ctrl1_bg, ctrl1_bg_n);

% TRULI
truli_ini_or = calculate_odds_ratio(pattern, truli_ini, truli_ini_n, truli_bg, truli_bg_n);
truli_res_or = calculate_odds_ratio(pattern, truli_res, truli_res_n, truli_bg, truli_bg_n);

% KI696
ki696_ini_or = calculate_odds_ratio(pattern, ki696_ini, ki696_ini_n, ki696_bg, ki696_bg_n);
ki696_res_or = calculate_odds_ratio(pattern, ki696_res, ki696_res_n, ki696_bg, ki696_bg_n);

% Control 2
ctrl2_ini_or = calculate_odds_ratio(pattern, ctrl2_ini, ctrl2_ini_n, ctrl2_bg, ctrl2_bg_n);
ctrl2_res_or = calculate_odds_ratio(pattern, ctrl2_res, ctrl2_res_n, ctrl2_bg, ctrl2_bg_n);

% TRULI
my_ini_or = calculate_odds_ratio(pattern, my_ini, my_ini_n, my_bg, my_bg_n);
my_res_or = calculate_odds_ratio(pattern, my_res, my_res_n, my_bg, my_bg_n);

% KI696
gne_ini_or = calculate_odds_ratio(pattern, gne_ini, gne_ini_n, gne_bg, gne_bg_n);
gne_res_or = calculate_odds_ratio(pattern, gne_res, gne_res_n, gne_bg, gne_bg_n);

%% plot odds ratio and confience interval
savefigure = 0;

% or_list = {ctrl_ini_or, ctrl_res_or, ctrl1_ini_or, ctrl1_res_or, truli_ini_or, truli_res_or, ki696_ini_or, ki696_res_or, ...
%     ctrl2_ini_or, ctrl2_res_or, my_ini_or, my_res_or, gne_ini_or, gne_res_or};
% perturb_list = ["ctrl", "ctrl", "ctrl1", "ctrl1", "truli", "truli", "ki696", "ki696", ...
%     "ctrl2", "ctrl2", "my", "my", "gne", "gne",];
% ini_res = ["ini" "res"+res_rep "ini" "res"+res_rep "ini" "res"+res_rep "ini" "res"+res_rep ...
%     "ini" "res"+res_rep "ini" "res"+res_rep "ini" "res"+res_rep ];

% or_list = {ctrl_ini_or, ctrl_res_or};
or_list = {truli_ini_or, truli_res_or}
perturb_list = ["ctrl", "ctrl"];
ini_res = ["ini" "res"+res_rep "ini"];

myColor = ["red", "blue", "#7030A0", "#7030A0", "black"];
x_range = 1:length(pattern);

for i = 1:length(or_list)
    or = or_list{i};

    figure
    for pttn = 1:length(pattern)
        hold on
        errorbar(pttn, or(pttn,1), or(pttn,2), or(pttn,3), 'bo', 'Color', myColor(pttn), 'MarkerFaceColor', myColor(pttn), 'LineWidth', 1.5)
        plot(pttn, or(pttn,1), '.', 'MarkerSize', 30, 'Color', myColor(pttn))
        yline(1, '--', 'Color', '#808080')
        ylim([0 3.5])
        xlim([0 6])
        xticks(x_range)
    %     yticklabels(pattern(pttn))
        ylabel("Odds ratio of occurrence")
        hold off
    end
    
    if savefigure
        outname = perturb_list(i)+"_pattern_odds_"+ini_size+"_photo_"+ini_res(i)+"_bg"+bg_rep;
        saveas(gcf, outpath+outname+".jpg", 'jpg')
        saveas(gcf, outpath+outname+".svg", 'svg')
    end
end

%% top pattern vs maln
%
%   process topological patterns vs densities and alignment
%

pttn_ini = ctrl_ini(ctrl_ini.phalf_150 > 0 | ctrl_ini.mhalf_150 > 0 | ctrl_ini.splay_150 > 0 | ctrl_ini.bend_150 > 0, :);
pttn_bg  = ctrl_bg(ctrl_bg.phalf_150 > 0   | ctrl_bg.mhalf_150 > 0  | ctrl_bg.splay_150 > 0  | ctrl_bg.bend_150 > 0, :);
pttn_res = ctrl_res(ctrl_res.phalf_150 > 0 | ctrl_res.mhalf_150 > 0 | ctrl_res.splay_150 > 0 | ctrl_res.bend_150 > 0, :);

npttn_ini = ctrl_ini(ctrl_ini.pattern == 0 | ctrl_ini.pattern == 5, :);
npttn_bg  = ctrl_bg(ctrl_bg.pattern == 0   | ctrl_bg.pattern == 5, :);
npttn_res = ctrl_res(ctrl_res.pattern == 0 | ctrl_res.pattern == 5, :);

% for initiation patterns
sp_pttn_or  = calculate_odds_ratio_by_sp(pttn_ini,  ctrl_ini_n, pttn_bg,  ctrl_bg_n); % height(pttn_ini) height(pttn_bg)
sp_npttn_or = calculate_odds_ratio_by_sp(npttn_ini, ctrl_ini_n, npttn_bg, ctrl_bg_n); % height(npttn_ini) height(npttn_bg)

% for resistent patterns
sp_pttn_r_or  = calculate_odds_ratio_by_sp(pttn_res,  ctrl_res_n, pttn_bg,  ctrl_bg_n); % height(pttn_res) height(pttn_bg)
sp_npttn_r_or = calculate_odds_ratio_by_sp(npttn_res, ctrl_res_n, npttn_bg, ctrl_bg_n); % height(npttn_res) height(npttn_bg)

%% plot for top pattern
%
%    plot odds ratio for topological patterns and without patterns
%

savefigure = 0;

pttn_or = sp_pttn_r_or;
npttn_or = sp_npttn_r_or;

mycol1 = '#ED7D31'; % '#ED7D31'; % #000000
mycol2 = '#843C0C'; % '#843C0C'; % '#595959'

figure
x_range = 1:length(pttn_or)
x_range = x_range - 0.1;

hold on
errorbar(x_range, pttn_or(:,1), pttn_or(:,2), pttn_or(:,3), 'o', 'LineWidth', 1.5, 'Color', mycol1)
plot(x_range, pttn_or(:,1), '.', 'MarkerSize', 30, 'LineWidth', 1, 'Color', mycol1)

x_range = x_range + 0.2;
errorbar(x_range, npttn_or(:,1), npttn_or(:,2), npttn_or(:,3), 'o', 'LineWidth', 1.5, 'Color', mycol2)
plot(x_range, npttn_or(:,1), '.', 'MarkerSize', 30, 'LineWidth', 1, 'Color', mycol2)
hold off

yline(1, '--', 'Color', '#808080')
ylim([0 3.5])
xlim([1.5 length(pttn_or)+0.5])
% xticks(x_range)
ylabel("Odds ratio of occurrences")

if savefigure
    outname = "ctrl_pttn-npttn_odds_"+ini_size+"_photo_ini_bg"+bg_rep;
    saveas(gcf, outpath+outname+".jpg", 'jpg')
    saveas(gcf, outpath+outname+".svg", 'svg')
end

%% FUNCTIONS
function or = calculate_odds_ratio_by_sp(ini, n_all, bg, bg_n_all)
    ini = ini(~isnan(ini.ent), :);
    bg = bg(~isnan(bg.ent), :);

    % % 0 20 40 60 80 100
    % d_thres = [7.0248   53.5752   63.7038   73.2201   85.7291  177.8219];
    % e_thres = [1.1417   -1.2398   -1.6740   -2.0369   -2.4471   -4.9712];

    % 0 25 50 75 100
    % d_thres = [7.0248 56.3820 68.3350 81.9792 177.8219];
    % e_thres = [1.1704 -1.3658 -1.8562 -2.3319 -4.9712];
    d_thres = [7.0248 56.3820 81.9792 177.8219];
    e_thres = [1.1704 -1.3658 -2.3319 -4.9712];
    combinations = {'LL', 'LM', 'LH' 'ML', 'MM', 'MH', 'HL', 'HM', 'HH'};

    % % 0 33 66 100
    % e_thres = [1.1704 -1.5587 -2.1684 -4.9712];
    % d_thres = [7.0248 60.3984 76.4317 177.8219];

    % % min mean max
    % e_thres = [1.1704 -1.8294 -4.9712];
    % d_thres = [7.0248 70.7203 177.8219];
    % combinations = {'LL', 'LH', 'HL', 'HH'};

    len = length(d_thres)-1;

    or = nan(length(combinations), 5);

    for e = 1:length(e_thres)-1
        for d = 1:length(d_thres)-1
            n_prop = numel( ini.ent( (ini.den > d_thres(d)) & (ini.den <= d_thres(d+1)) & (ini.ent <= e_thres(e)) & (ini.ent > e_thres(e+1)) ) );

            bg_n_prop = numel( bg.ent( (bg.den > d_thres(d)) & (bg.den <= d_thres(d+1)) & (bg.ent <= e_thres(e)) & (bg.ent > e_thres(e+1)) ) );

            fisher_data = [n_prop bg_n_prop; (n_all-n_prop) (bg_n_all-bg_n_prop)];
            [~,~,stats] = fishertest(fisher_data);

            or((e-1)*len+1+d, 1) = stats.OddsRatio;
            or((e-1)*len+1+d, 2) = stats.OddsRatio-stats.ConfidenceInterval(1);
            or((e-1)*len+1+d, 3) = stats.ConfidenceInterval(2)-stats.OddsRatio;
            or((e-1)*len+1+d, 4) = n_prop;
            or((e-1)*len+1+d, 5) = bg_n_prop;

        end
    end
end

function or = calculate_odds_ratio_by_prop(prop, ini, n_all, bg, bg_n_all)
    if strcmp(prop, 'den')
        pcol = 16;
        thres = [11.9358   53.5752   63.7038   73.2201   85.7291  153.8894];
    elseif strcmp(prop, 'ent')
        pcol = 14;
        thres = [-4.5307   -2.4471   -2.0369   -1.6740   -1.2398    1.1417]; %
    end

    ini = ini(~isnan(ini.(pcol)), :);
    bg = bg(~isnan(bg.(pcol)), :);

    pmean = mean(bg.(pcol), 'omitnan');
    pstd  = std(bg.(pcol), 'omitnan');

    or = nan(length(thres)-1, 3);

    for i = 1:length(thres)-1
        n_prop = numel(ini.(pcol)((ini.(pcol) > thres(i)) & (ini.(pcol) <= thres(i+1))));
        bg_n_prop = numel(bg.(pcol)((bg.(pcol) > thres(i)) & (bg.(pcol) <= thres(i+1))));

        % disp(n_prop+" "+bg_n_prop+" "+(n_all-n_prop)+" "+(bg_n_all-bg_n_prop))

        fisher_data = [n_prop bg_n_prop; (n_all-n_prop) (bg_n_all-bg_n_prop)];
        [hyp,pval,stats] = fishertest(fisher_data);

        or(i, 1) = stats.OddsRatio;
        or(i, 2) = stats.OddsRatio-stats.ConfidenceInterval(1);
        or(i, 3) = stats.ConfidenceInterval(2)-stats.OddsRatio;
        or(i, 4) = n_prop;
        or(i, 5) = bg_n_prop;
    end
end

function or =  calculate_odds_ratio(pattern, ini, n_all, bg, bg_n_all)
    pttn_col = 8;

    or = nan(length(pattern), 3);
    
    for pttn = 1:length(pattern)
        ini_pttn = ini.(pttn_col+pttn-1);
        n_pttn = sum(ini_pttn)/pttn;
    
        bg_ini_pttn = bg.(pttn_col+pttn-1);
        bg_n_pttn = sum(bg_ini_pttn)/pttn;

        fisher_data = [n_pttn bg_n_pttn; (n_all-n_pttn) (bg_n_all-bg_n_pttn)];
        [hyp,pval,stats] = fishertest(fisher_data);

        or(pttn, 1) = stats.OddsRatio;
        or(pttn, 2) = stats.OddsRatio-stats.ConfidenceInterval(1);
        or(pttn, 3) = stats.ConfidenceInterval(2)-stats.OddsRatio;
        or(pttn, 4) = n_pttn;
        or(pttn, 5) = bg_n_pttn;
    end
end

function [p_ini, p_ini_n, p_bg, p_bg_n, p_res, p_res_n] = load_by_perturb(perturb, ini, bg, res, prefix, imsize, ini_size)
    if strcmp(perturb, '')
        [p_ini, p_ini_n] = preprocess_ini_file(ini, prefix, imsize, ini_size);
        [p_bg,  p_bg_n]  = preprocess_ini_file(bg,  prefix, imsize, ini_size);
        [p_res, p_res_n] = preprocess_ini_file(res, prefix, imsize, ini_size);
    else
        idx = contains(ini.experiment, {perturb});
        p_ini = ini(idx, :);
        [p_ini, p_ini_n] = preprocess_ini_file(p_ini, prefix, imsize, ini_size);
    
        idx = contains(bg.experiment, {perturb});
        p_bg = bg(idx, :);
        [p_bg, p_bg_n] = preprocess_ini_file(p_bg, prefix, imsize, ini_size);
    
        idx = contains(res.experiment, {perturb});
        p_res = res(idx, :);
        [p_res, p_res_n] = preprocess_ini_file(p_res, prefix, imsize, ini_size);
    end
end

function [ini, n] = preprocess_ini_file(ini, prefix, imsize, sz)
    % n = 0;
    % offset = 51;
    % for i = 1:height(ini)
    %     p = find(prefix==ini.id(i));
    % 
    %     w = ini.well(i);
    %     x = ceil(ini.x(i));
    %     y = ceil(ini.y(i));
    % 
    %     imsz = imsize{p};
    % 
    %     if x-sz/2 < offset || y-sz/2 < offset || x+sz/2-1 > imsz(1)-offset || y+sz/2-1 > imsz(2)-offset || isnan(ini.ent(i))
    %         % ini.ent(i) = NaN;
    %         continue;
    %     end
    %     n = n + 1;
    % end
    % ini(isnan(ini.ent), :) = [];
    
    ini(isnan(ini.ent), :) = [];
    n = height(ini);

    ini.id = string(ini.id);
    ini.group = strcat(ini.id, "_", ini.experiment, "_", compose("%02d", ini.well));
end
