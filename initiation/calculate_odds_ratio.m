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

pattern = ["phalf" "mhalf" "splay" "bend" "aligned" "misaligned"];

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

% load initiations
s_ini_file = "D:/Spatiotemporal_analysis/Initiation/spontaneous_initiation_details_pattern_per5_cov10_sp150.txt";
s_ini = readtable(s_ini_file, 'Delimiter', '\t');

ini = [ini; s_ini];

% load random selected non-initiation background file
bg_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix + "_randbg" + bg_rep + "_details_pattern_per5_cov10_sp" + ini_size + ".txt";
bg = readtable(bg_file, 'Delimiter', '\t');

% load resistant files
res_file = "D:/Spatiotemporal_analysis/Initiation/" + all_prefix + "_resist" + res_rep + "_details_pattern_per5_cov10_sp" + ini_size + ".txt";
res = readtable(res_file, 'Delimiter', '\t');

[ctrl_ini, ctrl_ini_n, ctrl_bg, ctrl_bg_n, ctrl_res, ctrl_res_n] = load_by_perturb('', ini, bg, res, prefix, imsize, ini_size);

%% calculate odds ratio
% Control
ctrl_ini_or = calculate_odds_ratio(pattern, ctrl_ini, ctrl_ini_n, ctrl_bg, ctrl_bg_n);
ctrl_res_or = calculate_odds_ratio(pattern, ctrl_res, ctrl_res_n, ctrl_bg, ctrl_bg_n);

%% calculate odds ratio maln
% Control
ctrl_ini_or = calculate_odds_ratio_maln(pattern, ctrl_ini, ctrl_ini_n, ctrl_bg, ctrl_bg_n);
ctrl_res_or = calculate_odds_ratio_maln(pattern, ctrl_res, ctrl_res_n, ctrl_bg, ctrl_bg_n);

%% calculate odds ratio spatial
% Control
ctrl_ini_or = calculate_odds_ratio_by_sp(ctrl_ini, ctrl_ini_n, ctrl_bg, ctrl_bg_n);
ctrl_res_or = calculate_odds_ratio_by_sp(ctrl_res, ctrl_res_n, ctrl_bg, ctrl_bg_n);

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

bg_rep = 2;
res_rep = 2;

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

bg_rep = 1;
res_rep = 1;

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

% Control
ctrl_ini = [ctrl1_ini; ctrl2_ini];
ctrl_bg  = [ctrl1_bg;  ctrl2_bg];
ctrl_res = [ctrl1_res; ctrl2_res];

ctrl_ini_n = ctrl1_ini_n + ctrl2_ini_n;
ctrl_bg_n  = ctrl1_bg_n  + ctrl2_bg_n;
ctrl_res_n = ctrl1_res_n + ctrl2_res_n;

% %% All Control
% % Combine photoinduction-1-4 and ctrls of perturbation experiments
% ctrl_ini = [ctrl_ini; ctrl1_ini; ctrl2_ini];
% ctrl_bg  = [ctrl_bg;  ctrl1_bg;  ctrl2_bg];
% ctrl_res = [ctrl_res; ctrl1_res; ctrl2_res];
% 
% ctrl_ini_n = ctrl_ini_n + ctrl1_ini_n + ctrl2_ini_n;
% ctrl_bg_n  = ctrl_bg_n +  ctrl1_bg_n  + ctrl2_bg_n;
% ctrl_res_n = ctrl_res_n + ctrl1_res_n + ctrl2_res_n;

%% calculate odds ratio

truli_ini_or = calculate_odds_ratio_by_sp(truli_ini, truli_ini_n, truli_bg, truli_bg_n);
truli_res_or = calculate_odds_ratio_by_sp(truli_ini, truli_ini_n, truli_bg, truli_bg_n);
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
ylim([0 3])
% xlim([1 10])
xticks(2:5)


ylabel("Odds ratio of occurrence")

% calculate_odds_ratio_by_sp(pttn_ini,  ctrl_ini_n, pttn_bg,  ctrl_bg_n); 

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

% MY
my_ini_or = calculate_odds_ratio(pattern, my_ini, my_ini_n, my_bg, my_bg_n);
my_res_or = calculate_odds_ratio(pattern, my_res, my_res_n, my_bg, my_bg_n);

% GNE
gne_ini_or = calculate_odds_ratio(pattern, gne_ini, gne_ini_n, gne_bg, gne_bg_n);
gne_res_or = calculate_odds_ratio(pattern, gne_res, gne_res_n, gne_bg, gne_bg_n);

%% calculate odds ratio maln

% Control
ctrl_ini_or = calculate_odds_ratio_maln(pattern, ctrl_ini, ctrl_ini_n, ctrl_bg, ctrl_bg_n);
ctrl_res_or = calculate_odds_ratio_maln(pattern, ctrl_res, ctrl_res_n, ctrl_bg, ctrl_bg_n);

% Control 1
ctrl1_ini_or = calculate_odds_ratio_maln(pattern, ctrl1_ini, ctrl1_ini_n, ctrl1_bg, ctrl1_bg_n);
ctrl1_res_or = calculate_odds_ratio_maln(pattern, ctrl1_res, ctrl1_res_n, ctrl1_bg, ctrl1_bg_n);

% TRULI
truli_ini_or = calculate_odds_ratio_maln(pattern, truli_ini, truli_ini_n, truli_bg, truli_bg_n);
truli_res_or = calculate_odds_ratio_maln(pattern, truli_res, truli_res_n, truli_bg, truli_bg_n);

% KI696
ki696_ini_or = calculate_odds_ratio_maln(pattern, ki696_ini, ki696_ini_n, ki696_bg, ki696_bg_n);
ki696_res_or = calculate_odds_ratio_maln(pattern, ki696_res, ki696_res_n, ki696_bg, ki696_bg_n);

% Control 2
ctrl2_ini_or = calculate_odds_ratio_maln(pattern, ctrl2_ini, ctrl2_ini_n, ctrl2_bg, ctrl2_bg_n);
ctrl2_res_or = calculate_odds_ratio_maln(pattern, ctrl2_res, ctrl2_res_n, ctrl2_bg, ctrl2_bg_n);

% MY
my_ini_or = calculate_odds_ratio_maln(pattern, my_ini, my_ini_n, my_bg, my_bg_n);
my_res_or = calculate_odds_ratio_maln(pattern, my_res, my_res_n, my_bg, my_bg_n);

% GNE
gne_ini_or = calculate_odds_ratio_maln(pattern, gne_ini, gne_ini_n, gne_bg, gne_bg_n);
gne_res_or = calculate_odds_ratio_maln(pattern, gne_res, gne_res_n, gne_bg, gne_bg_n);

%% plot odds ratio and confience interval
savefigure = 1;

% or_list = {ctrl_ini_or, ctrl_res_or, ctrl1_ini_or, ctrl1_res_or, truli_ini_or, truli_res_or, ki696_ini_or, ki696_res_or, ...
%     ctrl2_ini_or, ctrl2_res_or, my_ini_or, my_res_or, gne_ini_or, gne_res_or};
% perturb_list = ["ctrl", "ctrl", "ctrl1", "ctrl1", "truli", "truli", "ki696", "ki696", ...
%     "ctrl2", "ctrl2", "my", "my", "gne", "gne",];
% ini_res = ["ini" "res"+res_rep "ini" "res"+res_rep "ini" "res"+res_rep "ini" "res"+res_rep ...
%     "ini" "res"+res_rep "ini" "res"+res_rep "ini" "res"+res_rep ];

or_list = {ctrl_ini_or, ctrl_res_or};
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
        ylabel("Odds ratio of occurrence")
        hold off
    end
    
    if savefigure
        outname = perturb_list(i)+"_pattern_odds_"+ini_size+"_photo_"+ini_res(i)+"_bg"+bg_rep;
        saveas(gcf, outpath+outname+".jpg", 'jpg')
        saveas(gcf, outpath+outname+".svg", 'svg')
    end
end

%% FUNCTIONS
function or = calculate_odds_ratio_by_sp(ini, n_all, bg, bg_n_all)
    ini = ini(~isnan(ini.ent), :);
    bg = bg(~isnan(bg.ent), :);

    % % 0 25 50 75 100
    d_thres = [-Inf 56.3820 81.9792 Inf];
    e_thres = [Inf -1.3658 -2.3319 -Inf];
    combinations = {'LL', 'LM', 'LH' 'ML', 'MM', 'MH', 'HL', 'HM', 'HH'};

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

function or = calculate_odds_ratio_maln(pattern, ini, n_all, bg, bg_n_all)
    pttn_col = 8;

    e_thres = -1.3658;

    or = nan(length(pattern), 3);

    for pttn = 1:length(pattern)-2
        n_pttn    = numel(ini.ent( (ini.ent > e_thres) & contains(string(ini.pattern), string(pttn)) ));

        bg_n_pttn = numel(bg.ent(  (bg.ent > e_thres) & contains(string(bg.pattern), string(pttn)) ));

        fisher_data = [n_pttn bg_n_pttn; (n_all-n_pttn) (bg_n_all-bg_n_pttn)];
        [hyp,pval,stats] = fishertest(fisher_data);

        or(pttn, 1) = stats.OddsRatio;
        or(pttn, 2) = stats.OddsRatio-stats.ConfidenceInterval(1);
        or(pttn, 3) = stats.ConfidenceInterval(2)-stats.OddsRatio;
        or(pttn, 4) = n_pttn;
        or(pttn, 5) = bg_n_pttn;
    end

    n_prop    = numel(ini.ent( (ini.pattern == 5) ) );
    bg_n_prop = numel(bg.ent(  (bg.pattern == 5)  ) );

    fisher_data = [n_prop bg_n_prop; (n_all-n_prop) (bg_n_all-bg_n_prop)];
    [hyp,pval,stats] = fishertest(fisher_data);

    or(length(pattern)-1, 1) = stats.OddsRatio;
    or(length(pattern)-1, 2) = stats.OddsRatio-stats.ConfidenceInterval(1);
    or(length(pattern)-1, 3) = stats.ConfidenceInterval(2)-stats.OddsRatio;
    or(length(pattern)-1, 4) = n_pttn;
    or(length(pattern)-1, 5) = bg_n_pttn;

    n_prop    = numel(ini.ent( (ini.ent > e_thres) & (ini.pattern == 0) ) );
    bg_n_prop = numel(bg.ent(  (bg.ent > e_thres) & (bg.pattern == 0)   ) );

    fisher_data = [n_prop bg_n_prop; (n_all-n_prop) (bg_n_all-bg_n_prop)];
    [~,~,stats] = fishertest(fisher_data);

    or(length(pattern), 1) = stats.OddsRatio;
    or(length(pattern), 2) = stats.OddsRatio-stats.ConfidenceInterval(1);
    or(length(pattern), 3) = stats.ConfidenceInterval(2)-stats.OddsRatio;
    or(length(pattern), 4) = n_prop;
    or(length(pattern), 5) = bg_n_prop;
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
   
    ini(isnan(ini.ent), :) = [];
    n = height(ini);

    ini.id = string(ini.id);
    ini.group = strcat(ini.id, "_", ini.experiment, "_", compose("%02d", ini.well));
end
