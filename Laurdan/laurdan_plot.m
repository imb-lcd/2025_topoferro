%
%   plot quantified length of lipid order
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
addpath(ROOTDIR + "/code");
addpath(ROOTDIR + "/code/util");
ROOTDIR = "E:/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% load the file for ridge of lipid order
file = ROOTDIR + "/Laurdan/laurdan_ridge_final.txt";
laurdan = readtable(file, "Delimiter", "\t");
laurdan(strcmp(laurdan.rep, "ctrl-aln-1"), :) = [];

laurdan.rep = categorical(laurdan.rep);

%% plot
% 
% reorder_exp = reordercats(laurdan.rep, {'ctrl-aln-2' 'ctrl-unaln-3' 'AA-aln-6' 'AA-unaln-2' ...
%                                     'AY-aln-1d-3' 'AY-unaln-1d'} );

figure
hold on
l = laurdan(contains(string(laurdan.rep), "ctrl-aln-2"), :);
boxchart(repelem(categorical(1), height(l)), l.length, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'BoxWidth', 0.2);
% swarmchart(repelem(categorical(1), height(l)), l.length, 'SizeData', 20, "MarkerEdgeColor", "k", "MarkerFaceColor", "k", "MarkerFaceAlpha", 0.5);
violinplot(repelem(categorical(1), height(l)), l.length);

l = laurdan(contains(string(laurdan.rep), "ctrl-unaln-3"), :);
boxchart(repelem(categorical(2), height(l)), l.length, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'BoxWidth', 0.2);
% swarmchart(repelem(categorical(2), height(l)), l.length, 'SizeData', 20, "MarkerEdgeColor", "k", "MarkerFaceColor", "k", "MarkerFaceAlpha", 0.5);
violinplot(repelem(categorical(2), height(l)), l.length);

l = laurdan(contains(string(laurdan.rep), "AA-aln-6"), :);
boxchart(repelem(categorical(3), height(l)), l.length, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'BoxWidth', 0.2);
% swarmchart(repelem(categorical(3), height(l)), l.length, 'SizeData', 20, "MarkerEdgeColor", "b", "MarkerFaceColor", "b", "MarkerFaceAlpha", 0.5);
violinplot(repelem(categorical(3), height(l)), l.length);

l = laurdan(contains(string(laurdan.rep), "AA-unaln-2"), :);
boxchart(repelem(categorical(4), height(l)), l.length, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'BoxWidth', 0.2);
% swarmchart(repelem(categorical(4), height(l)), l.length, 'SizeData', 20, "MarkerEdgeColor", "b", "MarkerFaceColor", "b", "MarkerFaceAlpha", 0.5);
violinplot(repelem(categorical(4), height(l)), l.length);

l = laurdan(contains(string(laurdan.rep), "AY-aln-1d-3"), :);
boxchart(repelem(categorical(5), height(l)), l.length, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'BoxWidth', 0.2);
% swarmchart(repelem(categorical(5), height(l)), l.length, 'SizeData', 20, "MarkerEdgeColor", "#77AC30", "MarkerFaceColor", "#77AC30", "MarkerFaceAlpha", 0.5);
violinplot(repelem(categorical(5), height(l)), l.length);

l = laurdan(contains(string(laurdan.rep), "AY-unaln-1d"), :);
boxchart(repelem(categorical(6), height(l)), l.length, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k', 'BoxWidth', 0.2);
% swarmchart(repelem(categorical(6), height(l)), l.length, 'SizeData', 20, "MarkerEdgeColor", "#77AC30", "MarkerFaceColor", "#77AC30", "MarkerFaceAlpha", 0.5);
violinplot(repelem(categorical(6), height(l)), l.length);

ylim([0 1500])
xticklabels(["Control Aln" "Control Unaln" "AA Aln" "AA Unaln" "AY Aln" "AY Unaln"])
xtickangle(90)

hold off

%%
pval1 = ranksum(laurdan(contains(string(laurdan.rep), "ctrl-aln-2"), :).length,  laurdan(contains(string(laurdan.rep), "ctrl-unaln-3"), :).length, "tail", "right");
pval2 = ranksum(laurdan(contains(string(laurdan.rep), "AA-aln-6"), :).length,    laurdan(contains(string(laurdan.rep), "AA-unaln-2"), :).length,   "tail", "right");
pval3 = ranksum(laurdan(contains(string(laurdan.rep), "AY-aln-1d-3"), :).length, laurdan(contains(string(laurdan.rep), "AY-unaln-1d"), :).length,  "tail", "right");
pval4 = ranksum(laurdan(contains(string(laurdan.rep), "ctrl-aln-2"), :).length,  laurdan(contains(string(laurdan.rep), "AA-aln-6"), :).length,     "tail", "right");
pval5 = ranksum(laurdan(contains(string(laurdan.rep), "ctrl-aln-2"), :).length,  laurdan(contains(string(laurdan.rep), "AY-aln-1d-3"), :).length,  "tail", "right");
pval6 = ranksum(laurdan(contains(string(laurdan.rep), "ctrl-unaln-3"), :).length,  laurdan(contains(string(laurdan.rep), "AA-unaln-2"), :).length,     "tail", "right");
pval7 = ranksum(laurdan(contains(string(laurdan.rep), "ctrl-unaln-3"), :).length,  laurdan(contains(string(laurdan.rep), "AY-unaln-1d"), :).length,  "tail", "right");

d1 = meanEffectSize(laurdan(contains(string(laurdan.rep), "ctrl-aln-2"), :).length,  laurdan(contains(string(laurdan.rep), "ctrl-unaln-3"), :).length, "Effect", "cohen");
d2 = meanEffectSize(laurdan(contains(string(laurdan.rep), "ctrl-aln-2"), :).length,  laurdan(contains(string(laurdan.rep), "AA-aln-6"), :).length, "Effect", "cohen");
d3 = meanEffectSize(laurdan(contains(string(laurdan.rep), "ctrl-aln-2"), :).length,  laurdan(contains(string(laurdan.rep), "AY-aln-1d-3"), :).length, "Effect", "cohen");

d4 = meanEffectSize(laurdan(contains(string(laurdan.rep), "ctrl-unaln-3"), :).length,  laurdan(contains(string(laurdan.rep), "AA-unaln-2"), :).length, "Effect", "cohen");
d5 = meanEffectSize(laurdan(contains(string(laurdan.rep), "ctrl-unaln-3"), :).length,  laurdan(contains(string(laurdan.rep), "AY-unaln-1d"), :).length, "Effect", "cohen");

disp(pval1 + " " + pval2 + " " + pval3  + " " + pval4  + " " + pval5 + " " + pval6  + " " + pval7)
disp(d1.Effect + " " + d2.Effect + " " + d3.Effect  + " " + d4.Effect  + " " + d5.Effect)