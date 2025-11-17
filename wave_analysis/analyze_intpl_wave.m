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

%% load wave information table

wvinfo_fname = "all_wave_info.txt";

wvinfo = readtable(PATH + wvinfo_fname);

% exclude waves that will not be analyzed
wvinfo(wvinfo.include == 0, :) = [];

%% load the streamline and boundary information
stream = cell(height(wvinfo),1);
bo = cell(height(wvinfo), 1);
nframe = cell(height(wvinfo), 1);

for ii = 1:height(wvinfo)
    % load streamline info
    fname = ROOTDIR + "wave_" + wvinfo.Set(ii) + "/" + ...
        "Well" + sprintf('%02d', wvinfo.well(ii)) + "/wave" + wvinfo.wave(ii) + "_analyses/" + ...
        set_filename(wvinfo.Set(ii), wvinfo.well(ii), wvinfo.start_frame(ii)+"-"+wvinfo.end_frame(ii), "") + "_streamline_table.mat";
    

    stream{ii} = struct2array(load(fname));

    % load wave boundary info
    fname = ROOTDIR + "wave_" + wvinfo.Set(ii) + "/" + ...
        "Well" + sprintf('%02d', wvinfo.well(ii)) + "/wave" + wvinfo.wave(ii) + "_analyses/" + ...
        set_filename(wvinfo.Set(ii), wvinfo.well(ii), wvinfo.start_frame(ii)+"-"+wvinfo.end_frame(ii), "") + "_bo.mat";

    bo{ii} = struct2array(load(fname));

    nframe{ii} = length(bo{ii});
end

%% load the smooth interpolated wave
disp("load the interpolated data")
reduced_sz = 50; % equivalent of windows size of 10 pixel (12.66 um)

all_data = [];

for ii = 1:height(wvinfo)
    % load the wave data information
    prefix = wvinfo.Set(ii);
    w = wvinfo.well(ii);
    wv = wvinfo.wave(ii);
    fr = wvinfo.spatial_frame(ii);
    frame_range = wvinfo.start_frame(ii)+"-"+wvinfo.end_frame(ii);

    path = ROOTDIR + "wave_" + prefix + "/" + "Well" + sprintf('%02d', w) + "/";

    ws = IMSIZE/reduced_sz;

    if ii == 15
        wv_name = path + "/wave" + wv + "_analyses/" + ...
            set_filename(prefix, w, frame_range, '') + "_wvdata_intpl" + ws + "_wCoh.mat";
    else
        wv_name = path + "/wave" + wv + "_analyses/" + ...
            set_filename(prefix, w, frame_range, '') + "_wvdata_intpl" + ws + ".mat";
    end

    wv_data = struct2array(load(wv_name));

    if sum(ismember(wv_data.Properties.VariableNames, 'angden'))
        wv_data.angden = [];
    end
    if sum(ismember(wv_data.Properties.VariableNames, 'coh'))
        wv_data.coh = [];
    end

    all_data = [all_data; wv_data];
end

% remove rows with NaNs to get only the wave data
wave = rmmissing(all_data);

%% plot scatter plot of directional discordance, density, and colored by speed
savefigure = 0;

figure
scatter(wave.angdiff, wave.density, 20, wave.speed, 'filled', 'MarkerFaceAlpha', 0.5)
xlim([0 0.85]);
ylim([30 130]);
xlabel('along alignment');
ylabel('density')
axis square

myparula = generate_black_parula(25);
colormap(myparula);


%% plot interpolated plot of directional discordance, density, and colored by speed
nbin = 50;

x = wave.angdiff;
y = wave.density;
v = wave.speed;

v_max = mean(v) + 4.5*std(v);

nbin = 40;
den_upr = 84.823; 
den_lwr = 29;

[xq, yq] = meshgrid(linspace(0, 1.6, nbin), linspace(den_lwr, den_upr, nbin));

[Xq,Yq,vq] = griddata(x, y, v, xq, yq, 'nearest');

vg = imgaussfilt(vq, 5);
%% 
figure
s = pcolor(xq, yq, vg);
set(s, 'FaceColor', 'interp', 'EdgeColor', 'none');

xlabel('along alignment');
ylabel('density');
axis square tight
myparula = generate_black_parula(25);
colormap(myparula);
colorbar

yticks([29 42.41150082 56.54866776 70.68583471 84.82300165])
yticklabels([1 1.5 2 2.5 3])

%% fit multiple regression model
a = normalize(wave.angdiff);
d = normalize(wave.density);
ad = a .* d;

corr(a, d, 'Rows', 'complete')

mdl = fitlm([a d ad], normalize(wave.speed))

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


%% FUNCTIONS
function [varargout] = block_nondead_cells(wave_bound, vector_per_side, IMSIZE, varargin)
    varargout = cell((nargin)-3, 1);
    for i = 1:(nargin-3)
        lg = imresize(varargin{i}, [IMSIZE IMSIZE], 'bilinear');
        lg(~wave_bound) = NaN;
        varargout{i} = imresize(lg, [vector_per_side vector_per_side], "bilinear");
    end
end

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

% idx = find(y>125 & x>1.2);
% v(idx) = v(idx)-20; 
% idx = find(y>90 & y<110 & x>1.3);
% v(idx) = v(idx)-5; 


















