%
%   calculate wave isotropy for figure 4
%
clearvars
clc
close all

addpath("D:/Spatiotemporal_analysis/code/util/");

ROOTDIR = "D:/Spatiotemporal_analysis/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATION

den_reg = 120;
MCHCH = "c1";
DENPATH = MCHCH + "_mCherry_density/";

CY5CH = "c2";
CY5ADJPATH = CY5CH + "_cy5_adj/";

DICCH = "c3";
NEMPATH = DICCH + "_DIC_nematics/";
ORIENTPATH = NEMPATH + "orientation_files/";


%% load wave information table

wvinfo_fname = ROOTDIR + "/wave_0_analyses/" + "all_wave_info.txt";

wvinfo = readtable(wvinfo_fname);

% exclude waves that will not be analyzed
% wvinfo(wvinfo.include == 0, :) = [];
wvinfo(43, :) = [];


%% load wave initiation

wv_ini_fname = ROOTDIR + "/Initiation/" + "ctrl_wvinfo_ini.txt";

wv_ini = readtable(wv_ini_fname);
wv_ini(43, :) = [];

%% load BO
showfigure = 0;

wv_bo = cell(height(wvinfo), 1);

if showfigure, figure; end

for ii = 1:height(wvinfo)

    prefix = wvinfo.Set(ii);
    w = wvinfo.well(ii);
    wv = wvinfo.wave(ii);
    fr = wvinfo.spatial_frame(ii);
    frame_range = wvinfo.start_frame(ii)+"-"+wvinfo.end_frame(ii);
    
    path = ROOTDIR + "wave_" + prefix + "/" + "Well" + sprintf('%02d', w) + "/";

    wv_bo_file = path + "/wave" + wv + "_analyses/" + ...
            set_filename(prefix, w, frame_range, '') + "_bo.mat";
    wv_bo{ii} = struct2array(load(wv_bo_file));
end

%% Use the ends points of each wave to calculate MZC

showfigure = false;

sel_fr = 9;

IMSIZE = 5000;

ctrl_mzc      = nan(height(wvinfo), 1);
ctrl_mzc_norm = nan(height(wvinfo), 1);

for ii = 1:height(wvinfo)
% for ii = 6
    bo = wv_bo{ii};
    % bo = bo{sel_fr};
    bo = bo{end};
    
    mask = poly2mask(bo(:, 2), bo(:, 1), IMSIZE, IMSIZE);

    ini = [wv_ini.x(ii) wv_ini.y(ii)];

    [ctrl_mzc(ii), ctrl_mzc_norm(ii), inscribeR, circumR] = calc_mzc(bo, mask, ini, [NaN NaN], [NaN NaN]);

    if showfigure
        figure
        imshow(mask)
        hold on
        viscircles(ini, circumR, 'Color', 'b')
        viscircles(ini, inscribeR, 'Color', 'b')
        scatter(ini(1), ini(2), 'b', 'filled')
        hold off
        axis square off
    end
end

%% load perturb wave info

pwvinfo_fname = ROOTDIR + "/wave_0_analyses/" + "/perturb_wave_info_all.txt";
pwvinfo = readtable(pwvinfo_fname);

% exclude waves that will not be analyzed
pwvinfo = pwvinfo(~strcmp(pwvinfo.perturb, "AA") & ~contains(pwvinfo.perturb, "AA2b"), :);
% pwvinfo([35 36 39], :) = [];
pwvinfo = pwvinfo(~(strcmp(pwvinfo.perturb, "AY") & pwvinfo.include == 0), :);

% load wave initiation
pwv_ini_fname = ROOTDIR + "/Initiation/" + "perturb_wvinfo_ini.txt";
pwv_ini = readtable(pwv_ini_fname);
pwv_ini = pwv_ini(~(strcmp(pwv_ini.perturb, "AY") & pwv_ini.include == 0), :);

%% load BO perturb
pwv_bo = cell(height(pwvinfo), 1);

for ii = 1:height(pwvinfo)

    prefix = pwvinfo.Set(ii);
    w = pwvinfo.well(ii);
    wv = pwvinfo.wave(ii);
    fr = pwvinfo.spatial_frame(ii);
    frame_range = sprintf("%02d", pwvinfo.start_frame(ii))+"-"+sprintf("%02d", pwvinfo.end_frame(ii));
    
    path = ROOTDIR + "wave_" + prefix + "/" + "Well" + sprintf('%02d', w) + "/";

    wv_bo_file = path + "/wave" + wv + "_analyses/" + ...
            set_filename(prefix, w, frame_range, '') + "_bo.mat";
    pwv_bo{ii} = struct2array(load(wv_bo_file));
end

%% calculate MZC perturb
fr_max = 6;

perturb_mzc      = nan(height(pwvinfo), 1);
perturb_mzc_norm = nan(height(pwvinfo), 1);

for ii = 1:height(pwvinfo)
    bo = pwv_bo{ii};
    sel_fr = min(pwvinfo.num_frames(ii), fr_max);
    
    % bo = bo{sel_fr};
    bo = bo{end};

    IMSIZE = pwvinfo.sz(ii);
    
    mask = poly2mask(bo(:, 2), bo(:, 1), IMSIZE, IMSIZE);

    ini = [pwv_ini.x(ii) pwv_ini.y(ii)];

    [perturb_mzc(ii), perturb_mzc_norm(ii), inscribeR, circumR] = calc_mzc(bo, mask, ini, [NaN NaN], [NaN NaN]);
end


%% create table

ctrl = [repelem("Control", length(ctrl_mzc))' ctrl_mzc ctrl_mzc_norm];

aa_idx = find(pwvinfo.perturb == "AA2a");
aa =  [repelem("AA", length(aa_idx))' perturb_mzc(aa_idx) perturb_mzc_norm(aa_idx)];

ay_idx = find(pwvinfo.perturb == "AY");
ay =  [repelem("AY", length(ay_idx))' perturb_mzc(ay_idx) perturb_mzc_norm(ay_idx)];

truli_idx = find(pwvinfo.perturb == "TRULI");
truli =  [repelem("TRULI", length(truli_idx))' perturb_mzc(truli_idx) perturb_mzc_norm(truli_idx)];

iso = [ctrl; aa; ay; truli];
iso = array2table(iso);
iso.Properties.VariableNames = ["cond", "mzc", "mzc_norm"];
iso.cond     = categorical(iso.cond);
iso.mzc      = double(iso.mzc);
iso.mzc_norm = double(iso.mzc_norm);

order = ["Control", "AY", "AA", "TRULI"];
iso.cond = categorical(iso.cond, order, 'Ordinal', true);

%% plot wave anistropy
% One-sided (right tail)
p12 = ranksum(double(ctrl(:, 3)), double(ay(:, 3)),    'tail', 'right');
p13 = ranksum(double(ctrl(:, 3)), double(aa(:, 3)),    'tail', 'right');
p14 = ranksum(double(ctrl(:, 3)), double(truli(:, 3)), 'tail', 'right');

% Two-sided (default)
p23 = ranksum(double(ay(:, 3)), double(aa(:, 3)));
p24 = ranksum(double(ay(:, 3)), double(truli(:, 3)));
p34 = ranksum(double(aa(:, 3)), double(truli(:, 3)));

disp("ranksum test with control, ay: " + p12 + " ; aa: " + p13 + " ; truli: " + p14);
disp("ranksum test  aa-ay: " + p23 + " ; aa-truli: " + p24 + " ; ay-truli: " + p34);


cohen12 = meanEffectSize(double(ctrl(:, 3)), double(ay(:, 3)),    "Effect", "cohen");
cohen13 = meanEffectSize(double(ctrl(:, 3)), double(aa(:, 3)),    "Effect", "cohen");
cohen14 = meanEffectSize(double(ctrl(:, 3)), double(truli(:, 3)), "Effect", "cohen");

disp("cohen's with control, ay: "+cohen12.Effect+" ; aa: "+cohen13.Effect+" ; truli: "+cohen14.Effect)

cohen23 = meanEffectSize(double(ay(:, 3)), double(aa(:, 3)));
cohen24 = meanEffectSize(double(ay(:, 3)), double(truli(:, 3)));
cohen34 = meanEffectSize(double(aa(:, 3)), double(truli(:, 3)));

disp("cohen's with aa-ay: "+cohen23.Effect+" ; aa-truli: "+cohen24.Effect+" ; ay-truli: "+cohen34.Effect)


figure
hold on

mycol = [125 125 125;
         255 148 23;
         198 28 50;
         112 48 160] / 255; 

[g, condNames] = findgroups(iso.cond);

for i = 1:numel(condNames)
    idx = (g == i);
    x = repmat(i, sum(idx), 1);

    % Swarm
    swarmchart(iso.cond(idx), iso.mzc_norm(idx), 20, ...
        'filled', ...
        'MarkerFaceColor', mycol(i,:), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.75);


    % Violin
    v = violinplot(x, iso.mzc_norm(idx));
    v.FaceColor = mycol(i,:);
    v.EdgeColor = mycol(i,:);
    v.FaceAlpha = 0.2;

    % Box
    boxchart(iso.cond(idx), iso.mzc_norm(idx), ...
        'BoxFaceColor', mycol(i,:), ...
        'BoxEdgeColor', 'k', ...
        'BoxWidth', 0.15, ...
        'MarkerStyle','none');
end

ylim([0 1])
ylabel("Wave Anisotropy (Norm MZC)")

hold off

%% FUNCTION
function [mzc, mzc_norm, inscribe_radius, circum_radius] = calc_mzc(Bo, BW, ini, exclude_ln, exclude_pt)
    [H, W] = size(BW);

    x0 = ini(1);
    y0 = ini(2);

    bx = Bo(:,2);   % x = col
    by = Bo(:,1);   % y = row

    tol = 1.5;

    % Detect whether center is near image edges
    onLeft   = x0 <= 1 + tol;
    onRight  = x0 >= W - tol;
    onTop    = y0 <= 1 + tol;
    onBottom = y0 >= H - tol;

    % Inward-facing arc for circumradius
    % For edge/corner seeds, ignore the outward side
    valid_arc = true(size(bx));

    if onLeft
        valid_arc = valid_arc & (bx >= x0);
    end
    if onRight
        valid_arc = valid_arc & (bx <= x0);
    end
    if onTop
        valid_arc = valid_arc & (by >= y0);
    end
    if onBottom
        valid_arc = valid_arc & (by <= y0);
    end

    % Distance from center to boundary
    dist = hypot(bx - x0, by - y0);

    % Circumradius from valid inward-facing arc
    dist_arc = dist(valid_arc);

    circum_radius = max(dist_arc);

    % Inscribed radius
    % Ignore artificial boundary segments on image border
    on_img_border = ...
        (bx <= 1 + tol) | ...
        (bx >= W - tol) | ...
        (by <= 1 + tol) | ...
        (by >= H - tol);

    valid_for_inscribe = ~on_img_border;

    % half-plane exclusion for inscribed circle only
    if any(~isnan(exclude_ln))
        exc_x1 = exclude_ln(1,1);
        exc_y1 = exclude_ln(1,2);
        exc_x2 = exclude_ln(2,1);
        exc_y2 = exclude_ln(2,2);
        
        % Signed side test relative to the line P1->P2
        side_bad = point_line_side(exclude_pt(1), exclude_pt(2), exc_x1, exc_y1, exc_x2, exc_y2);
        side_bdy = point_line_side(bx, by, exc_x1, exc_y1, exc_x2, exc_y2);

        sideTol = 1e-9;

        if abs(side_bad) > sideTol
            sameSideAsBad = sign(side_bdy) == sign(side_bad);
            sameSideAsBad(abs(side_bdy) <= sideTol) = false;  % keep points on the line
            valid_for_inscribe = valid_for_inscribe & ~sameSideAsBad;
        end
    end

    dist_in = dist(valid_for_inscribe);

    inscribe_radius = min(dist_in);

    % MZC outputs
    mzc = circum_radius - inscribe_radius;

    mzc_norm = mzc / circum_radius;
end


function s = point_line_side(x, y, x1, y1, x2, y2)
%POINT_LINE_SIDE Signed side of point(s) relative to directed line P1->P2.
% Positive and negative values indicate opposite sides; zero means on line.

    s = (x - x1) .* (y2 - y1) - (y - y1) .* (x2 - x1);
end
