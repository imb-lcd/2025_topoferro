%
%   calculate wave circularity
%
clearvars
clc
close all

addpath("D:/Spatiotemporal_analysis/code/util/");

ROOTDIR = "D:/Spatiotemporal_analysis/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% load wave information
ini_file = ROOTDIR + "\Initiation\diffday_initiation_details_adj.txt";
ini = readtable(ini_file);


fname_12h = "E:\movies_diff_day\Figure1\12h_diffday-1023-1d_s24\" + ...
    "diffday-1023-1d_s24t05c2_ORG_seg_bw_perwave.tif";

bw_12h = imread(fname_12h);


fname_24h = "E:\movies_diff_day\Figure1\24h_diffday-1118-1d_s14\" + ...
    "diffday-1118-1d_s14tAllc2_ORG_seg_bw_perwave.tif";

bw_24h = imread(fname_24h);


fname_48h = "E:\movies_diff_day\Figure1\48h_mCh300FBS10-1_s10\" + ...
    "mCh300FBS10-1_s10t52c2_ORG_seg_bw_perwave.tif";

bw_48h = imread(fname_48h);

%% calculate isotropy
IMSIZE = 5076;
[idx_12h, wv_area_12h, circ_12h, aniso_12h, aniso_norm_12h, aniso_area_12h, eccen_12h, pa_12h] = calc_isotropy(bw_12h, IMSIZE, ini, "12h");

[idx_24h, wv_area_24h, circ_24h, aniso_24h, aniso_norm_24h, aniso_area_24h, eccen_24h, pa_24h] = calc_isotropy(bw_24h, IMSIZE, ini, "24h");

IMSIZE = 4977;
[idx_48h, wv_area_48h, circ_48h, aniso_48h, aniso_norm_48h, aniso_area_48h, eccen_48h, pa_48h] = calc_isotropy(bw_48h, IMSIZE, ini, "48h");

%%

header = ["cond", "wv_area", "circ", "aniso", "aniso_norm", "aniso_area", "eccen", "pa_ratio"];

cond = categorical(repelem("12h", length(wv_area_12h))');
iso_12h = table(cond, wv_area_12h, circ_12h, aniso_12h, aniso_norm_12h, aniso_area_12h, eccen_12h, pa_12h);
iso_12h.Properties.VariableNames = header;

cond = categorical(repelem("24h", length(wv_area_24h))');
iso_24h = table(cond, wv_area_24h, circ_24h, aniso_24h, aniso_norm_24h, aniso_area_24h, eccen_24h, pa_24h);
iso_24h.Properties.VariableNames = header;

cond = categorical(repelem("48h", length(wv_area_48h))');
iso_48h = table(cond, wv_area_48h, circ_48h, aniso_48h, aniso_norm_48h, aniso_area_48h, eccen_48h, pa_48h);
iso_48h.Properties.VariableNames = header;

iso = [iso_12h; iso_24h; iso_48h];

%%
[~, p12] = ttest2(iso(iso.cond == "12h", :).aniso_norm, iso(iso.cond == "24h", :).aniso_norm, "Tail", "left");
[~, p23] = ttest2(iso(iso.cond == "24h", :).aniso_norm, iso(iso.cond == "48h", :).aniso_norm, "Tail", "left");
[~, p13] = ttest2(iso(iso.cond == "12h", :).aniso_norm, iso(iso.cond == "48h", :).aniso_norm, "Tail", "left");

disp("ranksum test for 12v24: " + p12 + " ; 24v48: " + p13 + " ; 12v48: " + p13);

cohen12 = meanEffectSize(iso(iso.cond == "12h", :).aniso_norm, iso(iso.cond == "24h", :).aniso_norm, "Effect", "cohen");
cohen23 = meanEffectSize(iso(iso.cond == "24h", :).aniso_norm, iso(iso.cond == "48h", :).aniso_norm, "Effect", "cohen");
cohen13 = meanEffectSize(iso(iso.cond == "12h", :).aniso_norm, iso(iso.cond == "48h", :).aniso_norm, "Effect", "cohen");

disp("cohen's for 12v24: "+cohen12.Effect+" ; 24vs48: "+cohen23.Effect+" ; 12v48"+cohen13.Effect)

figure
hold on
mycol = linspace(0.2, 0.8, 4)' * [1 1 1];

[g, condNames] = findgroups(iso.cond);

for i = 1:numel(condNames)
    idx = (g == i);
    x = repmat(i, sum(idx), 1);

    % % Swarm
    swarmchart(iso.cond(idx), iso.aniso_norm(idx), 20, ...
        'filled', ...
        'MarkerFaceColor', mycol(i,:), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.75);

    % Box
    boxchart(iso.cond(idx), iso.aniso_norm(idx), ...
        'BoxFaceColor', mycol(i,:), ...
        'BoxEdgeColor', 'k', ...
        'BoxWidth', 0.15, ...
        'MarkerStyle','none');
end

ylim([0 1])
ylabel("Wave Anisotropy")

hold off


%% FUNCTION
function [idx, wv_area, circ, aniso, aniso_norm, aniso_area, eccen, pa] = calc_isotropy(bw, IMSIZE, ini, seed_type)
    ini = ini(ini.part == seed_type, :);
    do_aniso = 1;

    sigma = 15;
    binarize_threshold = 0.75; 
    [bdy, Bo] = wave_boundary_gaussfilt(bw, 1, sigma, binarize_threshold);
    Bo = Bo{1};

    %%
    idx      = 1:length(Bo);
    wv_area  = nan(length(Bo), 1);
    circ     = nan(length(Bo), 1);
    aniso      = nan(length(Bo), 1);
    aniso_norm = nan(length(Bo), 1);
    aniso_area = nan(length(Bo), 1);
    eccen    = nan(length(Bo), 1);
    pa       = nan(length(Bo), 1);
    
    figure
    imshow(zeros(IMSIZE)+255)
    hold on

    for i = 1:length(Bo)
        bo = Bo{i};

        if isempty(bo)
            disp(i)
            continue;
        end
    
        x = bo(:, 2);
        y = bo(:, 1);
    
        plot(x, y, "k-");
        
    
        mask = poly2mask(x, y, IMSIZE, IMSIZE);
    
        stats = regionprops(mask, 'Area', 'Circularity', 'Eccentricity', 'Perimeter');
        
        circ(i)    = stats.Circularity;
        eccen(i)   = stats.Eccentricity;
        pa(i)      = stats.Perimeter ./ stats.Area;
        wv_area(i) = stats.Area;
        
        
        wv_ini_info = ini(ini.index == i, :);
        wv_ini = [wv_ini_info.x wv_ini_info.y];
        idx = find(ini.index == i);


        exclude_ln = [ini.mask1_x(idx) ini.mask1_y(idx); ini.mask2_x(idx) ini.mask2_y(idx)];
        exclude_pt = [ini.x_bad(idx) ini.y_bad(idx)];

        if do_aniso
            [aniso(i), aniso_norm(i), inscribeR, circumR] = calc_aniso(bo, mask, wv_ini, exclude_ln, exclude_pt);
   
    
            drawcircle("Center", wv_ini, "Radius", inscribeR);
            drawcircle("Center", wv_ini, "Radius", circumR, 'Color', 'k');
            scatter(wv_ini(1)+1, wv_ini(2)+1, 'b', 'LineWidth', 2)
    
            aniso_area(i) = wv_area(i) ./ (pi.*circumR.^2);

            text(mean(x), mean(y), string(i)+": "+sprintf("%0.2f", aniso_norm(i)))
        else
            text(mean(x), mean(y), string(i))
        end
    end
    hold off

end

function [aniso, aniso_norm, inscribe_radius, circum_radius] = calc_aniso(Bo, BW, ini, exclude_ln, exclude_pt)
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
    aniso = circum_radius - inscribe_radius;

    aniso_norm = aniso / circum_radius;
end


function s = point_line_side(x, y, x1, y1, x2, y2)
%POINT_LINE_SIDE Signed side of point(s) relative to directed line P1->P2.
% Positive and negative values indicate opposite sides; zero means on line.

    s = (x - x1) .* (y2 - y1) - (y - y1) .* (x2 - x1);
end


function [bdy, Bo] = wave_boundary_gaussfilt(bw, nframes, sigma, binarize_threshold)
    sm_area_thres = 3000;
    sm_hole_thres = 10000;

    bdy = cell(nframes, 1);
    Bo = cell(nframes, 1);

    run_medfilt = 1;
    if ~exist('medfilt_sz', 'var')
        run_medfilt = 0;
    end

    for i = 1:nframes
        B = imgaussfilt(double(bw), sigma);

        BW = imbinarize(B, binarize_threshold);    
        BW = bwareaopen(BW, sm_area_thres);         % Remove small objects and speckles
        BW = ~bwareaopen(~BW, sm_hole_thres); 

        bdy{i} = sparse(BW);
       
        bo = bwboundaries(BW);
        for j = 1:length(bo)
            if length(bo{j}) < 300
                bo{j} = [];
                continue;
            end

            b = bo{j};
            b(:, 2) = smoothdata(b(:, 2), 'sgolay', 50)';
            b(:, 1) = smoothdata(b(:, 1), 'sgolay', 50)';
            bo{j} = b;
        end
        


        Bo{i} = bo;
    end
end