clc; clear; close all;

cd("../../Jen-Hao/Laurdan/scripts_ccwu19/")

% ── Constants ──────────────────────────────────────────────────────────────
RIDGE_DIR   = fullfile(".", "ridges");
CHECK_DIR   = fullfile(RIDGE_DIR, "for_checking");
GP_RGB_PATH = fullfile(".", "cell_counts", "rgb_gp.tif");
LUT_PATH    = fullfile(".", "lut", "ICA.mat");

GP_IMIN     = 0;
GP_IMAX     = 195 / 255;
MONTAGE_ROWS = 4;

% ── Load data ──────────────────────────────────────────────────────────────
% tbl = readtable("fname_slices.csv");
load("cell_counts.mat");
load(LUT_PATH);

% ══════════════════════════════════════════════════════════════════════════
%  Main loop
% ══════════════════════════════════════════════════════════════════════════
for i = 1  % 1:height(tbl)

    % -- Row metadata --
    time  = tbl.time(i);
    well  = tbl.well(i);
    slice = tbl.slice(i);

    % ── Load ROIs ─────────────────────────────────────────────────────────
    pat     = sprintf("%02dh_s%02dz%02d", time, well, slice);
    entries = dir(fullfile(RIDGE_DIR, "*.zip"));
    entries = entries(contains({entries.name}, pat));
    ridgePath = fullfile(entries.folder, entries.name);

    roi      = ReadImageJROI(ridgePath);
    isRidge  = cellfun(@(s) s.strType == "FreeLine", roi);
    nRidges  = nnz(isRidge);

    % ── Extract ridge curves ───────────────────────────────────────────────
    curves = cell(nRidges, 1);
    idx    = 0;
    for j = 1:numel(roi)
        if ~isRidge(j), continue; end
        idx = idx + 1;
        x = roi{j}.mnCoordinates(:, 1);
        y = roi{j}.mnCoordinates(:, 2);
        curves{idx} = [x, y];
    end

    tbl.summary{i} = summarize_curve_image(curves, []);

    % ── Visual check ──────────────────────────────────────────────────────
    gp_rgb = tiffreadVolume(GP_RGB_PATH, "PixelRegion", {[1 Inf], [1 Inf], [i i]});
    gp_rgb = squeeze(gp_rgb);

    I_adj = im2double(gp_rgb);
    I_adj = (I_adj - GP_IMIN) / (GP_IMAX - GP_IMIN);
    I_adj = min(max(I_adj, 0), 1);           % clamp to [0, 1]

    % Capture axis limits and figure size from a hidden reference render
    close all;
    fig_ref = figure("Visible", "off");
    imshow(I_adj);
    xlimit = xlim;
    ylimit = ylim;
    fig_pos = get(fig_ref, "Position");

    % Draw ridge orientations over a blank axes
    figure;
    hold on;
    for j = 1:numel(roi)
        if string(roi{j}.strType) ~= "FreeLine", continue; end

        x = roi{j}.mnCoordinates(:, 1)';
        y = roi{j}.mnCoordinates(:, 2)';
        [~, theta] = calc_curvature(x, y);

        % Color-code by local orientation angle
        surface( ...
            [x; x], [y; y], ...
            zeros(2, numel(x)), ...
            [theta; theta], ...
            "FaceColor", "none", ...
            "EdgeColor", "interp", ...
            "LineWidth", 5);
    end
    hold off;

    colormap("hsv");  clim([0, pi]);
    axis square ij off;  box on;
    xlim(xlimit);  ylim(ylimit);
    set(gcf, "Position", fig_pos);

    copygraphics(gcf, "ContentType", "image", "Resolution", 300);

    % -- Optionally save to disk --
    % img_out = sprintf("%02dh_s%02d_z%02d.png", time, well, slice);
    % exportgraphics(gcf, fullfile(CHECK_DIR, img_out));

end
% save("ridge_length.mat", "tbl");

% ══════════════════════════════════════════════════════════════════════════
%  Montage of saved check images
% ══════════════════════════════════════════════════════════════════════════
entries = dir(CHECK_DIR);
entries = entries(startsWith({entries.name}, "48"));
imgList = arrayfun(@(s) string(fullfile(s.folder, s.name)), entries);

close all;
montage(imageDatastore(imgList), "Size", [MONTAGE_ROWS, ceil(numel(imgList) / MONTAGE_ROWS)]);
copygraphics(gcf, "ContentType", "image", "Resolution", 300);