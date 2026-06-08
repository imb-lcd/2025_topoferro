clc; clear; close all;

%% ── Configuration ────────────────────────────────────────────────────────────
DIR_SEGMENTED = "../img/2_segmented_nuc";
DIR_DENSITY   = "../img/4_density";
DIR_ORIENT    = "../img/5_orientation";

wells        = [49,50,51,52,53,54,55,56,57,58,59,60,85,86,87,88,89,90,91,92,93,94,95,96];
frameDensity = ones(size(wells));
frameDeath   = repelem({1:36}, 1, numel(wells));

SIGNAL  = "mch";   % "mch" | "cy5"
W       = 75;      % half-window size (pixels)
N_SAMPLE = 1e4;    % number of random windows per well
PAD     = 24;      % extra margin beyond W
MARGIN  = W + PAD;

SAVE_OUTPUT = false;
outFile = sprintf("data_sample_windows_%s.mat", SIGNAL);

%% ── Storage ──────────────────────────────────────────────────────────────────
windowTblByWell = cell(numel(wells), 1);
countTblByWell  = cell(numel(wells), 1);

%% ── Main loop (parallel over wells) ─────────────────────────────────────────
tic
for i = 1:numel(wells)
    wellID = wells(i);

    %── Load static maps ─────────────────────────────────────────────────────
    D = load(fullfile(DIR_DENSITY, sprintf("s%02d_mch_t%02d.mat", wellID, frameDensity(i))), "densityMap").densityMap;
    H = load(fullfile(DIR_ORIENT,  sprintf("s%02d_dic_t%02d.mat", wellID, frameDensity(i))), "H").H;
    [nRows, nCols] = size(D);

    %── Sample random window centres ─────────────────────────────────────────
    xValid = (1 + MARGIN):(nCols - MARGIN);
    yValid = (1 + MARGIN):(nRows - MARGIN);

    if isempty(xValid) || isempty(yValid)
        warning("Well %d skipped: image too small for chosen margin.", wellID);
        windowTblByWell{i} = table();
        countTblByWell{i}  = table();
        continue
    end

    rng(0);   % reproducible per well
    xCenter = xValid(randi(numel(xValid), N_SAMPLE, 1));
    yCenter = yValid(randi(numel(yValid), N_SAMPLE, 1));
    xCenter = xCenter(:);
    yCenter = yCenter(:);

    left   = xCenter - W;
    right  = xCenter + W;
    top    = yCenter - W;
    bottom = yCenter + W;
    nWin   = numel(xCenter);

    %── Per-window mean H and mean D ─────────────────────────────────────────
    meanH = NaN(nWin, 1);
    meanD = NaN(nWin, 1);
    for k = 1:nWin
        meanH(k) = mean(H(top(k):bottom(k), left(k):right(k)), "all", "omitmissing");
        meanD(k) = mean(D(top(k):bottom(k), left(k):right(k)), "all", "omitmissing");
    end

    %── Static window table ───────────────────────────────────────────────────
    windowID = (1:nWin)';
    windowTblByWell{i} = table( ...
        repmat(wellID, nWin, 1), windowID, xCenter, yCenter, meanH, meanD, ...
        'VariableNames', ["well","windowID","xCenter","yCenter","meanH","meanD"]);

    %── Frame-varying cell / death counts ────────────────────────────────────
    frames    = frameDeath{i};
    countRows = cell(numel(frames), 1);

    for j = 1:numel(frames)
        frame = frames(j);

        % Load live cells (mch, frame 1) and dead-cell marker
        cellsMch = load(fullfile(DIR_SEGMENTED, sprintf("s%02d_mch_t%02d.mat", wellID, 1)), "data").data;

        if SIGNAL == "mch"
            cellsSyt = load(fullfile(DIR_SEGMENTED, sprintf("s%02d_syt_t%02d.mat", wellID, frame)), "data").data;
        else   % cy5: dead cells are rows where isDead == 1
            raw      = load(fullfile(DIR_SEGMENTED, sprintf("s%02d_cy5_t%02d.mat", wellID, frame)), "data").data;
            cellsSyt = raw(raw.isDead == 1, 1:end-1);
        end

        % Count cells and dead cells falling inside each window
        nCell = zeros(nWin, 1);
        nDead = zeros(nWin, 1);
        for k = 1:nWin
            nCell(k) = nnz(cellsMch.x >= left(k)-0.5 & cellsMch.x <= right(k)+0.5 & ...
                           cellsMch.y >= top(k) -0.5 & cellsMch.y <= bottom(k)+0.5);
            nDead(k) = nnz(cellsSyt.x >= left(k)-0.5 & cellsSyt.x <= right(k)+0.5 & ...
                           cellsSyt.y >= top(k) -0.5 & cellsSyt.y <= bottom(k)+0.5);
        end

        pDead          = nDead ./ nCell;
        pDead(nCell==0) = NaN;

        countRows{j} = table( ...
            repmat(wellID, nWin, 1), repmat(frame, nWin, 1), windowID, ...
            nDead, nCell, pDead, ...
            'VariableNames', ["well","frame","windowID","nDead","nCell","pDead"]);
    end

    countTblByWell{i} = vertcat(countRows{:});
end
toc

%% ── Assemble and save output ─────────────────────────────────────────────────
windowTbl = vertcat(windowTblByWell{:});
countTbl  = vertcat(countTblByWell{:});
fullTbl   = innerjoin(countTbl, windowTbl, 'Keys', ["well","windowID"]);

if SAVE_OUTPUT
    save(outFile, "fullTbl");
end