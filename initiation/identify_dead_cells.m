clc; clear; close all;

%% ── Configuration ────────────────────────────────────────────────────────────
NEIGHBOR_RADIUS  = 120;
NBR_VAR          = "nNbr_" + NEIGHBOR_RADIUS;
THRESH_DEATH_SYT = 2400;
THRESH_DEATH_CY5 = 900;
FNAME_PATTERN    = 's(?<pos>\d+)_(?<channel>[A-Za-z0-9]+)_t(?<frame>\d+)\.tif$';

DIR_SEGMENTED = "../img/2_segmented_nuc";
DIR_BASIC     = "../img/1_BaSiC";

%% ── Collect unprocessed TIF files ────────────────────────────────────────────
listing  = dir(fullfile(DIR_SEGMENTED, "*.tif"));
allPaths = fullfile(string({listing.folder})', string({listing.name})');
allNames = string({listing.name})';

% % Skip files that already have a corresponding .mat
% isProcessed = arrayfun(@(p) isfile(replace(p, ".tif", ".mat")), allPaths);
% allPaths(isProcessed) = [];
% allNames(isProcessed) = [];

% Keep only files whose names match the expected pattern
tokens   = regexp(allNames, FNAME_PATTERN, "names");
isParsed = ~cellfun(@isempty, tokens);
paths    = allPaths(isParsed);
tokens   = tokens(isParsed);

channels = string(cellfun(@(t) t.channel, tokens, 'UniformOutput', false))';

%% ── Main processing loop ─────────────────────────────────────────────────────
t0 = tic;
for i = 1:numel(paths)

    fname = paths(i);
    seg   = imread(fname);
    tok   = tokens{i};
    ch    = channels(i);

    switch ch

        case "mch"
            data = processSegmentedFrame(seg, NEIGHBOR_RADIUS, NBR_VAR);

        case "syt"
            slice     = str2double(tok.frame);
            basicPath = fullfile(DIR_BASIC, "s"+tok.pos+"_"+tok.channel+".tif");
            V = tiffreadVolume(basicPath, 'PixelRegion', {[1 Inf],[1 Inf],[slice slice]});
            if ndims(V) == 3, V = V(:,:,1); end

            data = processSegmentedFrame(seg, NEIGHBOR_RADIUS, NBR_VAR, V, THRESH_DEATH_SYT);

        case "cy5"
            slice     = str2double(tok.frame);
            basicPath = fullfile(DIR_BASIC, "s"+tok.pos+"_"+tok.channel+".tif");
            V = tiffreadVolume(basicPath, 'PixelRegion', {[1 Inf],[1 Inf],[slice slice]});
            if ndims(V) == 3, V = V(:,:,1); end

            data = processSegmentedFrame(seg, NEIGHBOR_RADIUS, NBR_VAR, V, THRESH_DEATH_CY5);
            
        otherwise
            continue
    end

    % save(replace(fname, ".tif", ".mat"), "data");
    progress_bar(i, numel(paths), t0);
end

%% ════════════════════════════════════════════════════════════════════════════
%                            LOCAL FUNCTIONS
%% ════════════════════════════════════════════════════════════════════════════

function data = processSegmentedFrame(seg, radius, nbrVar, intensityImg, intensityThresh)
% PROCESSSEGMENTEDFRAME  Extract per-cell features from a labelled mask.
%
%   Required
%     seg      : labelled segmentation image
%     radius   : neighbour-search radius in pixels
%     nbrVar   : name for the neighbour-count column (e.g. "nNbr_120")
%
%   Optional (both must be supplied together)
%     intensityImg    : grayscale image for MeanIntensity measurement
%     intensityThresh : minimum MeanIntensity to retain a cell

    useIntensity = nargin >= 5;

    %── Compute region properties ────────────────────────────────────────────
    if useIntensity
        stats = regionprops(seg, intensityImg, ["EquivDiameter","MeanIntensity","Centroid"]);
    else
        stats = regionprops(seg, ["EquivDiameter","Centroid"]);
    end

    if isempty(stats)
        data = emptyResultTable(nbrVar);
        return
    end

    %── Filter cells by diameter (μ ± 2σ) ───────────────────────────────────
    diameters   = [stats.EquivDiameter];
    inDiamRange = diameters > (mean(diameters) - 2*std(diameters)) & ...
                  diameters < (mean(diameters) + 2*std(diameters));

    %── Optionally filter by fluorescence intensity ──────────────────────────
    if useIntensity
        isValid = inDiamRange & [stats.MeanIntensity] > intensityThresh;
    else
        isValid = inDiamRange;
    end

    if ~any(isValid)
        data = emptyResultTable(nbrVar);
        return
    end

    %── Count neighbours within radius ──────────────────────────────────────
    centroids = vertcat(stats(isValid).Centroid);   % N×2 [x, y]

    if size(centroids, 1) == 1
        nNeighbors = 1;   % single cell: counts itself, matching rangesearch behaviour
    else
        nNeighbors = cellfun('length', rangesearch(centroids, centroids, radius));
    end

    data = array2table([centroids, nNeighbors(:)], ...
        'VariableNames', ["x","y",nbrVar]);
end

function T = emptyResultTable(nbrVar)
    T = table(zeros(0,1), zeros(0,1), zeros(0,1), ...
              'VariableNames', ["x","y",nbrVar]);
end