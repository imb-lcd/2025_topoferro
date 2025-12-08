%
%   autocorr_selforg.m
%   quantify the autocorrelation for orientation and density
%
clearvars
clc
close all

addpath("D:/Spatiotemporal_analysis/code/util/");

ROOTDIR = "E:/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATION
prefix = ["diffday-1118-0.5d" "diffday-1118-1d" "diffday-1118-2d" "diffday-1118-2.5d"];

welllist = {[6], [6], [6], [6]};
n_wells = 24;

% select time points
% framelist = {[12], [7 19 30], [2 14 26 38], [9 21]}; % run for 10 time points
framelist = {[12], [19], [2 26], [9]};
nframes = 5;

% load file paths
sp = ["phi", "seg", "den"];
[ch, subpath] = load_file_paths();

%% load orientations phi and density
phi_all = cell(nframes, 1);
den_all = cell(nframes, 1);

cnt = 1;
for p = 1:length(prefix)
    for w = 1:length(welllist{p})
        well = sprintf("%02d", welllist{p}(w));

        PATH = ROOTDIR + "/wave_" + prefix(p) + "/" + "Well" + well + "/";
        phi_path = PATH + "/" + subpath("phi") + "/";
        den_path = PATH + "/" + subpath("den") + "/";
        
        for fr = framelist{p}
            tic
            fprintf("prefix " + prefix(p) + " well %02d, fr %02d\n", welllist{p}(w), fr);
           
            phi_file = phi_path + prefix(p)+"_s" + well + "t"+sprintf("%02d",fr)+ch("phi")+"_ORG_phi.mat";
            load(phi_file);

            den_file = den_path + prefix(p)+"_s" + well + "t"+sprintf("%02d",fr)+ch("den")+"_ORG_den120.mat";
            load(den_file);

            phi_all{cnt} = phi;
            den_all{cnt} = den;
            cnt = cnt+1;
        end
    end
end

%% calculate autocorrelation of nematic correlation

nframes = numel(phi_all);
[xi, radial, r] = compute_spatial_autocorr(phi_all, 'nematic');

%% plot the autocorrealtion for nematics

mycol = parula(nframes);

figure; 
hold on
for t = 1:nframes
    rad = radial{t};
    rad = rad(1:length(r));
    plot(r, rad, '-', 'LineWidth', 1.5, 'Color', mycol(t, :));
end
hold off
um = [0 100 200 300 400 500 600 700];
pixels = um/1.272;

xticks(pixels)
xticklabels(string(um))
xlim([0 um(end)/1.272])

xlabel('Distance (um)');
ylabel('Autocorrelation');

%% plot correlation length for nematics
figure;
plot(1:nframes-1, xi_t(1:nframes-1), 'r-o', 'LineWidth',2);

xticks([1:nframes])
xticklabels({'12h', '24h' '36h' '48h'})

yticks(pixels)
yticklabels(string(um))
ylim([pixels(1) pixels(end)])

xlabel('Time');
ylabel('Correlation Length \xi (um)');

%% calculate autocorrelation of nematic correlation

nframes = numel(den_all);
[xi, radial, r] = compute_spatial_autocorr(den_all, 'density');

%% plot the autocorrealtion for density

mycol = parula(nframes);

figure; 
hold on
for t = 1:nframes
    rad = radial{t};
    rad = rad(1:length(r));
    plot(r, rad, '-', 'LineWidth', 1.5, 'Color', mycol(t, :));
end
hold off
um = [0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300];
% um = [0 150 300 450 600];
pixels = um/1.272;
xticks(pixels)
xticklabels(string(um))
xlim([0 um(end)/1.272])
ylim([0 1])

xlabel('Distance (um)');
ylabel('Autocorrelation');

%% plot correlation length for density

um = [0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300];
pixels = um/1.272;
figure;
plot(1:nframes-1, xi(1:nframes-1), 'r-o', 'LineWidth',2);

xticks([1:nframes-1]);
xticklabels({'12h', '24h' '36h' '48h'})

yticks(pixels)
yticklabels(string(um))
ylim([pixels(1) pixels(end)])

xlabel('Time');
ylabel('Correlation Length \xi (um)');

%% FUNCTION
function [xi_t, all_rad, r] = compute_spatial_autocorr(field, type, threshold)
    if ~exist('threshold', 'var')
        threshold = exp(-1);
    end

    nframes = numel(field);

    xi_t    = nan(nframes, 1);
    all_rad = cell(nframes, 1);
    r       = [];

    for t = 1:nframes
        disp(t)
        % get the field for this frame
        A = field{t};

        % convert to 'f' depending on type
        switch lower(type)
            case 'nematic'
                phi = A;
                f = exp(2i * phi);
            
            case 'density'
                den = A;
                den(isnan(den)) = 0;
                denMean = mean(den(:), 'omitnan');
                f = den - denMean;
            otherwise
                error('Unknown type. Only nematics or density is allowed');
        end

        [M, N] = size(f);

        % Compute autocorrelation via FFT
        F = fft2(f);            % FFT of field
        powerSpec = abs(F).^2;  % power spectrum
        acf = ifft2(powerSpec); % complex number
        acf = fftshift(acf);    % move zero-lag to center

        % Edge normalization
        % Normalization factor: numebr of overlapping pixels at each lag
        w = ones(M, N);
        W = fft2(w);
        normFactor = ifft2(abs(W).^2);
        normFactor = fftshift(normFactor);

        normFactor(normFactor == 0) = NaN; % avoid dividion by zero

        % Normalized autocorrelation
        acf_norm = acf ./ normFactor;
        acf_norm = real(acf_norm); % take the real part of the acf

        % Make center = 1
        [rows, cols] = size(acf_norm);
        centerY = floor(rows/2) + 1;
        centerX = floor(cols/2) + 1;

        centerVal = acf_norm(centerY, centerX);
        acf_norm  = acf_norm / centerVal;

        % Radial averaging
        [X, Y] = meshgrid(1:cols, 1:rows);
        R = sqrt( (X - centerX).^2 + (Y - centerY).^2 ); % distance

        % Define radial bins with 1-pixel width
        maxR = floor(min(centerX, centerY));
        rEdges = 0:1:maxR;
        r_centers = 0.5 * (rEdges(1:end-1) + rEdges(2:end)); % bin centers
        radProfile = nan(length(r_centers), 1);

        for i = 1:length(r_centers)
            mask = (R >= rEdges(i)) & (R < rEdges(i+1));
            vals = acf_norm(mask);
            radProfile(i) = mean(vals, 'omitnan');
        end

        all_rad{t} = radProfile;
        if isempty(r)
            r = r_centers;
        end

        % Extract correlation length at threshold
        idx = find(radProfile <= threshold, 1, 'first');
        if ~isempty(idx)
            xi_t(t) = r(idx);
        else
            xi_t(t) = NaN;
        end
    end
end
