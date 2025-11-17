%
%   autocorr_selforg.m
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
%% quantify autocorrelation for phi

mycol = parula(nframes);

% Preallocate an array to store correlation lengths at each time point
xi_t = zeros(nframes, 1);

threshold = exp(-1);

figure
hold on
% Loop over each time point
for t = 1:(nframes-1)
    disp(t)
    % Get the orientation matrix for the current time
    phi = phi_all{t};  % orientations in radians
    [M, N] = size(phi);
    
    % Compute the complex order parameter (doubled angle for nematic symmetry)
    f = exp(2i * phi);
    
    % Compute the 2D autocorrelation via FFT
    F = fft2(f);
    % The autocorrelation is obtained by inverse FFT of |F|^2
    C = ifft2(abs(F).^2);
    % Shift the zero-displacement to the center of the matrix
    C = fftshift(C);
    
    % Compute the normalization factor
    % For each displacement, the number of overlapping pixels is given by
    % the convolution of a ones matrix with itself.
    w = ones(M, N);
    W = fft2(w);
    normFactor = ifft2(abs(W).^2);
    normFactor = fftshift(normFactor);
    
    % Normalize the autocorrelation function
    C_norm = C ./ normFactor;
    
    % Radially average the 2D autocorrelation
    % Create a grid of distances from the center
    [X, Y] = meshgrid(-(N-1)/2:(N-1)/2, -(M-1)/2:(M-1)/2);
    R = sqrt(X.^2 + Y.^2);
    % Bin distances to the nearest integer (in pixel units)
    R_bin = round(R);
    maxR = max(R_bin(:));
    C_radial = zeros(maxR+1,1);
    
    for r_val = 0:maxR
        mask = (R_bin == r_val);
        % Average over all pixels that fall into the same radial bin
        C_radial(r_val+1) = mean(real(C_norm(mask)));
    end
    
    % Determine the correlation length
    % The correlation length ξ is defined as the smallest r for which C_radial drops below exp(-1)
    idx = find(C_radial <= threshold, 1, 'first');
    if ~isempty(idx)
        xi_t(t) = idx - 1; % subtract one because bin indexing starts at 0
    else
        xi_t(t) = NaN; % if not found, assign NaN
    end
    
    % Optional: Plot the radially averaged correlation function for this time point
    plot(0:maxR, C_radial, '-', 'LineWidth', 1.5, 'Color', mycol(t, :));
end
hold off

um = [0 100 200 300 400 500 600 700];
pixels = um/1.272;

xticks(pixels)
xticklabels(string(um))
xlim([0 um(end)/1.272])

xlabel('Distance (um)');
ylabel('Autocorrelation');

%% plot correlation length vs time

figure;
plot(1:nframes, xi_t(1:nframes), 'r-o', 'LineWidth',2);

xticks([1:nframes])
xticklabels({'12h', '24h' '36h' '48h' '60h'})

yticks(pixels)
yticklabels(string(um))
ylim([pixels(2) pixels(end-2)])

xlabel('Time');
ylabel('Correlation Length \xi (um)');

%% quantify autocorrelation for density
L_max = 500;
dr = 5;

all_rad = cell(nframes, 1);
all_r = cell(nframes, 1);
cnt = 1;

for t = 1:nframes
    disp(t);
    den = den_all{t};

    den(isnan(den)) = 0;
    densityFluct = den - mean(den(:), 'omitnan');

    % Compute the FFT of the fluctuation image.
    F = fft2(densityFluct);

    % Compute the power spectrum (squared magnitude).
    powerSpec = abs(F).^2;

    % Compute the autocorrelation by inverse FFT (take the real part).
    acf = real(ifft2(powerSpec));

    % Shift the zero-lag term to the center of the image.
    acf = fftshift(acf);

    % Normalize the autocorrelation so that the center (zero lag) is 1.
    % ACF = acf / max(acf(:));
    ACF = acf / max(acf(:));

    % ---------------------
    % Radial Averaging:
    % ---------------------
    [rows, cols] = size(ACF);
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Determine the center of the image.
    centerX = ceil(cols / 2);
    centerY = ceil(rows / 2);
    
    % Compute distance of each pixel from the center.
    R = sqrt((X - centerX).^2 + (Y - centerY).^2);
    
    % Define radial bins (you can adjust the bin width here; we use 1 pixel here).
    maxR = min(centerX, centerY); % limit to the smallest half-dimension.
    rEdges = 0:1:maxR;
    r = rEdges(1:end-1) + 0.5;  % bin centers
    
    % Initialize the radial profile.
    radProfile = zeros(length(r), 1);
    
    % Average the ACF values in each annular ring.
    for i = 1:length(r)
        mask = (R >= rEdges(i)) & (R < rEdges(i+1));
        radProfile(i) = mean(ACF(mask));
    end

    all_rad{cnt} = radProfile;
    all_r{cnt} = r;
    cnt = cnt + 1;
end
%% plot the autocorrelation plot
mycol = parula(nframes);

figure
hold on
for i = 1:nframes
    rp = all_rad{i};
    r = all_r{i};
    plot(r, rp, 'Color', mycol(i, :), 'LineWidth', 2)
end
hold off

um = [0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300];
% um = [0 150 300 450 600];
pixels = um/1.272;


ylim([0 1])

xticks(pixels)
xticklabels(string(um))
xlim([0 um(end)/1.272])

xlabel('Distance (um)');
ylabel('Autocorrelation');

%%
um = [0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300];
% um = [0 150 300 450 600];
pixels = um/1.272;
xticks(pixels)
xticklabels(string(um))
xlim([0 um(end)/1.272])
%% calculate correlation length
xi_t = zeros(nframes,1);
threshold = exp(-1);

for i = 1:nframes
    rp = all_rad{i};
    r = all_r{i};
    idx = find(rp <= threshold, 1, 'first');
   
    if ~isempty(idx)
        xi_t(i) = idx - 1; % subtract one because bin indexing starts at 0
    else
        xi_t(i) = NaN; % if not found, assign NaN
    end
end

figure;
sm_xi_t = smoothdata(xi_t, 'movmean', 2);
plot(1:nframes, sm_xi_t(1:nframes), 'r-o', 'LineWidth',2);

xticks([1:nframes])
xticklabels({'12h', '24h' '36h' '48h' '60h'})
 
% yticks(pixels)
% yticklabels(string(um))
% ylim([pixels(4) pixels(end)])
ylim([300 1000])

xlabel('Time');
ylabel('Correlation Length \xi (um)');

