%
%   pdms_analysis_alnang_circCtrl.m
%
clearvars
clc
close all

addpath("D:/Spatiotemporal_analysis/code/util/");

ROOTDIR = "D:/Spatiotemporal_analysis/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATION
all_prefix = ["pdms-0404", "pdms-0404", "pdms-0418", "pdms-0418"];

welllist = [10 12 20 23];
fr = 1;
imsize = 5848;

n_prefix = length(all_prefix);
n_wells  = length(welllist);

%% load pdms positions
pdms_xy = [1580 1755;
           1950 1940
           1750 1830
           1910 2110];
pdms_r  = 2075/2;

%% load orientation and density
den_all  = cell(n_prefix, n_wells);
phi_all  = cell(n_prefix, n_wells);

for p = 1:n_prefix
    prefix = all_prefix(p);

    well = sprintf("%02d", welllist(p));
    disp(prefix + " " + well)

    PATH = ROOTDIR + "/wave_" + prefix + "/Well" + well + "/";
    den_path = PATH + "/c1_mCherry_density/density_reg120/";
    phi_path = PATH + "/c3_DIC_nematics/orientation_files/";
        
    den_file = den_path + prefix+"_s"+well+"t"+sprintf("%02d",fr)+"c1_ORG_den120.mat";
    load(den_file);
    den = imresize(den, [imsize imsize]);

    den_all{p} = den;

    phi_file = phi_path + prefix+"_s"+well+"t"+sprintf("%02d",fr)+"c3_ORG_phi.mat";
    load(phi_file)
    phi = imresize(phi, [imsize imsize], 'nearest');

    phi_all{p} = phi;
end

%% calculate angle diff

ctrl_ang_mean = nan(n_prefix, num_images);
ctrl_ang_std  = nan(n_prefix, num_images);
ctrl_ang_circstd = nan(n_prefix, num_images);
ctrl_ang_circvar = nan(n_prefix, num_images);

% precompute coordinate grid once
[Xgrid, Ygrid] = meshgrid(1:imsize, 1:imsize);

for p = 1:n_prefix

    prefix = all_prefix(p);

    well = sprintf("%02d", welllist(p));
    disp(prefix + " " + well)
    
    fname = ROOTDIR + "/wave_" + prefix + "/Well" + well + "/c2_cy5_adj/" + ...
        prefix + "_s"+well+"tAllc2_ORG_bw_adj.tif";

    info = imfinfo(fname);
    
    curr_phi0 = phi_all{p};
    
    prev_mask      = false(imsize);
    prev_ctrl_mask = false(imsize);

    cx = pdms_xy(p, 1) + pdms_r;
    cy = pdms_xy(p, 2) + pdms_r;

    % radial direction from initiation center to each pixel
    % Xgrid = column/x, Ygrid = row/y
    radial_angle = atan2(Ygrid - cy, Xgrid - cx);
    
    % angle difference between local cell orientation and radial direction
    % wrapped to [-pi/2, pi/2] because cell orientation is nematic
    radial_angdiff_map = mod(curr_phi0 - radial_angle + pi/2, pi) - pi/2;
    radial_angdiff_map = abs(radial_angdiff_map);

    for k = 3:num_images

        disp(k)
        bw = imread(fname, k);
      
    
        [~, Bo] = wave_boundary_gaussfilt(bw, 1, sigma, binarize_threshold);
        bo = Bo{1};
    
        mask = poly2mask(bo(:, 2), bo(:, 1), imsize, imsize);
        
        front_mask = mask & ~prev_mask;

        n_front_pix = nnz(front_mask);

        % same-area circle control
        wave_area = sum(mask(:), 'omitnan');

            r_ctrl = sqrt(wave_area / pi);

        
        ctrl_mask = ((Xgrid - cx).^2 + (Ygrid - cy).^2) <= r_ctrl^2;

        % control wavefront = newly added circular annulus
        ctrl_front_mask = ctrl_mask & ~prev_ctrl_mask;
        n_ctrl_front_pix = nnz(ctrl_front_mask);

        % circular-control angle difference
        
        ctrl_ad = radial_angdiff_map(ctrl_front_mask);
        ctrl_ad(isnan(ctrl_ad)) = [];
        
        if ~isempty(ctrl_ad)
        
            % ordinary mean and SD of wrapped angle difference
            ctrl_ang_mean(p, k) = mean(ctrl_ad, 'omitnan');
            ctrl_ang_std(p, k)  = std(ctrl_ad, 0, 'omitnan');
        
            % circular mean and circular SD, recommended for angle data
            z_ctrl = mean(exp(1i * 2 * ctrl_ad), 'omitnan');
        
            ctrl_ang_mean(p, k) = 0.5 * angle(z_ctrl);
        
            R_ctrl = abs(z_ctrl);
            R_ctrl = max(R_ctrl, eps);
        
            ctrl_ang_circstd(p, k) = 0.5 * sqrt(-2 * log(R_ctrl));
            ctrl_ang_circvar(p, k) = 1 - R_ctrl;
        
        end

        prev_mask = mask;
        prev_ctrl_mask = ctrl_mask;
    end
end


%% ctrl iso angular difference

savefile = true;

xplot = 3:num_images;

ctrl_ang_avg = abs(mean(ctrl_ang_mean(1:3, xplot), 1, 'omitnan'));
ctrl_ang_sd  = std(ctrl_ang_mean(1:3, xplot), 0, 1, 'omitnan');

sm_factor = 6;
ctrl_ang_avg = movmean(ctrl_ang_avg, sm_factor, 'Endpoints', pi/4);
ctrl_ang_sd = smoothdata(ctrl_ang_sd, 'movmean', 4);

figure
hold on


fill([xplot fliplr(xplot)], ...
     [ctrl_ang_avg - ctrl_ang_sd, fliplr(ctrl_ang_avg + ctrl_ang_sd)], ...
     'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');

plot(xplot, ctrl_ang_avg, 'b-', 'LineWidth', 2)

yline(0, '--')
xline([7 10], '--')

ylim([0 pi/2])
yticks(0:pi/4:pi/2)
yticklabels({'0', '\pi/4', '\pi/2'})

xlabel("frame")
ylabel("cell orientation - radial direction")
title("Radial angle difference: circle control")

hold off

if savefile
    outfile = "D:\Spatiotemporal_analysis\wave_0_analyses\wave_pdms\" + ...
        "pdms_angdiff_time circCtrl";
    exportgraphics(gcf, outfile + ".jpg");
    exportgraphics(gcf, outfile + ".pdf");
    savefig(outfile + ".fig");
end


%% FUNCTION
function [bdy, Bo] = wave_boundary_gaussfilt(bw, nframes, sigma, binarize_threshold, medfilt_sz)
    sm_area_thres = 3000;
    sm_hole_thres = 10000; %50000;

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
        BW = imclose(BW, strel('disk', 10)); 

        
        if run_medfilt
            
            BW = medfilt2(BW, [medfilt_sz medfilt_sz]);
        end

        bdy{i} = sparse(BW);
       
        bo = bwboundaries(BW);
        areas = cellfun(@(b) polyarea(b(:,2), b(:,1)), bo);
        
        % Find index of largest area
        [~, idx] = max(areas);
        
        % Keep only the largest boundary
        bo = bo{idx};

        bo(:, 2) = smoothdata(bo(:, 2), 'sgolay', 75)';
        bo(:, 1) = smoothdata(bo(:, 1), 'sgolay', 75)';

        Bo{i} = bo;
    end
end