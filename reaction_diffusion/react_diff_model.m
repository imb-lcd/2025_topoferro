% 
% mathematical model for simulating wave propagation using different
% distirbution of PUFA in the cells
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
addpath(ROOTDIR + "/code");
addpath(ROOTDIR + "/code/util");
ROOTDIR = "E:/reaction_diffusion/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATIONS
showfigure  = 1;
% savefile = 0;

scale = 10;
ori_sz = 5000;
sz = 1;

nuc_pos_exist = 1;

%% desginate the shape and distrbution of a cell

% % most polarized
cell_design1 = [1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1];     % 1 x 16
cell_design2 = [1 1 0 0 0 0 1 1; 1 1 0 0 0 0 1 1];    % 2 x 8
cell_design3 = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]; % 4 X 4

% cell_design1(:) = 1;
% cell_design2(:) = 1;
% cell_design3 = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]; % 4 X 4

% cell_design1 = cell_design3;
% cell_design2 = cell_design3;
% cell_design3 = cell_design3;

cell_design_map = {cell_design1, cell_design2, cell_design3};
cell_pufa_px = [nnz(cell_design1(:) ~= 0) nnz(cell_design2(:) ~= 0) nnz(cell_design3(:) ~= 0)];

%% Load wave information
wv = 1;

wave_info_file = ROOTDIR + "/simulation_wave_info.txt";
wave_info = readtable(wave_info_file, 'Delimiter', '\t');

crop_corner = [wave_info.ylim1(wv) wave_info.xlim1(wv)];
crop_sz = wave_info.sz(wv);

bo_file = ROOTDIR + wave_info.bo_file(wv);
phi_file = ROOTDIR + string(wave_info.phi_file(wv));
coh_file = ROOTDIR + wave_info.coh_file(wv);
den_file = ROOTDIR + wave_info.den_file(wv);
seg_file = ROOTDIR + wave_info.seg_file(wv);
nuc_pos_file = ROOTDIR + wave_info.nuc_pos_file(wv);

ini_cell_id = wave_info.ini_cell_id(wv);

imsize = fix([crop_sz/scale crop_sz/scale]);

bo_sel = str2num(wave_info.bo_sel{wv});

%% load borders of ground truth wave
all_bo = struct2array(load(bo_file));

bo_ref = length(all_bo);
bo_cnt = length(all_bo);

bo = cell(bo_cnt, 1);
bo_mask = cell(bo_cnt,1);

% figure
% imshow(zeros(171, 171)+255)
% hold on
for i = 1:bo_cnt
    curr_bo = all_bo{i};

    % create mask of the boundary and scale
    bo_mask{i} = poly2mask(curr_bo(:, 2), curr_bo(:, 1), ori_sz, ori_sz);
    bo_mask{i} = imcrop(bo_mask{i}, [fliplr(crop_corner) crop_sz crop_sz]);
    bo_mask{i} = imresize(bo_mask{i}, [fix(crop_sz/scale) fix(crop_sz/scale)], "nearest");
    bo_mask{i} = bo_mask{i}(1+sz:end-sz, 1+sz:end-sz);

    % scale the boundary
    curr_bo(:, 2) = (curr_bo(:, 2)-crop_corner(2))/scale;
    curr_bo(:, 1) = (curr_bo(:, 1)-crop_corner(1))/scale;

    bo{i} = curr_bo;
    % plot(curr_bo(:, 2), curr_bo(:, 1), 'k-')
end
% hold off
% axis ij off

%% load phi and coh

% load phi, coh, and density
phi = struct2array(load(phi_file));
coh = struct2array(load(coh_file));
den = struct2array(load(den_file));

% crop to the specific wave
phi = imcrop(phi, [fliplr(crop_corner) crop_sz crop_sz]);
coh = imcrop(coh, [fliplr(crop_corner) crop_sz crop_sz]);
den = den(1+sz:end-sz, 1+sz:end-sz);
den = imcrop(den, [fliplr(crop_corner) crop_sz crop_sz]);
den = imgaussfilt(den, 10);
ori_den = den;

% resize the image to fit the scale
phi = imresize(phi, [fix(crop_sz/scale) fix(crop_sz/scale)], "nearest");
coh = imresize(coh, [fix(crop_sz/scale) fix(crop_sz/scale)], "nearest");
den = imresize(den, [fix(crop_sz/scale) fix(crop_sz/scale)], "nearest");

% show the images
if showfigure
    figure
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
    
    nexttile
    imshow(phi, [-pi/2 pi/2], 'Border', 'tight', colormap=hsv);
    hold on
    plot(bo{bo_ref}(:,2), bo{bo_ref}(:,1), '-', 'Color', 'k', 'LineWidth', 3)
    hold off
    
    nexttile
    imshow(coh, 'Border', 'tight', colormap=parula)
    hold on
    plot(bo{bo_ref}(:,2), bo{bo_ref}(:,1), '-', 'Color', 'w', 'LineWidth', 3)
    hold off
       
    nexttile
    lwr = 20;
    upr = 140;
    imshow(den, [lwr upr], 'Border', 'tight', colormap=jet)
    hold on
    plot(bo{bo_ref}(:,2), bo{bo_ref}(:,1), '-', 'Color', 'k', 'LineWidth', 3)
    hold off
end

% calculate and show the nematic directors
[X, Y, U, V] = set_alignment(phi, imsize, showfigure, coh, bo{bo_ref});


%% load position of the cells

if nuc_pos_exist
    load(ROOTDIR + "\sim_wv" + wv + "_optim\sim_wv" + wv + "_nuc_pos.mat");
else
    % load segmented nucle file, and crop to the wave
    seg = imread(seg_file);
    seg = seg(1+sz:end-sz, 1+sz:end-sz);
    seg = imcrop(seg, [fliplr(crop_corner) crop_sz crop_sz]);
    
    % obtain centroid of each nucleus
    seg_pos = struct2table(regionprops(seg, 'Centroid'));
    seg_pos = rmmissing(seg_pos);
    
    nuc_pos = seg_pos.Centroid;
    nuc_pos = nuc_pos/scale; % scale the positions to scale;
    
    % show in the images
    if showfigure
        figure
        nexttile
        imshow(den, [lwr upr], 'Border', 'tight', colormap=jet)
        hold on
        scatter(nuc_pos(:,1), nuc_pos(:,2), 20, 'k', 'filled');
        hold off
        axis square off;
    end
    
    % adding more nucleus based on density to infer unlabeled nucleus
    rng(130);
    den_pos = ori_den;
    nuc_idx = sub2ind(size(den_pos), fix(seg_pos.Centroid(:,2)), fix(seg_pos.Centroid(:,1)));
    den_pos(nuc_idx) = 0;
    % cell_cnt = round(height(seg_pos) / 0.40 * 0.6);
    cell_cnt = round(height(seg_pos) / 0.40);
    den_pos(isnan(den_pos)) = 0;
    prob_weights = den_pos(:) / sum(den_pos(:));
    idx = randsample(length(prob_weights), fix(cell_cnt)-height(seg_pos), true, prob_weights);
    [row, col] = ind2sub(size(den_pos), idx);
    row = row/scale;
    col = col/scale;

    if showfigure
        nexttile
        imshow(den, [lwr upr], 'Border', 'tight', colormap=jet)
        hold on
        scatter(nuc_pos(:, 1), nuc_pos(:, 2), 20, 'k', 'filled')
        scatter(col, row, 20, 'k');
        hold off
        axis square off;
    end

    nuc_pos = [nuc_pos; col row];

    % mkdir(ROOTDIR + "\sim_wv" + wv + "_optim\");
    % nuc_pos_file = ROOTDIR + "\sim_wv" + wv + "_optim\sim_wv" + wv + "_nuc_pos.mat";
    % save(nuc_pos_file, 'nuc_pos');
end


%% layout the cells based on segment positions
cell_cnt = length(nuc_pos);

% for presentation
all_cell_id = zeros(imsize);
all_lipid_dist = zeros(imsize);

% for storing the positions and lipid distirbutions of the cell
lipid_dist = cell(cell_cnt, 1);
cell_id = cell(cell_cnt, 1);
cell_design = zeros(cell_cnt, 1);

% set the thresholds for coherency
coh_cut1 = 0.075; 
coh_cut2 = 0.25;

% loop through each cell and assign respective cell design to the position
for i = 1:cell_cnt
    cell_pos = false(imsize);
    lipid_pos = false(imsize);

    curr_x = fix(nuc_pos(i, 1));
    curr_x(curr_x==0) = 1;
    curr_x(curr_x>size(phi,1)) = size(phi, 1);
    curr_y = fix(nuc_pos(i, 2));
    curr_y(curr_y==0) = 1;
    curr_y(curr_y>size(phi,2)) = size(phi, 2);

    angle = phi(curr_y, curr_x);

    % calculate direction for the long side of the cell
    dir_x = cos(angle);
    dir_y = sin(angle);

    % calculate direction for the short side (perpendicular) of the cell
    perp_dir_x = -sin(angle);
    perp_dir_y = cos(angle);

    % set cell design based on coherency
    cell_design(i) = 1 * (coh(curr_y, curr_x) > coh_cut2) + ...
                     2 * (coh(curr_y, curr_x) <= coh_cut2 & coh(curr_y, curr_x) > coh_cut1) + ...
                     3 * (coh(curr_y, curr_x) <= coh_cut1);

    curr_design = cell_design_map{cell_design(i)};

    % lay out the cell design and pufa positions
    curr_design_sz = size(curr_design);
    for dy = 1:curr_design_sz(1)
        for dx = 1:curr_design_sz(2)
            % Calculate the position for the lement in the current cell
            % nucleus_position + horizontal offset + vertical offset
            x = round(nuc_pos(i, 1) + (dx-ceil(curr_design_sz(2)/2))*dir_x + (dy-ceil(curr_design_sz(1)/2))*perp_dir_x);
            y = round(nuc_pos(i, 2) + (dx-ceil(curr_design_sz(2)/2))*dir_y + (dy-ceil(curr_design_sz(1)/2))*perp_dir_y)-1;

            % Check bounds and assigns the cell if the cell position is
            % within the matrix
            if x >= 1 && x <= imsize(2) && y >= 1 && y <= imsize(1)
                cell_pos(y, x) = true;
                lipid_pos(y, x) = curr_design(dy, dx);

                all_cell_id(y, x) = i;
                all_lipid_dist(y, x) = all_lipid_dist(y, x) + curr_design(dy, dx);
            end
        end
    end

    % Store the positions
    cell_id{i} = imclose(cell_pos, strel('disk', 1));
    cell_id{i} = cell_id{i}(1+sz:end-sz, 1+sz:end-sz);

    lipid_dist{i} = lipid_pos(1+sz:end-sz, 1+sz:end-sz);

end

% calculate degree of overlap for each cell
overlap_cells = ones(imsize-2);
overlap_per_cell = cell(cell_cnt, 1);

for i = 1:cell_cnt
    overlap_cells = overlap_cells + cell_id{i}; %(i, :, :);
end

for i = 1:cell_cnt
    curr_cell_pos = cell_id{i};
    curr_lipid = lipid_dist{i};

    overlap_per_cell{i} = overlap_cells(curr_cell_pos & curr_lipid);
end

if showfigure
    figure
    % nexttile
    imshow(all_cell_id, [1 cell_cnt], 'Border', 'tight', colormap=[0 0 0; lines]);
    axis square off
    hold on
    plot(bo{bo_ref}(:,2), bo{bo_ref}(:,1), '-', 'Color', 'w', 'LineWidth', 3)
    hold off
    
    % nexttile
    % imshow(all_lipid_dist, [0 1], colormap=flipud(gray))
    % axis square off
    % hold on
    % plot(bo{bo_ref}(:,2), bo{bo_ref}(:,1), '-', 'Color', 'y', 'LineWidth', 3)
    % hold off
end

%% parameter of model
d_lwr      = 0.01;
c_lwr      = 0.01;

D          = 0.16;
a          = 1.12;
d_upr      = 0.82;
c_upr      = 0.22;
degrade    = 0.19;
death_time = 17;  

%% run the simulation
showfigure = 1;
savefile = 0;

%%% five different conditions:
% polar_sens, homo_sens, polar_insens, homo_insens, nonelong

run_condition = "polar_sens";

L = [1 1 0.75];

disp(run_condition);

if contains(run_condition, "death")
    tstart = 1;
    tend = 1500;
    t_interval = 10;
else
    tstart = 1;
    tend = 500;
    t_interval = 100;
end

disp(["pufa_weight:" L])

if savefile
    outpath = ROOTDIR + "/sim_wv" + wv + "_optim/";
    mkdir(outpath);
    
    date_str = string(datetime('today', 'Format', 'MMdd'));

    data_outfile = outpath + "sim_wv" + wv + "_" + date_str + "_" + run_condition + ".txt";
    if exist(data_outfile)
        delete(data_outfile);
    end

    if showfigure
        outpath_figure = outpath + date_str + "_" + run_condition + "_figure/";
        fig_file_prefix = "sim_wv" + wv + "_" + date_str + "_" + run_condition + "_t";
        mkdir(outpath_figure);
    end

    % write parameters to file
    data = ["a" a; ...
            "D" D; ...
            "d_upr" d_upr; ...
            "c_upr" c_upr; ...
            "degrade" degrade; ...
            "death" death_time];
    writematrix(data, data_outfile, 'WriteMode', 'append', 'Delimiter', '\t');
end

cnt = 1;

max_jac   = zeros(length(bo_sel), 1);
plot_jac  = cell(length(bo_sel), 1);
max_jac_t = nan(length(bo_sel), 1);
max_jac_b = nan(length(bo_sel), 1);
max_circ  = zeros(length(bo_sel), 1);

% record death
death_over_time = nan(length(tend), 1);

% simulate the model
% declare ros environment, initiator cells which is already at higher steady state
r = zeros(imsize) + padarray(lipid_dist{ini_cell_id}, [1 1], 0, 'both');

cell_ros = zeros(cell_cnt, 1);
cell_hss = false(cell_cnt, 1);

% set ROS and higher stready state in the initiator cell
cell_ros(ini_cell_id) = 1;
cell_hss(ini_cell_id) = 1;

% set up death count
death_cnt = ones(cell_cnt, 1) * death_time;

% set up the cell thresholds from density and coherency
d = den(1+sz:end-sz, 1+sz:end-sz);
d_upr_std = mean(d(:)) + 2*std(d(:));
d(d>d_upr_std) = d_upr_std;
d = rescale(d, d_lwr, d_upr);

c = coh(1+sz:end-sz, 1+sz:end-sz);
c = rescale(c, c_lwr, c_upr);

% ==== per-cell thresholds ====
thres = nan(cell_cnt, 1);

for i = 1:cell_cnt
    mask = cell_id{i} > 0;
    if any(mask, "all")
        thres(i) = max(d(mask)) .* max(c(mask));
    else
        thres(i) = NaN;  % unused, cell has no pixels
    end
end
if contains(run_condition, "nonelong")
    thres = ones(size(thres)) .* 0.01;
elseif contains(run_condition, "insens")
    thres = ones(size(thres)) .* 0.06;
end

% ==== precompute per-cell linear indices and inv(overlap) ====
cell_lipid_idx       = cell(cell_cnt, 1);
inv_overlap_per_cell = cell(cell_cnt, 1);

for i = 1:cell_cnt
    curr_lipid = lipid_dist{i};

    if ~isempty(curr_lipid) && any(curr_lipid(:))
        % linear indices of PUFA pixels for cell i
        cell_lipid_idx{i}       = find(curr_lipid);

        inv_overlap_per_cell{i} = 1 ./ overlap_per_cell{i};  % same length
    else
        cell_lipid_idx{i}       = [];
        inv_overlap_per_cell{i} = [];
    end
end

% ====
if savefile
    data = ["threshold", "CohDen"]; %["threshold", curr_thres];
    writematrix(data, data_outfile, 'WriteMode', 'append', 'Delimiter', '\t');
end

% label current run
curr_run = "a " + a + ...
           "; D " + D + ...
           "; den " + d_lwr + "-" + d_upr + ...
           "; coh " + c_lwr + "-" + c_upr + ...
           "; degrade " + degrade + ...
           "; death " + death_time;
disp(curr_run);  

% start simulation
if showfigure
    f = figure;

    plt_cnt = 1;
    plt_row = floor(tend/t_interval)+1;

    bo_curr = bo{1};

    nexttile
    plot_simulation_dead(0, r(1+sz:end-sz, 1+sz:end-sz), cell_cnt, cell_id, cell_hss, X, Y, U, V, bo_curr, 0, "", plt_cnt, plt_row);
    sgtitle(curr_run);
    
    plt_cnt = plt_cnt + 1;

    if savefile
        figure_outfile = outpath_figure + "/" + fig_file_prefix + "0000.jpg";
        exportgraphics(gcf, figure_outfile);
        close all
    end
end

% %%
% tstart = tend+1;
% tend = 1000;

lap_kernel = [0 1 0; 1 -4 1; 0 1 0] / 4;
tic
for t = tstart:tend
    R = r(2:end-1, 2:end-1);

    %%% diffusion term
    D_iso = D * conv2(r, lap_kernel, 'valid');

    %%% reaction term
    % count down for death
    death_cnt(cell_hss == 1) = max(0, death_cnt(cell_hss == 1)-1);

    % quantify reaction term for each cell
    for i = 1:cell_cnt
        if death_cnt(i) == 0
            continue;
        end

        idx = cell_lipid_idx{i};
        if isempty(idx)
            continue;
        end

        % calculate ROS in the current cell
        react = cell_ros(i) * (1 - cell_ros(i)) * (cell_ros(i) - thres(i));

        contrib = sum(D_iso(idx) * L(cell_design(i)) .* inv_overlap_per_cell{i});
        
        cell_ros(i) = cell_ros(i) + a * react + contrib;
        
                
        % if cell is in higher steady state, then add the
        % ROS of the cell to the environment
        if cell_hss(i) == 1 || cell_ros(i) > 0.9 
            cell_hss(i) = 1;

            R(lipid_dist{i}) = R(lipid_dist{i}) + cell_ros(i)/cell_pufa_px(cell_design(i));
        end
    end

    % incoprorate diffusion term in the environment with
    % degradation term
    R = max(0, R + D_iso - R * degrade);

    if showfigure && mod(t, t_interval) == 0
        % figure
        nexttile
        % disp(t);
        if t == 100
            bo_curr = bo{2};
        elseif t == 200
            bo_curr = bo{4};
        elseif t == 300
            bo_curr = bo{7};
        elseif t == 400
            bo_curr = bo{10};
        elseif t == 500
            bo_curr = bo{13};
        end
        plot_simulation_dead(t, R, cell_cnt, cell_id, cell_hss, X, Y, U, V, bo_curr, 0, "", plt_cnt, plt_row);
        plt_cnt = plt_cnt + 1;
        % title(t)

        if savefile
            figure_outfile = outpath_figure + "/" + fig_file_prefix + sprintf("%04d", t) + ".jpg";

            exportgraphics(gcf, figure_outfile);
            close all
        end
    end

    r = padarray(R, [1 1], 0, 'both');

    % if ismember(t, t_sel)
    %     % find the jaccard index at tsel
    %     idx = find(t_sel == t);
    % 
    %     b = bo_sel(idx);
    %     [jac, plot_dead] = calculate_jaccard_similarity(bo_mask{b}, R, cell_cnt, cell_id, cell_hss);
    % 
    %     max_jac(idx)   = jac;
    %     plot_jac{idx}  = plot_dead;
    %     max_jac_t(idx) = t;
    %     max_jac_b(idx) = b;
    % end

    % % for finding the best jaccard index at boundary
    % if mod(t, t_interval) == 0
    %     if ~contains(run_condition, "death")
    %         for b_idx = 1:length(bo_sel)
    %             b = bo_sel(b_idx);
    % 
    %         % for b_idx = 1:length(bo_cnt)
    %         %     b = bo_cnt;
    % 
    %             [jac, plot_dead] = calculate_jaccard_similarity(bo_mask{b}, R, cell_cnt, cell_id, cell_hss);
    % 
    %             bw = imbinarize(plot_dead, 0.1);
    %             stats = regionprops(bw, 'Circularity');
    % 
    %             if jac > max_jac(b_idx)
    %                 max_jac(b_idx)   = jac;
    %                 plot_jac{b_idx}  = plot_dead;
    %                 max_jac_t(b_idx) = t;
    %                 max_jac_b(b_idx) = b;
    % 
    %                 max_circ(b_idx) = stats.Circularity;
    %             end
    %         end
    %     else
    %         death_over_time(t) = sum(death_cnt==0);
    %     end
    % end
    % close all;
end
toc

if savefile
    data = [a D d_lwr d_upr c_lwr c_upr degrade death_time]; 

    writematrix(data, data_outfile, 'WriteMode', 'append', 'Delimiter', '\t');

    for ttt = 1:length(bo_sel)
        data = ["idx" ttt "max_jac_b" max_jac_b(ttt) "max_jac_t" max_jac_t(ttt) "max_jac" max_jac(ttt) "circ_at_max_jac" max_circ(ttt)];
        writematrix(data, data_outfile, 'WriteMode', 'append', 'Delimiter', '\t');
    end

    if contains(run_condition, "death")
        death_over_time = death_over_time(death_over_time ~= 0);
        death_over_time(1) = death_over_time(2);

        save(outpath + "/" + date_str + "_" + run_condition + "_death_over_time.mat", 'death_over_time');
    end
end
% close all


%% plot cell death and isotropy over time

polar_sens          = struct2array(load('E:\reaction_diffusion\sim_wv1_optim\1202_polar_sens_c2l_death_death_over_time.mat'));
homo_sens           = struct2array(load('E:\reaction_diffusion\sim_wv1_optim\1203_homo_sens_weightnew5_c2l_death_death_over_time.mat'));
polar_insen         = struct2array(load('E:\reaction_diffusion\sim_wv1_optim\1202_polar_insens_c2l_death_death_over_time.mat'));
homo_insen          = struct2array(load('E:\reaction_diffusion\sim_wv1_optim\1202_homo_insens_weightnew4_c2l_death_death_over_time.mat'));
homo_insen_nonelong = struct2array(load('E:\reaction_diffusion\sim_wv1_optim\1202_homo_insens_nonelong_c2l_death_death_over_time.mat'));

mycmap = gray(6);

tend = 1500;
total_cell = cell_cnt;

figure
hold on
plot(1:tend/10, polar_sens/total_cell,          'Color', mycmap(1,:), 'LineWidth', 3)
plot(1:tend/10, homo_sens/total_cell,           'Color', mycmap(2,:), 'LineWidth', 3)
plot(1:tend/10, polar_insen/total_cell,         'Color', mycmap(3,:), 'LineWidth', 3)
plot(1:tend/10, homo_insen/total_cell,          'Color', mycmap(4,:), 'LineWidth', 3)
plot(1:tend/10, homo_insen_nonelong/total_cell, 'Color', mycmap(5,:), 'LineWidth', 3)
hold off
xticklabels({0, 500, 1000, 1500})
ylim([0 1.1])
ylabel("Cell death (%)")
xlabel("Simulation time (t)")

%% FUNCTION
function [jac, plot_dead] = calculate_jaccard_similarity(bo_mask, R, cell_cnt, cell_id, cell_hss)
    plot_dead = R;
    for i = 1:cell_cnt
        if cell_hss(i)
            curr_cell_pos = cell_id{i};
            plot_dead(curr_cell_pos) = 1;
        end
    end

    plot_dead = imclose(plot_dead, strel('disk', 2));
    plot_dead = imgaussfilt(plot_dead, 1);

    sim_mask = imbinarize(double(plot_dead), 0.001);
    jac = jaccard(bo_mask, sim_mask);
end

function [plot_r_dead] = plot_simulation_dead(t, R, cell_cnt, cell_id, cell_hss, X, Y, U, V, bo, savefigure, outfile, plt_cnt, plt_row)
    plot_r_dead = R;
    for i = 1:cell_cnt

        curr_cell_pos = cell_id{i};

        if cell_hss(i)
            plot_r_dead(curr_cell_pos) = 1;
        end
    end
    plot_dead = plot_r_dead;
    plot_dead = imclose(plot_dead, strel('disk', 2));
    plot_dead = imgaussfilt(plot_dead, 1);

    imagesc(plot_dead)
    clim([0 0.25])
    axis square off

    title("t"+t)
    hold on
    % quiver(X, Y, U, V, 0.5, 'color', '#7F7F7F', 'ShowArrowHead', 'off', 'LineWidth', 1)
    plot(bo(:,2), bo(:,1), '-', 'Color', 'k', 'LineWidth', 2);
    hold off

    if savefigure
        exportgraphics(gcf, outfile);
    end
end

function [X, Y, U, V] = set_alignment(phi, imsize, showfigure, coh, bo)
    nem_ws = 10;
    
    phi_nem = imresize(phi, [length(phi)/nem_ws length(phi)/nem_ws], 'nearest');

    [X, Y] = meshgrid(1:pi/2:pi/2*length(phi_nem), 1:pi/2:pi/2*length(phi_nem));

    X = X - 1/2 * cos(phi_nem);
    Y = Y - 1/2 * sin(phi_nem);
    
    X = rescale(X, 1, imsize(2));
    Y = rescale(Y, 1, imsize(1));
    
    U = cos(phi_nem);
    V = sin(phi_nem);

    if showfigure
        nexttile
        imshow(zeros(size(phi))+1, [0, 1], colormap=gray)
        % imshow(coh)
        axis square off
        hold on
        quiver(X, Y, U, V, 0.5, 'color', 'k', 'ShowArrowHead', 'off', 'LineWidth', 1)
        plot(bo(:,2), bo(:,1), '-', 'Color', 'k', 'LineWidth', 3)
        hold off
    end
end