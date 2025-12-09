% 
% Mathematical model for simulating wave propagation by polarized PUFA
% Toy model
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
addpath(ROOTDIR + "/code");
addpath(ROOTDIR + "/code/util");

ROOTDIR = "E:/lipid_perox_model/sim_toy_run/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% designate shape and PUFA location in a cell

cell_design1 = [1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1]; 

cell_design_map = {cell_design1};
cell_pufa_px = [nnz(cell_design1)];

L = [1];

%% set up toy environment
imsize = [200 200];
sz = 1;

rng(135);

aln_type = categorical("random");

% set environment
phi = zeros(imsize);

test_angle = 4*pi/8;
angle_pos = 125; % 50

phi(1:end, angle_pos:end) = test_angle;
pivot_type = categorical("center");

start_pos = 40;

% for propagation along and against alignment
phi(start_pos:150, 80:190) = imgaussfilt(phi(start_pos:150, 80:190), 70);
phi(start_pos:200, 50:end) = imgaussfilt(phi(start_pos:200, 50:end), 50);

% set density and coherency
den = ones(imsize) * 0.001; % high_den = 0.01; high = 0.15
den(1:end, angle_pos:end) = 0.001;
coh = ones(imsize);

figure(1)
clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
[X, Y, U, V] = set_alignment(phi, imsize, 1, 8);

figure(2)
clf
hold on
imshow(zeros(imsize(1)-2, imsize(2))-2+255, 'Border', 'tight')
hold on
plot(1:imsize(2)-1, repelem(start_pos, imsize(2)-1), 'Color', 'b', 'LineWidth', 2)

lineobj = streamline(X, Y, U, V, 1, start_pos, 'Color', 'b', 'LineWidth', 2);
hold off

%% layout cells
cell_cnt = ceil((imsize(2)*1.25)^2/length(cell_design1)); % maximum possible number of cells allow 1/4 overlap

cell_id        = cell(cell_cnt,1);    % per-cell occupancy masks
lipid_dist     = cell(cell_cnt,1);    % per-cell lipid masks
all_cell_id    = zeros(imsize);        % final “which cell” label map
all_lipid_dist = zeros(imsize);        % final “which lipid” label map

cell_shape  = cell_design1; 
cell_design = repelem(1, cell_cnt)';

[H, W] = size(cell_shape);

% center of the cell
cx0 = ceil(W/2);
cy0 = ceil(H/2);

allowed_overlap = floor(H*W/4); % at most 1/4 of the cell may lie on existing cells
margin = 14;

cell_start_pos = nan(cell_cnt, 2);

placed = 0; % how many cells have been placed
max_attempts = 1000;
for k = 1:cell_cnt
    success = false;
    attempts = 0;

    while ~success && attempts < max_attempts
        attempts = attempts + 1;

        % set initiator cells
        if k == 1
            if pivot_type == "center"
                xc = round(W/2)+3;
                yc = start_pos;
            else
                xc = 3;
                yc = start_pos;
            end
        elseif k == 2 || k == 4
            if pivot_type == "center"
                xc = round(W/2);
                yc = start_pos-k/2;
            else
                xc = 1;
                yc = start_pos-k/2;
            end
        elseif k == 3 || k == 5
            if pivot_type == "center"
                xc = round(W/2);
                yc = start_pos+floor(k/2);
            else
                xc = 1;
                yc = start_pos+floor(k/2);
            end
        else
            if aln_type == "perfect"
                [yc, xc] = find(all_cell_id==0, 1, 'first');
                if isempty(yc)
                    break;
                else
                    success = true;
                end
            elseif aln_type == "random"
                % pick a random point as corner
                if pivot_type == "center"
                    xc = randi([1+cx0-margin, imsize(2)-cx0+margin]);
                    yc = randi([1+cy0-margin, imsize(1)-cy0+margin]);
                else
                    xc = randi(imsize(2));
                    yc = randi(imsize(1));
                end
            end
        end

        cell_mask  = false(imsize);
        lipid_mask = false(imsize);

        pivot_x = xc;
        pivot_y = yc;
        prev_x = 1;
        prev_y = 1;

        for dy = 1:H
            % snake order traversal
            if mod(dy, 2) ==1
                dx_range = 1:W;
            else
                dx_range = W:-1:1;
            end

            for dx = dx_range
                if pivot_type == "center"
                    dx0 = dx - cx0;
                    dy0 = dy - cy0;

                    pivot_x = xc + dx0;
                    pivot_y = yc + dy0;
                else
                    dy0 = dy - prev_y;
                    dx0 = dx - prev_x;
                end
                
                % get the local angle
                theta = phi(max(1,min(imsize(1),round(pivot_y))), ...
                            max(1,min(imsize(2),round(pivot_x))));

                % rotate
                if pivot_type == "center"
                    rotate_x = dx0*cos(theta) - dy0*sin(theta);
                    rotate_y = dx0*sin(theta) + dy0*cos(theta);

                    xi = round(xc + rotate_x);
                    yi = round(yc + rotate_y);
                else
                    pivot_x = pivot_x + dx0*cos(theta) - dy0*sin(theta);
                    pivot_y = pivot_y + dx0*sin(theta) + dy0*cos(theta);

                    xi = round(pivot_x);
                    yi = round(pivot_y);
                end

                % only if in‐bounds
                if xi>=1 && xi<=imsize(2) && yi>=1 && yi<=imsize(1)
                    cell_mask(yi, xi) = true;      % cell occupies here
                    if cell_shape(dy, dx) == 1
                        lipid_mask(yi, xi) = true; % lipid sits here
                    end
                end
                prev_x = dx;
                prev_y = dy;
            end
        end

        existing_cells = all_cell_id > 0;
        if nnz(cell_mask & existing_cells) <= allowed_overlap % num of non-zero element and exisiting cell
            success = true;
        end
    end

    cell_start_pos(k, :) = [xc yc];

    if ~success || isempty(yc)
        % no more cells can be placed
        break;         
    end

    % commit this cell
    placed = placed + 1;

    if test_angle ~= pi/2 % super hack!
        cell_mask = imclose(cell_mask, strel('disk',1));
    end

    cell_id{placed} = cell_mask;
    lipid_dist{placed} = lipid_mask;       % lipid positions

    cell_id{k} = cell_id{k}(1+sz:end-sz, 1+sz:end-sz);
    lipid_dist{k} = lipid_dist{k}(1+sz:end-sz, 1+sz:end-sz);

    all_cell_id(cell_mask) = placed;
    all_lipid_dist(lipid_mask) = placed;
end

% trim unused entries
cell_id    = cell_id(1:placed);
lipid_dist = lipid_dist(1:placed);

cell_cnt = placed;
% plot the toy environment
figure(3)
clf
nexttile
imagesc(all_cell_id);
axis square off;
colormap([0 0 0; lines]);
% colormap lines

hold on
plot(1:imsize(1)-2, repelem(start_pos, imsize(1)-2), 'Color', 'k', 'LineWidth', 5)
plot(1:imsize(1)-2, repelem(start_pos, imsize(1)-2), 'Color', 'w', 'LineWidth', 3)

lineobj = streamline(X, Y, U, V, 1, start_pos, 'Color', 'k', 'LineWidth', 5);
lineobj = streamline(X, Y, U, V, 1, start_pos, 'Color', 'w', 'LineWidth', 3);
lineobj = stream2(X, Y, U, V, 1, start_pos);
lineobj = lineobj{1};
hold off
%% calculate overlap
overlap_cells = zeros(imsize-2);
overlap_per_cell = cell(cell_cnt, 1);

for i = 1:cell_cnt
    overlap_cells = overlap_cells + cell_id{i}; %(i, :, :);
end

for i = 1:cell_cnt
    curr_cell_pos = cell_id{i};
    curr_lipid = lipid_dist{i};

    overlap_per_cell{i} = overlap_cells(curr_cell_pos & curr_lipid);
end

figure(4)
imshow(overlap_cells, [1 max(overlap_cells(:))])

%% FUNCTION
function [X, Y, U, V] = set_alignment(phi, imsize, showfigure, nem_ws)
    if ~exist('nem_ws', 'var')
        if imsize(1) < 50
            nem_ws = 2;
        else
            nem_ws = 15;
        end
    end

    phi_nem = imresize(phi, [length(phi)/nem_ws length(phi)/nem_ws], 'nearest');

    [X, Y] = meshgrid(1:pi/2:pi/2*length(phi_nem), 1:pi/2:pi/2*length(phi_nem));

    X = X - 1/2 * cos(phi_nem);
    Y = Y - 1/2 * sin(phi_nem);
    
    X = rescale(X, 1, imsize(2));
    Y = rescale(Y, 1, imsize(1));
    
    U = cos(phi_nem);
    V = sin(phi_nem);

    if showfigure
        % figure
        nexttile
        imshow(zeros(size(phi))+1, [0, 1], colormap=gray)
        axis square off
        hold on
        quiver(X, Y, U, V, 0.5, 'color', 'black', 'ShowArrowHead', 'off', 'LineWidth', 2)
        hold off
    end
end
