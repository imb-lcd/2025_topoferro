%
% Extract initiation images
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath D:\Spatiotemporal_analysis\code
addpath D:\Spatiotemporal_analysis\code\util
addpath D:\Spatiotemporal_analysis\code\nematics

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATIONS
% IMSIZE = 5000;
ini_size = 150;
offset = 51; % decrease the boundary of the image by this offset

ini_path = "D:/Spatiotemporal_analysis/Initiation/";

% photoinduct
% ini_type = 'ini';
% ini_file = ini_path + "photoinduct-all_initiation_details_pattern_per5_cov10.txt";
% outpath = "D:/Spatiotemporal_analysis/Initiation_clustering/photoinduct_ini_" + ini_size + "_rotated_10_fliplrud/";

% spontaneous
% ini_type = 'ini';
% ini_file = ini_path + "spontaneous_initiation_details_pattern_per5_cov10.txt";
% outpath = "D:/Spatiotemporal_analysis/Initiation_clustering/photoinduct_ini_" + ini_size + "_rotated_10_spon_fliplrud/";

% % background
ini_type = 'bg';
rep = 8;
ini_file = ini_path + "photoinduct-all_randbg" + rep + "_details_pattern_per5_cov10_sp150.txt";
outpath = "D:\Spatiotemporal_analysis\Initiation_clustering\photoinduct_randbg" + rep + "_" + ini_size + "_rotated_10\";

% % resistance
% ini_type = 'res';
% rep = 13;
% ini_file = ini_path + "photoinduct-all_resist" + rep + "_details_pattern_per5_cov10.txt";
% outpath = "D:\Spatiotemporal_analysis\Initiation_clustering\photoinduct_resist" + rep + "_" + ini_size + "_rotated_10\";

mkdir(outpath);

ini = readtable(ini_file, 'Delimiter', '\t');
ini.id = string(ini.id);

% prefix = ["mCh300FBS10-1", "mCh300FBS10-2" "overlap5-1" "overlap5-2"];
% n_wells = 24;
% welllist = {[1:24], [1:23], [1:20], [1:11 13:20]};
% framelist = load_framelist(prefix, welllist, n_wells);
% n_frame = 65;

prefix = ["photoinduct-1" "photoinduct-2" "photoinduct-3" "photoinduct-4"];
n_wells = 12;
welllist = {[1 4 5 8 9 12], [1:3 5 7:10 12], 1:12, 1:12};
framelist = num2cell(ones(length(prefix), n_wells));
n_frame = 1;

imsize = load_imsize(prefix, n_wells);

% im_type = ["nem" "mch" "ent" "den"];
im_type = ["nem" "dota"];

[ch, subpath] = load_file_paths();


%% load images
phi = cell(length(prefix), n_wells);
nem = cell(length(prefix), n_wells);
dota = cell(length(prefix), n_wells);


for p = 1:length(prefix)
% for p = 3
    disp(prefix(p));
    path = ROOTDIR + "/wave_" + prefix(p) + "/";
    for w = welllist{p}
        for fr = framelist{p,w}
            phi{p,w} = load_images(path, prefix(p), w, fr, "phi", imsize{p,w});
            nem{p,w} = load_images(path, prefix(p), w, fr, "nem", imsize{p,w});
            dota{p,w} = load_images(path, prefix(p), w, fr, "dota", imsize{p,w});
        end
    end
end

%% calculate the angle to rotate

rotate_sz = 150/2;
flip = 0;
flip_suffix = "";

for i = 1:size(ini, 1)
    p = find(prefix==ini.id(i));
    
    w = ini.well(i);
    x = ceil(ini.x(i));
    y = ceil(ini.y(i));

    fr = framelist{p,w};

    imsz = imsize{p, w};

    if x-ini_size/2 < offset || y-ini_size/2 < offset || x+ini_size/2-1 > imsz(1)-offset || y+ini_size/2-1 > imsz(2)-offset
        continue;
    end

    if isempty(phi{p,w})
        continue;
    end
   
    % calculate angle to be rotated, calculate from the mean of a half of
    % the initiation size window at initiation
    J = imcrop(phi{p,w}, [x-rotate_sz y-rotate_sz rotate_sz*2 rotate_sz*2]);
    mean_angle = circular_angle_mean(J(:));
    ang = 90-rad2deg(mean_angle);
    % ang = 0;

    % rotate the image according to the angle in phi
    phi_im = rotate_image(phi{p,w}, x, y, ang, ini_size);
    phi_im = circular_angle_difference(deg2rad(ang), phi_im, 'signed');
    dota_im = rotate_image(dota{p,w}, x, y, ang, ini_size);

    % replot nematic directors and save the image
    if flip == 1
        dota_im = fliplr(dota_im);

        flip_suffix = "_fliplr"; 
    elseif flip == 2
        dota_im = flipud(dota_im);
        
         flip_suffix = "_flipud";

    elseif flip == 3
        dota_im = flipud(fliplr(dota_im));
        
        flip_suffix = "_fliplrud";
    end

    ini_outname = outpath + set_filename(prefix(p), w, fr, ch('nem')) + "_" + ini_type + ini.index(i) + "_nem" + flip_suffix + ".tif";
    draw_rotated_nematic_director(phi_im, ini_size, ini_outname, flip);

    ini_outname = outpath + set_filename(prefix(p), w, fr, ch('dota')) + "_" + ini_type + ini.index(i) + "_dota" + flip_suffix + ".tif";
    imwrite(dota_im, ini_outname, 'tif');
end

%% FUNCTION
function J = rotate_image(I, x, y, angle, ini_size)
    % calculate the how much to enlarge the image using the hypotenuse of the image. 
    hyp = ceil(hypot(ini_size, ini_size))+10;
    large_I = imcrop(I, [x-hyp/2 y-hyp/2 hyp hyp]);

    % rotate the image by angle
    R = imrotate(large_I, angle, 'nearest');

    rx = length(R)/2;
    ry = width(R)/2;
    
    %recrop the image from rotated image, centered at the original center
    J = imcrop(R, [rx-ini_size/2 ry-ini_size/2 ini_size-1 ini_size-1]);
end

function [I, fname] = load_images(path, prefix, well, frame, type, imsize)
    if ~exist('imsize', 'var')
        imsize = [5000 5000];
    end

    [ch, subpath] = load_file_paths();

    well = num2str(well, '%02.f');
    frame = num2str(frame, '%02.f');
    
    den_reg = 60;

    % set the image file path
    fname = set_filename(prefix, well, frame, ch(type));
    if strcmp(type, "mch")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + ".tif";
    elseif strcmp(type, "seg")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + "_stardist.tif";
    elseif strcmp(type, "den")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + "_den" + den_reg + "_5000.mat";
        % im_file = path + "Well" + well + "/" + subpath(type) + fname + "_den" + den_reg + ".mat";
        % im_file = path + "Well" + well + "/" + subpath(type) + fname + "_den" + den_reg + "_dot.tif";
    elseif strcmp(type, "dot")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + "_den_dot.tif";
    elseif strcmp(type, "dota")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + "_den_dot.tif";
    elseif strcmp(type, "dic")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + "_dic.tif";
    elseif strcmp(type, "phi")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + "_phi.mat";
    elseif strcmp(type, "nem")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + "_nem.tif";
    elseif strcmp(type, "coh")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + "_coh.mat";
    elseif strcmp(type, "ent")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + "_ent.mat";
    elseif strcmp(type, "nemColor")
        im_file = path + "Well" + well + "/" + subpath(type) + fname + "nemColor.tif";
    else
        warning("image type: " + type + " is not recognized");
    end
    
    if contains(im_file, ".mat")
        I = struct2array(load(im_file));
    elseif contains(im_file, ".tif")
        I = imread(im_file);
    end

    % process for specific image type
    if contains(type, 'nem') % make nematics IMSIZE and 2D
        I = imresize(I, fliplr(imsize), 'nearest');
        I = rgb2gray(I);
    elseif contains(type, 'den')
        I = imresize(I, fliplr(imsize), 'nearest');
        % I = rgb2gray(I);
    elseif contains(type, 'dot')
        I = imresize(I, fliplr(imsize), 'nearest');
        I = rgb2gray(I);
    elseif strcmp(type, 'phi')
        I = rescale(I, -pi/2, pi/2); 
    elseif strcmp(type, 'seg')
        I(I > 0) = 1;
    end
    
    % pad the borders if the file is less than the image size
    sz = size(I);
    imsize = fliplr(imsize);
    if sz ~= imsize
        new_I = NaN(imsize);

        offset_x = ceil((imsize(1)-sz(1))/2);
        offset_y = ceil((imsize(2)-sz(2))/2);
    
        new_I(offset_x+1:(offset_x+sz(1)), offset_y+1:(offset_y+sz(2))) = [I];
    
        I = new_I;
    end
end

% function I = normalize_image(I, prefix, welllist)
%     min_I = min(min([I{:}]));
%     max_I = max(max([I{:}]));
% 
%     I = cellfun(@(x) uint16((x-min_I)/(max_I-min_I)*2^16), I, 'UniformOutput', false);
% end

function ang = circular_angle_difference(ang1, ang2, sign)
    % ang = min(pi/2-abs(ang1)+pi/2-abs(ang2), abs(ang1)+abs(ang2));
    if strcmp(sign, 'signed')
        ang2 = -ang2;
        ang1 = mod((ang1 + pi/2), pi) - pi/2;
        ang2 = mod((ang2 + pi/2), pi) - pi/2;
        
        % Calculate difference
        ang = ang1 - ang2;
        
        % Convert back to original range
        ang = mod((ang + pi/2), pi) - pi/2;
        if ang < -pi/2
            ang = ang + pi;
        elseif ang > pi/2
            ang = ang - pi;
        end

    elseif strcmp(sign, 'unsigned')
        ang = pi/2 - abs(abs(ang1-ang2) - pi/2);
    else
        error("Can only be signed or unsigned: "+sign);
    end
end

function mean_angle = circular_angle_mean(ang)
    % Convert angles to Cartesian coordinates
    x = cos(2*ang);
    y = sin(2*ang);
    
    % Calculate the mean of the Cartesian coordinates
    sum_x = sum(x, 'omitnan');
    sum_y = sum(y, 'omitnan');
    
    % Convert mean Cartesian coordinates back to an angle in degrees
    mean_angle = atan2(sum_y, sum_x)/2;
    
    % Ensure the angle is within the range [0, 360]
    if mean_angle <= 0
        mean_angle = mean_angle + pi;
    end
end

function draw_rotated_nematic_director(phi_im, ini_size, nem_outname, flip)
    scale_factor = 1;
    nem_color = [0 0 0];
    nem_ws = 18.75; % in cluster: 18.75; ini example: 75
    nem_width  = 5; % ini cluster: 5: ini example: 2

    % for redrawing nematics image
    phi_nematic = -phi_im;
    phi_nematic = rescale(phi_nematic, -pi/2, pi/2);

%     nexttile
%     imshow(phi_nematic, [min(phi_nematic(:)) max(phi_nematic(:))], colormap=hsv)
%     title('phi nematics')
% 
%     nexttile

    phi_nematic = imresize(phi_nematic, [length(phi_im)/nem_ws length(phi_im)/nem_ws], 'nearest');

    if flip == 1
        phi_nematic = fliplr(-phi_nematic);
    elseif flip == 2
        phi_nematic = flipud(-phi_nematic);
    elseif flip == 3
        phi_nematic = flipud(fliplr(-phi_nematic));
    end

    nem_sz = length(phi_nematic);
    [X, Y] = meshgrid(1:pi/2:pi/2*nem_sz, 1:pi/2:pi/2*nem_sz);
    X = X - scale_factor/2*cos(phi_nematic);
    Y = Y - scale_factor/2*sin(phi_nematic);

    % scale quiver accordingly
    U = cos(phi_nematic)*scale_factor;
    V = sin(phi_nematic)*scale_factor;

    % plot quiver
    quiver(X, Y, U, V, 0.5, 'color', nem_color, 'LineWidth', nem_width, 'ShowArrowHead','off'); % 'Autoscale', 'on'
    axis('image', 'off', 'tight', 'ij') %  set the borders to tight and without axis
    set(gcf, 'Color', 'w'); % set the background to white
    
    % Capture the content of the current axis as an image
    frame = getframe(gca);
    
    % Convert the frame to an image
    image_data = frame.cdata;
    
    % Resize the image to 150x150 pixels
    resized_image = imresize(image_data, [ini_size, ini_size]);
    
    imwrite(resized_image, nem_outname);
end
