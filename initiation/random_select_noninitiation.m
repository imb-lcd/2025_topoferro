%
% randomly select non-initiation background
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
%ROOTDIR = "/home/N417/Jen-Hao/Spatiotemporal/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath D:\Spatiotemporal_analysis\code
addpath D:\Spatiotemporal_analysis\code\util
addpath D:\Spatiotemporal_analysis\code\nematics

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATIONs
rep = 3;
ini_size = 150;

% prefix = ["photoinduct-1" "photoinduct-2" "photoinduct-3" "photoinduct-4"];
% welllist = {[1 4 5 8 9 12], [1:3 5 7:10 12], 1:12, 1:12};
% n_wells = 12;

% ini_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_initiation_details.txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_randbg" + rep + "_details.txt";
% outpath = "D:/Spatiotemporal_analysis/Initiation_clustering/photoinduct_back" + ini_size + "_rep" + rep + "/";


% ["photo-truli-KI-1" "photo-truli-KI-4" "photo-truli-ki-0730-1" "photo-truli-ki-0730-4"];
prefix = "photo-truli-KI-1";
welllist = {[1:12]};
% welllist = {[2:5 8:11]};
% welllist = {[1 6:7 12]};
n_wells = 12;

% prefix = "yapi-0827z4";
% welllist = {[1:24]};
% n_wells = 12;

ini_file = "D:/Spatiotemporal_analysis/Initiation/" + prefix + "_initiation_details.txt";
outfile = "D:/Spatiotemporal_analysis/Initiation/" + prefix + "_randbg" + rep + "_details.txt";
% outpath = "D:/Spatiotemporal_analysis/Initiation_clustering/" + prefix + "_back" + ini_size + "_rep" + rep + "/";
% mkdir(outpath);

IMSIZE = 4800; %5120; %4800; %5000;

%% load initiation file and data
ini = readtable(ini_file);
ini.id = string(ini.id);

%% Set unselectable initiation region
close all;
showfigure = 1;

all_mask = cell(length(prefix), n_wells);
ini_exp = cell(length(prefix), n_wells);

for p = 1:length(prefix)
    disp(prefix(p));

    % fill in noselect matrix
    noselect = cell(n_wells);
    for w = welllist{p}
        noselect{w} = zeros(IMSIZE);
    end
    
    % get unselectable region as a mask of ROIs

    for w = welllist{p}
        % select for the specific prefix and well for each iteration
        ini_pw = ini(strcmp(ini.id, prefix(p)) & ini.well == w, :);

        % save the experiment notes
        ini_exp{p, w} = string(cell2mat(unique(ini_pw.experiment)));
    
        % obtain upper left corner of the unselectable region (extend ini)
        x_ul_corner = ini_pw.x - ini_size;
        y_ul_corner = ini_pw.y - ini_size;
        
        % create a mask for each initiation
        mask = false(IMSIZE);
        for j = 1:length(x_ul_corner)
            h = images.roi.Rectangle('Position', [x_ul_corner(j), y_ul_corner(j), ini_size*2, ini_size*2]);
            mask = mask | createMask(h, noselect{w});
        end
        all_mask{p, w} = ~mask;

        if showfigure
            figure;
            imshow(mask)
        end
    end
end

%% randomly select 150 regions that are not initiations for each well
close all
showfigure = 1;

n = 150;

bg = cell(length(prefix), n_wells);

for p = 1:length(prefix)
    disp(prefix(p));

    for w = welllist{p}
        rand_x{w} = NaN(1, n);
        rand_y{w} = NaN(1, n);

        for r = 1:n
            while true % loop until it finds a region not in the mask
                rect = randomWindow2d([IMSIZE IMSIZE], [ini_size ini_size]);
               
                % note: need to flip the xycoord when comparing to mask
                if all_mask{p, w}(rect.YLimits(1), rect.XLimits(1)) ~= 0
                    rand_x{w}(r) = rect.XLimits(1);
                    rand_y{w}(r) = rect.YLimits(1);
                    break;
                end
            end

            % save the randomly selected non-ini background in the format
            % id, experiment, well, frame, index, x, y
            % bg{p, w} = [repelem(prefix(p), n)' repelem(ini_exp{p, w}, n)' repelem(w, n)' repelem(1, n)' (1:n)' rand_x{w}' rand_y{w}'];
            bg{p, w} = [repelem(prefix(p), n)' repelem("KI696", n)' repelem(w, n)' repelem(1, n)' (1:n)' rand_x{w}' rand_y{w}'];
        end
        
        if showfigure
            figure
            imshow(all_mask{p, w})
            hold on
            for k = 1:n
                drawrectangle('Position', [rand_x{w}(k) rand_y{w}(k) ini_size ini_size], 'Color', [1 0 0]);
            end
            hold off
        end
    end
end

%% write th non-initiation background into file
bg = vertcat(bg{:});
bg = array2table(bg, 'VariableNames', {'id', 'experiment', 'well', 'frame', 'index', 'x', 'y'});

writetable(bg, outfile, 'Delimiter', '\t');

%% load and crop the images
im_types = ["mch", "nem", "ent"];

for p = prefix
    disp(prefix(p))
    for w = welllist{p}
        for type = im_types
            % determine which frame to crop the background
            fr = 1; % ini(w);
    
            [I, fname] = load_images(w, fr, type, ini_size);
    %         figure
    %         imshow(I, [min(I(:)) max(I(:))])
    %         %imshow(all_mask{w})
    %         hold on
    %         for k = 1:ini_num{w}
    %             drawrectangle('Position', [rand_x{w}(k) rand_y{w}(k) ini_ws ini_ws], 'Color', [1 0 0]);
    %         end
    %         hold off
            
            for k = 1:length(rand_x{w})
                if ~isnan(rand_x{w}(k))
                    outfile = outpath + fname + "_back" + k + "_" + type + ".tif";
                    crop_background_image(outfile, I, rand_x{w}(k), rand_y{w}(k), type{1});
                end
            end
    
        end
    end
end


%% FUNCTIONS

% crop the images and save to designated location
function CI = crop_background_image(outfile, I, x, y, type, ini_size)
    for i = 1:length(x)
        CI = imcrop(I, [x(i) y(i) ini_size-1 ini_size-1]);
        if strcmp(type, "coh")
            CI = uint8(CI * 255);
        end
        imwrite(CI, outfile);
    end
end

% load the images
function [I, fname] = load_images(well, frame, type, IMSIZE)
    if ~exist('IMSIZE', 'var')
        IMSIZE = [5000 5000];
    end

    [ch, subpath] = load_file_paths();

    well = num2str(well, '%02.f');
    frame = num2str(frame, '%02.f');

    % set the image file path
    fname = set_filename(prefix, well, frame, ch(type));
    if strcmp(type, "mch")
        im_file = PATH + "Well" + well + "/" + subpath(type) + fname + ".tif";
    elseif strcmp(type, "den")
        im_file = PATH + "Well" + well + "/" + subpath(type) + fname + "_den" + den_reg + "_5000.mat";
    elseif strcmp(type, "dic")
        im_file = PATH + "Well" + well + "/" + subpath(type) + fname + "_dic.tif";
    elseif strcmp(type, "nem")
        im_file = PATH + "Well" + well + "/" + subpath(type) + fname + "_nem.tif";
    elseif strcmp(type, "coh")
        im_file = PATH + "Well" + well + "/" + subpath(type) + fname + "_coh.mat";
    elseif strcmp(type, "ent")
        im_file = PATH + "Well" + well + "/" + subpath(type) + fname + "_ent.tif";
    elseif strcmp(type, "nemColor")
        im_file = PATH + "Well" + well + "/" + subpath(type) + fname + "nemColor.tif";
    else
        warning("image type: " + type + " is not recognized");
    end

    if contains(im_file, ".mat")
        I = struct2array(load(im_file));
    elseif contains(im_file, ".tif")
        I = imread(im_file);
    end

    if contains(type, 'nem')
        I = imresize(temp_I, IMSIZE, 'nearest');
    end

    % pad the borders if the file is less than the image size
    sz = size(I);
    if sz ~= IMSIZE
    
        new_I = NaN(IMSIZE);
        
        offset_x = ceil((IMSIZE(1)-sz(1))/2);
        offset_y = ceil((IMSIZE(2)-sz(2))/2);
    
        new_I(offset_x+1:(offset_x+sz(1)), offset_y+1:(offset_y+sz(2))) = [I];
    
        I = new_I;
    end
end

