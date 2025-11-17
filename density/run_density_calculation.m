%
%   run density quantification
%   quantify density by given area
%
clearvars
clc
close all

MATLAB_FX = "D:/Matlab_FileExchange/";

addpath D:/code/util/
addpath D:/code/density/

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATIONS

den_reg = 120; %220;

prefix = ["mCh300FBS10-1"];
welllist = {[1:24]};
framelist = {[1]};

ch = "c1"; %"c3";

im_type = "jpg"; % tif;

is_plot_den = 1;
is_plot_den_by_circ = 0;
is_plot_den_by_nuc = 0;

lwr = 20;
upr = 120;

%% run density for a time series density

savemat = 1;

% set density path 
if strcmp(ch, "c1")
    DENPATH = ch + "_mCherry_density/";
elseif strcmp(ch, "c2")
    DENPATH = ch + "_cy5_density";
elseif strcmp(ch, "c3")
    DENPATH = ch + "_dapi";
else
    error("channel is not recognized: "+ch);
end

% run density calculation for time series data
for p = 1:length(prefix)
    PATH = "D:/wave_" + prefix(p) + "/";

    for w = welllist{p}
        well = sprintf("%02d", w);

        % disp("prefix " + p + ", well " + well)
        
        inpath  = PATH + "Well" + well + "/" + DENPATH + "/segment/";
        outpath = PATH + "Well" + well + "/" + DENPATH + "/density_reg" + den_reg + "/";
        mkdir(outpath);

        for fr = framelist{p}
            frame = sprintf('%02d', fr);
            
            fname = inpath + set_filename(prefix(p), w, fr, ch) + "_stardist.tif";
            
            
            den_file = outpath + set_filename(prefix(p), w, fr, ch) + "_den" + den_reg + ".mat";

            if ~exist(den_file)
                [den, x, y, n_cell, imsize, valid_idx] = calculate_density(fname, den_reg);
            else
                den = struct2array(load(den_file));
            end

            if savemat == 1
                save(den_file, 'den');
            end

            if is_plot_den
                
                plot_den(den, 0, lwr, upr); % 0 or 1 for normalize to current image

                outfile = outpath + set_filename(prefix(p), w, fr, ch) + "_den" + den_reg + "." + im_type;
                exportgraphics(gcf, outfile);
            end

            if is_plot_den_by_circ
                plot_den_by_circ(den, x, y, n_cell, imsize);

                outfile = outpath + set_filename(prefix(p), w, fr, ch) + "_den" + den_reg + "_dot." + im_type;
                exportgraphics(gcf, outfile, 'Resolution', 314);
            end

            if is_plot_den_by_nuc
                plot_den_by_nuc(fname, n_cell, imsize, valid_idx);

                outfile = outpath + set_filename(prefix(p), w, fr, ch) + "_den" + den_reg + "_nuc." + im_type;
                exportgraphics(gcf, outfile);
            end
            close all
        end
    end
end

%% FUNCTIONS
function plot_den(den, norm, lwr, upr)
    if norm
        im_den = (den-min(den(:))) / (max(den(:))-min(den(:)));
    else
        im_den = den;
    end

    if ~exist('lwr', 'var')
        lwr = 0;
    end
    if ~exist('upr', 'var')
        upr = mean(im_den(:),'omitnan') + 2*std(im_den(:),'omitnan');
    end

    figure
    imshow(im_den,[lwr upr], 'Border', 'tight')
    colormap jet
    axis square
end

function plot_den_by_circ(den, x, y, n_cell, imsize)
    [f, order] = sort(n_cell);
    
    figure
    scatter(x(order), y(order), 10, f, 'filled');
    xticks([])
    yticks([])
    colormap jet
    axis square
    axis ij
    axis tight

    set(gcf,'Units','pixels','Position',[0 0 imsize(1) imsize(2)]);
end

function plot_den_by_nuc(seg, n_cell, imsize, valid_idx) % plot the density by nucleus shape
    seg = imread(seg);

    den_nuc = zeros(imsize);

    label_map = zeros(max(seg(:)), 1); % preallocate
    label_map(valid_idx) = n_cell; % only valid labels get values 

    for i = 1:length(label_map)
        if label_map(i) > 0
            den_nuc(seg==i) = label_map(i);
        end
    end

    figure
    % imshow(den_nuc, [1 max(den_nuc(:))], 'Border', 'tight')
    imshow(den_nuc, [], 'Border', 'tight')
    colormap([0 0 0; jet(256)])
    axis square
end
