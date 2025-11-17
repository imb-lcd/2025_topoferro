%
%   initation analysis
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath D:\Spatiotemporal_analysis\code
addpath D:\Spatiotemporal_analysis\code\util
addpath D:\Spatiotemporal_analysis\code\nematics

addpath D:\Matlab_FileExchange\scatter_kde

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATIONS
% ini_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_initiation_details_pattern_per5_cov10.txt";
ini_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-1-1_initiation_details.txt"

rep = 10;
bg_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_randbg"+rep+"_details_pattern_per5_cov10.txt";
outpath = "D:\Spatiotemporal_analysis\Initiation_clustering\figure\";


prefix = ["photoinduct-1" "photoinduct-2" "photoinduct-3" "photoinduct-4"];
% prefix = ["photoinduct-1-1"];

sp = ["ent" "coh" "den"];

[ch, subpath] = load_file_paths();

welllist = {[1 4 5 8 9 12], [1:3 5 7:10 12], 1:12, 1:12};
% welllist = {[2:4 6 9:11]}
n_wells = 12;
fr = 1;

IMSIZE = 5000;
ini_size = 150;
offset = 51; % decrease the boundary of the image by this offset

%% load initiation data frame
ini = readtable(ini_file, 'Delimiter', '\t');

ini.id = string(ini.id);
ini.experiment = string(ini.experiment);

%% load random selected non-initiation background
bg = readtable(bg_file, 'Delimiter', '\t');

bg.id = string(bg.id);
bg.experiment = string(bg.experiment);

%% load spatial property files

ent = cell(length(prefix), n_wells);
coh = cell(length(prefix), n_wells);
den = cell(length(prefix), n_wells);

for p = 1:length(prefix)
    disp(prefix(p));

    % load entire image and initiations
    for w = welllist{p}
        % ent{p, w} = load_image('ent', '.mat', prefix(p), subpath('ent'), w, fr, ch('ent'));
        ent{p, w} = load_image('ent', '.mat', prefix(p), "c3_DIC_nematics", w, fr, ch('ent'));
        coh{p, w} = load_image('coh', '.mat', prefix(p), subpath('coh'), w, fr, ch('coh'));
        den{p, w} = load_image('den120', '_5000.mat', prefix(p), subpath('den'), w, fr, ch('den'));
    end
end
%% load spatial properties for each initiation
ini_ent = cell(height(ini),1);
ini_coh = cell(height(ini),1);
ini_den = cell(height(ini),1);

for i = 1:height(ini)
    p = convertStringsToChars(ini.id(i));
    p = str2num(p(end));

    w = ini.well(i);
    x = ceil(ini.x(i));
    y = ceil(ini.y(i));

    % check if initiation occurs near edge. If so, skip the initiation
    if x-ini_size/2 < offset || y-ini_size/2 < offset || x+ini_size/2-1 > IMSIZE-offset || y+ini_size/2-1 > IMSIZE-offset
        continue;
    end

    ini_ent{i} = crop_circle(ent{p,w}, x, y, ini_size/2);
    ini_coh{i} = crop_circle(coh{p,w}, x, y, ini_size/2);
    ini_den{i} = crop_circle(den{p,w}, x, y, ini_size/2);
end

%% load spatial properties for randomly selected non-initiation backgrounds

bg_ent = cell(height(bg),1);
bg_coh = cell(height(bg),1);
bg_den = cell(height(bg),1);

for i = 1:height(bg)
    p = convertStringsToChars(bg.id(i));
    p = str2num(p(end));

    w = bg.well(i);
    x = ceil(bg.x(i));
    y = ceil(bg.y(i));

    % check if initiation occurs near edge. If so, skip the initiation
    if x-ini_size/2 < offset || y-ini_size/2 < offset || x+ini_size/2-1 > IMSIZE-offset || y+ini_size/2-1 > IMSIZE-offset
        continue;
    end

    bg_ent{i} = crop_circle(ent{p,w}, x, y, ini_size/2);
    bg_coh{i} = crop_circle(coh{p,w}, x, y, ini_size/2);
    bg_den{i} = crop_circle(den{p,w}, x, y, ini_size/2);
end

%% plot initiation vs. non-initiation background

outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\ini_dist_"+ini_size+"_rand"+rep+"_ent.jpg";
plot_distr_spatialprop(bg_ent, ini_ent, [-4 1], 'all', 1, outname);

outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\ini_dist_"+ini_size+"_rand"+rep+"_coh.jpg";
plot_distr_spatialprop(bg_coh, ini_coh, [0 0.7], 'all', 1, outname);

outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\ini_dist_"+ini_size+"_rand"+rep+"_den.jpg";
plot_distr_spatialprop(bg_den, ini_den, [20 160], 'all', 1, outname);



%% plot 2d scatter contour plots for all
outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\" + ...
        "ini"+ini_size+"_rand"+rep;
plot_scatter_contour_spatialprop(ini_ent, ini_den, bg_ent, bg_den, [-4 1], [20 160], 'ent', 'den', 'all', 1, outname);
plot_scatter_contour_spatialprop(ini_coh, ini_den, bg_coh, bg_den, [0 0.7], [20 160], 'coh', 'den', 'all', 1, outname);


%% obtain different patterns

pattern = ["phalf", "mhalf", "splay", "bend"];
pttn_col = 8:11;
pttn_idx = cell(length(pattern));

nan_idx = true(height(ini), length(pattern));
for i = 1:height(ini)
    p = convertStringsToChars(ini.id(i));
    p = str2num(p(end));

    x = ceil(ini.x(i));
    y = ceil(ini.y(i));

    % check if initiation occurs near edge. If so, skip the initiation
    if x-ini_size/2 < offset || y-ini_size/2 < offset || x+ini_size/2-1 > IMSIZE-offset || y+ini_size/2-1 > IMSIZE-offset
        nan_idx(i,:) = 0;
        continue
    end

    for pttn = 1:length(pttn_col)
        if ini{i, pttn_col(pttn)} == 0
            nan_idx(i, pttn) = 0;
        end
    end
end
%% plot distributions for each cellular pattern
for pttn = 1:length(pattern)
    pttn_ini_ent = ini_ent(nan_idx(:,pttn));
    pttn_ini_coh = ini_coh(nan_idx(:,pttn));
    pttn_ini_den = ini_den(nan_idx(:,pttn));

    outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\ini_dist_"+ini_size+"_rand"+rep+"_ent_"+pattern(pttn)+".jpg";
    plot_distr_spatialprop(bg_ent, pttn_ini_ent, [-4 1], pattern(pttn), 1, outname);
    
    outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\ini_dist_"+ini_size+"_rand"+rep+"_coh_"+pattern(pttn)+".jpg";
    plot_distr_spatialprop(bg_coh, pttn_ini_coh, [0 0.7], pattern(pttn), 1, outname);
    
    outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\ini_dist_"+ini_size+"_rand"+rep+"_den_"+pattern(pttn)+".jpg";
    plot_distr_spatialprop(bg_den, pttn_ini_den, [20 160], pattern(pttn), 1, outname);
    
end

%% plot 2d scatter contour plots for each cellular pattern
for pttn = 1:length(pattern)
    pttn_ini_ent = ini_ent(nan_idx(:,pttn));
    pttn_ini_coh = ini_coh(nan_idx(:,pttn));
    pttn_ini_den = ini_den(nan_idx(:,pttn));

    outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\" + ...
        "ini"+ini_size+"_rand"+rep;
    plot_scatter_contour_spatialprop(pttn_ini_ent, pttn_ini_den, bg_ent, bg_den, [-4 1], [20 160], 'ent', 'den', pattern(pttn), 1, outname);
    plot_scatter_contour_spatialprop(pttn_ini_coh, pttn_ini_den, bg_coh, bg_den, [0 0.7], [20 160], 'coh', 'den', pattern(pttn), 1, outname);

end


%% FUNCTIONS

function I = load_image(spatial_type, suffix, prefix, sp_path, w, fr, ch)
    im_file = "D:/Spatiotemporal_analysis/wave_" + prefix + "/Well" + sprintf("%02d", w) + "/" + sp_path + "/" + ...
          set_filename(prefix, w, fr, ch) + "_" + spatial_type + suffix;

    if contains(suffix, "mat")
        I = struct2array(load(im_file));
    else
        temp_I = imread(im_file);
        if length(temp_I) ~= 5000 && width(temp_I) ~= 5000
            temp_I = imresize(temp_I, [5000 5000], 'nearest');
        end
        I = temp_I;
    end
end

function I = crop_circle(I, x0, y0, r)
    th = linspace(0,2*pi);
    x = r*cos(th) + x0;
    y = r*sin(th) + y0;

    [size_x, size_y] = size(I);
    [X, Y] = meshgrid(1:size_y, 1:size_x);
    idx = inpolygon(X(:), Y(:), x, y);

    I(~idx) = NaN;
    I = imcrop(I, [x0-r, y0-r, r*2-1 r*2-1]);
end

function plot_distr_spatialprop(bg, ini, xlimits, pattern, savefigure, outname)
    pttn_color = containers.Map({'all', 'phalf', 'mhalf', 'splay', 'bend'}, ["#ED7D31", "red", "blue", "#FFC000", "#2CA02C"]);

    figure
    % plot rand select non-ini background
    bg = vertcat(bg{:});
    
    [values, edges] = histcounts(bg(:), 'Normalization', 'pdf');
    centers = (edges(1:end-1)+edges(2:end))/2;
    values = (values);
    
    area(centers, values, 'LineStyle', 'none', 'FaceColor','#7F7F7F', 'FaceAlpha', 0.2)
    
    % plot initiation
    ini = vertcat(ini{:});
    
    [values, edges] = histcounts(ini(:), 'Normalization', 'pdf');
    centers = (edges(1:end-1)+edges(2:end))/2;
    values = (values);

    hold on
    area(centers, values, 'LineStyle', 'none', 'FaceColor',pttn_color(pattern), 'FaceAlpha', 0.2)
    hold off
    
    % set(gca,'box','off','xcolor', 'w', 'ycolor', 'w', 'XTick',[], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [])
    set(gca,'box','off', 'ycolor', 'w', 'YTick', [], 'YTickLabel', [])
    xlim(xlimits)
    
    if savefigure
        exportgraphics(gcf, outname);
    end
end

function plot_scatter_contour_spatialprop(x_ini, y_ini, x_bg, y_bg, xlimits, ylimits, x_sprop, y_sprop, pttn, savefigure, outname)
    % hashes for color and label for plotting
    pttn_rgb = containers.Map({'all', 'phalf', 'mhalf', 'splay', 'bend'}, ...
        {[237/255 125/255 49/255], [1 0 0], [0 0 1], [1 192/255 0] [44/255 160/255 44/255]});
    lab = containers.Map({'ent', 'den', 'coh'}, {'entropy', 'density', 'coherency'});

    % calculate the mean of the spatial properties
    x_mean = cell2mat(cellfun(@(x) mean(x(:), "omitnan"), x_ini, 'UniformOutput', false));
    y_mean = cell2mat(cellfun(@(x) mean(x(:), "omitnan"), y_ini, 'UniformOutput', false));
    
    % plot density plots.
    figure
    h = scatter_kde(x_mean, y_mean,  'filled', 'MarkerSize', 10);
    v = h.CData;
    set(gca, 'XLim', xlimits, 'YLim', ylimits, 'XLabel', text('String', lab(x_sprop)), 'YLabel', text('String', lab(y_sprop)), 'PlotBoxAspectRatio', [1,1,1])

    % remove nan
    idx = isnan(x_mean);
    x_mean(idx)=[];
    y_mean(idx)=[];
    
    % set up vq
    v(idx) = [];
    x_ws = (xlimits(2) - xlimits(1))/100; % (0.5 - -3.5)/100;
    y_ws = (ylimits(2) - ylimits(1))/100; % (160-20)/100;
    [Xq, Yq] = ndgrid(xlimits(1):x_ws:xlimits(2), ylimits(1):y_ws:ylimits(2));
    vq = griddata(x_mean',y_mean',v,Xq,Yq, 'natural');
    
    % non-ini background: calculate mean and contour height
    x_bg_mean = cell2mat(cellfun(@(x) mean(x(:), "omitnan"), x_bg, 'UniformOutput', false));
    y_bg_mean = cell2mat(cellfun(@(x) mean(x(:), "omitnan"), y_bg, 'UniformOutput', false));
    
    figure
    bg_h = scatter_kde(x_bg_mean, y_bg_mean,  'filled', 'MarkerSize', 10);
    bg_v = bg_h.CData;
    set(gca, 'XLim', xlimits, 'YLim', ylimits, 'XLabel', text('String', lab(x_sprop)), 'YLabel', text('String', lab(y_sprop)), 'PlotBoxAspectRatio', [1,1,1])
    
    % remove nan
    idx = isnan(x_bg_mean);
    x_bg_mean(idx)=[];
    y_bg_mean(idx)=[];
    
    % set up bg_vq
    bg_v(idx) = [];
    bg_vq = griddata(x_bg_mean',y_bg_mean',bg_v,Xq,Yq, 'natural');
    
    % plot 2d figure
    figure
    hold on
    scatter(x_bg_mean, y_bg_mean, 5, 'MarkerFaceColor', '#696969', 'MarkerEdgeColor','#696969', 'MarkerFaceAlpha', 0.1,'MarkerEdgeAlpha',.1)
    scatter(x_mean, y_mean, 5, 'MarkerFaceColor', pttn_rgb(pttn), 'MarkerEdgeColor',pttn_rgb(pttn), 'MarkerFaceAlpha', 0.1,'MarkerEdgeAlpha',.2)
   set(gca, 'XLim', xlimits, 'YLim', ylimits, 'XLabel', text('String', lab(x_sprop)), 'YLabel', text('String', lab(y_sprop)), 'PlotBoxAspectRatio', [1,1,1])

    hold off
    
    if savefigure
        outfile = outname + "_2dscatter_" + x_sprop + "_" + y_sprop + "_" + pttn + ".jpg"; 
        exportgraphics(gcf, outfile);
    end

    % plot contour for initiation
    figure
    contourf(xlimits(1):x_ws:xlimits(2), ylimits(1):y_ws:ylimits(2), vq', 1,'LineColor', 'none')
    colormap([1 1 1; pttn_rgb(pttn)])
    set(gca,'XTick',[], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'XLim', xlimits, 'YLim', ylimits, 'PlotBoxAspectRatio', [1,1,1])
    if savefigure
        outfile = outname + "_2dcontourIni_" + x_sprop + "_" + y_sprop + "_" + pttn + ".jpg"; 
        exportgraphics(gcf, outfile);
    end
    
    % plot contour for rand select non-ini background
    figure
    contourf(xlimits(1):x_ws:xlimits(2), ylimits(1):y_ws:ylimits(2), bg_vq', 1, 'LineColor', 'none')
    colormap([1 1 1; 105/255 105/255 105/255])
    set(gca,'XTick',[], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'XLim', xlimits, 'YLim', ylimits, 'PlotBoxAspectRatio', [1,1,1])
    if savefigure
        outfile = outname + "_2dcontourBg_" + x_sprop + "_" + y_sprop + "_" + pttn + ".jpg"; 
        exportgraphics(gcf, outfile);
    end
end