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

prefix = ["mCh300FBS10-1" "mCh300FBS10-2" "overlap5-1" "overlap5-2"]; %["photo-truli-KI-4"]; 

sp = ["ent" "coh" "den"];

[ch, subpath] = load_file_paths();

ini_size = 150;
offset = 51; % decrease the boundary of the image by this offset


%% load initiation from initiation
welllist = {[1:12]};
n_wells = 12;
framelist = 1;
n_frame = 1;

% ini_file = "D:/Spatiotemporal_analysis/Initiation/" + prefix + inhibitor-1_initiation_details_pos_pattern_per5_cov10.txt";
% % ini_file = "D:\Spatiotemporal_analysis\Initiation\" + prefix + "_resist1_details_pattern_per5_cov10.txt";
ini_file = "D:\Spatiotemporal_analysis\Initiation\" + prefix + "_initiation_details_pattern_per5_cov10.txt";
ini = readtable(ini_file, 'Delimiter', '\t');
ini.id = string(ini.id);

imsize = [5210 5210]; % [4800 4800]; [3072 3072];

%% load initiation from spontaneous initiation
welllist = {[1:24], [1:23], [1:20], [1:11 13:20]};
n_wells = 24;

framelist = load_framelist(prefix, welllist, n_wells);
n_frame = 65;

imsize = load_imsize(prefix, n_wells);

s_ini_file = "D:/Spatiotemporal_analysis/Initiation/spontaneous_initiation_details_pattern_per5_cov10.txt";
s_ini = readtable(s_ini_file, 'Delimiter', '\t');
s_ini.id = string(s_ini.id);

%% load initiation from photoinduction
welllist = {[1 4 5 8 9 12], [1:3 5 7:10 12], 1:12, 1:12};
n_wells = 12;

framelist = num2cell(ones(length(prefix), n_wells));
n_frame = 1;

imsize = load_imsize(prefix, n_wells);

p_ini_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_initiation_details_pattern_per5_cov10.txt";
p_ini = readtable(p_ini_file, 'Delimiter', '\t');
p_ini.id = string(p_ini.id);

%% load random selected non-initiation background
bg_rep = 1;
bg_file = "D:/Spatiotemporal_analysis/Initiation/" + prefix + "_randbg"+bg_rep+"_details_pattern_per5_cov10.txt";

bg = readtable(bg_file, 'Delimiter', '\t');
bg.id = string(bg.id);

%% load random selected non-initiation background
res_rep = 1;
res_file = "D:/Spatiotemporal_analysis/Initiation/" + prefix + "_resist"+res_rep+"_details_pattern_per5_cov10.txt";

res = readtable(res_file, 'Delimiter', '\t');
res.id = string(res.id);

%% FOLDERS
% CY5CH = "c1";
% ADJPATH = CY5CH + "_cy5_adj/";
% DENPATH = CY5CH + "_cy5_density/density_reg75/";
% 
% DICCH = "c2"; % channel of the DIC file
% DICPATH = DICCH + "_DIC/";
% NEMPATH = DICCH + "_DIC_nematics/";
% ORIENTPATH = NEMPATH + "/orientation_files/";

%% load spatial property files

ent = cell(length(prefix), n_wells);
coh = cell(length(prefix), n_wells);
den = cell(length(prefix), n_wells);

for p = 1:length(prefix)
    disp(prefix(p));
    path = ROOTDIR + "/wave_" + prefix(p) + "/";
    % load entire image and initiations
    for w = welllist{p}
        for fr = framelist{p,w}
        % for fr = 1
            ent{p, w} = load_image('ent', '.mat', prefix(p), subpath('ent'), w, fr, ch('ent'));
            coh{p, w} = load_image('coh', '.mat', prefix(p), subpath('coh'), w, fr, ch('coh'));
            den{p, w} = load_image('den120', '.mat', prefix(p), subpath('den'), w, fr, ch('den'));
        end
    end
end
%% load spatial properties for each initiation

% ini may be spontaneous, photoinduction, background or resistant
sp = res;

sp_ent = cell(height(sp),1);
sp_coh = cell(height(sp),1);
sp_den = cell(height(sp),1);

for i = 1:height(sp)
    p = find(prefix==sp.id(i));

    w = sp.well(i);
    x = ceil(sp.x(i));
    y = ceil(sp.y(i));

    % imsz = imsize{p, w};
    imsz = imsize;

    % check if initiation occurs near edge. If so, skip the initiation
    if x-ini_size/2 < offset || y-ini_size/2 < offset || x+ini_size/2-1 > imsz(1)-offset || y+ini_size/2-1 > imsz(2)-offset
        continue;
    end

    if isempty(den{p,w})
        continue;
    end

    sp_ent{i} = crop_circle(ent{p,w}, x, y, ini_size/2);
    sp_coh{i} = crop_circle(coh{p,w}, x, y, ini_size/2);
    sp_den{i} = crop_circle(den{p,w}, x, y, ini_size/2);
end

%% save spatial properties as mat file
save_path = "D:\Spatiotemporal_analysis\initiation_statisitcs\" + prefix + "\sp_mat_files\";
mkdir(save_path);

% % save ini
% ini_ent = sp_ent;
% ini_coh = sp_coh;
% ini_den = sp_den;
% 
% save(save_path+"ini_ent.mat", 'ini_ent', '-v7.3');
% save(save_path+"ini_coh.mat", 'ini_coh', '-v7.3');
% save(save_path+"ini_den.mat", 'ini_den', '-v7.3');

% % save bg
% bg_ent = sp_ent; 
% bg_coh = sp_coh;
% bg_den = sp_den;
% 
% save(save_path+"bg_ent_rep"+bg_rep+".mat", 'bg_ent', '-v7.3');
% save(save_path+"bg_coh_rep"+bg_rep+".mat", 'bg_coh', '-v7.3');
% save(save_path+"bg_den_rep"+bg_rep+".mat", 'bg_den', '-v7.3');

% save res
res_ent = sp_ent; 
res_coh = sp_coh; 
res_den = sp_den; 

save(save_path+"resist_ent_rep"+res_rep+".mat", 'res_ent', '-v7.3');
save(save_path+"resist_coh_rep"+res_rep+".mat", 'res_coh', '-v7.3');
save(save_path+"resist_den_rep"+res_rep+".mat", 'res_den', '-v7.3');

%% all
all_ini_ent = [ini_ent; s_ini_ent];
all_ini_coh = [ini_coh; s_ini_coh];
all_ini_den = [ini_den; s_ini_den];

%% store the mean spatial properties into the table
df = ini; %bg; %p_ini;
df_ent = sp_ent; %bg_ent; %ini_ent;
df_coh = sp_coh; %bg_coh; %ini_coh;
df_den = sp_den; %bg_den; %p_ini_den;

mean_ent = nan(height(df), 1);
mean_coh = nan(height(df), 1);
mean_den = nan(height(df), 1);

for i = 1:height(df)
    mean_ent(i) = mean(df_ent{i}(:), 'omitnan');
    mean_coh(i) = mean(df_coh{i}(:), 'omitnan');
    mean_den(i) = mean(df_den{i}(:), 'omitnan');
end

df.ent = mean_ent;
df.coh = mean_coh;
df.den = mean_den;

%% Write to file
disp("writing to file")
outfile = "D:/Spatiotemporal_analysis/Initiation/spontaneous_initiation_details_pattern_per5_cov10_sp"+ini_size+".txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_initiation_details_pattern_per5_cov10_sp"+ini_size+".txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_randbg"+bg_rep+"_details_pattern_per5_cov10_sp"+ini_size+".txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_resist"+res_rep+"_details_pattern_per5_cov10_sp"+ini_size+".txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/photoinduct-1-1_initiation_details_pattern_per5_cov10_sp"+ini_size+".txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/inhibitor-1_initiation_details_pos_pattern_per5_cov10_sp"+ini_size+".txt";
% outfile = "D:\Spatiotemporal_analysis\Initiation\" + prefix + "_resist1_details_pattern_per5_cov10_sp"+ini_size+".txt";
writetable(df, outfile, 'Delimiter', '\t')

%% create outpath
outpath = "D:\Spatiotemporal_analysis\initiation_patterns\spatial_distribution_bg"+bg_rep+"_res"+res_rep+"\";
mkdir(outpath);

%% plot ini vs bg, res vs bg, and ini vs res
savefigure = 1;

outname = outpath + "\distr_" + ini_size + "_photospon_bg" + bg_rep + "_ent.svg";
% outname = outpath + "\distr_" + ini_size + "_res" + res_rep + "_bg" + bg_rep + "_ent.svg";
% outname = outpath + "\distr_" + ini_size + "_photospon_res" + res_rep + "_ent.svg";
plot_distr_spatialprop(bg_ent, ini_ent, [-4 1], 'all', savefigure, outname);

outname = outpath + "\distr_" + ini_size + "_photospon_bg" + bg_rep + "_coh.svg";
% outname = outpath + "\distr_" + ini_size + "_res" + res_rep + "_bg" + bg_rep + "_coh.svg";
% outname = outpath + "\distr_" + ini_size + "_photospon_res" + res_rep + "_coh.svg";
plot_distr_spatialprop(bg_coh, ini_coh, [0 0.7], 'all', savefigure, outname);

outname = outpath + "\distr_" + ini_size + "_photospon_bg" + bg_rep + "_den.svg";
% outname = outpath + "\distr_" + ini_size + "_res" + res_rep + "_bg" + bg_rep + "_den.svg";
% outname = outpath + "\distr_" + ini_size + "_photospon_res" + res_rep + "_den.svg";
plot_distr_spatialprop(bg_den, ini_den, [20 160], 'all', savefigure, outname);

%% plot initiation vs. non-initiation background vs. resistant region
savefigure = 1;

outname = outpath + "distr_"+ini_size+"_photo_rand"+bg_rep+"_resist"+res_rep+"_ent.jpg";
plot_distr_spatialprop_3(bg_ent, ini_ent, res_ent, [-4 1], 'all', savefigure, outname);

outname = outpath + "distr_"+ini_size+"_photo_rand"+bg_rep+"_resist"+res_rep+"_coh.jpg";
plot_distr_spatialprop_3(bg_coh, ini_coh, res_coh, [0 0.7], 'all', savefigure, outname);

outname = outpath + "distr_"+ini_size+"_photo_rand"+bg_rep+"_resist"+res_rep+"_den.jpg";
plot_distr_spatialprop_3(bg_den, ini_den, res_den, [20 160], 'all', savefigure, outname);

%% plot spontaneous initiation vs. non-initiation background vs. resistant region
savefigure = 1;

outname = outpath + "dist_"+ini_size+"_spon_rand"+bg_rep+"_resist"+res_rep+"_ent.jpg";
plot_distr_spatialprop_3(bg_ent, s_ini_ent, res_ent, [-4 1], 'all', savefigure, outname);

outname = outpath + "dist_"+ini_size+"_spon_rand"+bg_rep+"_resist"+res_rep+"_coh.jpg";
plot_distr_spatialprop_3(bg_coh, s_ini_coh, res_coh, [0 0.7], 'all', savefigure, outname);

outname = outpath + "dist_"+ini_size+"_spon_rand"+bg_rep+"_resist"+res_rep+"_den.jpg";
plot_distr_spatialprop_3(bg_den, s_ini_den, res_den, [20 160], 'all', savefigure, outname);

%% plot all initiation vs. non-initiation background vs. resistant region
savefigure = 1;

outname = outpath + "dist_"+ini_size+"_allini_rand"+bg_rep+"_resist"+res_rep+"_ent.jpg";
plot_distr_spatialprop_3(bg_ent, all_ini_ent, res_ent, [-4 1], 'all', savefigure, outname);

outname = outpath + "dist_"+ini_size+"_allini_rand"+bg_rep+"_resist"+res_rep+"_coh.jpg";
plot_distr_spatialprop_3(bg_coh, all_ini_coh, res_coh, [0 0.7], 'all', savefigure, outname);

outname = outpath + "dist_"+ini_size+"_allini_rand"+bg_rep+"_resist"+res_rep+"_den.jpg";
plot_distr_spatialprop_3(bg_den, all_ini_den, res_den, [20 160], 'all', savefigure, outname);

%% plot spontaneous initiation vs. initiation vs. non-initiation background

savefigure = 1;

outname = outpath + "dist_"+ini_size+"_photo_spon_rand"+bg_rep+"_ent.jpg";
plot_distr_spatialprop_3(bg_ent, s_ini_ent, ini_ent, [-4 1], 'all', savefigure, outname);

outname = outpath + "dist_"+ini_size+"_photo_spon_rand"+bg_rep+"_coh.jpg";
plot_distr_spatialprop_3(bg_coh, s_ini_coh, ini_coh, [0 0.7], 'all', savefigure, outname);

outname = outpath + "dist_"+ini_size+"_photo_spon_rand"+bg_rep+"_den.jpg";
plot_distr_spatialprop_3(bg_den, s_ini_den, ini_den, [20 160], 'all', savefigure, outname);


%% plot 2d scatter contour plots for all
outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\" + ...
        "ini"+ini_size+"_rand"+rep;
plot_scatter_contour_spatialprop(ini_ent, ini_den, bg_ent, bg_den, [-4 1], [20 160], 'ent', 'den', 'all', 1, outname);
plot_scatter_contour_spatialprop(ini_coh, ini_den, bg_coh, bg_den, [0 0.7], [20 160], 'coh', 'den', 'all', 1, outname);

%% plot 2d scatter contour plots for ini and resist regardless of pattern
savefigure = 1;
outname = outpath + "ini"+ini_size+"_photo_rand"+bg_rep+"_resist"+res_rep;
plot_scatter_contour_spatialprop_3(ini_ent, ini_den, bg_ent, bg_den, res_ent, res_den, [-4 1], [20 160], 'ent', 'den', 'all', savefigure, outname, 'jpg');
plot_scatter_contour_spatialprop_3(ini_coh, ini_den, bg_coh, bg_den, res_coh, res_den, [0 0.7], [20 160], 'coh', 'den', 'all', savefigure, outname, 'jpg');

%% plot 2d scatter contour plots for spontaneous ini and resist regardless of pattern
savefigure = 1;
outname = outpath + "ini"+ini_size+"_spon_rand"+bg_rep+"_resist"+res_rep;
plot_scatter_contour_spatialprop_3(s_ini_ent, s_ini_den, bg_ent, bg_den, res_ent, res_den, [-4 1], [20 160], 'ent', 'den', 'all', savefigure, outname, 'jpg');
plot_scatter_contour_spatialprop_3(s_ini_coh, s_ini_den, bg_coh, bg_den, res_coh, res_den, [0 0.7], [20 160], 'coh', 'den', 'all', savefigure, outname, 'jpg');


%% plot 2d scatter contour plots for spontaneous ini and resist regardless of pattern
savefigure = 1;
outname = outpath + "ini"+ini_size+"_allini_rand"+bg_rep+"_resist"+res_rep;
plot_scatter_contour_spatialprop_3(all_ini_ent, all_ini_den, bg_ent, bg_den, res_ent, res_den, [-4 1], [20 160], 'ent', 'den', 'all', savefigure, outname, 'jpg');
plot_scatter_contour_spatialprop_3(all_ini_coh, all_ini_den, bg_coh, bg_den, res_coh, res_den, [0 0.7], [20 160], 'coh', 'den', 'all', savefigure, outname, 'jpg');


%% obtain different patterns

pattern = ["phalf", "mhalf", "splay", "bend", "aligned"];
pttn_col = 8:12;
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

    outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\ini_dist_"+ini_size+"_rand"+rep+"_ent_"+pattern(pttn)+".svg";
    plot_distr_spatialprop(bg_ent, pttn_ini_ent, [-4 1], pattern(pttn), 1, outname);
    
    outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\ini_dist_"+ini_size+"_rand"+rep+"_coh_"+pattern(pttn)+".svg";
    plot_distr_spatialprop(bg_coh, pttn_ini_coh, [0 0.7], pattern(pttn), 1, outname);
    
    outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\ini_dist_"+ini_size+"_rand"+rep+"_den_"+pattern(pttn)+".svg";
    plot_distr_spatialprop(bg_den, pttn_ini_den, [20 160], pattern(pttn), 1, outname);
    
end

%% plot 2d scatter contour plots for each cellular pattern
for pttn = 1:length(pattern)
% for pttn = length(pattern)
    pttn_ini_ent = ini_ent(nan_idx(:,pttn));
    pttn_ini_coh = ini_coh(nan_idx(:,pttn));
    pttn_ini_den = ini_den(nan_idx(:,pttn));

    outname = "D:\Spatiotemporal_analysis\initiation_patterns\pattern_distribution\" + ...
        "ini"+ini_size+"_rand"+rep;
    plot_scatter_contour_spatialprop(pttn_ini_ent, pttn_ini_den, bg_ent, bg_den, [-4 1], [20 160], 'ent', 'den', pattern(pttn), 1, outname, 'svg');
    plot_scatter_contour_spatialprop(pttn_ini_coh, pttn_ini_den, bg_coh, bg_den, [0 0.7], [20 160], 'coh', 'den', pattern(pttn), 1, outname, 'svg');

end
close all;

%% FUNCTIONS

function I = load_image(spatial_type, suffix, prefix, sp_path, w, fr, ch)
    if contains(spatial_type, "den")
        try
            im_file = "D:/Spatiotemporal_analysis/wave_" + prefix + "/Well" + sprintf("%02d", w) + "/" + sp_path + "/" + ...
                set_filename(prefix, w, fr, ch) + "_" + spatial_type + suffix;
        catch
            im_file = "D:/Spatiotemporal_analysis/wave_" + prefix + "/Well" + sprintf("%02d", w) + "/" + sp_path + "/" + ...
                set_filename(prefix, w, fr, ch) + "_" + spatial_type + ".mat";
        end
    else
        im_file = "D:/Spatiotemporal_analysis/wave_" + prefix + "/Well" + sprintf("%02d", w) + "/" + sp_path + "/" + ...
                set_filename(prefix, w, fr, ch) + "_" + spatial_type + suffix;
    end


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
    pttn_color = containers.Map({'all', 'phalf', 'mhalf', 'splay', 'bend', 'aligned'}, ["#ED7D31", "red", "blue", "#FFC000", "#2CA02C", "black"]);

    figure
    % plot rand select non-ini background
    if isa(bg, 'cell')
        bg = vertcat(bg{:});
    else
        bg = vertcat(bg(:));
    end
    
    [values, edges] = histcounts(bg(:), 'Normalization', 'pdf');
    centers = (edges(1:end-1)+edges(2:end))/2;
    values = (values);
    
    area(centers, values, 'LineStyle', 'none', 'FaceColor','black', 'FaceAlpha', 0.2) % '#7F7F7F'
    
    % plot initiation
    
    if isa(ini, 'cell')
        ini = vertcat(ini{:});
    else
        ini = vertcat(ini(:));
    end
    
    [values, edges] = histcounts(ini(:), 'Normalization', 'pdf');
    centers = (edges(1:end-1)+edges(2:end))/2;
    values = (values);

    hold on
    area(centers, values, 'LineStyle', 'none', 'FaceColor',pttn_color(pattern), 'FaceAlpha', 0.2)
    hold off
    
    set(gca,'box','off','xcolor', 'w', 'ycolor', 'w', 'XTick',[], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [])
    xlim(xlimits)
    
    if savefigure
        if contains(outname, 'svg')
            saveas(gcf, outname, 'svg');
        elseif contains(outname, 'pdf')
            exportgraphics(gcf, outname, 'ContentType', 'vector');
        else
            exportgraphics(gcf, outname);
        end
    end
end

function plot_distr_spatialprop_3(bg, ini, res, xlimits, pattern, savefigure, outname)
    pttn_color = containers.Map({'all', 'phalf', 'mhalf', 'splay', 'bend', 'aligned'}, ["#ED7D31", "red", "blue", "#FFC000", "#2CA02C", "black"]);

    figure
    % plot rand select non-ini background
    bg = vertcat(bg{:});
    
    [values, edges] = histcounts(bg(:), 'Normalization', 'pdf');
    centers = (edges(1:end-1)+edges(2:end))/2;
    values = (values);
    
    area(centers, values, 'LineStyle', 'none', 'FaceColor','#7F7F7F', 'FaceAlpha', 0.2)
    
    % plot resistant region
    res = vertcat(res{:});
    
    [values, edges] = histcounts(res(:), 'Normalization', 'pdf');
    centers = (edges(1:end-1)+edges(2:end))/2;
    values = (values);

    hold on
    area(centers, values, 'LineStyle', 'none', 'FaceColor','black', 'FaceAlpha', 0.2)
    hold off

    % plot initiation
    ini = vertcat(ini{:});
    
    [values, edges] = histcounts(ini(:), 'Normalization', 'pdf');
    centers = (edges(1:end-1)+edges(2:end))/2;
    values = (values);

    hold on
    area(centers, values, 'LineStyle', 'none', 'FaceColor',pttn_color(pattern), 'FaceAlpha', 0.3)
    hold off
    
    set(gca,'box','off','xcolor', 'w', 'ycolor', 'w', 'XTick',[], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [])
    xlim(xlimits)
    
    if savefigure
        if contains(outname, 'svg')
            saveas(gcf, outname, 'svg');
        elseif contains(outname, 'pdf')
            exportgraphics(gcf, outname, 'ContentType', 'vector');
        else
            exportgraphics(gcf, outname);
        end
    end
end

function plot_scatter_contour_spatialprop(x_ini, y_ini, x_bg, y_bg, xlimits, ylimits, x_sprop, y_sprop, pttn, savefigure, outname, outtype)
    % hashes for color and label for plotting
    pttn_rgb = containers.Map({'all', 'phalf', 'mhalf', 'splay', 'bend', 'aligned'}, ...
        {[237/255 125/255 49/255], [1 0 0], [0 0 1], [1 192/255 0] [44/255 160/255 44/255], [0 0 0]});
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
    scatter(x_bg_mean, y_bg_mean, 5, 'MarkerFaceColor', '#696969', 'MarkerEdgeColor','#696969', 'MarkerFaceAlpha', 0.2,'MarkerEdgeAlpha',0.3)
    scatter(x_mean, y_mean, 5, 'MarkerFaceColor', pttn_rgb(pttn), 'MarkerEdgeColor',pttn_rgb(pttn), 'MarkerFaceAlpha', 0.4,'MarkerEdgeAlpha',0.4)
    set(gca, 'XLim', xlimits, 'YLim', ylimits, 'XLabel', text('String', lab(x_sprop)), 'YLabel', text('String', lab(y_sprop)), 'PlotBoxAspectRatio', [1,1,1])

    hold off
    
    if savefigure
        outfile = outname + "_2dscatter_" + x_sprop + "_" + y_sprop + "_" + pttn + "." + outtype; 
        if contains(outfile, 'svg')
            saveas(gcf, outfile, 'svg');
        elseif contains(outfile, 'pdf')
            exportgraphics(gcf, outfile, 'ContentType', 'vector');
        else
            exportgraphics(gcf, outfile);
        end
    end

    % plot contour for initiation
    figure
    contourf(xlimits(1):x_ws:xlimits(2), ylimits(1):y_ws:ylimits(2), vq', 1,'LineColor', 'none')
    colormap([1 1 1; pttn_rgb(pttn)])
    set(gca,'XTick',[], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'XLim', xlimits, 'YLim', ylimits, 'PlotBoxAspectRatio', [1,1,1])
    if savefigure
        outfile = outname + "_2dcontourIni_" + x_sprop + "_" + y_sprop + "_" + pttn + "." + outtype; 
        if contains(outfile, 'svg')
            saveas(gcf, outfile, 'svg');
        elseif contains(outfile, 'pdf')
            exportgraphics(gcf, outfile, 'ContentType', 'vector');
        else
            exportgraphics(gcf, outfile);
        end
    end
    
    % plot contour for rand select non-ini background
    figure
    contourf(xlimits(1):x_ws:xlimits(2), ylimits(1):y_ws:ylimits(2), bg_vq', 1, 'LineColor', 'none')
    colormap([1 1 1; 105/255 105/255 105/255])
%     colormap([1 1 1; 0 0 0])
    set(gca,'XTick',[], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'XLim', xlimits, 'YLim', ylimits, 'PlotBoxAspectRatio', [1,1,1])
    if savefigure
        outfile = outname + "_2dcontourBg_" + x_sprop + "_" + y_sprop + "_" + pttn + "." + outtype; 
        if contains(outfile, 'svg')
            saveas(gcf, outfile, 'svg');
        elseif contains(outfile, 'pdf')
            exportgraphics(gcf, outfile, 'ContentType', 'vector');
        else
            exportgraphics(gcf, outfile);
        end
    end
end

function plot_scatter_contour_spatialprop_3(x_ini, y_ini, x_bg, y_bg, x_res, y_res, xlimits, ylimits, x_sprop, y_sprop, pttn, savefigure, outname, outtype)
    % hashes for color and label for plotting
    pttn_rgb = containers.Map({'all', 'phalf', 'mhalf', 'splay', 'bend', 'aligned'}, ...
        {[237/255 125/255 49/255], [1 0 0], [0 0 1], [1 192/255 0] [44/255 160/255 44/255], [0 0 0]});
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

    % resistant region: calculate mean and contour height
    x_res_mean = cell2mat(cellfun(@(x) mean(x(:), "omitnan"), x_res, 'UniformOutput', false));
    y_res_mean = cell2mat(cellfun(@(x) mean(x(:), "omitnan"), y_res, 'UniformOutput', false));
    
    figure
    res_h = scatter_kde(x_res_mean, y_res_mean,  'filled', 'MarkerSize', 10);
    res_v = res_h.CData;
    set(gca, 'XLim', xlimits, 'YLim', ylimits, 'XLabel', text('String', lab(x_sprop)), 'YLabel', text('String', lab(y_sprop)), 'PlotBoxAspectRatio', [1,1,1])
    
    % remove nan
    idx = isnan(x_res_mean);
    x_res_mean(idx)=[];
    y_res_mean(idx)=[];
    
    % set up bg_vq
    res_v(idx) = [];
    res_vq = griddata(x_res_mean',y_res_mean',res_v,Xq,Yq, 'natural');
    
    % plot 2d figure
    figure
    hold on
    
    scatter(x_bg_mean, y_bg_mean, 5, 'MarkerFaceColor', 'black', 'MarkerEdgeColor','black', 'MarkerFaceAlpha', 0.2,'MarkerEdgeAlpha', 0.3)
    
    scatter(x_res_mean, y_res_mean, 5, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor','blue', 'MarkerFaceAlpha', 0.3,'MarkerEdgeAlpha', 0.4)
    scatter(x_mean, y_mean, 5, 'MarkerFaceColor', pttn_rgb(pttn), 'MarkerEdgeColor',pttn_rgb(pttn), 'MarkerFaceAlpha', 0.3,'MarkerEdgeAlpha', 0.4)
    
    set(gca, 'XLim', xlimits, 'YLim', ylimits, 'XLabel', text('String', lab(x_sprop)), 'YLabel', text('String', lab(y_sprop)), 'PlotBoxAspectRatio', [1,1,1])

    hold off
    
    if savefigure
        outfile = outname + "_2dscatter_" + x_sprop + "_" + y_sprop + "_" + pttn + "." + outtype; 
        if contains(outfile, 'svg')
            saveas(gcf, outfile, 'svg');
        elseif contains(outfile, 'pdf')
            exportgraphics(gcf, outfile, 'ContentType', 'vector');
        else
            exportgraphics(gcf, outfile);
        end
    end

    % plot contour for initiation
    figure
    contourf(xlimits(1):x_ws:xlimits(2), ylimits(1):y_ws:ylimits(2), vq', 1,'LineColor', 'none')
    colormap([1 1 1; pttn_rgb(pttn)])
    set(gca,'XTick',[], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'XLim', xlimits, 'YLim', ylimits, 'PlotBoxAspectRatio', [1,1,1])
    if savefigure
        outfile = outname + "_2dcontourIni_" + x_sprop + "_" + y_sprop + "_" + pttn + "." + outtype; 
        if contains(outfile, 'svg')
            saveas(gcf, outfile, 'svg');
        elseif contains(outfile, 'pdf')
            exportgraphics(gcf, outfile, 'ContentType', 'vector');
        else
            exportgraphics(gcf, outfile);
        end
    end
    
    % plot contour for rand select non-ini background
    figure
    contourf(xlimits(1):x_ws:xlimits(2), ylimits(1):y_ws:ylimits(2), bg_vq', 1, 'LineColor', 'none')
    colormap([1 1 1; 105/255 105/255 105/255])
    set(gca,'XTick',[], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'XLim', xlimits, 'YLim', ylimits, 'PlotBoxAspectRatio', [1,1,1])
    if savefigure
        outfile = outname + "_2dcontourBg_" + x_sprop + "_" + y_sprop + "_" + pttn + "." + outtype; 
        if contains(outfile, 'svg')
            saveas(gcf, outfile, 'svg');
        elseif contains(outfile, 'pdf')
            exportgraphics(gcf, outfile, 'ContentType', 'vector');
        else
            exportgraphics(gcf, outfile);
        end
    end

    % plot contour for rand select resistant region
    figure
    contourf(xlimits(1):x_ws:xlimits(2), ylimits(1):y_ws:ylimits(2), res_vq', 1, 'LineColor', 'none')
    colormap([1 1 1; 0 0 1])
    set(gca,'XTick',[], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'XLim', xlimits, 'YLim', ylimits, 'PlotBoxAspectRatio', [1,1,1])
    if savefigure
        outfile = outname + "_2dcontourRes_" + x_sprop + "_" + y_sprop + "_" + pttn + "." + outtype; 
        if contains(outfile, 'svg')
            saveas(gcf, outfile, 'svg');
        elseif contains(outfile, 'pdf')
            exportgraphics(gcf, outfile, 'ContentType', 'vector');
        else
            exportgraphics(gcf, outfile);
        end
    end
end
