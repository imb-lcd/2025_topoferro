%
%   assign spatial properties to initiation/background/resistant files
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
% prefix = ["photo-truli-KI-1" "photo-truli-KI-4" "photo-truli-ki-0730-1" "photo-truli-ki-0730-4"];
% prefix = ["yapi-0816" "yapi-0819" "yapi-0827z4"]; 
prefix = "photo-truli-ki-0730-4";

sp = ["ent" "coh" "den"];

[ch, subpath] = load_file_paths();

imsize = {[5120 5120]}; %, [5120 5120], [4800 4800], [5120 5120]};
ini_size = 150;
offset = 51; % decrease the boundary of the image by this offset

% welllist = {[1 3 4 9 10 12 13 15 16 21 22 24], [1 2 3 10 11 12 13 14 15 22 23 24], [1:24]};
welllist = {[1:12]}; %, [1:12], [1:12], [1:12]};
n_wells = 12;
framelist = 1;
n_frame = 1;

% save_path = "D:\Spatiotemporal_analysis\initiation_statistics\" + prefix + "\sp_mat_files\";
save_path = "D:\Spatiotemporal_analysis\initiation_statistics\perturb-photo\sp_mat_files\";
mkdir(save_path);

%% load spatial property files

ent = cell(length(prefix), n_wells);
coh = cell(length(prefix), n_wells);
den = cell(length(prefix), n_wells);

for p = 1:length(prefix)
    disp(prefix(p));
    path = ROOTDIR + "/wave_" + prefix(p) + "/";
    % load entire image and initiations
    for w = welllist{p}
        for fr = 1
            ent{p, w} = load_image('ent', '.mat', prefix(p), subpath('ent'), w, fr, ch('ent'));
            coh{p, w} = load_image('coh', '.mat', prefix(p), subpath('coh'), w, fr, ch('coh'));
            den{p, w} = load_image('den120', '_5000.mat', prefix(p), subpath('den'), w, fr, ch('den'));
        end
    end
end


%% load initiation from initiation
savefile = 1;

% sp_prefix = "yapi-0816";
all_ftype = ["randbg2" "randbg3" "resist2" "resist3"];
for ftype = all_ftype
    % ftype = "resist1";
    
    ini_file = "D:\Spatiotemporal_analysis\Initiation\" + prefix + "_"+ftype+"_details_pattern_per5_cov10.txt";
    outfile = "D:\Spatiotemporal_analysis\Initiation\" + prefix + "_"+ftype+"_details_pattern_per5_cov10_sp"+ini_size+".txt";
    
    
    ini = readtable(ini_file, 'Delimiter', '\t');
    ini.id = string(ini.id);
    
    [ini_ent, ini_coh, ini_den] = get_sp_for_ini(prefix, ini, imsize, ini_size, ent, coh, den, offset);
    ini = store_sp_for_ini(ini, ini_ent, ini_coh, ini_den);
    
    if savefile
        disp("saving");
        save(save_path+prefix+"_"+ftype+"_ent.mat", 'ini_ent', '-v7.3');
        save(save_path+prefix+"_"+ftype+"_coh.mat", 'ini_coh', '-v7.3');
        save(save_path+prefix+"_"+ftype+"_den.mat", 'ini_den', '-v7.3');
        
        writetable(ini, outfile, 'Delimiter', '\t')
    end
end
%% load initiation from spontaneous initiation
% welllist = {[1:24], [1:23], [1:20], [1:11 13:20]};
% n_wells = 24;
% 
% framelist = load_framelist(prefix, welllist, n_wells);
% n_frame = 65;
% 
% imsize = load_imsize(prefix, n_wells);
% 
% s_ini_file = "D:/Spatiotemporal_analysis/Initiation/spontaneous_initiation_details_pattern_per5_cov10.txt";
% s_ini = readtable(s_ini_file, 'Delimiter', '\t');
% s_ini.id = string(s_ini.id);

% %% load initiation from photoinduction
% welllist = {[1 4 5 8 9 12], [1:3 5 7:10 12], 1:12, 1:12};
% n_wells = 12;
% 
% framelist = num2cell(ones(length(prefix), n_wells));
% n_frame = 1;
% 
% imsize = load_imsize(prefix, n_wells);
% 
% p_ini_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_initiation_details_pattern_per5_cov10.txt";
% p_ini = readtable(p_ini_file, 'Delimiter', '\t');
% p_ini.id = string(p_ini.id);
%
% % all initiations
% all_ini_ent = [ini_ent; s_ini_ent];
% all_ini_coh = [ini_coh; s_ini_coh];
% all_ini_den = [ini_den; s_ini_den];

%% FUNCTIONS
function df = store_sp_for_ini(df, sp_ent, sp_coh, sp_den)
    % store the mean spatial properties into the table
    mean_ent = nan(height(df), 1);
    mean_coh = nan(height(df), 1);
    mean_den = nan(height(df), 1);
    
    for i = 1:height(df)
        mean_ent(i) = mean(sp_ent{i}(:), 'omitnan');
        mean_coh(i) = mean(sp_coh{i}(:), 'omitnan');
        mean_den(i) = mean(sp_den{i}(:), 'omitnan');
    end
    
    df.ent = mean_ent;
    df.coh = mean_coh;
    df.den = mean_den;
end
                                                
function [sp_ent, sp_coh, sp_den] = get_sp_for_ini(prefix, sp, imsize, ini_size, ent, coh, den, offset)
    % load spatial properties for each initiation
    % sp may be spontaneous, photoinduction, background or resistant
    sp_ent = cell(height(sp),1);
    sp_coh = cell(height(sp),1);
    sp_den = cell(height(sp),1);
    
    for i = 1:height(sp)
        p = find(prefix==sp.id(i));
        w = sp.well(i);
        x = ceil(sp.x(i));
        y = ceil(sp.y(i));

        sz = imsize{p};

        % check if initiation occurs near edge. If so, skip the initiation
        if x-ini_size/2 < offset || y-ini_size/2 < offset || x+ini_size/2-1 > sz(1)-offset || y+ini_size/2-1 > sz(2)-offset
            continue;
        end
    
        if isempty(den{p,w})
            continue;
        end
    
        % sp_ent{i} = crop_circle(ent{p,w}, x, y, ini_size/2);
        % sp_coh{i} = crop_circle(coh{p,w}, x, y, ini_size/2);
        % sp_den{i} = crop_circle(den{p,w}, x, y, ini_size/2);
    
        sp_ent{i} = imcrop(ent{p,w}, [x y ini_size/2 ini_size/2]);
        sp_coh{i} = imcrop(coh{p,w}, [x y ini_size/2 ini_size/2]);
        sp_den{i} = imcrop(den{p,w}, [x y ini_size/2 ini_size/2]);
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