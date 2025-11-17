%
%   Match cellular pattern topologies with labeled initiation
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
pattern = ["phalf" "mhalf" "splay" "bend" "aligned"];
sp = ["phalf" "mhalf" "splay" "bend" "ent"];
[ch, subpath] = load_file_paths();

% spontaneous initiations
% prefix = ["mCh300FBS10-1", "mCh300FBS10-2" "overlap5-1" "overlap5-2"];
% 
% welllist = {1:24, 1:23, 1:20, [1:11 13:20]};
% n_wells = 24;
% 
% framelist = load_framelist(prefix, welllist, n_wells);
% n_frame = 65;

% photoinduction
% prefix = ["photoinduct-1" "photoinduct-2" "photoinduct-3" "photoinduct-4"];
%
% welllist = {[1 4 5 8 9 12], [1:3 5 7:10 12], 1:12, 1:12};
% n_wells = 12;
%
% framelist = num2cell(ones(length(prefix), n_wells));
% n_frame = 1;

% prefix = ["spon-perturb-4"];  % "yapi-0819" "yapi-0827z4"];
% prefix = ["yapi-0816" "yapi-0819" "yapi-0827z4"]; 
prefix = ["photo-truli-ki-0730-4"]; % "photo-truli-KI-4" "photo-truli-ki-0730-1" "photo-truli-ki-0730-4"];

% welllist = {[1 3 4 9 10 12 13 15 16 21 22 24], [1 2 3 10 11 12 13 14 15 22 23 24], [1:24]};
% n_wells = 24;
welllist = {[1:12]}; %, [1:12], [1:12] [1:12]};
n_wells = 12;

% framelist = num2cell(ones(length(prefix), n_wells));
framelist = 1;
n_frame = 1;

imsize = {[5120 5120]}; %, [5120 5120], [4800 4800] [5120 5120]}; %[4800 4800]; %[3072 3072]; %[4949 4949];

%%
% CY5CH = "c1";
% ADJPATH = CY5CH + "_cy5_adj/";
% DENPATH = CY5CH + "_cy5_density/den_regn120/";
% 
% DICCH = "c2"; % channel of the DIC file
% DICPATH = DICCH + "_DIC/";
% NEMPATH = DICCH + "_DIC_nematics/";
% ORIENTPATH = NEMPATH + "/orientation_files/";

%% load defects, splay and bend files
phalf_coord = cell(length(prefix), n_wells);
mhalf_coord = cell(length(prefix), n_wells);
splay_coord = cell(length(prefix), n_wells);
bend_coord = cell(length(prefix), n_wells);

for p = 1:length(prefix)
    disp(prefix(p));
    for w = welllist{p}
        % for fr = framelist{p,w}
        for fr = 1
            phalf_coord{p, w} = load_pattern_coord('phalf', '.txt', prefix(p), subpath("phalf"), w, fr, ch("phalf"));
            mhalf_coord{p, w} = load_pattern_coord('mhalf', '.txt', prefix(p), subpath("mhalf"), w, fr, ch("mhalf"));

            splay_coord{p, w} = load_pattern_coord('splay', '_5.mat', prefix(p), subpath("splay"), w, fr, ch("splay"));
            bend_coord{p, w} = load_pattern_coord('bend', '_5.mat', prefix(p), subpath("bend"), w, fr, ch("bend"));
        end
    end
end

%% load entropy file

ent = cell(length(prefix), n_wells);
ent_std = nan(length(prefix)*n_wells, 1);

for p = 1:length(prefix)
    disp(prefix(p));

    for w = welllist{p}
        % for fr = framelist{p,w}
        for fr = 1
           ent_file = "D:/Spatiotemporal_analysis/wave_" + prefix(p) + "/Well" + sprintf("%02d", w) + "/" +  subpath("aligned") + "/" +...
              set_filename(prefix(p), w, fr, ch("aligned")) + "_ent.mat"; 

           e = struct2array(load(ent_file));
    
           ent_std((p-1)*n_wells + w) = mean(e(:), 'omitnan') - 1*std(e(:), 'omitnan');
           ent{p,w} = e;
        end
    end
end
ent_cutoff = mean(ent_std, 'omitnan');
disp(mean(ent_std, 'omitnan'));

%% load initiation file

% for initiation file
% ini_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_initiation_details.txt";
% ini_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_initiation_details_pattern_per5_cov10.txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/photoinduct-all_initiation_details_pattern_per5_cov10.txt";

% for photoinduct-1-1 initiation file
% ini_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-1-1_initiation_details.txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/photoinduct-1-1_initiation_details_pattern_per5_cov10.txt";

% % for photoinduct-0721 initiation file
% ini_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-0721_initiation_details.txt";
% ini_file = "D:/Spatiotemporal_analysis/Initiation/photoinduct-0721_initiation_details_pattern_per5_cov10_sp150.txt"
% outfile = "D:/Spatiotemporal_analysis/Initiation/photoinduct-0721_initiation_details_pattern_per5_cov10.txt";

% % for spontaneous initiation file
% ini_file = "D:/Spatiotemporal_analysis/Initiation/spontaneous_initiation_details.txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/spontaneous_initiation_details_pattern_per5_cov10.txt";

% % for random non-initiation backgrounds
% rep = 10;
% disp("running random "+rep);
% ini_file = "D:/Spatiotemporal_analysis/Initiation/" + "photoinduct-all_randbg" + rep + "_details.txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/" + "photoinduct-all_randbg" + rep + "_details_pattern_per5_cov10.txt";

% for random non-initiation backgrounds
% rep = 15;
% disp("running resistant "+rep);
% ini_file = "D:/Spatiotemporal_analysis/Initiation/" + "photoinduct-all_resist" + rep + "_details.txt";
% outfile = "D:/Spatiotemporal_analysis/Initiation/" + "photoinduct-all_resist" + rep + "_details_pattern_per5_cov10.txt";

% % for cell migration initiation data
% ini_file = "D:\Spatiotemporal_analysis\Initiation\cellmigration_initiation_details.txt";
% outfile = "D:\Spatiotemporal_analysis\Initiation\cellmigration_initiation_details_pattern_per5_cov10.txt";

% % for random background
% ini_file = "D:\Spatiotemporal_analysis\Initiation\cellmigration_randbg2_details.txt";
% outfile = "D:\Spatiotemporal_analysis\Initiation\cellmigration_randbg2_details_pattern_per5_cov10.txt";

% for inhibitors
% ini_file = "D:\Spatiotemporal_analysis\Initiation\inhibitor-1_initiation_details_pos.txt";
% outfile = "D:\Spatiotemporal_analysis\Initiation\inhibitor-1_initiation_details_pos_pattern_per5_cov10.txt";

% YAP NRF2 perturbation
% ini_file = "D:\Spatiotemporal_analysis\Initiation\YAP-NRF2-perturb_initiation_details.txt";
% outfile = "D:\Spatiotemporal_analysis\Initiation\YAP-NRF2-perturb_initiation_details_pattern_per5_cov10.txt";

% ini_file = "D:\Spatiotemporal_analysis\Initiation\Z1-YAP-NRF2-perturb_initiation_details.txt";
% outfile = "D:\Spatiotemporal_analysis\Initiation\Z1-YAP-NRF2-perturb_initiation_details_pattern_per5_cov10.txt";

% ["photo-truli-KI-1" "photo-truli-KI-4" "photo-truli-ki-0730-1" "photo-truli-ki-0730-4"];

all_ftype = ["randbg2" "randbg3" "resist2" "resist3"];

for ftype = all_ftype
% ftype = "resist3";
ini_file = "D:\Spatiotemporal_analysis\Initiation\" + prefix + "_" + ftype + "_details.txt";
outfile  = "D:\Spatiotemporal_analysis\Initiation\" + prefix + "_" + ftype + "_details_pattern_per5_cov10.txt";

disp(prefix + "_" + ftype)

% ini_file_type = "initiations";
% ini_file = "D:\Spatiotemporal_analysis\Initiation\" + prefix + "_" + ini_file_type + "_details.txt";
% outfile = "D:\Spatiotemporal_analysis\Initiation\" + prefix + "_" + ini_file_type + "_details_pattern_per5_cov10.txt";

% sp_prefix = "yapi-0827z4";

% ini_file = "D:\Spatiotemporal_analysis\Initiation\" + sp_prefix + "_initiation_details.txt";
% outfile =  "D:\Spatiotemporal_analysis\Initiation\" + sp_prefix + "_initiation_details_pattern_per5_cov10.txt";

% ini_file = "D:\Spatiotemporal_analysis\Initiation\" + sp_prefix + "_resist1_details.txt";
% outfile =  "D:\Spatiotemporal_analysis\Initiation\" + sp_prefix + "_resist1_details_pattern_per5_cov10.txt";

ini = readtable(ini_file, 'Delimiter', '\t');
ini.id = string(ini.id);

%% match initiations to splay and bend
area_cutoff = 0.1;
% ent_cutoff = -2.5617; %-2.7731 % YAP perturb-wave -2.5721
ent_coutoff = -2.5617;


ini_size = 150;
offset = 51;

all_coord = {phalf_coord, mhalf_coord, splay_coord, bend_coord};

for i = 1:length(pattern)
% for i = length(pattern)
    disp(pattern(i));

    if i ~= length(pattern)
        sp_coord = all_coord{i};
    end
    
    colname = pattern(i) + "_" + string(ini_size);

    is_pattern = zeros(height(ini),1);
    
    for j = 1:height(ini)
        p = find(prefix==ini.id(i));
    
        w = ini.well(j);
    
        % skip empty well data
        if isempty(sp_coord{p,w})
            continue;
        end

        % define initation center
        x = ceil(ini.x(j));
        y = ceil(ini.y(j));

        sz = imsize{p};

        if x-ini_size/2 < offset || y-ini_size/2 < offset || x+ini_size/2-1 > sz(1)-offset || y+ini_size/2-1 > sz(2)-offset
            continue;
        end

        % if the initiation window contains the spatial patterns
        % within an area for splay/bend, and a single point for defect
        % else if it has the small mean entropy for aligned regions
        if contains(pattern(i), "half")
            in = in_circle(x, y, ini_size, sp_coord{p,w}(:,1), sp_coord{p,w}(:,2));
            if in > 0
                is_pattern(j) = i;
            end
        elseif strcmp(pattern(i), 'splay') || strcmp(pattern(i), 'bend')
            in = in_circle(x, y, ini_size, sp_coord{p,w}(:,1), sp_coord{p,w}(:,2));
            if in/(pi*(ini_size/2)^2) > area_cutoff
                is_pattern(j) = i;
            end
        elseif strcmp(pattern(i), 'aligned')
            % e = crop_circle(ent{p,w}, x, y, ini_size/2);
            e = imcrop(ent{p,w}, [x y ini_size/2 ini_size/2]);
            if mean(e(:), 'omitnan') <= ent_cutoff
                is_pattern(j) = i;
            end
        end
    end
    
    ini.(colname) = is_pattern;
end


%%
ini.pattern = strcat(string(ini.phalf_150), string(ini.mhalf_150), string(ini.splay_150), string(ini.bend_150));
ini.pattern((ini.aligned_150 == 5) & (ini.splay_150 == 0) & (ini.bend_150 == 0)) = "5";
ini.pattern = regexprep(ini.pattern, '0', '');
ini.pattern(strcmp(ini.pattern, "")) = "0";
ini.pattern = str2double(ini.pattern);
ini = movevars(ini, 'pattern', 'Before', 13);
%% Write to file
disp("writing to file...")
writetable(ini, outfile, 'Delimiter', '\t');
end

%% FUNCTIONS
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

function in = in_circle(x0, y0, sz, xq, yq)
    % draw circle
    angle = 0:0.01:2*pi;
    xv = sz/2*cos(angle) + x0;
    yv = sz/2*sin(angle) + y0;

    % crop the values to a squared based on input size to boost up speed
    xcoord = (xq>=(x0-sz/2)) & (xq<=(x0+sz/2));
    ycoord = (yq>=(y0-sz/2)) & (yq<=(y0+sz/2)); 
    xycoord = xcoord & ycoord;

    % select only the values in the square
    xq = xq(xycoord);
    yq = yq(xycoord);

    % calculate how many spatial properties are in the circle
    in = sum(inpolygon(xq, yq, xv, yv));
end

function pat_coord = load_pattern_coord(spatial_type, suffix, prefix, subpath, w, fr, ch)
    pfile = "D:/Spatiotemporal_analysis/wave_" + prefix + "/Well" + sprintf("%02d", w) + "/" + subpath + "/" + ...
          set_filename(prefix, w, fr, ch) + "_" + spatial_type + suffix;

    if contains(suffix, "mat")
        pat_coord = struct2array(load(pfile));

        [ycoord, xcoord] = find(pat_coord == 1);
        pat_coord = [xcoord ycoord];
    elseif contains(suffix, "txt")
        pat_coord = struct2array(tdfread(pfile));
    end
end
