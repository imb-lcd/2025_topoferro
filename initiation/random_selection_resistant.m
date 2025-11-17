%
% randomly select resistant region
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
%ROOTDIR = "/home/N417/Jen-Hao/Spatiotemporal/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath(ROOTDIR + "/code");
addpath(ROOTDIR + "/code/util");

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATIONs
rep = 3;
ini_size = 150;
IMSIZE = [5120 5120]; %[5120 5120]; %[4800 4800]; % [5000 5000];

prefix = "photo-truli-ki-0730-4"; %"yapi-0827z4";

ini_file = "D:/Spatiotemporal_analysis/Initiation/" + prefix + "_initiation_details.txt";
outfile = "D:/Spatiotemporal_analysis/Initiation/" + prefix + "_resist" + rep + "_details.txt";

welllist = {[1:12]};
n_wells = 12;
fr = 45;

lg_wv_idx = {[2], [2], [2], [2], [2], [2], [2], [2], [2], [2], [2], [2]};

%% Load the images

cy5 = cell(length(prefix), n_wells);

for p = 1:length(prefix)
    disp(prefix(p));
    for w = welllist{p}
        well = sprintf("%02d", w);

        path = "D:\Spatiotemporal_analysis\wave_" + prefix(p) + "\Well" + well + "\c2_cy5_adj\";
        
%         fname = path + prefix + "_s" + we04t25c2_ORG_bw.tif";
        fname = path + set_filename(prefix(p), w, fr, 'c2') + "_bw.tif";
        
        % get the actual filename without paths
        [~,name,~] = fileparts(fname);
        
        % import multi-page tif file into a cell array
        info = imfinfo(fname);
        num_images = numel(info); % number of array elements
        
        cy5{p, w} = imread(fname, 1, 'Info', info); %store as array of objects
    end
end


%% Outline dead cells population
%close all
sm_nucl_thres = 15; % 20
lg_nucl_thres = 1000;
sm_area_thres = 3500; %50000;
sm_hole_thres = 50000;
sm_bo = 25000; %50000;

dilate_sz = 20;
close_sz = 20; %15;
erode_sz = dilate_sz - 15;

% bdy keeps the mask for dead cell population
bdy = cell(length(prefix), n_wells);
bdy_unsel = cell(length(prefix), n_wells);

% for each image, outline the wave
for p = 1:length(prefix)
    disp(prefix(p));
    for w = welllist{p}

        I = cy5{p, w};
    
        % remove any labels with area > lg_nucl_thres, which is assumed to be debris
        L = bwlabel(I);
        obj = regionprops('table', L, 'Area'); % calculate the area for each component
        if any(obj.Area > lg_nucl_thres)
            I(ismember(L,find(obj.Area > lg_nucl_thres))) = 0; % rm by setting to 0
        end
    
        % remove any labels with area < sm_nucl_thres, which is assumed to be debris
        L = bwlabel(I);
        obj = regionprops('table', L, 'Area'); % calculate the area for each component
        if any(obj.Area < sm_nucl_thres)
            I(ismember(L,find(obj.Area < sm_nucl_thres))) = 0; % rm by setting to 0
        end
    
        % Dilate, Close, Rm sm obj, Rm holes, Erode, Rm sm obj
        BW = imdilate(I, strel('disk', dilate_sz));   % Dilate to connect the dead cells
        BW = bwareaopen(BW, sm_area_thres);         % Remove small objects and speckles
        BW = imerode(BW, strel('disk', erode_sz));  % Erode the image back
    
        % remove any wave area < sm_bo, which is assumed to be debris
        L = bwlabel(BW);
        obj = regionprops('table', L, 'Area'); % calculate the area for each component
        if any(obj.Area < sm_bo)
            BW(ismember(L,find(obj.Area < sm_bo))) = 0; % rm by setting to 0
        end
    
        % store the outline in bdy
        bdy{p, w} = sparse(BW);
    
        BW = imdilate(BW, strel('disk', 75)); 
        bdy_unsel{p, w} = sparse(BW);
    end
end

%% Determine unselectable initiation region
close all
showfigure = 0;

offset = 51;
border = ini_size/2+offset;

mask = cell(length(prefix), n_wells);

[xgrid, ygrid] = meshgrid(1:IMSIZE(2), 1:IMSIZE(1));

for p = 1:length(prefix)
    for w = welllist{p}
        if w==2 || w==5 || w==8 || w==11
            continue
        end

        mask{p,w} = false(IMSIZE);

        Bo_unsel = bwboundaries(full(bdy_unsel{p, w}));
        % continue;
        % run boundary selection
        if isempty(Bo_unsel)
            mask{p, w} = false(IMSIZE);
        else
            for i = 1:(lg_wv_idx{p,w}-1)
                mask{p,w} = mask{p,w} | poly2mask(Bo_unsel{i}(:,2), Bo_unsel{i}(:,1), IMSIZE(2), IMSIZE(1));
            end
            for j = lg_wv_idx{p,w}:length(Bo_unsel)
                mask{p,w} = mask{p,w} & ~poly2mask(Bo_unsel{j}(:,2), Bo_unsel{j}(:,1), IMSIZE(2), IMSIZE(1));
            end
        end
    
        % add a border of unselectable region to the mask
        mask{p,w}(1:border, :) = 1;
        mask{p,w}(end-border:end, :) = 1;
        mask{p,w}(:, 1:border) = 1;
        mask{p,w}(:, end-border:end) = 1;

        if showfigure
            path = "D:\Spatiotemporal_analysis\wave_" + prefix(p) + "\Well" + sprintf("%02d", w) + "\c2_cy5_adj\";
            fname = path + set_filename(prefix(p), w, fr, 'c2') + "_bw.tif";
            I = imread(fname);
            figure
            nexttile
            imshow(cy5{p,w}, 'Border', 'tight')
            title(prefix(p) + " Well" + sprintf("%02d", w))
            nexttile
            imshow(mask{p,w})    
        end
    end
end


%%
ini = readtable(ini_file);
ini.id = string(ini.id);

%% randomly select 150 resistant regions
close all
showfigure = 0;

ini_size = 150;

res = cell(length(prefix), n_wells);

for p = 1:length(prefix)
    for w = welllist{p}

        ini_pw = ini(strcmp(ini.id, prefix(p)) & ini.well == w, :);
        ini_exp = string(cell2mat(unique(ini_pw.experiment)));

        m = ~mask{p,w};

        n = ceil(sum(m(:))/15000);
        if n > 150
            n = 150;
        end
        
        rand_x{w} = NaN(1, n);
        rand_y{w} = NaN(1, n);
        
        for r = 1:n
            while true % loop until it finds a region not in the mask
                rect = randomWindow2d([IMSIZE(2) IMSIZE(1)], [ini_size ini_size]);
               
                % note: need to flip the xycoord when comparing to mask
                if m(rect.YLimits(2)-ini_size/2, rect.XLimits(2)-ini_size/2) ~= 0
                    rand_x{w}(r) = rect.XLimits(2)-ini_size/2;
                    rand_y{w}(r) = rect.YLimits(2)-ini_size/2;
                    
                    % generate additional mask per initation for lower
                    % degree of overlap
                    xlimits = [rect.XLimits(1)-25 rect.XLimits(1)-25 rect.XLimits(2)+25 rect.XLimits(2)+25];
                    ylimits = [rect.YLimits(1)-25 rect.YLimits(2)+25 rect.YLimits(2)+25 rect.YLimits(1)-25];

                    m = m & ~poly2mask(xlimits, ylimits, IMSIZE(2), IMSIZE(1));

                    mask{p,w} = ~m;

                    break;
                end
            end
        
            % save the randomly selected resistant region in the format
            % id, experiment, well, frame, index, x, y
            res{p, w} = [repelem(prefix(p), n)' repelem(ini_exp, n)' repelem(w, n)' repelem(1, n)' (1:n)' rand_x{w}' rand_y{w}'];
        end
        
        if showfigure
            % figure
            imshow(mask{p,w}, 'Border', 'tight')
            hold on
            for k = 1:n
                drawrectangle('Position', [rand_x{w}(k)-ini_size/2 rand_y{w}(k)-ini_size/2 ini_size ini_size], 'Color', [1 0 0]);
            end
            Bo = bwboundaries(full(bdy{p,w}));
            for j = 1:length(Bo)
                plot(Bo{j}(:,2), Bo{j}(:,1), 'b-', 'LineWidth', 1.5);
            end
        
            hold off
            close all
        end
        
    end
end

%% write the resistant coordinates to file
res = vertcat(res{:});
res = array2table(res, 'VariableNames', {'id', 'experiment', 'well', 'frame', 'index', 'x', 'y'});

writetable(res, outfile, 'Delimiter', '\t');