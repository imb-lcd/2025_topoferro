%% 
% calculate_streamline.m
% calculates the streamlines based on binarized data cy5 data
% requires the Mapping ToolBox
%
clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath(ROOTDIR + "/code");
addpath(ROOTDIR + "/code/util");

%% CONFIGURATIONS
prefix = "AA-0916";

% assign frame range
start_frame = 14;
end_frame = 20;
frame_range = sprintf("%02d", start_frame) + "-" + sprintf("%02d", end_frame);
num_frame = end_frame-start_frame+1;

% assign well
w = 93;
well = sprintf("%02d", w);

% designate the wave number used in this well
wv = "2";

IMSIZE = 2900; % assign image size. Assumed square
% IMSIZE = 2968; %4977; %3072;

% CONSTANTS
PATH = ROOTDIR + "wave_" + prefix + "/";

den_reg = 120;
MCHCH = "c1";
DENPATH = MCHCH + "_mCherry_density/";

CY5CH = "c2";
CY5ADJPATH = CY5CH + "_cy5_adj/";

DICCH = "c3";
NEMPATH = DICCH + "_DIC_nematics/";
ORIENTPATH = NEMPATH + "orientation_files/";

%% Load the image
fname = "D:/Spatiotemporal_analysis/wave_" + prefix + "/Well" + well + "/wave" + wv + "_analyses/" + ...
   set_filename(prefix, w, frame_range, CY5CH) + "_bw_wave" + wv + ".tif";

% get the actual filename without paths
[~,name,~] = fileparts(fname);

% import multi-page tif file into a cell array
info = imfinfo(fname);
num_images = numel(info); % number of array elements

cy5 = cell(num_images, 1);
prev_cy5 = zeros(IMSIZE, IMSIZE);
for i = 1:num_images
    cy5{i} = imread(fname, i, 'Info', info); %store as array of objects
    temp = cy5{i};

    temp = temp(1:IMSIZE, 1:IMSIZE);
    cy5{i} = double(temp) + prev_cy5;

    prev_cy5 = cy5{i};
end

imsize = size(cy5{1});
IMSIZE = imsize(1);

%% Outline dead cells population
sm_nucl_thres = 5; %15 20
lg_nucl_thres = 500;
sm_area_thres = 30000;
sm_hole_thres = 50000;

dilate_sz = 20;
close_sz = 20; %15;
erode_sz = dilate_sz - 15;
medfilt_sz = 30; %25;

% bdy keeps the mask for dead cell population
bdy = cell(num_images, 1);
bdy_smooth = cell(num_images, 1);

% for each image, outline the wave
for i = 1:num_images
    if i <= 2
        sm_nucl_thres = 1;
        sm_area_thres = 3000;
    else
        sm_nucl_thres = 5;
        sm_area_thres = 30000;
    end

    I = cy5{i};

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
    BW = imclose(BW, strel('disk', close_sz));  % Close the gaps in the dead cells
    BW = bwareaopen(BW, sm_area_thres);         % Remove small objects and speckles
    BW = ~bwareaopen(~BW, sm_hole_thres);       % Remove holes
    BW = medfilt2(BW, [medfilt_sz medfilt_sz]);
    
    % store the outline in bdy
    bdy{i} = sparse(BW);
end

%% construct direction components based on the all frames of interest
disp('constructing direction components')
tic
vector_per_side = 250;

foi = find(~cellfun(@(x) isempty(bwboundaries(full(x))), bdy)); % list the non-empty arrays

C = bdy(foi);
filteredImage = zeros(IMSIZE, IMSIZE, length(C));
for i = 1:1:length(C)
    filteredImage(:,:,i) = full(C{i});  % assign each non-empty array
end

% store the total number of frame after removing empty frames
nframe = length(C);

% Get the outermost boundary
wave_bound = sum(filteredImage, 3);
wave_bound(wave_bound > 0) = 1;

Z = length(C)+1-sum(filteredImage,3);   % construct contours (Z dim)

% Re-scale dimensions to vector_per_side x vector_per_side 
[X, Y] = meshgrid(1:1:IMSIZE);  % original image size

ws = IMSIZE/vector_per_side;    % scale factor for 200 vectors in one dimension
[Xq, Yq] = meshgrid(1:ws:IMSIZE);
% [Xq, Yq] = meshgrid(1:75:IMSIZE);

% interpolate vectors to new dimensions
Vq = interp2(X, Y, Z, Xq, Yq);  
Vq(Vq == max(Vq(:))) = NaN;  

% % generate the vector field
[Gmag, Gdir] = imgradient(Vq, 'central');
Gmag(Gmag==0) = NaN;
invGmag = 1./Gmag;  % inverse magnitude to become larger mag for steeper slope

% % obtain original direction component
[UG, VG] = obtain_direction_component(Gmag, Gdir, 3.5);
% [UG, VG] = obtain_direction_component(invGmag, Gdir, 1);

% set limits for plots later
[xlim1, ylim1, xlim2,  ylim2] = find_plot_xylim(Xq, Yq, Vq);

toc

% Draw vector field or streamlines
figure
imshow(Z, [0 length(C)+1], 'Border', 'tight')
axis square
xlim([xlim1 xlim2])
ylim([ylim1 ylim2])

% hold on
% % Draw vector fields
% quiver(Xq, Yq, UG, VG, 1, 'color', 'k')

% Draw streamlines
% lineobj = streamslice(Xq, Yq, UG, VG, 30, 'noarrows');
% lineobj = streamslice(Xq, Yq, UG, VG, streamline_density, 'noarrows');
% set(lineobj, 'Color', 'r', 'LineWidth', 0.75)
% hold off

%% visualize the wave outline
figure
tiledlayout('flow','Padding','none','TileSpacing','none')
for i = 1:num_images
    nexttile
    imshow(cy5{i},[0 255])
    % xlim([xlim1 xlim2])
    % ylim([ylim1 ylim2])
    
    hold on
    Bo = bwboundaries(full(bdy{i}));
    for j = 1:length(Bo)
        hold on
        if i == num_images
            plot(Bo{j}(:,2), Bo{j}(:,1), '-', 'Color', "#EDB120", 'LineWidth', 4);
        else
            plot(Bo{j}(:,2), Bo{j}(:,1), '-', 'Color', "#EDB120",  'LineWidth', 2);
        end
    end
    hold off
end

%% save wave boundaries as svg files
outpath = "D:/Spatiotemporal_analysis/wave_" + prefix + "/Well" + well + "/wave" + wv + "_analyses/";

figure
imshow(zeros(size(cy5{i}))+255, 'Border', 'tight')
hold on
for i = 1:num_images
    disp(i)
    Bo = bwboundaries(full(bdy{i}));
    for j = 1:length(Bo)
        plot(Bo{j}(:,2), Bo{j}(:,1), '-', 'Color', "k",  'LineWidth', 2);
    end
end
hold off
outname = set_filename(prefix, w, frame_range, CY5CH) + "_bw_wave" + wv + ".pdf";
saveas(gcf, outpath+outname);
%% erode the boundaries to ensure there streamlines would intersect with boundaries
disp("calcualte the wave boundaries");

tic
erode_c = cellfun(@(x) imerode(full(x), strel('disk', erode_sz)), C, 'un', false);

bo = cellfun(@(x) bwboundaries(x), erode_c, 'un', false);

% select for the largest area as the boundary
for i = 1:length(bo)
    area = 0;
    tmp_Bo = [];
    for j = 1:length(bo{i})
        new_area = polyarea(int16(bo{i}{j}(:,1)), int16(bo{i}{j}(:,2)));
        if area < new_area
            area = new_area;
            tmp_Bo = bo{i}{j};
        end
    end
    bo{i} = tmp_Bo;
end

% store the wave boundaries
bo_path = PATH + "Well" + well + "/wave" + wv + "_analyses/";
bo_name = bo_path + set_filename(prefix, w, frame_range, "") + "_bo.mat";
save(bo_name, 'bo');
toc

%% Obtain starting points for reversed streamlines
% Lines of streamslice are not continuous from initiation to boundary. 
% Obtain starting points as the intersection between forward streamlines
% wave boundaries to obtain start points for full streamlines

disp('obtain starting points for reversed streamline')

% Determine the density of streamlines for the wave. 
% Important step to get an evenly distributed density
perim = bwperim(wave_bound);
streamline_density = round(sum(perim(:))/100);
% streamline_density = round(sum(perim(:))/2000);

% streamlines calculated right up to the wave boundary
intersect_bound = imerode(wave_bound, strel('disk', erode_sz));
Bo_all = bwboundaries(intersect_bound);
Bo_all = Bo_all{1}; % Assume no holes in the boundary

% streamlines calculated to a dilated wave boundary (to extpl beyond bound)
extpl_bound = imdilate(wave_bound, strel('disk', dilate_sz*1));
Bo_extpl = bwboundaries(extpl_bound);
Bo_extpl = Bo_extpl{1}; % Assume no holes in the boundary

% calculate streamlines.
tic
[verts, averts] = streamslice(Xq, Yq, UG, VG, streamline_density, 'noarrows');
voi = find(~cellfun(@isempty, verts)); % list the non-empty arrays
verts = verts(voi);

% plot outermost boundary
% figure
% imshow(intersect_bound)
% axis square
% hold on 
% plot(Bo_all(:,2),Bo_all(:,1), 'r-', 'LineWidth', 1);
% lineobj = streamslice(Xq, Yq, UG, VG, streamline_density, 'noarrows');
% set(lineobj, 'Color', 'b', 'LineWidth', 0.75)
% hold off

% Obtain starting point by intersections of streamslince with outermost boundary
wv_extpl_intrsct = cell(length(verts),1); %     wv_bo_intrsct = cell(length(verts),1);
wv_bo_intrsct = cell(length(verts),1);
Bo_extpl_x = Bo_extpl(:,2); %   Bo_x = Bo_all(:,2);
Bo_extpl_y = Bo_extpl(:,1); %   Bo_y = Bo_all(:,1);
Bo_x = Bo_all(:,2);
Bo_y = Bo_all(:,1);
parfor i = 1:length(verts)
    XData = verts{i}(:,1);
    YData = verts{i}(:,2);

    [xi, yi] = polyxpoly(Bo_extpl_x, Bo_extpl_y, XData, YData); %   [xi, yi] = polyxpoly(Bo_x, Bo_y, XData, YData);
    wv_extpl_intrsct{i} = [xi yi]; %     wv_bo_intrsct{i} = [xi yi];

    [xi, yi] = polyxpoly(Bo_x, Bo_y, XData, YData);
    wv_bo_intrsct{i} = [xi yi];

    % for plotting intersects with a circle
    % mapshow(xi,yi,'DisplayType','point','Marker','o')
end
wv_extpl_intrsct = vertcat(wv_extpl_intrsct{:}); % wv_bo_intrsct = vertcat(wv_bo_intrsct{:});
wv_bo_intrsct = vertcat(wv_bo_intrsct{:});
toc

%% Calculate streamlines
disp('calculating streamlines')

% reverse the streamlines and block out non-dead cells
[UG, VG] = obtain_direction_component(Gmag, Gdir-180, 3.5);
% [UG, VG] = obtain_direction_component(Gmag, Gdir-180, 2);

% store the streamlines without plotting them
tic
XY = stream2(Xq, Yq, UG, VG, wv_bo_intrsct(:,1), wv_bo_intrsct(:,2)); % calculated for plotting

XYextpl = stream2(Xq, Yq, UG, VG, wv_extpl_intrsct(:,1), wv_extpl_intrsct(:,2)); 
toc
 
figure
imshow(zeros(IMSIZE)+255,[255 255], 'Border', 'tight')
xlim([xlim1 xlim2])
ylim([ylim1 ylim2])

hold on

plot(Bo_all(:,2),Bo_all(:,1), 'k-', 'LineWidth', 1);
for j = 1:length(bo) %:-2:1
    Bo = bo{j};
    plot(Bo(:,2), Bo(:,1), 'Color', 'k', 'LineWidth', 1);
end

lineobj = streamline(Xq, Yq, UG, VG, wv_extpl_intrsct(:,1), wv_extpl_intrsct(:,2)); % %lineobj = streamline(Xq, Yq, UG, VG, wv_bo_intrsct(:,1), wv_bo_intrsct(:,2));
set(lineobj, 'Color', 'r', 'LineWidth', 2)

for i = 1:length(XYextpl)
    text(XYextpl{i}(1,1), XYextpl{i}(1,2), num2str(i));
end

hold off

%% Process streamlines for extrapolated XY
disp('Processing streamlines')

tic
% % remove streamlines with NaN (also do not originate from center)
xy_no_nan = find(~cellfun(@(x) any(isnan(x(:))), XYextpl));
XYextpl = XYextpl(xy_no_nan);

% find initiation point. Initiation point is the mean of
% (minX, maxX, minY, maxY) of the unconverged region
ends = cellfun(@(x) x(end,:), XYextpl, 'un', false);
ends = reshape(cell2mat(ends), 2, [])';
[minX, maxX, minY, maxY]  = deal(min(ends(:,1)), max(ends(:,1)), min(ends(:,2)), max(ends(:,2)));
ini = [mean([minX maxX]), mean([minY maxY])];
r = max(maxX-ini(1), maxY-ini(2))+1.25;
xv = r*cos(0:0.1:2*pi)+ini(1);
yv = r*sin(0:0.1:2*pi)+ini(2);

% remove unconverged streamline ends
XYextpl = cellfun(@(x) x(~inpolygon(x(:,1), x(:,2), xv, yv),:), XYextpl, 'un', false);
toc

figure
imshow(zeros(IMSIZE)+255,[255 255], 'Border', 'tight')
xlim([xlim1 xlim2])
ylim([ylim1 ylim2])

hold on

plot(Bo_all(:,2),Bo_all(:,1), 'k-', 'LineWidth', 1);
for j = 1:length(bo) %:-2:1
    Bo = bo{j};
    plot(Bo(:,2), Bo(:,1), 'Color', 'k', 'LineWidth', 1);
end

for i = 1:length(XYextpl)
    plot(XYextpl{i}(:,1), XYextpl{i}(:,2), 'r', 'LineWidth', 2)
    text(XYextpl{i}(1,1), XYextpl{i}(1,2), num2str(i));
end

scatter(ini(1), ini(2))

hold off


%% finding intersection points between streamline and wave boundaries
disp('calculating intersection points')

% calculate intersection points
tic
intx = cell(length(C), 1);
inty = cell(length(C), 1);
parfor i = 1:length(bo)
    [intx{i}, inty{i}] = cellfun(@(x) polyxpoly(bo{i}(:,2), bo{i}(:,1), x(:,1), x(:,2)), XYextpl, 'un', false);
end

% change the data structure of intersection points to fit other data
stream_intrsct_extpl = cell(length(intx{1}), 1);
for j = 1:length(intx{1})
    stream_intrsct_extpl{j} = NaN(length(intx), 2);
    for k = 1:length(intx)
        intxy = horzcat(cell2mat(intx{k}(j)), cell2mat(inty{k}(j)));
        if ~isempty(intxy)
            stream_intrsct_extpl{j}(k, :) = intxy(1,:); % '1' is the furtherest intersect, 'end' is the nearest intersect
        end
    end
end
toc

%% Calculate stream lengths as speed and frame
disp('Calculating speed and frame as streamline lengths')
tic
frame_extpl = cell(length(XYextpl), 1);
stream_len_extpl = cell(length(XYextpl), 1);
for i = 1:length(XYextpl)
    % declare variables 
    frame_extpl{i} = NaN(length(XYextpl{i}), 1);
    stream_len_extpl{i} = NaN(length(XYextpl{i}), 1);

    % assign start to be end of the streamline (initiation)
    int_idx1 = length(XYextpl{i});
    for j = 1:length(stream_intrsct_extpl{i})

        % find the intersection point
        int_idx2 = dsearchn(XYextpl{i}, stream_intrsct_extpl{i}(j, :)); 
        %disp(i + " "+fr+" "+int_idx1+" "+int_idx2)

        % assign length between two intersections
        if int_idx1 == int_idx2
            stream_len_extpl{i}(int_idx2:int_idx1) = 0;
        else
            % calculate distance of the streamline
            d = diff(XYextpl{i}(int_idx2:int_idx1, :));
            stream_len_extpl{i}(int_idx2:int_idx1) = sum(sqrt(sum(d.*d,2)));
        end
        % assign corresponding frame
        frame_extpl{i}(int_idx2:int_idx1) = j;
    
        % if already at wave boundary
        if int_idx2 == 1 || (int_idx2 ~= 1 && j == length(stream_intrsct_extpl{i}))
            break;
        end

        % shift indices
        int_idx1 = int_idx2;
    end
end
toc

%% Plot the streamlines by frames
colorOrder = ['y', 'm', 'c', 'r', 'g', 'b','y', 'm', 'c', 'r', 'g', 'b','y', 'm', 'c', 'r', 'g', 'b'];

figure
imshow(zeros(IMSIZE)+255,[255 255], 'Border', 'tight')
xlim([xlim1 xlim2])
ylim([ylim1 ylim2])

hold on
for j = 1:length(bo)
    Bo = bo{j};
    plot(Bo(:,2),Bo(:,1), 'k', 'LineWidth', 2);
end

% plot streamline with per boundary segment
for i = 1:length(XYextpl)
    for j = 1:length(bo)
        idx = find(frame_extpl{i} == j);
        if isempty(idx)
            continue;
        end
        if idx(1) ~= 1
            idx = [idx(1)-1; idx];
        end
        plot(XYextpl{i}(idx,1), XYextpl{i}(idx,2), colorOrder(j), 'LineWidth', 3)
        text(XYextpl{i}(1,1), XYextpl{i}(1,2), num2str(i))
    end
%     text(XYextpl{i}(1,1), XYextpl{i}(1,2), num2str(i))
end
hold off


%% change streamline into table
disp('change data structure of streamlines into table');

stream = cell(length(XYextpl), 1);
for i = 1:length(XYextpl)
% for i = 1
    Idx = repelem(i, length(XYextpl{i}))';
    XData = round(XYextpl{i}(:, 1),3);
    YData = round(XYextpl{i}(:, 2),3);

    slope = diff(XYextpl{i}(1:end, 2)) ./ diff(XYextpl{i}(1:end, 1));
    rad = [nan; -1*atan(slope)];
    
    stream{i} = table(Idx, XData,YData, frame_extpl{i}, stream_len_extpl{i}, rad);
    stream{i}.Properties.VariableNames = [{'Index'} {'XData'} {'YData'} {'Frame'} {'Speed'} {'WvAngle'}];
end

%% smoothen speed
disp('smoothening the speed');

sm_ws = 130;

for i = 1:length(stream)
% for i = 1
    % add the missing frames and assign the speed to 0, with spatial
    % properties equalling to nearest value of the last frame
    missing_frame = find(ismember(1:nframe, stream{i}.Frame)==0);

    for m = 1:length(missing_frame)
        fr = missing_frame(m);

        prev_fr = fr - 1;

        s = stream{i};

        row_idx = find(s.Frame == prev_fr, 1, 'first');

        new_row = s(row_idx, :);
        if isempty(new_row)
            disp(m)
            continue;
        end
        new_row.Frame = fr;
        new_row.Speed = 0;

        s = [s(1:row_idx-1, :); new_row; s(row_idx:end, :)];

        stream{i} = s;
    end

    % smoothing out the speed
    % The devised method smoothes to the last value of the streamline
    % But the built in method smoothes both edges, while one 1 end is
    % desired. To fix this, break into two smoothing segments. 1)
    % smooth by shrinking, 2) smooth to last value. Use the smooth by
    % shrinking on the begining and smooth to last towards the end

    spd = stream{i}.Speed; % retrieve the speed and flip
    spd = spd(~isnan(spd));

    if length(spd) < sm_ws
        sm_ws = ((round(length(spd)/10)-1)-1)*10;
    end
    if sm_ws == 0
        sm_ws = 10;
    end

    shrink_spd = spd(find(~isnan(spd), 1));

    sm_0 = movmean(spd, sm_ws, 'omitnan', 'Endpoints', shrink_spd); % smooth to the last value to both ends
    sm_shrink = movmean(spd, sm_ws,'omitnan', 'Endpoints', 'shrink'); % smooth by shrinking

    shrink_idx = round(find(sm_0 == sm_shrink, 1, 'first')); % find the idx where sm_0 and sm_shrink matches
    shrink_idx = round(shrink_idx);

    % smooth to the endpoints up to the shrink index. Continue until
    % value is similar to original
    while abs(sm_0(1)-shrink_spd) > 0.01
        sm_0(1:shrink_idx) = movmean(sm_0(1:shrink_idx), shrink_idx-1, 'omitnan', 'Endpoints', shrink_spd);
    end

    sm_ini_shrink = spd(end);
    sm_ini = movmean(fliplr(spd), sm_ws, 'omitnan', 'Endpoints', sm_ini_shrink); % smooth to the last value to both ends

    shrink_idx_ini = round(find(sm_ini == fliplr(sm_shrink), 1, 'first')); % find the idx where sm_0 and sm_shrink matches
    
    shrink_idx_ini = round(shrink_idx_ini);

    while abs(sm_ini(1)-sm_ini_shrink) > 0.01
        sm_ini(1:shrink_idx_ini) = movmean(sm_ini(1:shrink_idx_ini), shrink_idx_ini-1, 'omitnan', 'Endpoints', sm_ini_shrink);
    end
    sm_ini = fliplr(sm_ini);


    % sm = [sm_0(1:shrink_idx); sm_shrink(shrink_idx+1:end)];
    sm = [sm_0(1:shrink_idx); sm_shrink(shrink_idx+1:shrink_idx_ini-1); sm_ini(shrink_idx_ini:end)];

    % reconstruct smoothed array

    sm_sgolay = smoothdata(sm, 'sgolay', sm_ws);
    sm_sgolay(sm_sgolay<0) = 0;

    try
        stream{i}.SmoothSpeed = [nan(length(stream{i}.Speed)-length(sm_sgolay), 1); sm_sgolay];
        stream{i}.SmoothSpeed(1) = nan;
    catch
        stream{i}.SmoothSpeed = repelem(nan, length(stream{i}.Speed))';
    end
end

%% store the stream
disp('store the streamline mat file')
stream_path = PATH + "Well" + well + "/wave" + wv + "_analyses/";
% stream_name = stream_path + set_filename(prefix, w, frame_range, "") + "_streamline_update.mat";
stream_name = stream_path + set_filename(prefix, w, frame_range, "") + "_streamline.mat";
save(stream_name, 'stream');

%% 
%
% add spatial information to the streamline
%

%% load spatial information
disp('load spatial info');

fr = 1;

% load phi
fname = PATH + "Well" + well + "/" + ORIENTPATH + "/" + set_filename(prefix, w, fr, DICCH) + "_phi.mat";

phi = struct2array(load(fname));
phi = imresize(phi, [IMSIZE IMSIZE]);
phi = rescale(phi, -pi/2, pi/2);

% load coh
fname = PATH + "Well" + well + "/" + ORIENTPATH + "/" + set_filename(prefix, w, fr, DICCH) + "_coh.mat";
coh = struct2array(load(fname));
coh = imresize(coh, [IMSIZE IMSIZE]);

% fr = 26;

% load den
fname = PATH + "Well" + well + "/" + DENPATH + "/density_reg" + den_reg + "/" + ...   
    set_filename(prefix, w, fr, MCHCH) + "_den" + den_reg + ".mat";   
den = struct2array(load(fname));

% load stream
stream_path = PATH + "Well" + well + "/wave" + wv + "_analyses/";
stream_name = stream_path + set_filename(prefix, w, frame_range, "") + "_streamline.mat";
load(stream_name);

%% add spatial information
disp('add spatial info to streamline');

for i = 1:length(stream)

    xpos = fix(stream{i}.XData);
    ypos = fix(stream{i}.YData);

    stream{i}.Nematics = phi(sub2ind(size(phi), ypos, xpos));
    
    norm_ang = mod(stream{i}.Nematics-stream{i}.WvAngle, pi);
    stream{i}.AngDiff = min(pi-norm_ang, norm_ang);
    stream{i}.AngDiff360 = angdiff(stream{i}.Nematics, stream{i}.WvAngle);

    stream{i}.Coh = coh(sub2ind(size(phi), ypos, xpos));

    stream{i}.Density = den(sub2ind(size(den), ypos, xpos));
end

%% store the stream
stream_path = PATH + "Well" + well + "/wave" + wv + "_analyses/";
stream_name = stream_path + set_filename(prefix, w, frame_range, "") + "_streamline_table.mat";
save(stream_name, 'stream');

%% interpolate wave speed on density and along cell orientation
nbin = 40; 

wave = vertcat(stream{:});

x = wave.AngDiff;
y = wave.Density;
v = wave.SmoothSpeed;

idx = isnan(x) | isnan(y) | isnan(v);
x = x(~idx);
y = y(~idx);
v = v(~idx);

y(y>130) = 130;
y(y<30) = 30;

v(v>300) = 300;

[xq, yq] = meshgrid(linspace(0, 1.6, nbin), linspace(30, 130, nbin));

[Xq,Yq,vq] = griddata(x, y, v, xq, yq, 'nearest');

vg = imgaussfilt(vq, 5);
%
figure
s = pcolor(xq, yq, vg);
set(s, 'FaceColor', 'interp', 'EdgeColor', 'none');
xlabel('along alignment');
ylabel('density');
axis square tight
myparula = generate_black_parula(25);
colormap(myparula);
colorbar

%%
% 
% create interpolated wave data
%

%% map the speed into grids
disp('interpolate the speed');

sp_smooth = 5;

allxy = vertcat(stream{:});

x = allxy.XData;
y = allxy.YData;
v = allxy.SmoothSpeed;

[xq, yq] = meshgrid(1:1:IMSIZE);

vq = griddata(x, y, v, xq, yq, 'nearest'); % 'natural');

% vq(vq>250) = 250;
vq(vq<0) = 0;

vg = imgaussfilt(vq, sp_smooth); %10

speed = block_nondead_cells(wave_bound, IMSIZE, IMSIZE, vg);

% plot speed
figure
pcolor(xq, yq, speed, 'FaceColor', 'interp', 'EdgeColor', 'none');
xlim([xlim1 xlim2])
ylim([ylim1 ylim2])
axis square off
colormap([1 1 1; parula])

%% store the intpl speed

fname = PATH + "Well" + well + "/wave" + wv + "_analyses/" ...
    + set_filename(prefix, w, frame_range, DICCH) + "_speed_" + IMSIZE + ".mat";
save(fname, 'speed');

%% map the directions as wave angles into grids
disp("calculate wave angle as grid")

tic
x = vertcat(stream{:}).XData';
y = vertcat(stream{:}).YData';

v = vertcat(stream{:}).WvAngle';

% interpolate the wave angle as a square
[Xq,Yq] = meshgrid(1:1:IMSIZE-2);

vq = griddata(x,y,v,Xq,Yq, 'nearest'); % 'natural');

wvangle = block_nondead_cells(wave_bound, IMSIZE-2, IMSIZE, vq);
toc

% % visualize the directions
figure
im = imshow(wvangle ,[-pi/2 pi/2]);
xlim([xlim1 xlim2])
ylim([ylim1 ylim2])
colormap HSV
set(im, 'AlphaData', ~isnan(wvangle))
axis square

%% store the interpolated wave angle

fname = PATH + "Well" + well + "/wave" + wv + "_analyses/" + ... 
    set_filename(prefix, w, frame_range, DICCH) + "_wvangle.mat";
save(fname, 'wvangle');

%% calculate the interpolated angle differece
disp('calculate the interpolated angle difference');

x = allxy.XData;
y = allxy.YData;
v = vertcat(allxy.AngDiff);

tic
[xq, yq] = meshgrid(1:1:IMSIZE);

vq = griddata(x, y, v, xq, yq, 'nearest'); % 'natural');

vg = imgaussfilt(vq, sp_smooth);
angdiff = block_nondead_cells(wave_bound, IMSIZE, IMSIZE, vg);
toc

% visualize the directions
figure
im = imshow(angdiff ,[-pi/2 pi/2]);
xlim([xlim1 xlim2])
ylim([ylim1 ylim2])
colormap parula
axis square

%% calculate the interpolated angle differece
disp('calculate the interpolated angle difference for 360');

wave_path = PATH + "Well" + well + "/wave" + wv + "_analyses/";
stream_name = wave_path + set_filename(prefix, w, frame_range, "") + "_streamline_table.mat";
stream = struct2array(load(stream_name));

sp_smooth = 5;
allxy = vertcat(stream{:});

x = allxy.XData;
y = allxy.YData;
v = vertcat(allxy.AngDiff360);
tic
[xq, yq] = meshgrid(1:1:IMSIZE);

vq = griddata(x, y, v, xq, yq, 'nearest'); % 'natural');

vg = imgaussfilt(vq, sp_smooth);
angdiff360 = block_nondead_cells(wave_bound, IMSIZE, IMSIZE, vg);
toc

% visualize the directions
figure
nexttile
im = imshow(angdiff360, [-pi pi]);
xlim([xlim1 xlim2])
ylim([ylim1 ylim2])
colormap parula
axis square

nexttile
polarhistogram(angdiff360)

%
fname = PATH + "Well" + well + "/wave" + wv + "_analyses/" + ... 
    set_filename(prefix, w, frame_range, DICCH) + "_angdiff360.mat";
save(fname, 'angdiff360');

%% store the information as waves
disp('store all information as interpolated wave');

intpl_sz = 100;

rspeed = imresize(speed, [intpl_sz intpl_sz], 'nearest');
rphi = imresize(phi, [intpl_sz intpl_sz], 'nearest');
rwvangle = imresize(wvangle, [intpl_sz intpl_sz], 'nearest');
rangdiff = imresize(angdiff, [intpl_sz intpl_sz], 'nearest');
rangdiff360 = imresize(angdiff360, [intpl_sz intpl_sz], 'nearest');
rcoh = imresize(coh, [intpl_sz intpl_sz], 'nearest');
rden = imresize(den, [intpl_sz intpl_sz], 'nearest');

intpl_ws = IMSIZE/intpl_sz;
[xq,yq] = ndgrid(1:intpl_ws:IMSIZE);

wave = table(xq(:), yq(:), round(rphi(:),3), round(rwvangle(:),3), round(rangdiff(:),3), round(rangdiff360(:),3), round(rcoh(:), 3), round(rden(:),3), round(rspeed(:),3));
wave.Properties.VariableNames = [{'xq'} {'yq'} {'phi'} {'wvangle'} {'angdiff'} {'angdiff360'} {'coh'} {'density'} {'speed'}];

%  store the interpolated wave 
wave_path = PATH + "Well" + well + "/wave" + wv + "_analyses/";
wave_name = wave_path + set_filename(prefix, w, frame_range, '') + "_wvdata_intpl" + intpl_sz + ".mat";
 
save(wave_name, 'wave');


%% FUNCTION


function [varargout] = block_nondead_cells(wave_bound, vector_per_side, IMSIZE, varargin)
    varargout = cell((nargin)-3, 1);
    for i = 1:(nargin-3)
        lg = imresize(varargin{i}, [IMSIZE IMSIZE], 'bilinear');
        lg(~wave_bound) = NaN;
        varargout{i} = imresize(lg, [vector_per_side vector_per_side], "bilinear");
    end
end

function [xlim1, ylim1, xlim2, ylim2] = find_plot_xylim(x, y, v)
    offset = 100;

    a = [x(:) y(:) v(:)];

    xmin = min(a(~isnan(a(:,3)),1));
    ymin = min(a(~isnan(a(:,3)),2));
    xmax = max(a(~isnan(a(:,3)),1));
    ymax = max(a(~isnan(a(:,3)),2));

    lim_sz = max([xmax-xmin ymax-ymin]);

    xmid = mean([xmin xmax]);
    ymid = mean([ymin ymax]);
    

    xlim1 = xmid - lim_sz/2 - offset;
    ylim1 = ymid - lim_sz/2 - offset;
    xlim2 = xmid + lim_sz/2 + offset;
    ylim2 = ymid + lim_sz/2 + offset;
end

function [UG, VG] = obtain_direction_component(mag, dir, sigma)
    % obtain component
    U = mag.*cosd(dir);
    V = -1.*mag.*sind(dir);
    
    U(isnan(U)) = 0;
    V(isnan(V)) = 0;
        
    % smoothen component
    UG = imgaussfilt(U, sigma);
    VG = imgaussfilt(V, sigma);
end

function myparula = generate_black_parula(steps_black)
    if ~exist('steps_black')
        steps_black = 25;
    end

    % Number of steps for each segment
    n_black_to_blue = steps_black;
    n_parula = 412;
    
    % Get the parula colormap
    parula_map = parula(n_parula);
    
    % Identify the dark blue of parula (first color)
    parula_start = parula_map(1, :);  % usually a dark blue
    
    % Interpolate from black to parula's blue
    black = [50/256 50/256 50/256];
    ramp = [linspace(black(1), parula_start(1), n_black_to_blue)', ...
            linspace(black(2), parula_start(2), n_black_to_blue)', ...
            linspace(black(3), parula_start(3), n_black_to_blue)'];
    
    % Concatenate the colormaps
    myparula = [ramp; parula_map];
end
