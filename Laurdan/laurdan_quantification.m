%
% crop laurdan image
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
addpath(ROOTDIR + "/code");
addpath(ROOTDIR + "/code/util");

warning('off', 'MATLAB:MKDIR:DirectoryExists');


%% load laurdan image
showfigure = 0;
use_hsb = 1;

prefix = "laurdan-0925-control-02"; %"laurdan-0706-AA"; % "laurdan-0707-control"; %"laurdan-0803";

w_start = 1;
w_end = 1;
welllist = w_start:w_end;

z_start = 1;
z_end = 36;

imsize = [1347 1347]; %[1536 1536]; %[1693 1693]; % [1536 1536];

hsb = cell(w_end-w_start+1, z_end);
I = cell(w_end-w_start+11, z_end);

yellow_H_low = 0.105;
yellow_H_high = 0.2;

for w = welllist
    well = sprintf("%02d", w);
    for i = z_start:z_end
        try
            file = ROOTDIR + "/Laurdan/" + prefix + "/Well" + well + "/HSB_images/" + ...
                prefix + "_s" + well + "z" + sprintf("%02d", i) + "c1_ORG_HSB.tif";
            hsb{w, i} = imread(file);
        catch
            file = ROOTDIR + "/Laurdan/" + prefix + "/Well" + well + "/HSB_images/" + ...
                prefix + "_s" + well + "t" + sprintf("%02d", i) + "c1_ORG_HSB.tif";
            hsb{w, i} = imread(file);
        end
    
   
        if use_hsb
            im = rgb2hsv(hsb{w, i});
            H = im(:, :, 1);
            S = im(:, :, 2);
            V = im(:, :, 3);
        
            yellow_mask = (H >= yellow_H_low) & (H <= yellow_H_high); 
        
            yellow_intensity = V .* yellow_mask;
        
            I{w, i} = yellow_intensity;
        else
            % Uses GP image for selecting ordered-ness
            file = ROOTDIR + "/Laurdan/" + prefix + "/Well" + well + "/GP_images/" + ...
                prefix + "_s" + well + "t" + sprintf("%02d", i) + "c1_ORG_GP.tif";

            gp = imread(file);
            gp_mask = imbinarize(gp_mask, 0);
        
            im = rgb2hsv(hsb{w, i});
            V = im(:, :, 3);
            order_im = V .* gp_mask;
        
            I{w, i} = order_im;
        end
    
        if showfigure
            figure
            tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact')
            nexttile
            imshow(hsv2rgb(im))
    
            if use_hsb
                nexttile
                imshow(yellow_mask)
        
                nexttile
                imshow(yellow_intensity, [0 0.4], colormap=copper)
            else
                nexttile
                imshow(gp_mask, [0 1])
        
                nexttile
                imshow(order_im, [0 0.2], colormap=copper)
            end
        end   
    end
end


%% load position
file = ROOTDIR + "/Laurdan/" + prefix + "_nuc_positions.txt";
data = readtable(file, 'Delimiter', '\t');
data.exp = string(data.exp);

%% 
% close all

outpath = ROOTDIR + "/Laurdan/" + prefix + "_crop_cell_smooth/";
mkdir(outpath);

summary = [];

overlay = zeros([300 300]);

periph_ratio = [];
nuc_ratio = [];

curr_well = 1;

for i = 1:height(data) %168:188 %81:86 %:95 %96:125 %3 %1:6 %117  %1:22 %135:137 %96:125
    x = data.x(i);
    y = data.y(i);
    w = data.well(i);
    s = data.slice(i);

    well = sprintf("%02d", w);
    idx = sprintf("%02d", data.index(i));

    mkdir(outpath + "Well" + well);

    if w ~= curr_well
        continue;
    end

    if isnan(x)
        disp("skips " + i + ", x is nan");
        continue;
    end

    im = I{w, s}; % load the masked image
    im_hsb = hsb{w, s}; % load the hsb image

    % get the xy coordinates of the corner of semi-major/minor axes
    len_x = data.length_x(i);
    len_y = data.length_y(i);

    wid_x = data.width_x(i);
    wid_y = data.width_y(i);

    side_x = max( [abs(x-len_x) abs(x-wid_x)] )*2; % determine the long side
    side_y = max( [abs(y-len_y) abs(y-wid_y)] )*2;
    % side = max([side_x side_y]) + max([side_x side_y]); % add half of its nucleus length in the cropped im
    
    % old code
    side = max([side_x side_y]) + 125; % add half of its nucleus length in the cropped im
    side_vertical = side_x+100;
    side_horizontal = side_y+10;
    % end of old code
    
    hyp = ceil(hypot(side, side));

    corner_x = x - hyp/2;
    corner_y = y - hyp/2;

    if ~any([corner_x corner_y]) || any([corner_x+hyp corner_y+hyp] > imsize(1))
        disp("skips " + i + ", corner larger than image size");
        continue
    end

    angle = atan2(len_y-y, len_x-x) - pi/2;

    figure
    % tiledlayout(1, 4, "TileSpacing", "compact", 'Padding', 'compact')
    tiledlayout(1, 4, "TileSpacing", "compact", 'Padding', 'compact')
    sgtitle( data.exp(i)+" well"+data.well(i)+" cell"+data.index(i))

    % original image
    nexttile
    hold on
    imshow(hsb{w, s});
    scatter(x, y, 'yellow')
    scatter(len_x, len_y, 'r')
    scatter(wid_x, wid_y, 'b')
    drawrectangle('Position', [corner_x corner_y hyp hyp]);
    hold off

    % focus in on the selected nucleus
    large_im = imcrop(im, [corner_x corner_y hyp hyp]);
    large_hsb = imcrop(hsb{w, s}, [corner_x corner_y hyp hyp]);
    
    nexttile
    imshow(large_hsb)

    % nexttile
    % imshow(crop_R_hsb)

    % rotate the image to be vertical
    R_im = imrotate(large_im, rad2deg(angle), 'nearest');
    R_hsb = imrotate(large_hsb, rad2deg(angle), 'nearest');

    nexttile
    imshow(R_hsb)
    center_x = size(R_hsb, 1)/2;
    center_y = size(R_hsb, 2)/2;
    semiaxes = [pdist([wid_x wid_y; x y]) pdist([len_x len_y; x y])];

    % inner_ellipse = drawellipse('Center', [center_x, center_y], 'SemiAxes', semiaxes);

    %%% new code
    nuc_length = semiaxes(2)*2;
    nuc_width = semiaxes(1)*2;

    crop_length = nuc_length*2;
    crop_width = nuc_width*2.5;

    crop_corner_x = center_x - nuc_width*1.25;
    crop_corner_y = center_y - nuc_length;

    drawrectangle('Position', [crop_corner_x crop_corner_y crop_width crop_length]);

    % crop the rotated image
    crop_R_im = imcrop(R_im, [crop_corner_x crop_corner_y crop_width crop_length]);
    crop_R_hsb = imcrop(R_hsb, [crop_corner_x crop_corner_y crop_width crop_length]);

    nexttile
    imshow(crop_R_hsb)
    %%% end of new code


    outfile = outpath + "Well" + well + "/" + data.exp(i)+"_well"+well+"_cell"+idx+"_crop.tif";
    imwrite(crop_R_im, outfile);

    close all

end


%% Create a custom colormap bwr
yellow = [1, 1, 0];
black = [0, 0, 0];
blue = [0, 0, 1];

nSteps = 25;

mycol = [linspace(yellow(1), black(1), nSteps)' linspace(yellow(2), black(2), nSteps)' linspace(yellow(3), black(3), nSteps)'; ...
        linspace(black(1), blue(1), nSteps)' linspace(black(2), blue(2), nSteps)' linspace(black(3), blue(3), nSteps)'];
mycol = flipud(mycol);