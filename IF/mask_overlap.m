%
%   Calculate boundaries and overlaps
%

clearvars
clc
close all

ROOTDIR = "E:/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATIONS

% all_prefix = ["YAP-1018-day1-tile" "YAP-1018-day2-tile" "YAP-1018-day3-tile"];
% all_prefix = ["NRF2-1125-100k-day1-tile", "NRF2-1125-100k-day2-tile", "NRF2-1125-100k-day3-tile"];

all_prefix = ["YAP-0613-TRULI", "NRF2-0613-KI696-200k"];

all_end_frame = [15, 16];

has_well = 0;
imsize = [5000 5000];

% for p = 1:length(all_prefix)
for p = [2]
    prefix = all_prefix(p);
    disp(prefix);

    [w, w_name, well] = get_well_info(has_well);

    outpath = ROOTDIR + "IF/" + prefix + "/" + well + "/c1_sig_adj/";
    mkdir(outpath);
    
    outfile_prefix = outpath + prefix + "_ratio";

    % get frame/z information
    z_start = 1;
    z_end = all_end_frame(p);

    fnames = strings(z_end);
    for z = z_start:z_end
        if z_end > 9
            fnames(z) = prefix + "_z" + sprintf("%02d", z);
        else
            fnames(z) = prefix + "_z" + sprintf("%01d", z);
        end
    end

    % load files
    seg_file = ROOTDIR + "/IF/" + prefix + "/c3_dapi_density/segment/" + ...
        prefix + "_zMaxc3_ORG_stardist.tif";
    seg = imread(seg_file);
    
    focus_dapi_file = ROOTDIR + "/IF/" + prefix + "/c3_dapi_density/"  + ...
        prefix + "_zMaxc3_ORG_stardist_table_focus_dapi.txt";
    focus_dapi = readtable(focus_dapi_file, 'Delimiter', '\t');
    
    focus_sig_file = ROOTDIR + "/IF/" + prefix + "/c3_dapi_density/"  + ...
        prefix + "_zMaxc3_ORG_stardist_table_focus_sig.txt";
    focus_sig = readtable(focus_sig_file, 'Delimiter', '\t');

    % initialize variables
    nuc_bo = cell(height(focus_sig), 1);
    outer_bo = cell(height(focus_sig), 1);

    all_bo = cell(z_end, height(focus_sig));
    all_bo0 = cell(z_end, height(focus_sig));
    all_bo1 = cell(z_end, height(focus_sig));
    all_bo2 = cell(z_end, height(focus_sig));

    overlap_mask = cell(length(z_end), 1);

    for z = z_start:z_end
        overlap_mask{z} = zeros(imsize);
    end

    for i = 1:height(focus_sig)
        mod(i, 1000) == 0 && fprintf('%d\n', i); % prints i every 1000 iterations

        % calculate the centroid
        centroidx = focus_sig.X(i)+focus_sig.Width(i)/2;
        centroidy = focus_sig.Y(i)+focus_sig.Height(i)/2;
    
        % skip nucelus with position or centroids that are outside of the image
        if skip_outofbounds(focus_sig.X(i), focus_sig.Y(i), 0, imsize(1)), continue; end
        if skip_outofbounds(centroidx, centroidy, 0, imsize(1)), continue; end

        % load the most focused z from IF and dapi images
        focus_sig_z  = focus_sig.Focus(i);
        focus_dapi_z = focus_dapi.Focus(i);
    
        % skip nucleus that do not have a proper z
        if skip_outofbounds(focus_sig_z, focus_dapi_z, z_start-1, z_end), continue; end

        % skip nucleus that are not on the same level or having 1 higher z of nucleus over cytoplasm
        if focus_dapi_z - focus_sig_z > 1 || focus_dapi_z - focus_sig_z < 0, continue; end
    
        % calculate nucleus signal
        seg_value = seg(fix(centroidy), fix(centroidx));
        J = seg==seg_value;
        
        se = strel('disk', 2);
        J0 = imerode(J, se);

        se = strel('disk',6);
        J2 = imdilate(J, se);
        se = strel('disk',2);
        J1 = imdilate(J,se);

        B = bwboundaries(J);
        B = B{1};
        B0 = bwboundaries(J0);
        B0 = B0{1};
        B1 = bwboundaries(J1);
        B1 = B1{1};
        B2 = bwboundaries(J2);
        B2 = B2{1};
        
        overlap_mask{focus_sig_z} = overlap_mask{focus_sig_z} + poly2mask(B2(:, 2), B2(:, 1), imsize(1), imsize(2));
        
        all_bo{focus_sig_z, i} = B;
        all_bo0{focus_sig_z, i} = B0;
        all_bo1{focus_sig_z, i} = B1;
        all_bo2{focus_sig_z, i} = B2;

        % figure;
        % imshow(overlap_mask{focus_sig_z})
    end

    % determine overlapped regions and store binary mask
    for z = z_start:z_end
        overlap_mask{z} = imbinarize(overlap_mask{z}, 1.5);
                
        figure;
        
        imshow(overlap_mask{z})
        
        hold on
        for k = 1:height(focus_sig)
            B = all_bo2{z, k};
            if ~isempty(B)
                plot(B(:, 2), B(:, 1), 'w');
            end
        end
        title(gca, prefix + " z" + z)
        drawnow;
        pause(0.1);
        hold off
    end

    mask_outfile = ROOTDIR + "/IF/" + prefix + "/c3_dapi_density/"  + ...
        prefix + "_zMaxc3_ORG_overlap_mask.mat";
    save(mask_outfile, 'overlap_mask');
end

%% FUNCTIONS
function result = skip_outofbounds(x, y, low, high)
    result = x <= low || y <= low || x > high || y > high || isnan(x) || isnan(y);
end


function [w, w_name, well] = get_well_info(has_well, w)
    if has_well
        w_name = "s" + w;
        well = "Well" + sprintf("%01d", w) + "/";
    else
	    w = "";
	    w_name = "";
	    well = "";
    end
end