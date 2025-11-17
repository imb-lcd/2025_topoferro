%
% Calculate YAP nucleus to cytosol ratio
% Select the most focused Z
%
clearvars
clc
close all

ROOTDIR = "/lab_home/N417/Jen-Hao/Spatiotemporal/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% run all prefixes
% all_prefix = ["YAP-1018-150k-day1", "YAP-1018-150k-day2", "YAP-1018-150k-day3", "YAP-1208-100k"];
% all_prefix = ["YAP-0625-Control", "YAP-0625-TRULI", "NRF2-0625-Control", "NRF2-0625-KI696"];
% 
% all_end_frame = [19, 16, 18, 14];
% all_end_frame = [13 19 16 7]; %, 10 10 10];
% all_end_frame = [8, 12, 14];
all_prefix = ["YAP-0613-TRULI", "NRF2-0613-KI696-200k"];

all_end_frame = [15, 16];

hasWell = 0;
imsize = [5000 5000];
cmap = [0 0 0; parula(256)];

has_overlap = 1;
overlap_threshold = 0.1;

showfigure = 0;

%% load images (cy5, dapi, and segmented nucleus-labelled)
% for p = 1:length(all_prefix)
for p = 1:length(all_prefix)
    prefix = all_prefix(p); %%"YAP-1018-day1-tile";
    disp(prefix);

    if hasWell
	    w = 1;
	    w_name = "s" + w;
	    well = "Well" + sprintf("%01d", w) + "/";
    else
	    w = "";
	    w_name = "";
	    well = "";
    end
    
    % load IF and segementation files
    z_start = 1;
    z_end = all_end_frame(p);
    
    fnames = strings(z_end);
    for z = z_start:z_end
        if z_end > 0
            fnames(z) = prefix + "_z" + sprintf("%02d", z);
        else
            fnames(z) = prefix + "_z" + sprintf("%01d", z);
        end
    end
    
    if has_overlap
        mask_path = ROOTDIR + "IF/" + prefix + "/" + well + "/c3_dapi_density/" + ...
            prefix + "_zMaxc3_ORG_overlap_mask.mat";
        load(mask_path);
    else
        overlap_mask = cell(z_end, 1); % to remove
    end

    IF = cell(z_end, 1);
    for z = z_start:z_end
        path = ROOTDIR + "IF/" + prefix + "/" + well + "/c1_sig_mod/";
        IF{z} = imread(path+fnames(z)+"c1_ORG.tif");
        IF{z} = medfilt2(IF{z});
        
        if has_overlap
            overlap_mask{z} = double(overlap_mask{z});
            overlap_mask{z}(overlap_mask{z} == 1) = -1;
            overlap_mask{z}(overlap_mask{z} == 0) = 1;
    
            IF{z} = single(IF{z}) .* overlap_mask{z};
            IF{z}(IF{z} < 0) = NaN;
        end
    end
    
    seg_file = ROOTDIR + "/IF/" + prefix + "/c3_dapi_density/segment/" + ...
        prefix + "_zMaxc3_ORG_stardist.tif";
    seg = imread(seg_file);
    
    focus_dapi_file = ROOTDIR + "/IF/" + prefix + "/c3_dapi_density/"  + ...
        prefix + "_zMaxc3_ORG_stardist_table_focus_dapi.txt";
    focus_dapi = readtable(focus_dapi_file, 'Delimiter', '\t');
    
    focus_sig_file = ROOTDIR + "/IF/" + prefix + "/c3_dapi_density/"  + ...
        prefix + "_zMaxc3_ORG_stardist_table_focus_sig.txt";
    focus_sig = readtable(focus_sig_file, 'Delimiter', '\t');
    
    %
    outpath = ROOTDIR + "IF/" + prefix + "/" + well + "/c1_sig_adj/";
    mkdir(outpath);
    
    outfile_prefix = outpath + prefix + "_ratio";
    
    % initialize variables
    seg_value = nan(height(focus_sig), 1);

    nuc_bo = cell(height(focus_sig), 1);
    outer_bo = cell(height(focus_sig), 1);

    nuc_area = nan(height(focus_sig), 1);
    nuc_circ = nan(height(focus_sig), 1);
    nuc_ecc = nan(height(focus_sig), 1);
    nuc_solid = nan(height(focus_sig), 1);
    nuc_aspect = nan(height(focus_sig), 1);

    px_sig_n = cell(height(focus_sig), 1);
    px_sig_c = cell(height(focus_sig), 1);

    overlap_n = nan(height(focus_sig), 1);
    overlap_c = nan(height(focus_sig), 1);
    
    [sig_n_sd, sig_n_mean, sig_n_top25, sig_n_top50, sig_n_topBttm] = deal(nan(height(focus_sig), 1), nan(height(focus_sig), 1), nan(height(focus_sig), 1), nan(height(focus_sig), 1), nan(height(focus_sig), 1));
    [sig_c_sd, sig_c_mean, sig_c_top25, sig_c_top50, sig_c_topBttm] = deal(nan(height(focus_sig), 1), nan(height(focus_sig), 1), nan(height(focus_sig), 1), nan(height(focus_sig), 1), nan(height(focus_sig), 1));

    [ratio, ratio_top25, ratio_top50, ratio_topBttm] = deal(nan(height(focus_sig), 1), nan(height(focus_sig), 1), nan(height(focus_sig), 1), nan(height(focus_sig), 1));
    
    if showfigure
        [plt_ratio, plt_r_top25, plt_r_top50, plt_r_topBttm] = deal(cell(z_end, 1), cell(z_end, 1), cell(z_end, 1), cell(z_end, 1));
        [plt_n, plt_n_top25, plt_n_top50, plt_n_topBttm] = deal(cell(z_end, 1), cell(z_end, 1), cell(z_end, 1), cell(z_end, 1));
        [plt_c, plt_c_top25, plt_c_top50, plt_c_topBttm] = deal(cell(z_end, 1), cell(z_end, 1), cell(z_end, 1), cell(z_end, 1));
    
        % initialize plt images for each z
        for z = z_start:z_end
            [plt_ratio{z}, plt_r_top25{z}, plt_r_top50{z}, plt_r_topBttm{z}] = deal(zeros(imsize), zeros(imsize), zeros(imsize), zeros(imsize));
    
            [plt_n{z}, plt_n_top25{z}, plt_n_top50{z}, plt_n_topBttm{z}] = deal(zeros(imsize), zeros(imsize), zeros(imsize), zeros(imsize));
            [plt_c{z}, plt_c_top25{z}, plt_c_top50{z}, plt_c_topBttm{z}] = deal(zeros(imsize), zeros(imsize), zeros(imsize), zeros(imsize));
        end
    end

    %
    for i = 1:height(focus_sig)
        mod(i, 500) == 0 && fprintf('%d\n', i); % prints i every 1000 iterations
   
        % calculate the centroid
        centroidx = focus_sig.X(i)+focus_sig.Width(i)/2;
        centroidy = focus_sig.Y(i)+focus_sig.Height(i)/2;
    
        % skip nucelus with position or centroids that are outside of the image
        if skip_outofbounds(focus_sig.X(i), focus_sig.Y(i), 0, imsize(1)), continue; end
        if skip_outofbounds(centroidx, centroidy, 0, imsize(1)), continue; end
    
        % load the most focused z from IF and dapi images
        focus_sig_z  = focus_sig.Focus(i);
        focus_dapi_z = focus_dapi.Focus(i);
    
        if isnan(focus_sig_z)
            continue;
        end

        % skip nucleus that do not have a proper z
        if skip_outofbounds(focus_sig_z, focus_dapi_z, z_start-1, z_end), continue; end
 
        % calculate nucleus signal
        seg_value(i) = seg(fix(centroidy), fix(centroidx));
        J = seg==seg_value(i);

        stat = regionprops(J, IF{focus_sig_z}, {'Area', 'Circularity', 'Eccentricity', 'Solidity', 'MajorAxisLength', 'MinorAxisLength'});
        nuc_area(i) = double(stat(1).Area);
        nuc_circ(i) = double(stat(1).Circularity);
        nuc_ecc(i) = double(stat(1).Eccentricity);
        nuc_solid(i) = double(stat(1).Solidity);
        nuc_aspect(i) = double(stat(1).MinorAxisLength) / double(stat(1).MajorAxisLength);

        % erode to get the inner nucleus
        se = strel('disk', 2);
        J0 = imerode(J, se);

        stat = regionprops(J0, IF{focus_sig_z}, {'PixelValues'});

        [px_sig_n{i}, sig_n_sd(i), sig_n_mean(i), sig_n_top25(i), sig_n_top50(i), sig_n_topBttm(i)] = get_sig_mean_prctile(stat);

        overlap_n(i) = sum(isnan(px_sig_n{i})) / numel(px_sig_n{i});

        % calculate outer (cytoplasm) signal
        se = strel('disk',6);
        J2 = imdilate(J,se);
        se = strel('disk',2);
        J1 = imdilate(J,se);

        J12 = J1==0 & J2==1;

        stat = regionprops(J12, IF{focus_sig_z},{'PixelValues'});

        [px_sig_c{i}, sig_c_sd(i), sig_c_mean(i), sig_c_top25(i), sig_c_top50(i), sig_c_topBttm(i)] = get_sig_mean_prctile(stat);

        overlap_c(i) = sum(isnan(px_sig_c{i})) / numel(px_sig_c{i});

        % calculate ratio
        ratio(i) = sig_n_mean(i) / sig_c_mean(i);
        ratio_top25(i) = sig_n_top25(i) / sig_c_top25(i);
        ratio_top50(i) = sig_n_top50(i) / sig_c_top50(i);
        ratio_topBttm(i) = sig_n_topBttm(i) / sig_c_topBttm(i);

        % store boundaries
        bo = bwboundaries(J);
        nuc_bo{i} = bo(1);

        bo = bwboundaries(J2);
        outer_bo{i} = bo(1);

        % calculate figure
        if showfigure
            % figure % figure(focus_sig_z);
            % imshow(IF{focus_sig_z}, [1 30])
            % imshow(seg, [0 1])
            % % plot boundary
            % hold on
            % plot_boundary(J0, 'c') % inner
            % plot_boundary(J1, 'm') % outer 1
            % plot_boundary(J2, 'm') % outer 2
            % hold off

            % store nuc values
            plt_ratio{focus_sig_z}(J0 == 1) = ratio(i);
            plt_r_top25{focus_sig_z}(J0 == 1) = ratio_top25(i);
            plt_r_top50{focus_sig_z}(J0 == 1) = ratio_top50(i);
            plt_r_topBttm{focus_sig_z}(J0 == 1) = ratio_topBttm(i);
    
            plt_n{focus_sig_z}(J0 == 1) = sig_n_mean(i);
            plt_n_top25{focus_sig_z}(J0 == 1) = sig_n_top25(i);
            plt_n_top50{focus_sig_z}(J0 == 1) = sig_n_top50(i);
            plt_n_topBttm{focus_sig_z}(J0 == 1) = sig_n_topBttm(i);
            
            plt_c{focus_sig_z}(J0 == 1) = sig_c_mean(i);
            plt_c_top25{focus_sig_z}(J0 == 1) = sig_c_top25(i);
            plt_c_top50{focus_sig_z}(J0 == 1) = sig_c_top50(i);
            plt_c_topBttm{focus_sig_z}(J0 == 1) = sig_c_topBttm(i);
        end
    end
    
    if showfigure
        % sum the plt variables
        all_ratio = sum(cat(3, plt_ratio{:}), 3);
        all_r_top25 = sum(cat(3, plt_r_top25{:}), 3);
        all_r_top50 = sum(cat(3, plt_r_top50{:}), 3);
        all_r_topBttm = sum(cat(3, plt_r_topBttm{:}), 3);

        all_n = sum(cat(3, plt_n{:}), 3);
        all_n_top25 = sum(cat(3, plt_n_top25{:}), 3);
        all_n_top50 = sum(cat(3, plt_n_top50{:}), 3);
        all_n_topBttm = sum(cat(3, plt_n_topBttm{:}), 3);

        all_c = sum(cat(3, plt_c{:}), 3);
        all_c_top25 = sum(cat(3, plt_c_top25{:}), 3);
        all_c_top50 = sum(cat(3, plt_c_top50{:}), 3);
        all_c_topBttm = sum(cat(3, plt_c_topBttm{:}), 3);

        all_ratio(all_ratio == 0) = nan;
        all_r_top25(all_r_top25 == 0) = nan;
        all_r_top50(all_r_top50 == 0) = nan;
        all_r_topBttm(all_r_topBttm == 0) = nan;

        all_n(all_n == 0) = nan;
        all_n_top25(all_n_top25 == 0) = nan;
        all_n_top50(all_n_top50 == 0) = nan;
        all_n_topBttm(all_n_topBttm == 0) = nan;

        all_c(all_c == 0) = nan;
        all_c_top25(all_c_top25 == 0) = nan;
        all_c_top50(all_c_top50 == 0) = nan;
        all_c_topBttm(all_c_topBttm == 0) = nan;

        %
        data = {all_ratio, all_r_top25, all_r_top50, all_r_topBttm, all_n, all_n_top25, all_n_top50, all_n_topBttm, all_c, all_c_top25, all_c_top50, all_c_topBttm};
        all_title = ["ratio" "top25" "top50" "top-bttom", "nuc-int" "nuc-top25" "nuc-top50" "nuc-top-bttom", "cyto-int" "cyto-top25" "cyto-top50" "cyto-top-bttom"];

        % for i = 1:length(data)
        for i = 1
            figure
            imagesc(data{i}); %, 'AlphaData', ~isnan(nuc_int{z})); %, 'AlphaDataMapping', 'none')
            axis square ij off tight
            colormap(cmap)
            if i <= 4
                clim([0.5 1.5])
            else
                clim([5 15])
            end
            title(all_title(i));

            fig_file = outfile_prefix + "_" + all_title(i) + ".jpg";
            exportgraphics(gcf, fig_file);
        end
        % close all
    end
    
    diff_focus = focus_dapi.Focus - focus_sig.Focus;

    if_ratio = [focus_sig.index focus_sig.X focus_sig.Y focus_sig.Width focus_sig.Height nuc_area nuc_circ nuc_ecc nuc_solid nuc_aspect ...
        overlap_n, sig_n_mean sig_n_top25 sig_n_top50 sig_n_topBttm sig_n_sd ...
        overlap_c, sig_c_mean sig_c_top25 sig_c_top50 sig_c_topBttm sig_c_sd ...
        ratio ratio_top25 ratio_top50 ratio_topBttm ...
        focus_dapi.Focus focus_sig.Focus diff_focus];

    if_ratio = array2table(if_ratio, 'VariableNames',{'index', 'x','y', 'Width', 'Height', 'Area', 'Circularity', 'Eccentricity', 'Solidity', 'Aspect_ratio' ...
        'nuc_overlap_percent', 'nuc_intensity', 'nuc_top25', 'nuc_top50', 'nuc_topBttm', 'nuc_std', ...
        'cyto_overlap_percent', 'cyto_intensity', 'cyto_top50', 'cyto_top25', 'cyto_topBttm', 'cyto_std', ...
        'ratio', 'ratio_top25', 'ratio_top50', 'ratio_topBttm', ...
        'Dapi_Focus' 'IF_Focus' 'Focus_diff'});

    writetable(if_ratio, outfile_prefix+"_summary.txt", 'Delimiter', '\t');

    % store boundaries and pixel_intensties
    px_n_outfile = outfile_prefix + "_px_nuc.mat";
    save(px_n_outfile, 'px_sig_n');

    px_c_outfile = outfile_prefix + "_px_cyto.mat";
    save(px_c_outfile, 'px_sig_c');
end

%% FUNCTION
function result = skip_outofbounds(x, y, low, high)
    result = x <= low || y <= low || x > high || y > high || isnan(x) || isnan(y);
end


function [px_sig, sig_sd, sig_mean, sig_top25, sig_top50, sig_topBttm] = get_sig_mean_prctile(stat)

    [~, idx] = max(cellfun(@numel, {stat.PixelValues})); % find the largest area in the nucleus

    px_sig = double(stat(idx).PixelValues);
    % px_sig = px_sig(~isnan(px_sig));

    sig_mean = mean(px_sig, 'omitnan');
    sig_sd = std(px_sig, 'omitnan');

    % sig_mean = 

    % discard bottom 25%
    tmp = px_sig(px_sig >= prctile(px_sig, 25));
    sig_top25 = mean(tmp);

    % discard bottm 50%
    tmp = px_sig(px_sig > median(px_sig, 'omitnan'));
    sig_top50 = mean(tmp);

    % discard top and bottom 25%
    tmp = px_sig(px_sig >= prctile(px_sig, 25) & px_sig <= prctile(px_sig, 75));
    sig_topBttm = mean(tmp);
end

function plot_boundary(J, mycol, linewidth)
    if ~exist('linewidth', 'var')
        linewidth = 1;
    end

    B = bwboundaries(J);
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), mycol, 'LineWidth', linewidth)
    end
end
