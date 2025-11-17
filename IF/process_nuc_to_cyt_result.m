%
%   Process file of nucleus-cytoplasm intensities from IF results
%

clearvars
clc
close all

ROOTDIR = "E:/IF/";

addpath D:/Spatiotemporal_analysis/code/
addpath D:/Spatiotemporal_analysis/code/util/

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATION

all_prefix = ["YAP-0613-TRULI" "NRF2-0613-KI696-200k"];


z_start = [1 1];  % [1 1 1];    % z_start = [1 1 1, 1 1 1];
z_end   = [15 16]; %[10 9]; % [10 11 11]; % z_end   = [13 10 13, 10 10 11];

imsize = [5000 5000];

%% match ilastik output to ratio file

for p = 1:length(all_prefix)
% for p = 4
    prefix = all_prefix(p);
    
    ilastik_file = ROOTDIR + "/" + prefix + "/ilastik/" + prefix + "_c1_sig_mod_pred_table.csv";
    ilastik = readtable(ilastik_file, "Delimiter", ",");

    seg_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/segment/" + prefix + "_zMaxc3_ORG_stardist.tif";
    seg = imread(seg_file);
    
    ratio_file = ROOTDIR + "/" + prefix + "/c1_sig_adj/" + prefix + "_ratio_summary.txt";
    ratio = readtable(ratio_file, "Delimiter", "\t");
    
    outfile = ROOTDIR + "/" + prefix + "/c1_sig_adj/" + prefix + "_ratio_summary_ilastik.txt";

    % load and match ilastick labels to ratio file
    label = strings(height(ratio), 1);
    
    missing_idx = nan(height(ilastik), 1);
    
    for i = 1:height(ilastik)
        % for i = 1756
        center_x = fix(ilastik.CenterOfTheObject_0(i));
        center_y = fix(ilastik.CenterOfTheObject_1(i));

        center_x = max(1, center_x);
        center_y = max(1, center_y);
    
        index = double(seg(center_y, center_x));
        if index == 0
            missing_idx(i) = i;
        else
            label(index) = ilastik.PredictedClass(i);
        end
    end
    missing_idx(isnan(missing_idx)) = [];
    disp(prefix + " " + length(missing_idx));
    
    % write ilastik label to file
    ratio.ilastik = strings(height(ratio), 1);
    for i = 1:height(ratio) %265 %1768
        if ~isempty(label{i})
            ratio.ilastik(i) = label(i);
        end
    end
    writetable(ratio, outfile, "Delimiter", "\t");
end

%%  label dividing cell, off-focused rows, and ilastik removed rows
disp("run rm dividing cell, off-focused rows, ilastik");

jac = [0.07 0.09 0.07 0.09 0.07 0.06];

for p = 1:length(all_prefix)
    prefix = all_prefix(p);
    disp(prefix)

    seg_table_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/" + prefix + "_zMaxc3_ORG_stardist_table.txt";
    seg = readtable(seg_table_file, "Delimiter", "\t");
    
    div_table_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/" + prefix + "_zMaxc3_ORG_stardist_divcell_table.txt";
    div = readtable(div_table_file, "Delimiter", "\t");
    
    ratio_file = ROOTDIR + "/" + prefix + "/c1_sig_adj/" + prefix + "_ratio_summary_ilastik.txt";
    ratio = readtable(ratio_file, "Delimiter", "\t");

    outfile = ROOTDIR + "/" + prefix + "/c1_sig_adj/" + prefix + "_ratio_summary_ilastik_filter.txt";

    % get the row of dividing cell
    divcell_row = get_row_dividing_cell(div, seg, jac(p));

    off_focus_row = ratio.IF_Focus < z_start(p) | ratio.IF_Focus > z_end(p) | isnan(ratio.IF_Focus) | isnan(ratio.Dapi_Focus);
    off_focus_row = ~(off_focus_row | (ratio.Dapi_Focus-ratio.IF_Focus <= 1 & ratio.Dapi_Focus-ratio.IF_Focus >= 0));

    ratio.divcell = divcell_row;
    ratio.off_focus = off_focus_row;
    ratio.filter = ratio.divcell | ratio.off_focus | double(ismember(ratio.ilastik, ["Label 2", ""]));
    
    writetable(ratio, outfile, "Delimiter", "\t")
end

%% add spatial properties to the IF files 
disp("adding spatial properties");

sp_sz = 220;
den_reg = 220;

for p = 1:length(all_prefix)
% for p = 4
    prefix = all_prefix(p);
    disp(prefix);
    
    ratio_file = ROOTDIR + "/" + prefix + "/c1_sig_adj/" + prefix + "_ratio_summary_ilastik_filter.txt";
    outfile    = ROOTDIR + "/" + prefix + "/c1_sig_adj/" + prefix + "_ratio_summary_ilastik_filter_sp.txt";

    ratio = readtable(ratio_file, "Delimiter", "\t");

    % load density file
    den_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/" + prefix + "_zMaxc3_ORG_den" + den_reg + ".mat";
    den = struct2array(load(den_file));
    den = den'; % necessary step to fit with xy coordinates of the ratio file
    den = imresize(den, imsize);

    % load coherency file
    coh = cell(z_end(p), 1);
    % for i = z_start:10 %z_end
    %     coh_file = ROOTDIR + "/" + prefix + "/c3_dapi_nematics/orientation_files/" + prefix + "_z" + sprintf("%02d", i) + "c3_ORG_coh.mat";
    %     coh{i} = struct2array(load(coh_file));
    %     coh{i} = imresize(coh{i}, imsize);
    % end

    nuc_sp = nan(height(ratio), 2);
    for i = 1:height(ratio)
        focus = ratio.IF_Focus(i);

        centroidx = fix(ratio.x(i) + ratio.Width(i)/2);
        centroidy = fix(ratio.y(i) + ratio.Height(i)/2);
    
        % skip nucelus centroids that are outside of the image
        if centroidx <= 0 || centroidy <= 0 || centroidx > imsize(1) || centroidy > imsize(2)
            continue;
        end
        if focus > 10 || isnan(focus)
            continue;
        end

        start_x = max(centroidx-sp_sz/2, 1);
        end_x   = min(centroidx+sp_sz/2, imsize(1));
    
        start_y = max(centroidy-sp_sz/2, 1);
        end_y   = min(centroidy+sp_sz/2, imsize(2));
    
        den_region = den(start_x:end_x, start_y:end_y);
        % coh_region = coh{focus}(start_x:end_x, start_y:end_y);

        % nuc_sp(i, :) = [mean(den_region(:), 'omitnan') mean(coh_region(:), 'omitnan')];
        nuc_sp(i, :) = [mean(den_region(:), 'omitnan') nan];
    end

    ratio.den = nuc_sp(:, 1);
    ratio.coh = nuc_sp(:, 2);
    
    writetable(ratio, outfile, 'Delimiter', '\t');
end


%%  Determine threshold for dividing cell

% all_prefix = ["YAP-1018-day2-tile", "YAP-1018-day3-tile"];
all_prefix = ["NRF2-1125-100k-day1-tile" "NRF2-1125-100k-day2-tile" "NRF2-1125-100k-day3-tile"];

for p = 1:length(all_prefix)
    prefix = all_prefix(p);
    disp(prefix);
    
    seg_table_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/" + prefix + "_zMaxc3_ORG_stardist_table.txt";
    seg = readtable(seg_table_file, "Delimiter", "\t");

    % div_im_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/segment/" + prefix + "_zMaxc3_ORG_stardist_divcell.tif";
    % div_im = imread(div_im_file);
    % div_im(div_im > 0) = 1;
    
    div_table_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/" + prefix + "_zMaxc3_ORG_stardist_divcell_table.txt";
    div = readtable(div_table_file, "Delimiter", "\t");

    jac_range = [0.01:0.01:0.12];

    headers = {'jac', 'TP', 'FP', 'FN', 'precision', 'recall', 'F1'};
    out = array2table(nan(length(jac_range), numel(headers)), 'VariableNames', headers);

    for i = 1:length(jac_range)
        jac = jac_range(i);
        % figure
        % imshow(div_im, [0 1], 'Border', 'tight')
        % hold on

        jac_cnt = 0;
        div_cnt = 0;

        for dv = 1:height(div)
            x1 = div.X(dv);
            y1 = div.Y(dv);
            w1 = div.Width(dv);
            h1 = div.Height(dv);
            
            counted = 0;

            for s = 1:height(seg)
                x2 = seg.X(s);
                y2 = seg.Y(s);
                w2 = seg.Width(s);
                h2 = seg.Height(s);

                [jac_index, ~, ~] = calculate_jaccard(x1, y1, w1, h1, x2, y2, w2, h2);

                if jac_index > jac
                    jac_cnt = jac_cnt + 1;
                    counted = 1;
                    % drawrectangle("Position", [x2 y2 w2 h2], "MarkerSize", 1, "Color", "y", "LineWidth", 1);
                end
            end
            if counted
                div_cnt = div_cnt + 1;
            end
        end
        % hold off

        % calculate TP, FN, FP. Note this is only a rough estimate
        out.jac(i) = jac;
        out.TP(i) = height(div);
        out.FN(i) = max(out.TP(i)-div_cnt, 0);
        out.FP(i) = max(jac_cnt-out.TP(i), 0);

        out.precision(i) = out.TP(i) / (out.TP(i) + out.FP(i));
        out.recall(i) = out.TP(i) / (out.TP(i) + out.FN(i));
        out.F1(i) = 2 * out.precision(i) * out.recall(i) / (out.precision(i) + out.recall(i));
    end

    outfile = ROOTDIR + "/" + prefix + "/" + prefix + "_jaccard_dividing_cells.txt";
    writetable(out, outfile, 'Delimiter', '\t');
end


%%  FUNCTION
function divcell_row = get_row_dividing_cell(div, seg, jac)
    divcell_row = false(height(seg), 1);
    for dv = 1:height(div)
        x1 = div.X(dv);
        y1 = div.Y(dv);
        w1 = div.Width(dv);
        h1 = div.Height(dv);
    
        for s = 1:height(seg)
            % record row if overlap with dividing cellls
            x2 = seg.X(s);
            y2 = seg.Y(s);
            w2 = seg.Width(s);
            h2 = seg.Height(s);
    
            [jac_index, ~, ~] = calculate_jaccard(x1, y1, w1, h1, x2, y2, w2, h2);
    
            if jac_index > jac
                divcell_row(s) = 1 | divcell_row(s);
            end
        end
    end
end

function [jac_index, union_area, intersect_area] = calculate_jaccard(x1, y1, w1, h1, x2, y2, w2, h2)
    % Calculate the coordinates of the intersection rectangle
    intersect_x = max(x1, x2);
    intersect_y = max(y1, y2);

    intersect_width  = min(x1+w1, x2+w2) - intersect_x;
    intersect_height = min(y1+h1, y2+h2) - intersect_y;

    % Check for overlap
    if intersect_width <= 0 || intersect_height <= 0
        % No overlap, Jaccard Index is 0
        jac_index = 0;
        union_area = 0;
        intersect_area = 0;
    else
        % Calculate intersection and union area
        intersect_area = intersect_width * intersect_height;

        union_area = (w1*h1) + (w2*h2) - intersect_area;

        % return Jaccard Index
        jac_index = intersect_area / union_area;
    end
end