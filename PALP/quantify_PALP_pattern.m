%
%   quantify_PALP_pattern.m
%
clearvars
clc
close all

ROOTDIR = "E:/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath D:\Spatiotemporal_analysis\code
addpath D:\Spatiotemporal_analysis\code\util

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATIONS

prefix = "PALP-0224-";
path = ROOTDIR + "PALP/PALP-pattern-0224-2/";
type_list = ["HDaln-2-", "HDmaln-2-", "LDaln-2-", "LDmaln-2-", "negdef-2-", "posdef-2-"];
fend_list = [5, 7, 5, 2, 4, 1];
 

fstart = 1;

%%
summary = table('Size', [0 12], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double','double','double','double', 'double'}, ...
    'VariableNames', {'pattern', 'intensity%', 'o', 'r', 'or', 'o_mask1', 'r_mask1', 'or_mask1', 'o_mask2', 'r_mask2', 'or_mask2', 'median_o_mask1'});

all_data = cell(sum(fend_list), 1);

cnt = 1;

% quantify PALP signals
for k = 1:length(type_list)
    ftype = type_list(k);
    fend = fend_list(k);


    for ii = fstart:fend
        disp([ftype ii])
        rep = sprintf('%02d', ii);
    
        % obtain the coordinate of induction sites and create mask on the sites
        t = 1;
        try
            fname = path + "/" + prefix + ftype + ii + "/" + prefix + ftype + ii + "_t"+t+"c2_labeled.tif";
            I = load_image_file(t, fname);
        catch
            fname = path + "/" + prefix + ftype + ii + "/" + prefix + ftype + ii + "_t"+t+"c2.tif";
            I = load_image_file(t, fname);
        end
        
    
        % find a cirlce 
        [centers, radii] = imfindcircles(I(:, :, 1), [180 210], 'Sensitivity', 0.95);
        radii = radii - 5; % choose a slightly smaller circle to be conservative
        
        imsize = size(I(:, :, 1));
        center_x = centers(:, 1);
        center_y = centers(:, 2);
        
        [X, Y] = meshgrid(1:imsize(2), 1:imsize(1));
        mask = false(imsize);
        mask = mask | hypot(X - center_x(1), Y - center_y(1)) <= radii(1);
    
        % load time 2 after induction
        t = 2;
        oxi_im_file = path + "/" + prefix + ftype + ii + "/" + prefix + ftype + ii + "_t"+t+"c1_ORG.tif";
        oxi_im = load_image_file(t, oxi_im_file);
    
        red_im_file = path + "/" + prefix + ftype + ii + "/" + prefix + ftype + ii + "_t"+t+"c2_ORG.tif";
        red_im = load_image_file(t, red_im_file);
    
        norm_oxi_im = oxi_im;
        norm_red_im = red_im;
    
        % get the pixel values
        stats = regionprops(mask, norm_oxi_im, 'PixelValues');
        o = vertcat(stats.PixelValues);
    
        stats = regionprops(mask, norm_red_im, 'PixelValues');
        r = vertcat(stats.PixelValues);;
    

        % calculate ( oxi / (oxi + red) )
        or = o ./ ( o + r );
    
        % load background masks
        % background mask 1
        mask_file1 = path + "/" + prefix + ftype + ii + "/" + prefix + ftype + ii + "_t1c2_ORG_1bitmask.tif";
        m1 = imread(mask_file1);
        mask1 = mask & logical(m1);
   
        stats = regionprops(mask1, norm_oxi_im, 'PixelValues');
        o_m1 = vertcat(stats.PixelValues);
        
        stats = regionprops(mask1, norm_red_im, 'PixelValues');
        r_m1 = vertcat(stats.PixelValues);
        
        
        med_o_m1 = median(double(o_m1), 'omitnan');

        or_m1 = o_m1 ./ ( o_m1 + r_m1 );

        % background mask 1 medfilt
        mask_file2 = path + "/" + prefix + ftype + ii + "/" + prefix + ftype + ii + "_t2c2_ORG_1bitmask_medfilt.tif";
        m2 = imread(mask_file2);
        mask2 = mask & logical(m2);
    
        stats = regionprops(mask2, norm_oxi_im, 'PixelValues');
        o_m2 = vertcat(stats.PixelValues);
        
        stats = regionprops(mask2, norm_red_im, 'PixelValues');
        r_m2 = vertcat(stats.PixelValues);
        
        or_m2 = o_m2 ./ ( o_m2 + r_m2 );
    
        parts = strsplit(ftype, '-');
        newline = {string(parts{1}), string(parts{2}), ...
            mean(o, 'omitnan'), mean(r, 'omitnan'), mean(or, 'omitnan'), ...
            mean(o_m1, 'omitnan'), mean(r_m1, 'omitnan'), mean(or_m1, 'omitnan'), ...
            mean(o_m2, 'omitnan'), mean(r_m2, 'omitnan'), mean(or_m2, 'omitnan'), ...
            med_o_m1};
    
        summary(end+1, :) = newline;

        all_data{cnt} = struct;

        all_data{cnt}.o  = o;
        all_data{cnt}.r  = r;
        all_data{cnt}.or = or;
        all_data{cnt}.o_m1  = o_m1;
        all_data{cnt}.r_m1  = r_m1;
        all_data{cnt}.or_m1 = or_m1;
        all_data{cnt}.med_o_m1 = med_o_m1;

    end
    close all;

end

summary


%% plot the summary

order = ["HDaln", "HDmaln", "LDaln", "LDmaln", "negdef", "posdef", "MDaln"];
ylim_range = [0 1];

summary.pattern = categorical(summary.pattern, order, 'Ordinal', true);

group = summary.pattern;
mask = group == "LDmaln" | group == "negdef" | group == "posdef";
group(mask) = "LDmaln";
group = categorical(group);
group = removecats(group);


figure
tiledlayout(1,2)
nexttile
hold on
swarmchart(group, summary.or, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', 0.5)
boxchart(group, summary.or, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k');
xticklabels(["HDaln", "HDmaln", "LDaln", "LDmaln"])
ylim(ylim_range)
ylabel("o/(o+r)")
title("no masking background")
hold off

nexttile
hold on
swarmchart(group, summary.or_mask1, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', 0.5)
boxchart(group, summary.or_mask1, 'MarkerStyle','none', 'BoxFaceColor', 'none', 'BoxEdgeColor', 'k');

xticklabels(["HDaln", "HDmaln", "LDaln", "LDmaln"])
 ylim(ylim_range)
ylabel("o/(o+r)")
title("mask background")
hold off


%% FUNCTIONS
function I = load_image_file(t, file_path)
    try
        I = imread(file_path);
    catch
        file_path = strrep(file_path, "t"+t, "t0"+t);                                                                                                       
        I = imread(file_path);
    end
end

function I = min_max_normalization(I, Imin, Imax, a, b)
    I = (I - Imin) / (Imax - Imin) * (b-a) + a;
end

function quantile_normalization(images)
    flattened_im = cellfun(@im)
end

function norm_I = percentile_normalization(I, lower_p, upper_p)
    p_low = prctile(I(:), lower_p);
    p_high = prctile(I(:), upper_p);

    norm_I = (I - p_low) ./ (p_high - p_low+1e-20);
    norm_I = reshape(norm_I, size(I));
end