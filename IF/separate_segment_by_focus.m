%
%   Create segmented nuclues files for their most focused nucleus
%

clearvars
clc
close all

addpath D:/Spatiotemporal_analysis/code/
addpath D:/Spatiotemporal_analysis/code/util/


ROOTDIR = "E:/IF/";

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% load focused Z file

all_prefix = ["YAP-0613-TRULI", "NRF2-0613-KI696-200k"];
all_end_frame = [15, 16];

imsize = [5000 5000];

% for p = 1:length(all_prefix)
for p = 2
    prefix = all_prefix(p);
    disp(prefix);

    zseg = repmat({zeros(imsize)}, all_end_frame(p), 1);

    ratio_file = ROOTDIR + "/" + prefix + "/c1_sig_adj/" + prefix + "_ratio_summary.txt";
    ratio = readtable(ratio_file, "Delimiter", "\t");

    seg_file = ROOTDIR + "/" + prefix + "/c3_dapi_density/segment/" + prefix + "_zMaxc3_ORG_stardist.tif";
    seg = imread(seg_file);

    outpath = ROOTDIR + "/" + prefix + "/c3_dapi_density/segment_focus/";
    mkdir(outpath);

    for i = 1:height(ratio)
        focus = ratio.IF_Focus(i);
        
        % if ratio.Focus_diff(i) < 0 || ratio.Focus_diff(i) > 1 || isnan(focus)
        %     continue;
        % end
        if isnan(focus)
            continue
        end

        J = seg==(ratio.index(i)+1);
        J = J * (ratio.index(i)+1);

        % figure
        % imshow(J)
        
        zseg{focus} = zseg{focus} + J;
    end

    for z = 1:all_end_frame(p)
        if max(zseg{z}) > max(seg(:))
            disp("error " + z)
        end
        
        zseg{z} = uint16(zseg{z});

        outfile = prefix + "_z" + sprintf("%02d", z) + "c3_ORG_stardist_focus.tif";
        imwrite(zseg{z}, outpath+outfile)
    end
end