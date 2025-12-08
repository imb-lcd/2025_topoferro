% run_orientation_analysis.m
% pipeline for orientation analysis for calculation of nematics, coherency, and defects
clearvars
clc
close all

ROOTDIR = "D:/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath(ROOTDIR + "/code/util/");
addpath(ROOTDIR + "/code/orientation/");
addpath(MATLAB_FX + "differential_entropy");

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% CONFIGURATIONS
% orientation
ws_px = 2;
sigma = 45; % defects and entropy 45.0; phi is 20

% generating nematics
nem_ws = 50; % super large presentation: 100 or 75; large presenation: 50; for nematics image: 20
nem_color = 'k'; %[0 0 0];
scale_factor = 1;

% threshold for defect
defthres = 0.3;

% threshold for splay and bend cutoff
pttn_cutoff = 5;

% entropy window size
ent_ws = 50;

% load file paths
sp = ["dic", "phi", "coh", "ent", "phalf", "mhalf", "splay", "bend", "aligned"];
[ch, subpath] = load_file_paths();

%% INPUT IMAGES
prefix = ["mCh300FBS10-1"];
welllist = {[1:24]};
framelist = {[1]};


%% Calculate phi and coh
showfigure = 1;

for p = 1:length(prefix)
    PATH = ROOTDIR + "/wave_" + prefix(p) + "/";

    for w = welllist{p}
        outpath = PATH + "Well" + sprintf("%01d", w) + "/" + subpath("phi") + "/";
        
        mkdir(outpath);
        
        for fr = framelist{p}
            fprintf("prefix: " + prefix(p) + " well %02d, f %02d\n", w, fr)
            
            imfile = ROOTDIR + "/wave_" + prefix(p) + "/Well" + sprintf("%02d", w) + "/" + subpath("dic") + "/" + ...
                prefix(p)+"_s"+sprintf("%02d", w)+"t"+sprintf("%02d",fr)+ch("dic")+"_ORG.tif";

            im = imread(imfile);
            imsize = size(im);

            [coh, phi] = calculate_tensor(imfile);

            coh_outfile = outpath + "/" + ...
                prefix(p)+"_s"+sprintf("%02d", w)+"t"+sprintf("%02d",fr)+ch("dic")+"_ORG_coh.mat";
            save(coh_outfile, 'coh');

            phi_outfile = outpath + "/" + ...
                prefix(p)+"_s"+sprintf("%02d", w)+"t"+sprintf("%02d",fr)+ch("dic")+"_ORG_phi.mat";
            save(phi_outfile, 'phi');

            if showfigure
                figure
                tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
                imshow(phi, [-pi/2 pi/2], 'Border', 'tight');
                colormap hsv
                axis square off

                phi_outfile = outpath + "/" + ...
                    prefix(p)+"_s"+sprintf("%01d", w)+"t"+sprintf("%03d",fr)+ch("dic")+"_ORG_phi.jpg";
                exportgraphics(gcf, phi_outfile)

                close all
            end

            % draw nematic directors
            outfile = outpath + prefix(p)+"_s" + sprintf("%02d", w) + "t"+sprintf("%02d",fr)+ch("phi")+"_ORG_nem"+nem_ws+".svg";
            draw_nematic_directors(phi, "", size(phi), nem_ws, 2, nem_color, 0, 1, outfile);

        end
    end
end

%% run entropy
showfigure = 1;

for p = 1:length(prefix)
    for w = welllist{p}
        PATH = ROOTDIR + "/wave_" + prefix(p) + "/" + "Well" + sprintf("%02d", w) + "/";
        outpath = PATH + subpath("phi") + "/";
        mkdir(outpath);
        
        for fr = 1
            fprintf("prefix " + prefix(p) + " well %02d, fr %02d\n", w, fr);
            outpath = PATH + "/" + subpath("phi") + "/";
            mkdir(outpath);
            
            phi_file = outpath + prefix(p)+"_s" + sprintf("%02d", w) + "t"+sprintf("%02d",fr)+ch("ent")+"_ORG_phi.mat";

            H = calculate_entropy(phi_file);
            
            if showfigure
                ent = uint16(normalize(H(:), 'range')*2^16);
                ent = reshape(ent, size(H));
                figure
                imshow(ent, 'Border', 'tight')
                colormap jet
        
                exportgraphics(gcf, outpath + prefix(p)+"_s" + sprintf("%02d",w) + "t"+sprintf("%02d",fr)+ch("ent")+"_ORG_ent.jpg");
                close all
            end

            ent_outfile = outpath + prefix(p)+"_s" + sprintf("%02d",w) + "t"+sprintf("%02d",fr)+ch("ent")+"_ORG_ent.mat";
            save(ent_outfile, 'H');
        end
    end
end

%% calculate defects

for p = 1:length(prefix)
    for w = welllist{p}
        PATH = ROOTDIR + "/wave_" + prefix(p) + "/" + "Well" + sprintf("%02d", w) + "/";
        outpath = PATH + "/" + subpath("phi") + "/";
        mkdir(outpath);
        
        for fr = 1
            fprintf("prefix " + prefix(p) + " well %02d, fr %02d\n", w, fr);
            outpath = PATH + "/" + subpath("phi") + "/";
            mkdir(outpath);
            
            phi_file = outpath + prefix(p)+"_s" + sprintf("%02d", w) + "t"+sprintf("%02d",fr)+ch("phi")+"_ORG_phi.mat";

            [mhalf, phalf] = calculate_defects(phi_file);
            
            % wite half defects to file
            mhalf_outfile = outpath + set_filename(prefix(p), w, fr, ch("mhalf")) + "_mhalf.txt";
            writematrix(mhalf, mhalf_outfile);
    
            phalf_outfile = outpath + set_filename(prefix(p), w, fr, ch("phalf")) + "_phalf.txt";
            writematrix(phalf, phalf_outfile);
        end
    end
end

%% show defects on DIC on a figure
for p = 1:length(prefix)
    for w = welllist{p}
        PATH = ROOTDIR + "/wave_" + prefix(p) + "/" + "Well" + sprintf("%02d", w) + "/";

        for fr = 1
            dic_file = PATH + subpath("dic") + "/" + ...
                prefix(p)+"_s" + sprintf("%02d", w) + "t"+sprintf("%02d",fr)+ch("dic")+"_ORG.tif";
            dic = imread(dic_file);

            phalf_file = PATH + subpath("phalf") + "/" + ...
                prefix(p)+"_s" + sprintf("%02d", w) + "t"+sprintf("%02d",fr)+ch("phalf")+"_ORG_phalf.txt";
            phalf = readmatrix(phalf_file);

            mhalf_file = PATH + subpath("mhalf") + "/" + ...
                prefix(p)+"_s" + sprintf("%02d", w) + "t"+sprintf("%02d",fr)+ch("mhalf")+"_ORG_mhalf.txt";
            mhalf = readmatrix(mhalf_file);

            figure
            imshow(dic, [min(dic(:)) max(dic(:))])
            hold on
            plot(phalf(:,1), phalf(:,2), 'r.', 'MarkerSize', 30)
            plot(mhalf(:,1), mhalf(:,2), 'b.', 'MarkerSize', 30)
            hold off
        end
    end
end

%% calculate splay and dend
showfigure = 0;

phi_diff = cell(length(prefix), all_well);
splay = cell(length(prefix), all_well);
bend = cell(length(prefix), all_well);

for p = 1:length(prefix)
    for w = welllist{p}
        PATH = ROOTDIR + "/wave_" + prefix(p) + "/" + "Well" + sprintf("%02d", w) + "/";
        outpath = PATH + "Well" + sprintf("%02d", w) + "/" + subpath("phi") + "/";
        mkdir(outpath);
        
        for fr = 1
            fprintf("prefix " + prefix(p) + " well %02d, fr %02d\n", w, fr);
            outpath = PATH + "/" + subpath("phi") + "/";
            mkdir(outpath);
            
            phi_file = outpath + prefix(p)+"_s" + sprintf("%02d", w) + "t"+sprintf("%02d",fr)+ch("phi")+"_ORG_phi.mat";
            
            [b, s, pdiff] = calculate_top_pattern(phi_file);

            if showfigure
                figure
                imshow(pdiff, [min(pdiff(:)) max(pdiff(:))], 'Border', 'tight')
                colormap jet
                
                figure
                imshow(s, [min(s(:)) max(s(:))], 'Border', 'tight')
                colormap jet
                
                figure
                imshow(b, [min(b(:)) max(b(:))], 'Border', 'tight')
                colormap jet
            end

            % store nematic orderings
            pdiff_outfile = outpath + set_filename(prefix(p), w, fr, ch("dic")) + "_phidiff_"+pttn_cutoff+".mat";
            save(pdiff_outfile, 'pdiff');

            splay_outfile = outpath + set_filename(prefix(p), w, fr, ch("splay")) + "_splay_"+pttn_cutoff+".mat";
            save(splay_outfile, 's');

            bend_outfile = outpath + set_filename(prefix(p), w, fr, ch("bend")) + "_bend_"+pttn_cutoff+".mat";
            save(bend_outfile, 'b');
        end
    end
end


