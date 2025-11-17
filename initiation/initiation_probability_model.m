%
%   probablity model for initiation
%

clearvars
clc
close all

ROOTDIR = "/home/jenhao/";

addpath /home/jenhao/code/util/
addpath /home/jenhao/code/Matlab_FileExchange/swtest/

warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% Configuration
welllist = {[1 4 5 8 9 12], [1:3 5 7:10 12], 1:12, 1:12}; % [2:4 6 9:11]};

imsize = [5000 5000];
ini_size = 150;
offset = 51; 

%% load computed circlar spatial properties for initiation and non-ini background
load(ROOTDIR + 'initiation/initiation_sp/ini_coh.mat');
load(ROOTDIR + 'initiation/initiation_sp/ini_den.mat');
load(ROOTDIR + 'initiation/initiation_sp/ini_ent.mat');

load(ROOTDIR + 'initiation/initiation_sp/bg_coh_rep10.mat');
load(ROOTDIR + 'initiation/initiation_sp/bg_den_rep10.mat');
load(ROOTDIR + 'initiation/initiation_sp/bg_ent_rep10.mat');

load(ROOTDIR + 'initiation/initiation_sp/resist_coh_rep11.mat');
load(ROOTDIR + 'initiation/initiation_sp/resist_den_rep11.mat');
load(ROOTDIR + 'initiation/initiation_sp/resist_ent_rep11.mat');

ini_coh_data = vertcat(ini_coh{:});
ini_den_data = vertcat(ini_den{:});
ini_ent_data = vertcat(ini_ent{:});
ini_coh_data = ini_coh_data(:);
ini_den_data = ini_den_data(:);
ini_ent_data = ini_ent_data(:);
ini_coh_data(isnan(ini_coh_data)) = [];
ini_den_data(isnan(ini_den_data)) = [];
ini_ent_data(isnan(ini_ent_data)) = [];

nini_coh_data = vertcat(bg_coh{:});
nini_den_data = vertcat(bg_den{:});
nini_ent_data = vertcat(bg_ent{:});
nini_coh_data = nini_coh_data(:);
nini_den_data = nini_den_data(:);
nini_ent_data = nini_ent_data(:);
nini_coh_data(isnan(nini_coh_data)) = [];
nini_den_data(isnan(nini_den_data)) = [];
nini_ent_data(isnan(nini_ent_data)) = [];

% % load spontaneous initiations
% load('D:\Spatiotemporal_analysis\initiation_patterns\ini_rand_spatial_patterns\ini_spon_coh.mat')
% load('D:\Spatiotemporal_analysis\initiation_patterns\ini_rand_spatial_patterns\ini_spon_den.mat')
% load('D:\Spatiotemporal_analysis\initiation_patterns\ini_rand_spatial_patterns\ini_spon_ent.mat')

% load(ROOTDIR + 'initiation/initiation_sp/resist_ent_rep11.mat');
% load(ROOTDIR + 'initiation/initiation_sp/resist_coh_rep11.mat');
% load(ROOTDIR + 'initiation/initiation_sp/resist_den_rep11.mat');
% 
% % concate photoinduction and spontaneous initiations
% ini_coh = [ini_coh; s_ini_coh];
% ini_den = [ini_den; s_ini_den];
% ini_ent = [ini_ent; s_ini_ent];

%% test normality for the spatial properties
% results, null hypothesis is rejected using ks test, but ploting cdfplot
% shows they are very similar. 

[h, p] = kstest(ini_coh_data); % null hypothesis rejected
[h, p] = kstest(ini_den_data); % null hypothesis rejected
[h, p] = kstest(ini_ent_data); % null hypothesis rejected
[h, p] = kstest(nini_coh_data);  % null hypothesis rejected   
[h, p] = kstest(nini_den_data);  % null hypothesis rejected   
[h, p] = kstest(nini_ent_data);  % null hypothesis rejected   
     
% plot cdf to check if data satisfies normality assumption
figure
cdfplot(ini_ent_data)
hold on
x_values = linspace(min(ini_ent_data),max(ini_ent_data));
plot(x_values,normcdf(x_values,mean(ini_ent_data),std(ini_ent_data)),'r-')

%% calculate p(ini) and p(not_ini)
num_ini = 3172 - 120; % (total ini - those in the edges), 3172 is photoinduction
ini_radius = fix(ini_size/2);
total_area = imsize(1) - offset;
num_wells = sum(cellfun(@numel, welllist));

p_ini = num_ini*(pi*ini_radius^2) / (total_area^2*num_wells);
p_nini = 1 - p_ini;

log_p_ini = log(p_ini);
log_p_nini = log(p_nini);

%% Use normal distribution
% 
% contruct pdf for sp for coh/ent and den assuming they fit normal distribution
%

aln_type = "coh";

if strcmp(aln_type, "ent")
    ini_aln_data = ini_ent_data;
    nini_aln_data = nini_ent_data;
    
    aln_interval = 0.01; %0.0001;
    min_aln = -6.5;
    max_aln = 1.5;
elseif strcmp(aln_type, "coh")
    ini_aln_data = ini_coh_data;
    nini_aln_data = nini_coh_data;

    aln_interval = 0.001; %0.0001;
    min_aln = 0;
    max_aln = 1;
end

den_interval = 0.1;

min_den = 0;
max_den = 180;

mu_ini = [mean(ini_aln_data) mean(ini_den_data)];
mu_nini = [mean(nini_aln_data) mean(nini_den_data)];

sigma_ini = cov(ini_aln_data, ini_den_data);
sigma_nini = cov(nini_aln_data, nini_den_data);

%% use log to calculate the probability
prob_map_outfile = ROOTDIR + "/initiation/initiation_probability/" + "norm_probmap_" + aln_type + "den.mat";

prob_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
tic
for curr_aln = min_aln:aln_interval:max_aln
    for curr_den = min_den:den_interval:max_den
        key = sprintf('%04f_%01f', curr_aln, curr_den);
        
        % calculate law of total prob
        log_p_xy_given_i = log_bivariate_normal_pdf(curr_aln, curr_den, mu_ini, sigma_ini);
        log_p_xy_given_ni = log_bivariate_normal_pdf(curr_aln, curr_den, mu_nini, sigma_nini);
        % p_xy_given_i = bivariate_normal_pdf(curr_aln, curr_den, mu_ini, sigma_ini);
        % p_xy_given_ni = bivariate_normal_pdf(curr_aln, curr_den, mu_nini, sigma_nini);

        % Compute log P(X = new_x, Y = new_y)
        log_p_xy = logsumexp([log_p_xy_given_i + log_p_ini, log_p_xy_given_ni + log_p_nini]);
        % p_xy = p_xy_given_i * p_ini + p_xy_given_ni * p_nini;

        log_prob = log_p_xy_given_i + log_p_ini - log_p_xy;
        % prob = (p_xy_given_i * p_ini) ./ p_xy;

        prob_map(key) = exp(log_prob);
        % prob_map(key) = (prob);
    end
end
toc

save(prob_map_outfile, 'prob_map');

%% Empirical distribution
%
% construct pdf for sp for coh, for empirical distirbution
% P(e, d) = P(e, d|ini)(P(ini) + P(e, d|not_ini)P(not_ini);
%

aln_type = "coh"; % = "ent";
pdf_ini_outfile  = ROOTDIR + "/initiation/initiation_probability/emp_pdf_ini_" + aln_type + "den.mat";
pdf_nini_outfile = ROOTDIR + "/initiation/initiation_probability/emp_pdf_nini_" + aln_type + "den.mat";
prob_map_outfile = ROOTDIR + "/initiation/initiation_probability/" + "emp_probmap_" + aln_type + "den_bw0.05.mat";

min_den = 0;
max_den = 180;
den_interval = 0.1;

if strcmp(aln_type, "ent")
    ini_aln_data = ini_ent_data;
    nini_aln_data = nini_ent_data;

    min_aln = -6.5;
    max_aln = 1.5;
    aln_interval = 0.01;
elseif strcmp(aln_type, "coh")
    ini_aln_data = ini_coh_data;
    nini_aln_data = nini_coh_data;

    min_aln = 0;
    max_aln = 1;
    aln_interval = 0.001;
end

ini_data = [ini_aln_data ini_den_data];
nini_data = [nini_aln_data nini_den_data(:)];

%% calcualte pdf
% define space to evaluate the kde
grid_aln = linspace(min([ini_aln_data; nini_aln_data]), max([ini_aln_data; nini_aln_data]), 1000);
grid_den = linspace(min([ini_den_data; nini_den_data]), max([ini_den_data; nini_den_data]), 1000);

[gridX, gridY] = meshgrid(grid_aln, grid_den);
gridPoints = [gridX(:), gridY(:)];
disp("done calculating gridPoints");

tic
[pdf_ini, ~] = ksdensity(ini_data, gridPoints, 'Function', 'pdf', 'Bandwidth', 0.05);
toc
disp("done calculating pdf for ini_data");
tic
[pdf_nini, ~] = ksdensity(nini_data, gridPoints, 'Function', 'pdf', 'Bandwidth', 0.05);
toc
disp("done calculating pdf for not_ini_data");

pdf_ini = reshape(pdf_ini, size(gridX));
pdf_nini = reshape(pdf_nini, size(gridX));

save(pdf_ini_outfile, 'pdf_ini');
save(pdf_nini_outfile, 'pdf_nini');
% 
% load(ROOTDIR + "/initiation/initiation_probability/emp_pdf_ini_" + aln + "den.mat")
% load(ROOTDIR + "/initiation/initiation_probability/emp_pdf_nini_" + aln + "den.mat"')

%% code for calculating bayes theorem and law of total prob, using log

prob_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

tic
for curr_aln = min_aln:aln_interval:max_aln
    for curr_den = min_den:den_interval:max_den
        key = sprintf('%04f_%01f', curr_aln, curr_den);
        
        % calculate law of total prob
        p_xy_given_i  = interp2(grid_aln, grid_den, pdf_ini, curr_aln, curr_den, 'spline');
        p_xy_given_ni = interp2(grid_aln, grid_den, pdf_nini, curr_aln, curr_den, 'spline');
        
        epsilon = 1e-300;

        p_xy_given_i(p_xy_given_i<epsilon) = epsilon;
        p_xy_given_ni(p_xy_given_ni<epsilon) = epsilon;

        log_p_xy_given_i  = log(p_xy_given_i);
        log_p_xy_given_ni = log(p_xy_given_ni);

        % p_xy = p_xy_given_i * p_ini + p_xy_given_ni * p_nini;
        log_p_xy = logsumexp([log_p_xy_given_i + log_p_ini, log_p_xy_given_ni + log_p_nini]);

        % compute p(ini|a, d)
        % if p_xy == 0
        if log_p_xy == 0
            prob_map(key) = 0;
        else
            % p_ini_given_xy = (p_xy_given_i * p_ini) ./ p_xy;
            log_p_ini_given_xy = log_p_xy_given_i + log_p_ini - log_p_xy;

            % prob_map(key) = p_ini_given_xy;
            prob_map(key) = exp(log_p_ini_given_xy);
        end
    end
end
toc

save(prob_map_outfile, 'prob_map')

%%
% 
% visualize the probability model
% 

ent_file = ROOTDIR + "/lipid_perox/wave_mCh300FBS10-1/Well10/mCh300FBS10-1_s10t28c3_ORG_ent.mat";
coh_file = ROOTDIR + "/lipid_perox/wave_mCh300FBS10-1/Well10/mCh300FBS10-1_s10t28c3_ORG_coh.mat";
den_file = ROOTDIR + "/lipid_perox/wave_mCh300FBS10-1/Well10/mCh300FBS10-1_s10t28c1_ORG_den120_5000.mat";

aln_type = "coh";

% load('/home/jenhao/initiation/initiation_probability/emp_probmap_cohden.mat');
% load('/home/jenhao/initiation/initiation_probability/archive/prob_map_log.mat');
load('/home/jenhao/initiation/initiation_probability/norm_probmap_' + aln_type + 'den.mat');

if strcmp(aln_type, "ent")
    aln_file = ent_file;
    aln_interval = 0.01;
elseif strcmp(aln_type, "coh")
    aln_file = coh_file;
    aln_interval = 0.001;
end

aln = struct2array(load(aln_file));
den = struct2array(load(den_file));

den = den(2:end-1, 2:end-1);
den_interval = 0.1;

aln = ceil(aln./aln_interval) .* aln_interval;
den = ceil(den./den_interval) .* den_interval;

%%
disp("calculating probabilities for the selected image")

prob_size = [4998 4998];

x_start = 100;
x_end = 2200 + x_start;

y_start = 50;
y_end = 4105+y_start;

prob = nan(prob_size);

for xi = x_start:x_end
    for yi = y_start:y_end
        curr_aln = aln(xi, yi);
        curr_den = den(xi, yi);

        % % prob(xi, yi) = calc_probability(grid_aln, grid_den, pdf_ini, pdf_nini, curr_aln, curr_den, p_ini, p_nini);
        curr_den(curr_den<16.3) = 16.3;
        curr_den(curr_den>171) = 171;

        if strcmp(aln_type, "coh")
            curr_aln(curr_aln<1.0000e-04) = 1.0000e-04;
            curr_aln(curr_aln>0.687) = 0.687;
        elseif strcmp(aln_type, "ent")
            curr_aln(curr_aln<-4.64) = -4.64;
            curr_aln(curr_aln>1.15) = 1.15;
        end

        if curr_aln == 0
            curr_aln = 0;
        end
        key = sprintf('%04f_%01f', curr_aln, curr_den);
        prob(xi, yi) = prob_map(key);
        % if isnan(prob_map(key))
        %     disp(xi + "," + yi + ";"+curr_aln+","+curr_den+";"+prob_map(key)+";")
        %     break;
        % end
    end
    % if isnan(prob_map(key))
    %     break;
    % end
end
%%
prob_im = imcrop(prob, [y_start x_start 4105 2200]);
upr = mean(prob_im(:), 'omitnan') + 3*std(prob_im(:), 'omitnan');
prob_im = imgaussfilt(prob_im,20);
prob_im = imresize(prob_im, [550 1026]);
figure
% set(f, 'Position', [10, 1500, 1026, 550]);
imshow(prob_im, [0 upr], 'Border', 'tight')
% clim([0, 0.2])
colormap parula

%%
im = imread('/home/jenhao/initiation/initiation_probability/emp_probmap_log_entden.jpg');
% im = imgaussfilt(im, 10);
lwr = mean(im(:), 'omitnan') - 1*std(double(im(:)), 'omitnan');
upr = mean(im(:), 'omitnan') + 1*std(double(im(:)), 'omitnan');
figure;
imshow(im, [lwr upr], 'Border', 'tight')
colormap parula

%% FUNCTION

function p = bivariate_normal_pdf(x, y, mu, Sigma)
    % x and y are scalar values
    % mu is a 1x2 vector [mu_x, mu_y]
    % Sigma is a 2x2 covariance matrix
    diff = [x - mu(1), y - mu(2)];
    p = (1 / (2 * pi * sqrt(det(Sigma)))) * ...
        exp(-0.5 * diff * (Sigma \ diff'));
end

function log_p = log_bivariate_normal_pdf(x, y, mu, Sigma)
    % Computes the logarithm of the bivariate normal PDF
    diff = [x - mu(1), y - mu(2)];
    log_det_Sigma = log(det(Sigma));
    inv_Sigma = inv(Sigma);
    log_p = -0.5 * (log_det_Sigma + diff * inv_Sigma * diff' + 2 * log(2 * pi));
end

function s = logsumexp(a) %, dim)
    % Stable computation of log(sum(exp(a), dim))
    max_a = max(a); %, [], dim);
    s = max_a + log(sum(exp(a - max_a))); %, dim));
end
