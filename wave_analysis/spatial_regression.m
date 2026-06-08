%
%   Spatial regression for relative contribution
%

clearvars
clc
close all

ROOTDIR = "D:/Spatiotemporal_analysis/";
MATLAB_FX = "D:/Matlab_FileExchange/";

addpath(ROOTDIR + "/code");
addpath(ROOTDIR + "/code/util");

warning('off', 'MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:singularMatrix')

%% CONFIGURATION

PATH = ROOTDIR + "/wave_0_analyses/";

den_reg = 120;
MCHCH = "c1";
DENPATH = MCHCH + "_mCherry_density/";

CY5CH = "c2";
CY5ADJPATH = CY5CH + "_cy5_adj/";

DICCH = "c3";
NEMPATH = DICCH + "_DIC_nematics/";
ORIENTPATH = NEMPATH + "orientation_files/";

IMSIZE = 5000;

%% load wave information table

wvinfo_fname = "all_wave_info.txt";

wvinfo = readtable(PATH + wvinfo_fname);

% exclude waves that will not be analyzed
wvinfo(wvinfo.include == 0, :) = [];

%% load BO
showfigure = false;

wv_bo = cell(height(wvinfo), 1);

if showfigure, figure; end

for ii = 1:height(wvinfo)
    prefix = wvinfo.Set(ii);
    w = wvinfo.well(ii);
    wv = wvinfo.wave(ii);
    fr = wvinfo.spatial_frame(ii);
    frame_range = wvinfo.start_frame(ii)+"-"+wvinfo.end_frame(ii);
    
    path = ROOTDIR + "wave_" + prefix + "/" + "Well" + sprintf('%02d', w) + "/";

    wv_bo_file = path + "/wave" + wv + "_analyses/" + ...
            set_filename(prefix, w, frame_range, '') + "_bo.mat";
    wv_bo{ii} = struct2array(load(wv_bo_file));

    if showfigure
        nexttile
        imshow(ones(IMSIZE))
        hold on
        Bo = wv_bo{ii};
        color_bo = turbo(length(Bo));
        for b = 1:length(Bo)
            bo = Bo{b};
            plot(bo(:, 2), bo(:, 1), 'Color', color_bo(b,:))
        end
        title("wave "+ii,'color','w')
        hold off
    end
end

%% load interpolated wave information
showfigure = false;

ws = 59; % 75 um

wave = cell(height(wvinfo), 1);

all_spd = cell(height(wvinfo), 1);
all_den = cell(height(wvinfo), 1);
all_ang = cell(height(wvinfo), 1);

for ii = 1:height(wvinfo)
    prefix = wvinfo.Set(ii);
    w = wvinfo.well(ii);
    wv = wvinfo.wave(ii);
    fr = wvinfo.spatial_frame(ii);
    frame_range = wvinfo.start_frame(ii)+"-"+wvinfo.end_frame(ii);
    

    path = ROOTDIR + "wave_" + prefix + "/" + "Well" + sprintf('%02d', w) + "/";

    wv_name = path + "/wave" + wv + "_analyses/" + ...
        set_filename(prefix, w, frame_range, '') + "_wvdata_step" + ws + "_new.mat";


    wave{ii} = struct2array(load(wv_name));


    if sum(ismember(wave{ii}.Properties.VariableNames, 'angden'))
        wave{ii}.angden = [];
    end
    if sum(ismember(wave{ii}.Properties.VariableNames, 'coh'))
        wave{ii}.coh = [];
    end

    resz = ceil(IMSIZE/ws);

    spd = nan(resz, resz);
    den = nan(resz, resz);
    ang = nan(resz, resz);

    x = round((wave{ii}.xq - 1) / ws) + 1;
    y = round((wave{ii}.yq - 1) / ws) + 1;

    idx = sub2ind([resz, resz], x, y);

    spd(idx) = wave{ii}.speed;
    den(idx) = wave{ii}.density;
    ang(idx) = wave{ii}.angdiff;

    if showfigure
        figure
        nexttile
        imagesc([min(wave{ii}.xq) max(wave{ii}.xq)],[min(wave{ii}.yq) max(wave{ii}.yq)],spd)
        colormap(gca,"parula")
        clim([0 200])
        title("speed")
        xlim([0.5 5000.5])
        ylim([0.5 5000.5])
        axis square
    
        nexttile
        imagesc([min(wave{ii}.xq) max(wave{ii}.xq)],[min(wave{ii}.yq) max(wave{ii}.yq)],den)
        colormap(gca,"jet")
        clim([20 120])
        title("density")
        xlim([0.5 5000.5])
        ylim([0.5 5000.5])
        axis square

        nexttile
        imagesc([min(wave{ii}.xq) max(wave{ii}.xq)],[min(wave{ii}.yq) max(wave{ii}.yq)],ang)
        colormap(gca,"parula")
        clim([0, pi/2])
        title("angle difference")
        xlim([0.5 5000.5])
        ylim([0.5 5000.5])
        axis square
    end

    all_spd{ii} = spd;
    all_den{ii} = den;
    all_ang{ii} = ang;
end

%% load direction
showfigure = false;
all_Gdir = cell(height(wvinfo), 1);

for ii = 1:height(wvinfo)
    prefix = wvinfo.Set(ii);
    w = wvinfo.well(ii);
    wv = wvinfo.wave(ii);
    fr = wvinfo.spatial_frame(ii);
    frame_range = wvinfo.start_frame(ii)+"-"+wvinfo.end_frame(ii);

    path = ROOTDIR + "wave_" + prefix + "/" + "Well" + sprintf('%02d', w) + "/";
    gdir_file = path + "/wave" + wv + "_analyses/" + ...
            set_filename(prefix, w, frame_range, '') + "_ws" + ws + "_Gdir.mat";

    Gdir = struct2array(load(gdir_file));

    % convert angle to unit vectors first
    Ux = cosd(Gdir);
    Uy = -sind(Gdir);
    
    resz = size(all_spd{ii}, 1);
    
    % resize vector components, not angles
    Ux = imresize(Ux, [resz resz]);
    Uy = imresize(Uy, [resz resz]);
    
    % normalize after resizing
    mag = hypot(Ux, Uy);
    Ux = Ux ./ mag;
    Uy = Uy ./ mag;

    Gdir = atan2d(-Uy, Ux); 
    if showfigure
        step = ws;

        [Xq, Yq] = meshgrid( ...
                    linspace(ws/2, IMSIZE - ws/2, resz), ...
                    linspace(ws/2, IMSIZE - ws/2, resz));
        
        figure
        imshow(zeros(IMSIZE, IMSIZE)+255)
        hold on
        quiver(Xq, Yq, Ux, Uy, 1, 'color', 'k')
        hold off
        axis ij square off
    end

    all_Gdir{ii} = Gdir;
end


%% Run Spatial Durbin Model
% Spatial dependence in speed was constrained to upstream neighbors to 
% respect causal propagation, whereas spatial variables (density and angle) 
% were allowed to influence locally from all directions

debug = 0;

out_all = cell(height(wvinfo), 1);

offs = [ -1 0; 1 0; 0 -1; 0 1; -1 -1; -1 1; 1 -1; 1 1 ];

tolDeg = 1;  % keep neighbors within +/- 45 degrees of upstream direction

wrap180 = @(a) mod(a + 180, 360) - 180;  % wrap to (-180, 180]

for n = 1:height(wvinfo)
    spd = all_spd{n};
    den = all_den{n};
    ang = all_ang{n};
    dir = all_Gdir{n};

    % correct cartesian coordinate into image coordinates
    dir = 90-dir';

    % % directiona and others are the same coordinates
    if debug
        figure
        imagesc(imresize(spd, [IMSIZE, IMSIZE], 'nearest'))
        axis square
        hold on
    
        Ux = cosd(dir);
        Uy = sind(dir);
    
        [X, Y] = meshgrid(1:step:IMSIZE, 1:step:IMSIZE);
    
        quiver(X, Y, Ux, Uy, 1, 'color', 'k')
        hold off
    end
    
    [nr, nc] = size(spd);

    mask_wave  = true(nr, nc);

    % only wave-defined regions
    mask_wave  = mask_wave & ~isnan(spd) & ~isnan(ang) & ~isnan(dir);
    N = nnz(mask_wave);

    resz = nr;

    % create look-up table for 2D grid to 1D vector index
    idxMap = zeros(nr, nc, 'int32');
    idxMap(mask_wave) = int32(1:N);

    % build response variable
    y = spd(mask_wave);
    y = y(:);

    % build design matrix X
    x_den = den(mask_wave);
    x_den = x_den(:);
    x_ang = ang(mask_wave);
    x_ang = x_ang(:);

    % make angle into cosine (why?)
    x_ang = cos(2*x_ang);
    
    x_den = zscore(x_den);
    x_ang = zscore(x_ang);

    x_int = x_den .* x_ang;
    
    intercept = ones(N, 1);

    X = [intercept, x_den, x_ang, x_int];
    X_noInt = [intercept, x_den, x_ang];

    % construct spatial weights W on grid
    ii_x = []; jj_x = [];
    ii_y = []; jj_y = []; vv_y = [];

    % extract all valid pixels
    [rr, cc] = find(mask_wave);

    % loop over each valid piexel
    for k = 1:numel(rr)
        r = rr(k); c = cc(k);

        % convert grid locatoin to vector index
        i = idxMap(r,c);

        phi = dir(r, c);             % local direction at pixel i (deg, -180..180)
        phi_up = wrap180(phi + 180); % upstream direction (where wave came from)

        % loop over neighbor offset
        for o = 1:size(offs, 1) 
            % compute neighbor location
            r2 = r + offs(o, 1);
            c2 = c + offs(o,2);

            if r2>=1 && r2<=nr && c2>=1 && c2<=nc && mask_wave(r2,c2) % check if neighbor is valid
                % compute neighbor to vector index
                j = idxMap(r2,c2);

                % % build W_x: include all valide neighbors
                ii_x(end+1, 1) = i;
                jj_x(end+1, 1) = j;

                % build W_y: includ eonly upstream neighbors
                dx = c2 - c;
                dy = r2 - r;

                bearing = atan2d(dy, dx); % neighbor bearing in image coords (deg)

                d = abs(wrap180(bearing - phi_up));

                if d <= tolDeg
                    % use weights to determine the upstream angle contribution
                    w = cosd(d);

                    if w > 0
                        % store the neighbor
                        ii_y(end+1, 1) = i;
                        jj_y(end+1, 1) = j;
                        vv_y(end+1, 1) = w; % weights
                    end
                end
            end
        end
    end

    % compute density neighbor across ALL 8 neighbors (inside or outside wave)
    den_nb = nan(N,1);
    for k = 1:N
        vals = [];
        r = rr(k); c = cc(k);
    
        for o = 1:size(offs,1)
            r2 = r + offs(o,1);
            c2 = c + offs(o,2);
    
            if r2>=1 && r2<=nr && c2>=1 && c2<=nc
                if ~isnan(den(r2,c2))
                    vals(end+1,1) = den(r2,c2);
                end
            end
        end
    
        if isempty(vals)
            den_nb(k) = NaN;
        else
            den_nb(k) = mean(vals);
        end
    end
    
    % standardize neighbor density term
    mu = mean(den_nb, 'omitnan');
    sd = std(den_nb, [], 'omitnan');
    den_nb = (den_nb - mu) ./ sd;
    den_nb(isnan(den_nb)) = 0;
    WX_den = den_nb;

    % adjacency (binary weights)
    W_x = sparse(ii_x, jj_x, 1, N, N);
    W_y = sparse(ii_y, jj_y, vv_y, N, N);

    if debug
        figure;spy(W_x)
        title("C_x")
        copygraphics(gcf,'ContentType','image','Resolution',300)

        figure;spy(W_y)
        title("C_y")
        copygraphics(gcf,'ContentType','image','Resolution',300)
    end

    % remove self-loops if any
    W_x = W_x - spdiags(diag(W_x), 0, N, N);
    W_y = W_y - spdiags(diag(W_y), 0, N, N);

    % Row-standardize so each row sums to 1 (common in SDM/SLM/SAR)
    rs_x = sum(W_x, 2);
    rs_x(rs_x==0) = 1; % isolated pixels (edge of mask) -> keep row all zeros
    W_x = spdiags(1./rs_x, 0, N, N) * W_x;

    rs_y = sum(W_y, 2);
    rs_y(rs_y==0) = 1; % isolated pixels (edge of mask) -> keep row all zeros
    W_y = spdiags(1./rs_y, 0, N, N) * W_y;


    % fit SDM
    % precompute spatially lagged X
    WX = W_x * X(:, 2:end);

    % lag only angle and interaction within wave
    WX_ang = W_x * X(:,3);
    WX_int = W_x * X(:,4);
    
    % final regressor matrix:
    % [local terms, neighbor density(all 8), neighbor angle(in-wave), neighbor int(in-wave)]
    Z = [X, WX_den, WX_ang, WX_int];
    Z_noInt = [X_noInt, WX_den, WX_ang];

    if rank(full(Z)) < size(Z,2)
        warning("Wave %d: Z is rank deficient (rank %d < %d).", n, rank(full(Z)), size(Z,2));
    end

    % Negative log-likelihood function
    negLL = @(rho) sdm_negloglik(rho, y, Z, W_y);

    % Optimize rho
    % rho bounds (for row-standardized W, eigenvalues in [-1,1] typically)
    rho_hat = fminbnd(negLL, -0.99, 0.99);

    % Given rho_hat, estimate gamma = [beta; theta] by OLS on
    % (I - rho W) y = X beta + W X theta + e
    % A = I - rho W
    A = speye(N) - rho_hat * W_y;
    y_tilde = A * y;
    gamma_hat = Z \ y_tilde;  

    % Residuals in original SDM equation form:
    % (I - rho W) y = Z*gamma + e
    e_hat = y_tilde - Z * gamma_hat;
    sigma2_hat = (e_hat' * e_hat) / (N - size(Z,2));

    % fitted values and residuals on ORIGINAL y scale
    y_hat = A \ (Z * gamma_hat);     % fitted values on original response scale
    resid_y = y - y_hat;             % prediction residuals on original scale
    
    % RMSE
    RMSE = sqrt(mean(resid_y.^2));
    
    % NRMSE
    NRMSE = RMSE / (max(y) - min(y)); % normalize by range


    % Basic fit stats (pseudo-R2 on transformed equation)
    SST = sum((y_tilde - mean(y_tilde)).^2);
    SSE = sum(e_hat.^2);
    R2_tilde = 1 - SSE / SST;

    % AIC
    [LL_sdm, sigma2_sdm_LL, ~, logdetA] = sdm_loglik(rho_hat, y, Z, W_y);
    
    k_sdm = size(Z,2) ...   % gamma params (beta, theta, den_out)
            + 1 ...         % rho
            + 1;            % sigma2
    
    AIC_sdm = -2*LL_sdm + 2*k_sdm;
    
    % OLS baseline log-likelihood + AIC (no spatial lag, no WX terms)
    X_ols = [X];   % same non-spatial covariates as SDM
    [LL_ols, sigma2_ols, e_ols] = ols_loglik(y, X_ols);
    
    k_ols = size(X_ols,2) + 1;  % beta + sigma2
    AIC_ols = -2*LL_ols + 2*k_ols;

    % Reduced SDM
    negLL_noInt = @(rho) sdm_negloglik(rho, y, Z_noInt, W_y);
    rho_noInt = fminbnd(negLL_noInt, -0.99, 0.99);
    
    [LL_noInt, ~, ~, logdetA_noInt] = sdm_loglik(rho_noInt, y, Z_noInt, W_y);
    
    k_noInt = size(Z_noInt,2) + 1 + 1;
    AIC_noInt = -2*LL_noInt + 2*k_noInt;

    % ---- fitted values and residuals for NO-INTERACTION model on ORIGINAL y scale ----
    A_noInt = speye(N) - rho_noInt * W_y;
    y_tilde_noInt = A_noInt * y;
    
    gamma_noInt = Z_noInt \ y_tilde_noInt;
    
    y_hat_noInt = A_noInt \ (Z_noInt * gamma_noInt);
    resid_y_noInt = y - y_hat_noInt;
    
    RMSE_noInt = sqrt(mean(resid_y_noInt.^2));
    NRMSE_noInt = RMSE_noInt / (max(y) - min(y));
  

    % calculate delta AIC
    dAIC_int = AIC_noInt - AIC_sdm;   
    
    % Residual Moran’s I (show OLS vs SDM residual autocorrelation)
    % For SDM, e_hat is the innovation error from (I-rhoW)y = Zg + e
    % Use W_y (or whichever W you consider the spatial structure).
    [I_sdm, pI_sdm] = moransI_perm(e_hat, W_y, 999, n);    % seed with wave index
    
    % For comparison, do OLS residual Moran’s I too (optional but useful)
    [I_ols, pI_ols] = moransI_perm(e_ols, W_y, 999, n);

    % calculate partial R2
    SSE_full = SSE;
    p_full   = size(Z, 2);

    fit_reduced = @(keepCols) deal( ...
        Z(:,keepCols) \ y_tilde, ...
        y_tilde - Z(:,keepCols) * (Z(:,keepCols) \ y_tilde) ...
    );

    % % define groups to drop (by Z column indices)
    grp.den = [2 5];    % direct den + lag(den)
    grp.ang = [3 6];    % direct ang + lag(ang)
    grp.int = [4 7];    % direct interaction + lag(interaction)

    grp.den_int = [2 5 4 7];   % direct den & int + lag(den & int)
    grp.ang_int = [3 6 4 7];   % direct ang & int + lag(ang & int)

    grp.lag = [5 6 7];       % lag(den & ang & int)
    grp.den_local = [2 4];   % direct den & int
    grp.den_lag = [5 7];     % lag den
    grp.ang_local = [3 4];     % direct ang & int
    grp.ang_lag = [6 7];       % lag ang

    groupNames = fieldnames(grp);
    
    partialR2 = struct();
    deltaR2   = struct();
    SSE_red   = struct();
    
    for gi = 1:numel(groupNames)
        name = groupNames{gi};
        dropCols = grp.(name);
    
        keepCols = setdiff(1:p_full, dropCols);
    
        % fit reduced
        [gamma_red, e_red] = fit_reduced(keepCols);
        SSEr = sum(e_red.^2);
    
        SSE_red.(name) = SSEr;
    
        % partial R2 (given the rest)
        partialR2.(name) = (SSEr - SSE_full) / SSEr;
    
        % delta R2 in your pseudo-R2 sense (uses same SST as full)
        % compute R2_red with same SST definition (on y_tilde)
        R2_red = 1 - SSEr / SST;
        deltaR2.(name) = R2_tilde - R2_red;
    end

    % results
    pX = size(X,2);

    out = struct();
    out.N = N;
    out.rho = rho_hat;
    out.gamma = gamma_hat;
    out.beta = gamma_hat(1:pX);
    out.theta = gamma_hat(pX+1:end);
    out.sigma2 = sigma2_hat;
    out.R2_tilde = R2_tilde;
    out.Wx = W_x;
    out.Wy = W_y;
    out.mask = mask_wave;
    out.idxMap = idxMap;
    out.X = X;
    out.WX = WX;
    out.y = y;
    out.y_tilde = y_tilde;
    out.resid = e_hat;

    out.y_hat = y_hat;
    out.resid_y = resid_y;
    
    out.RMSE = RMSE;
    out.NRMSE = NRMSE;

    out.RMSE_noInt = RMSE_noInt;
    out.NRMSE_noInt = NRMSE_noInt;

    out.partialR2 = partialR2;   % fractions ~ 0..1
    out.deltaR2   = deltaR2;     % fractions ~ -1..1
    out.SSE_red   = SSE_red;     % raw SSE numbers (like 38627)
    out.SSE_full  = SSE_full;
    out.SST       = SST;

    out.LL_sdm   = LL_sdm;
    out.AIC_sdm  = AIC_sdm;
    out.logdetA  = logdetA;

    out.AIC_noInt = AIC_noInt;

    out.dAIC_int  = dAIC_int;

    out.LL_ols   = LL_ols;
    out.AIC_ols  = AIC_ols;
     
    out.MoranI_sdm = I_sdm;
    out.MoranP_sdm = pI_sdm;
    
    out.MoranI_ols = I_ols;
    out.MoranP_ols = pI_ols;

    out_all{n} = out;
end

%% Full model scalar stats across waves
nWaves = numel(out_all);

% Collect per-wave scalars
N_vec      = nan(nWaves,1);
rho_vec    = nan(nWaves,1);
sigma2_vec = nan(nWaves,1);
R2_vec     = nan(nWaves,1);
RMSE_vec   = nan(nWaves,1);
NRMSE_vec  = nan(nWaves,1);
RMSE_vec_noInt  = nan(nWaves,1);
NRMSE_vec_noInt = nan(nWaves,1);

for n = 1:nWaves
    out = out_all{n};

    if isfield(out,'N');           N_vec(n)      = out.N;        end
    if isfield(out,'rho');         rho_vec(n)    = out.rho;      end
    if isfield(out,'sigma2');      sigma2_vec(n) = out.sigma2;   end
    if isfield(out,'R2_tilde');    R2_vec(n)     = out.R2_tilde; end
    if isfield(out,'RMSE');        RMSE_vec(n)   = out.RMSE;     end
    if isfield(out,'NRMSE');       NRMSE_vec(n)  = out.NRMSE;    end
    if isfield(out,'RMSE_noInt');  RMSE_vec_noInt(n)  = out.RMSE_noInt;  end
    if isfield(out,'NRMSE_noInt'); NRMSE_vec_noInt(n) = out.NRMSE_noInt; end
end

% Summary helper
mean_ = @(x) mean(x, 'omitnan');
sd_   = @(x) std(x, 0, 'omitnan');
ci95_ = @(x) 1.96 * sd_(x) / sqrt(sum(~isnan(x)));

Terms = {'N'; 'rho'; 'sigma2'; 'R2_tilde'; 'RMSE'; 'NRMSE'; 'RMSE_noInt'; 'NRMSE_noInt'};
Mean  = [mean_(N_vec); mean_(rho_vec); mean_(sigma2_vec); mean_(R2_vec); mean_(RMSE_vec); mean_(NRMSE_vec); mean_(RMSE_vec_noInt); mean_(NRMSE_vec_noInt)];
SD    = [sd_(N_vec);   sd_(rho_vec);   sd_(sigma2_vec);   sd_(R2_vec);   sd_(RMSE_vec);   sd_(NRMSE_vec)  ; sd_(RMSE_vec_noInt);   sd_(NRMSE_vec_noInt)];
CI95  = [ci95_(N_vec); ci95_(rho_vec); ci95_(sigma2_vec); ci95_(R2_vec); ci95_(RMSE_vec); ci95_(NRMSE_vec); ci95_(RMSE_vec_noInt); ci95_(NRMSE_vec_noInt)];

T_full = table(Terms, Mean, SD, CI95, ...
    'VariableNames', {'Term','Mean','SD','CI95'});

T_full

%% Residual Moran's I
nWaves = numel(out_all);

MoranI_sdm_vec = nan(nWaves,1);
MoranP_sdm_vec = nan(nWaves,1);

MoranI_ols_vec = nan(nWaves,1);
MoranP_ols_vec = nan(nWaves,1);

for n = 1:nWaves
    out = out_all{n};

    if isfield(out,'MoranI_sdm');  MoranI_sdm_vec(n) = out.MoranI_sdm; end
    if isfield(out,'MoranP_sdm');  MoranP_sdm_vec(n) = out.MoranP_sdm; end

    if isfield(out,'MoranI_ols');  MoranI_ols_vec(n) = out.MoranI_ols; end
    if isfield(out,'MoranP_ols');  MoranP_ols_vec(n) = out.MoranP_ols; end
end

pct_  = @(x) 100 * mean(x, 'omitnan');

alpha = 0.05;

Moran_ols_sig_pct = pct_(MoranP_ols_vec < alpha);
Moran_sdm_sig_pct = pct_(MoranP_sdm_vec < alpha);

% Moran significance drop
pct_moran_drop = pct_((MoranP_ols_vec < alpha) & (MoranP_sdm_vec >= alpha));

T_Moran = table( ...
    mean_(MoranI_ols_vec), mean_(MoranI_sdm_vec), ...
    Moran_ols_sig_pct, Moran_sdm_sig_pct, ...
    pct_(MoranP_ols_vec<alpha), pct_(MoranP_sdm_vec<alpha), pct_moran_drop, ...
    'VariableNames', {...
        'Mean_MoranI_OLS','Mean_MoranI_SDM', ...
        'Pct_OLS_sig','Pct_SDM_sig' ...
        'MoranP_OLS_pct_sig','MoranP_SDM_pct_sig','Pct_OLSsig_to_SDMnonsig' ...
    });

% T_Moran

%% AIC
nWaves = numel(out_all);

% Preallocate
AIC_sdm_vec   = nan(nWaves,1);
AIC_ols_vec   = nan(nWaves,1);
dAIC_vec      = nan(nWaves,1);

AIC_noInt_vec = nan(nWaves, 1);

dAIC_int_vec   = nan(nWaves,1);   % = AIC_noInt - AIC_full

for n = 1:nWaves
    out = out_all{n};

    if isfield(out,'AIC_sdm');     AIC_sdm_vec(n) = out.AIC_sdm; end
    if isfield(out,'AIC_ols');     AIC_ols_vec(n) = out.AIC_ols; end

    if isfield(out,'AIC_sdm') && isfield(out,'AIC_ols')
        dAIC_vec(n) = out.AIC_sdm - out.AIC_ols;
    end

    if isfield(out,'AIC_noInt');  AIC_noInt_vec(n) = out.AIC_noInt; end

    if isfield(out,'dAIC_int');    dAIC_int_vec(n) = out.dAIC_int; end
end

pct_  = @(x) 100 * mean(x, 'omitnan');

T_dAIC = table( ...
    mean_(dAIC_vec), ...
    sd_(dAIC_vec), ...
    ci95_(dAIC_vec), ...
    pct_(dAIC_vec < 0), ...
    'VariableNames', {'Mean_dAIC','SD','CI95','Pct_SDM_better'});

alpha = 0.05;

% AIC evidence thresholds
pct_dAIC_lt_0   = pct_(dAIC_vec < 0);
pct_dAIC_le_m2  = pct_(dAIC_vec <= -2);
pct_dAIC_le_m10 = pct_(dAIC_vec <= -10);

% AIC evidence for interaction: dAIC_int > 0 favors interaction model
T_dAIC_int = table( ...
    mean_(AIC_sdm_vec), ci95_(AIC_sdm_vec), mean_(AIC_noInt_vec), ci95_(AIC_noInt_vec),...
    mean_(dAIC_int_vec), sd_(dAIC_int_vec), ci95_(dAIC_int_vec), ...
    pct_(dAIC_int_vec > 0), pct_(dAIC_int_vec >= 2), pct_(dAIC_int_vec >= 10), ...
    'VariableNames', { ...
        'AIC_sdm_mean', 'AIC_sdm_CI95', 'AIC_noInt_mean', 'AIC_noInt_CI95', ...
        'dAIC_int_mean','dAIC_int_sd','dAIC_int_CI95', ...
        'dAIC_int_pct_gt0','dAIC_int_pct_ge2','dAIC_int_pct_ge10', ...
    } );

T_dAIC_int

%% get the partial and delta R2
% get group names from first wave
groupNames = fieldnames(out_all{1}.partialR2);
nGroups = numel(groupNames);
nWaves  = numel(out_all);

% build matrics
partial_mat = nan(nWaves, nGroups);
delta_mat   = nan(nWaves, nGroups);

for n = 1:nWaves
    for g = 1:nGroups
        name = groupNames{g};
        
        if isfield(out_all{n}.partialR2, name)
            partial_mat(n,g) = out_all{n}.partialR2.(name);
        end
        
        if isfield(out_all{n}.deltaR2, name)
            delta_mat(n,g) = out_all{n}.deltaR2.(name);
        end
    end
end

% average across waves
mean_partial = mean(partial_mat, 1, 'omitnan');
std_partial  = std(partial_mat, 0, 1, 'omitnan');

mean_delta   = mean(delta_mat, 1, 'omitnan');
std_delta    = std(delta_mat, 0, 1, 'omitnan');

% 95% confidence interval
sem_partial = std_partial ./ sqrt(nWaves);
sem_delta   = std_delta   ./ sqrt(nWaves);

CI95_partial = 1.96 * sem_partial;
CI95_delta   = 1.96 * sem_delta;

% summarize
T_summary = table( ...
    groupNames, ...
    mean_partial', std_partial', CI95_partial', ...
    mean_delta',   std_delta',   CI95_delta', ...
    'VariableNames', { ...
        'Term', ...
        'Mean_PartialR2', 'SD_PartialR2', 'CI95_PartialR2', ...
        'Mean_DeltaR2',   'SD_DeltaR2',   'CI95_DeltaR2' ...
    });

% T_summary

%%
function nll = sdm_negloglik(rho, y, Z, W)
    % sdm is                 y = rho*W*y + Z*gamma + e
    % rearrange: (I - rho*W) y = Z*gamma + e
    % define:                A = I - rho*W 
    % get:                 A*y = Z*gamma + e

    % Concentrated (profile) negative log-likelihood for SDM
    N = length(y);
    A = speye(N) - rho * W;
    
    % log|A| using sparse LU (works reasonably up to ~1e5, fine for N<=1e4)
    % For N~1e4, you may want a faster approximation. For 100x100 fully valid, N<=1e4.
    try
        [L,U,P,Q] = lu(A);
        diagU = abs(diag(U));
        diagU(diagU==0) = realmin;
        logdetA = sum(log(diagU));
    catch
        % fallback: penalize invalid rho
        nll = inf;
        return
    end
    
    y_tilde = A * y;
    gamma = (Z' * Z) \ (Z' * y_tilde);
    e = y_tilde - Z * gamma;
    sigma2 = (e' * e) / N;
    
    % Negative concentrated log-likelihood (up to constant)
    nll = (N/2)*log(sigma2) - logdetA;
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

function [LL, sigma2, e, logdetA] = sdm_loglik(rho, y, Z, W)
% Log-likelihood of SDM/SAR-type model:
% (I - rho W) y = Z*gamma + e,  e ~ N(0, sigma2 I)
%
% LL = log|I - rho W| - (N/2)*(log(2*pi*sigma2) + 1)
%
% Uses concentrated (profile) likelihood in gamma, sigma2.

    N = length(y);
    A = speye(N) - rho * W;

    % Solve for gamma (conditional on rho)
    ytilde = A * y;
    gamma  = Z \ ytilde;
    e      = ytilde - Z * gamma;

    sigma2 = (e' * e) / N;

    % log|A| via sparse LU
    [~, U, ~] = lu(A);
    du = diag(U);
    % If A is singular/near-singular at some rho, du may contain zeros.
    if any(~isfinite(du)) || any(abs(du) < 1e-14)
        LL = -Inf; logdetA = -Inf;
        return;
    end
    logdetA = sum(log(abs(du)));

    LL = logdetA - (N/2) * (log(2*pi*sigma2) + 1);
end


function [LL, sigma2, e] = ols_loglik(y, X)
% OLS Gaussian log-likelihood
    N = length(y);
    b = X \ y;
    e = y - X*b;
    sigma2 = (e' * e) / N;
    LL = -(N/2) * (log(2*pi*sigma2) + 1);
end


function [I, p_perm] = moransI_perm(e, W, nperm, seed)
% Moran's I with permutation p-value (two-sided).
% W can be row-standardized or not; uses S0 = sum(W(:)).

    if nargin < 3 || isempty(nperm), nperm = 999; end
    if nargin < 4, seed = 0; end
    rng(seed);

    e = e(:);
    N = length(e);

    e0 = e - mean(e);                 % center
    denom = (e0' * e0);
    if denom <= 0
        I = NaN; p_perm = NaN; return;
    end

    S0 = full(sum(W(:)));
    num = e0' * (W * e0);
    I = (N / S0) * (num / denom);

    % permutation distribution
    Iperm = zeros(nperm,1);
    for b = 1:nperm
        ep = e0(randperm(N));
        Iperm(b) = (N / S0) * ((ep' * (W*ep)) / (ep' * ep));
    end

    % two-sided p-value
    p_perm = (sum(abs(Iperm) >= abs(I)) + 1) / (nperm + 1);
end