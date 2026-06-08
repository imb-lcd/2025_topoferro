function summary = summarize_curve_image(curves, A)
% curves: cell array, each cell = [Nx2] matrix of (x,y)
% A: image area (optional, can be [])
%
% Output:
%   summary.image  -> one-row table with image-level summaries
%   summary.curves -> one-row-per-curve table with individual summaries

    n = numel(curves);

    curve_id = (1:n).';
    n_points = zeros(n,1);

    % Geometry
    L = NaN(n,1);
    D = NaN(n,1);
    mean_abs_curv = NaN(n,1);
    rms_curvature = NaN(n,1);
    total_abs_curvature = NaN(n,1);
    max_abs_curvature = NaN(n,1);
    straightness = NaN(n,1);

    % Orientation
    mean_nematic_orientation = NaN(n,1);   % radians in [0, pi)
    orientation_order = NaN(n,1);          % in [0, 1]
    orientation_dispersion = NaN(n,1);     % = 1 - orientation_order

    for i = 1:n
        xy = curves{i};
        x = xy(:,1);
        y = xy(:,2);

        n_points(i) = numel(x);

        if numel(x) < 3 || numel(y) < 3
            continue
        end

        % Smooth coordinates
        xs = smoothdata(x(:), 'gaussian', 5);
        ys = smoothdata(y(:), 'gaussian', 5);

        % Derivatives
        dx  = gradient(xs);
        dy  = gradient(ys);
        ddx = gradient(dx);
        ddy = gradient(dy);

        speed = hypot(dx, dy);
        valid_dir = isfinite(speed) & (speed > eps);

        % Arc length element
        ds = hypot(diff(xs), diff(ys));
        if isempty(ds)
            continue
        end
        ds = [ds; ds(end)];

        % Geometry
        denom = (dx.^2 + dy.^2).^(3/2);
        kappa = (dx .* ddy - dy .* ddx) ./ denom;
        kappa(denom < eps) = NaN;

        L(i) = sum(ds, 'omitnan');
        D(i) = hypot(xs(end) - xs(1), ys(end) - ys(1));

        total_abs_curvature(i) = sum(abs(kappa) .* ds, 'omitnan');
        mean_abs_curv(i) = total_abs_curvature(i) / L(i);
        rms_curvature(i) = sqrt(sum((kappa.^2) .* ds, 'omitnan') / L(i));
        max_abs_curvature(i) = max(abs(kappa), [], 'omitnan');
        straightness(i) = D(i) / L(i);

        % Local nematic orientation along curve
        % theta in [0, pi)
        theta = 0.5 * atan2(2*dx.*dy, dx.^2 - dy.^2);
        theta = mod(theta, pi);
        theta(~valid_dir) = NaN;

        % Nematic mean and order, arc-length weighted
        w = ds;
        valid_theta = isfinite(theta) & isfinite(w) & (w > 0);

        if any(valid_theta)
            c = sum(w(valid_theta) .* cos(2 * theta(valid_theta)));
            s = sum(w(valid_theta) .* sin(2 * theta(valid_theta)));
            wsum = sum(w(valid_theta));

            mean_nematic_orientation(i) = 0.5 * atan2(s, c);
            mean_nematic_orientation(i) = mod(mean_nematic_orientation(i), pi);

            orientation_order(i) = sqrt(c^2 + s^2) / wsum;
            orientation_dispersion(i) = 1 - orientation_order(i);
        end
    end

    % Per-curve table
    curveTbl = table( ...
        curve_id, ...
        n_points, ...
        L, ...
        D, ...
        mean_abs_curv, ...
        rms_curvature, ...
        total_abs_curvature, ...
        max_abs_curvature, ...
        straightness, ...
        mean_nematic_orientation, ...
        orientation_order, ...
        orientation_dispersion, ...
        'VariableNames', { ...
            'curve_id', ...
            'n_points', ...
            'length', ...
            'end_to_end', ...
            'mean_abs_curvature', ...
            'rms_curvature', ...
            'total_abs_curvature', ...
            'max_abs_curvature', ...
            'straightness', ...
            'mean_nematic_orientation', ...
            'orientation_order', ...
            'orientation_dispersion'});

    % Valid curves for image-level summaries
    valid_geom = isfinite(L) & isfinite(mean_abs_curv) & isfinite(straightness) & (L > 0);
    Lv = L(valid_geom);
    mv = mean_abs_curv(valid_geom);
    sv = straightness(valid_geom);

    imageSummary = struct();

    % Geometry summaries
    imageSummary.n_curves        = numel(Lv);
    imageSummary.total_length    = sum(Lv);
    imageSummary.mean_curv_w     = sum(Lv .* mv) / sum(Lv);
    imageSummary.median_curv     = median(mv);
    imageSummary.iqr_curv        = iqr(mv);
    imageSummary.mean_straight_w = sum(Lv .* sv) / sum(Lv);
    imageSummary.median_straight = median(sv);
    imageSummary.iqr_straight    = iqr(sv);

    % Image-level nematic summary across curves
    valid_orient = isfinite(mean_nematic_orientation) & isfinite(orientation_order) & isfinite(L) & (L > 0);

    if any(valid_orient)
        % Weight by both curve length and within-curve orientational coherence
        w_img = L(valid_orient) .* orientation_order(valid_orient);
        th_img = mean_nematic_orientation(valid_orient);

        C = sum(w_img .* cos(2 * th_img));
        S = sum(w_img .* sin(2 * th_img));
        W = sum(w_img);

        imageSummary.global_nematic_orientation = 0.5 * atan2(S, C);
        imageSummary.global_nematic_orientation = mod(imageSummary.global_nematic_orientation, pi);
        imageSummary.global_orientation_order = sqrt(C^2 + S^2) / W;

        imageSummary.mean_curve_orientation_order = mean(orientation_order(valid_orient));
        imageSummary.median_curve_orientation_order = median(orientation_order(valid_orient));
        imageSummary.iqr_curve_orientation_order = iqr(orientation_order(valid_orient));
    else
        imageSummary.global_nematic_orientation = NaN;
        imageSummary.global_orientation_order = NaN;
        imageSummary.mean_curve_orientation_order = NaN;
        imageSummary.median_curve_orientation_order = NaN;
        imageSummary.iqr_curve_orientation_order = NaN;
    end

    if nargin > 1 && ~isempty(A)
        imageSummary.length_density = imageSummary.total_length / A;
        imageSummary.count_density  = imageSummary.n_curves / A;
    end

    summary = struct();
    summary.image = struct2table(imageSummary);
    summary.curves = curveTbl;
end