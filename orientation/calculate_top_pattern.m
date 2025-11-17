%
%   Calculate phi_diff, splay, and bend
%

function [b, s, pdiff] = calculate_top_pattern(phi_file, pttn_cutoff)
    if ~exist('pttn_cutoff', 'var') % pattern cutoff percentage
        pttn_cutoff = 0.5;
    end    

    load(phi_file);

    n_vec_x = sin(phi);
    n_vec_y = cos(phi);
    
    phi_diff_m_x = phi(ws_px+1:end, ws_px:end-1) - phi(1:end-2, ws_px:end-1);
    phi_diff_m_y = phi(ws_px:end-1, ws_px+1:end) - phi(ws_px:end-1, 1:end-2);
    
    div_n = diff_n_vec(n_vec_x(ws_px+1:end, ws_px:end-1), n_vec_x(1:end-2, ws_px:end-1), phi_diff_m_x) + ...
        diff_n_vec(n_vec_y(ws_px:end-1, ws_px+1:end), n_vec_y(ws_px:end-1, 1:end-2), phi_diff_m_y);
    
    rot_n_z = diff_n_vec(n_vec_x(ws_px:end-1, 1:end-2), n_vec_x(ws_px:end-1, ws_px+1:end), phi_diff_m_y) + ...
        diff_n_vec(n_vec_y(ws_px+1:end, ws_px:end-1), n_vec_y(1:end-2, ws_px:end-1), phi_diff_m_x);
    
    n_rot_n_x = n_vec_y(ws_px:end-1, ws_px:end-1) .* rot_n_z;
    n_rot_n_y = -n_vec_x(ws_px:end-1, ws_px:end-1) .* rot_n_z;
    
    % calculate phi difference
    pdiff = diff_n_vec(phi(1:end-2, ws_px:end-1), phi(ws_px+1:end, ws_px:end-1), phi_diff_m_x).^2 + ...
        diff_n_vec(phi(ws_px:end-1, 1:end-2), phi(ws_px:end-1, ws_px+1:end), phi_diff_m_y).^2;
    
    pdiff = normalize_nematic_ordering(pdiff, pttn_cutoff, 100-pttn_cutoff);
    
    % calculate splay
    s = div_n.^2;
    s = normalize_nematic_ordering(s, pttn_cutoff, 100-pttn_cutoff);
    
    % calculate bend
    b = n_rot_n_x.^2 + n_rot_n_y.^2;
    b = normalize_nematic_ordering(b, pttn_cutoff, 100-pttn_cutoff);
end

%% FUNCTIONS
function input = normalize_nematic_ordering(input, btm_cut, top_cut)
    if nargin == 3
        btm = prctile(input(:), btm_cut);
        top = prctile(input(:), top_cut);
    else
        btm = prctile(input(:), 1);
        top = prctile(input(:), 99);
    end
    
    input = 2*(input-median(input))./(top-btm);
    input(input>1) = 1;
    input(input<-1) = -1;
end

function diff = diff_n_vec(n1, n2, phi_diff)
    plus_reg = phi_diff > 0.45;
    min_reg = phi_diff < -0.45;
    diff = n1 - n2;

    diff(plus_reg) = n1(plus_reg) + n2(plus_reg);
    diff(min_reg) = -n1(min_reg) - n2(min_reg);
end