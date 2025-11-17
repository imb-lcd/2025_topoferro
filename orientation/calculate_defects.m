%
%   Calculate defects
%

function [mhalf, phalf] = calcualte_defects(phi_file, defthres)
    if ~exist('defthres', 'var') % entropy window size
        defthres = 0.3;
    end
    
    ws_px = 2;

    load(phi_file);
    
    % Reconstruct tensor with amplitude set to 1
    GGxx = sin(phi).^2 - 0.5;
    GGyy = cos(phi).^2 - 0.5;
    GGxy = sin(phi).*cos(phi);
    
    % Calculate charge deisty to identify the position of topological defect
    divXXX = GGxx(ws_px+1:end, ws_px:end-1) - GGxx(1:end-ws_px, ws_px:end-1);
    divXXY = GGxx(ws_px:end-1, ws_px+1:end) - GGxx(ws_px:end-1, 1:end-ws_px);
    divXYX = GGxy(ws_px+1:end, ws_px:end-1) - GGxy(1:end-ws_px, ws_px:end-1);
    divXYY = GGxy(ws_px:end-1, ws_px+1:end) - GGxy(ws_px:end-1, 1:end-ws_px);
    chargedensity = divXXX.*divXYY - divXYX.*divXXY;
    
    % find position of defect (mimic peak_local_max from python)
    % **centroid method may detect >1 coordinates for a single local maximum
    mhalf_max = chargedensity >= defthres;
    m_obj = regionprops('table', mhalf_max, 'Centroid');
    
    phalf_max = -chargedensity >= defthres;
    p_obj = regionprops('table', phalf_max, 'Centroid');
    
    mhalf = [m_obj.Centroid(:,1), m_obj.Centroid(:,2)];
    phalf = [p_obj.Centroid(:,1), p_obj.Centroid(:,2)];
end