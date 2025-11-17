%
%   calculate orientation based on phase-contrast image
%

function [coh, phi] = calculate_tensor(dic_file, ws_px, sigma)
    if ~exist('ws_px', 'var')
        ws_px = 2;
    end
    if ~exist('sigma', 'var')
        sigma = 45;
    end

    I = imread(dic_file);

    % normalize image
    I = double(I);
    I = (I-min(I(:))) / (max(I(:))-min(I(:)));
    
    % calcute the derivatives
    Ix = I(ws_px+1:end, ws_px:end-1) - I(1:end-ws_px, ws_px:end-1);
    Iy = I(ws_px:end-1, ws_px+1:end) - I(ws_px:end-1, 1:end-ws_px);
    
    % construct the tensor and blur with Gaussian filter
    Gxx = imgaussfilt(Ix.*Ix, sigma);   %   Gxx = imgaussfilt(-Ix.*-Ix, sigma);
    Gyy = imgaussfilt(Iy.*Iy, sigma);   %   Gyy = imgaussfilt(Iy.*Iy, sigma);
    Gxy = imgaussfilt(Ix.*Iy, sigma);   %   Gxy = imgaussfilt(-Ix.*Iy, sigma);
    
    % calculate coherency (the 'amplitude' of nematic order) and phase (angle)
    coh = ((Gxx-Gyy).*(Gxx-Gyy)+4.*Gxy.*Gxy) ./ ((Gxx+Gyy).*(Gxx+Gyy));
    phi = atan2(2.*Gxy, Gxx-Gyy) ./ 2;
end