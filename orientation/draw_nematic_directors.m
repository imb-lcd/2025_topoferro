%
% draw nematics directors
%

function draw_nematic_directors(phi, coh, imsize, nem_ws, nem_width, nem_color, show_coherency, savefigure, outfile)
    scale_factor = 1;

    phi_nematic = -phi; % phi;
           
    % rescale phi from (-0.5, 0.5) to (-pi/2, pi/2) for correct nematics
    % this is due to previous code that normalized the data by pi.
    phi_nematic = rescale(phi_nematic, -pi/2, pi/2);

    % subselect the angles to draw nematics 
    phi_nematic = imresize(phi_nematic, [length(phi)/nem_ws length(phi)/nem_ws], 'nearest'); % for 100*100 nematics 

    if show_coherency
        coh_nematic = imresize(coh, [length(coh)/nem_ws length(coh)/nem_ws], 'nearest'); % for 100*100 nematics 
    end

    % offset X and Y for nematic position
    sz = length(phi_nematic);
    [X, Y] = meshgrid(1:pi/2:pi/2*sz, 1:pi/2:pi/2*sz);  % for 100*100 nematics 

    X = X - scale_factor/2*cos(phi_nematic);
    Y = Y - scale_factor/2*sin(phi_nematic);

    % scale quiver accordingly
    U = cos(phi_nematic)*scale_factor;
    V = sin(phi_nematic)*scale_factor;
    
    % set the length based on coherency
    if show_coherency
        U = U .* coh_nematic;
        V = V .* coh_nematic;
    end

    f = figure;
    imshow(zeros(imsize)+255, 'Border', 'tight')
    
    rX = rescale(X, 25, imsize(2)-25);
    rY = rescale(Y, 25, imsize(1)-25);
    hold on
    quiver(rX, rY, U, V, 0.5, 'color', nem_color, 'LineWidth', nem_width, 'ShowArrowHead','off'); % 'Autoscale', 'on'
    
    set(gcf, 'Color', 'w'); % set the background to white
    
    hold off

    if savefigure
        set(gcf,'Units','pixels','Position',[0 0 imsize]);
        saveas(gcf, outfile, 'svg');
        close all;

    end
end
