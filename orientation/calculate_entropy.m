%
%   Calculate differential entropy
%

function [H] = calculate_entropy(phi_file, ent_ws)
    if ~exist('ent_ws', 'var') % entropy window size
        ent_ws = 50;
    end

    load(phi_file);
    
    H = zeros(size(phi)-ent_ws+1);

    [row_count, col_count] = size(H);

    % sliding window to select for entropy
    for r = 1:row_count
        for c = 1:col_count
            theta = phi(r:r+ent_ws-1, c:c+ent_ws-1);

            % differential entropy
            H(r, c) = differential_entropy(theta(:));
        end
        if mod(r, 1000) == 0
            disp(r);
        end
    end
end