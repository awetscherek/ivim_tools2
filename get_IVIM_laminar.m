function F = get_IVIM_laminar(b, T, l, v, profile)

    % data must be ordered such that the ratio T*v/l is constant along the 
    % first dimension, as it often is the case, e.g. for calculating a
    % signal attenuation curve for a given diffusion gradient profile.

    % b in s/mm²
    % v in mm/s
    % T in ms
    % l in mm

    % number of steps for parabolic velocity profile (= velocity
    % homogeneously distributed between 0 and 2v).
    steps = 25;
    shape = size(v);
    
    F = zeros(shape);
    
    % calling get_IVIM_signal for each different N = T*v/l, summing over
    % the velocity distribution:
    for k = 1:size(v, 2)
        vs = v(:,k) * 2 * (0.5:steps) / steps;
        ls = l(:,k) * ones(1, steps);
        Ts = T(:,k) * ones(1, steps);
        bs = b(:,k) * ones(1, steps);
        F(:, k) = reshape(sum(get_IVIM_signal(bs, Ts, ls, vs, profile), 2), size(v, 1), 1) / steps;
    end 

end