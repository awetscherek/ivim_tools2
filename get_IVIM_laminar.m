function F = get_IVIM_laminar(b, T, l, v, profile)

    steps = 25;
    shape = size(v);
    
    F = zeros(size(v));
    
    for k = 1:size(v, 2)    
        vs = v(:,k) * 2 * (0.5:steps) / steps;
        ls = l(:,k) * ones(1, steps);
        Ts = T(:,k) * ones(1, steps);
        bs = b(:,k) * ones(1, steps);
        F(:, k) = reshape(sum(get_IVIM_signal(bs, Ts, ls, vs, profile), 2), size(v, 1), 1) / steps;
    end 

end