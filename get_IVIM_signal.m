function F = get_IVIM_signal(b, T, l, v, profile)

    % data must be ordered such that the ratio T*v/l is constant along the 
    % first dimension, as it often is the case, e.g. for calculating a
    % signal attenuation curve for a given diffusion gradient profile.

    % b in s/mm²
    % v in mm/s
    % T in ms
    % l in mm

    global generic
    
    shape = size(v);
    phds = reshape(get_norm_phd(T(1,:) .* v(1,:) ./ l(1,:) / 1000, profile), size(v,2), numel(generic.('phis')));
    phds = reshape(ones(size(v,1), 1) * phds(:)', [numel(v) numel(generic.('phis'))]);
    F = reshape(sum(phds .* cos(v(:) .* sqrt(b(:) .* T(:) / 1000) * generic.('phis')), 2), shape);

end