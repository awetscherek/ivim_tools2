function F = get_IVIM_signal(b, T, l, v, profile)

    % data must be ordered such that the ratio T*v/l is constant along the 
    % first dimension, as it often is the case, e.g. for calculating a
    % signal attenuation curve for a given diffusion gradient profile.

    % b in s/mm�
    % v in mm/s
    % T in ms
    % l in mm

    % phase distribution data is assumed to be available via global variable
    % generic
    global generic
    
    % output dimensions:
    shape = size(v);
    
    % fetching phase distributions for different N = T*v/l
    phds = reshape(get_norm_phd(T(1,:) .* v(1,:) ./ l(1,:) / 1000, profile), size(v,2), numel(generic.('phis')));
    phds = reshape(ones(size(v,1), 1) * phds(:)', [numel(v) numel(generic.('phis'))]);
    
    % calculate signal attenuation
    F = reshape(sum(phds .* cos(v(:) .* sqrt(b(:) .* T(:) / 1000) * generic.('phis')), 2), shape);

end