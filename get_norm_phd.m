function phd = get_norm_phd(N, profile)

% fetches normalized phase distribution for the specified N and profile 
% out of global variable generic. Linear interpolation between the values
% of N for which phase distributions are available are performed.
global generic

xmax = max(generic.('phis'));
        
phd = zeros([size(N), numel(generic.('phis'))]); 

for i = 1:numel(N)

    iceil = find(generic.('N') > N(i), 1);

    % if N > max(generic.('N')) use approximation by Gaussian
    if (numel(iceil) == 0)
        phd(i + ((1:numel(generic.('phis'))) - 1) * numel(N)) ...
            = erf(min(xmax, generic.('phis') + generic.('phis')(2) / 2) * sqrt(3/2 * N(i))) ...
            - erf(max(   0, generic.('phis') - generic.('phis')(2) / 2) * sqrt(3/2 * N(i)));
    else
        % use linear interpolation (exact for 0 <= N <= 1)
        ifloor = find(generic.('N') <= N(i), 1, 'last');
        c = (N(i) - generic.('N')(ifloor)) / (generic.('N')(iceil) - generic.('N')(ifloor));
            phd(i + ((1:numel(generic.('phis'))) - 1) * numel(N)) ...
            = (1-c) * generic.(profile)(ifloor, :) + c * generic.(profile)(iceil, :);
    end

end

end