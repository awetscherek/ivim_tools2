function S = generic_model_laminar(x, b, T, betta)
% calculates IVIM-signal for the tissue parameters specified via x (see
% below). b, T and betta all need to have the same size. This function
% currently only supports flow-compensated and bipolar profiles, but can
% easily be generalized by just modifying the calls to get_IVIM_laminar...

D = abs(x(1)); % 10^-3 mm²/s
f = abs(x(2));
l = abs(x(3)) * ones(size(b)); % mm
v = abs(x(4)) * ones(size(b)); % mm/s

F = zeros(size(b));

% signal attenuation of blood compartment for bipolar gradients:
if (numel(betta(betta == 1)) > 0) 
    F(betta == 1) = get_IVIM_laminar(b(betta == 1)', T(betta == 1)', l(betta == 1)', v(betta == 1)', 'bip');
end

% signal attenuation of blood compartment for flow-compensated gradients
if (numel(betta(betta == 2)) > 0)
    F(betta == 2) = get_IVIM_laminar(b(betta == 2)', T(betta == 2)', l(betta == 2)', v(betta == 2)', 'fc1');
end

% IVIM signal assuming a diffusion coefficient of 1.6 * 10^-3 mm²/s for
% blood.
S = ((1 - f) * exp(-b * D / 1000) + f * exp(-b * 1.6 / 1000) .* F);
        
end
