
%% First part: Normalized phase distributions for different gradient profiles
global generic
load('generic.mat');

% choose N for which phase distributions are displayed
Ns = [0 1 3 10 30];
legendstr = 'legend(';
for i = 1:numel(Ns)
    if (i < numel(Ns))
        legendstr = [legendstr '''N = ' num2str(Ns(i)) ''', '];
    else
        legendstr = [legendstr '''N = ' num2str(Ns(i)) ''');'];
    end
end

% Figure 1: normalized phase distributions for bipolar gradient profile
phd = squeeze(get_norm_phd(Ns, 'bip'));
phd(:, 1) = phd(:, 1) * 2;
phd(:, size(phd, 2)) = phd(:, size(phd, 2)) * 2;

figure(1)
hold off
for i = 1:numel(Ns)
    plot(generic.phis, phd(i, :));
    if (i == 1)
        hold on
    end
end
ylim([0 0.01]);
title('Normalized phase distributions for bip profile');
xlabel('\vartheta');
ylabel('\rho(\vartheta)');
eval(legendstr);

% Figure 2: normalized phase distributions for flow-compensated gradient profile
phd = squeeze(get_norm_phd(Ns, 'fc1'));
phd(:, 1) = phd(:, 1) * 2;
phd(:, size(phd, 2)) = phd(:, size(phd, 2)) * 2;

figure(2)
hold off
for i = 1:numel(Ns)
    plot(generic.phis, phd(i, :));
    if (i == 1)
        hold on
    end
end
ylim([0 0.01]);
title('Normalized phase distributions for fc+ profile');
xlabel('\vartheta');
ylabel('\rho(\vartheta)');
eval(legendstr);

% Figure 3: normalized phase distributions for sine gradient profile
phd = squeeze(get_norm_phd(Ns, 'sin4'));
phd(:, 1) = phd(:, 1) * 2;
phd(:, size(phd, 2)) = phd(:, size(phd, 2)) * 2;

figure(3)
hold off
for i = 1:numel(Ns)
    plot(generic.phis, phd(i, :));
    if (i == 1)
        hold on
    end
end
ylim([0 0.01]);
title('Normalized phase distributions for sin4 profile');
xlabel('\vartheta');
ylabel('\rho(\vartheta)');
eval(legendstr);

% Figure 4: normalized phase distributions for cosine gradient profile
phd = squeeze(get_norm_phd(Ns, 'cos1_4'));
phd(:, 1) = phd(:, 1) * 2;
phd(:, size(phd, 2)) = phd(:, size(phd, 2)) * 2;

figure(4)
hold off
for i = 1:numel(Ns)
    plot(generic.phis, phd(i, :));
    if (i == 1)
        hold on
    end
end
ylim([0 0.01]);
title('Normalized phase distributions for cos4 profile');
xlabel('\vartheta');
ylabel('\rho(\vartheta)');
eval(legendstr);

% Figure 5: apparent pseudo-diffusion coefficients for bip/fc gradient profiles
figure(5)
Ns = [0 1./fliplr(generic.N(3:numel(generic.N)).') generic.N(2:numel(generic.N)).'];
hold off
ADC_bip = sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'bip')), 2) / 2;
ADC_fc1 = sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'fc0')), 2) / 2;
ADC_fc0 = sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'fc1')), 2) / 2;
loglog(Ns, ADC_bip);
xlim([0.1 100]);
ylim([0.00001 0.2]);
hold on
loglog(Ns, ADC_fc1);
loglog(Ns, ADC_fc0);
title('Apparent pseudo-diffusion coefficient bip/fc');
xlabel('N');
ylabel('AD^{*}C');
legend('bip', 'fc0', 'fc1', 'Location', 'southeast');

% Figure 6: apparent pseudo-diffusion coefficients for cosine gradient profiles
figure(6)
Ns = [0 1./fliplr(generic.N(3:numel(generic.N)).') generic.N(2:numel(generic.N)).'];
hold off
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'cos0_1')), 2) / 2);
xlim([0.1 100]);
ylim([0.00001 0.2]);
hold on
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'cos0_2')), 2) / 2);
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'cos0_4')), 2) / 2);
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'cos0_6')), 2) / 2);
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'cos0_10')), 2) / 2);
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'cos0_16')), 2) / 2);
title('Apparent pseudo-diffusion coefficient cosine');
xlabel('N');
ylabel('AD^{*}C');
legend('cos 1 osc', 'cos 2 osc', 'cos 4 osc', 'cos 6 osc', 'cos 10 osc', 'cos 16 osc', 'Location', 'southeast');

% Figure 7: apparent pseudo-diffusion coefficients for sine gradient profiles
figure(7)
Ns = [0 1./fliplr(generic.N(3:numel(generic.N)).') generic.N(2:numel(generic.N)).'];
hold off
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'sin1')), 2) / 2);
xlim([0.1 100]);
ylim([0.00001 0.2]);
hold on
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'sin2')), 2) / 2);
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'sin4')), 2) / 2);
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'sin6')), 2) / 2);
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'sin10')), 2) / 2);
loglog(Ns, sum((ones(numel(Ns), 1) * (generic.phis .* generic.phis)) .* squeeze(get_norm_phd(Ns.', 'sin16')), 2) / 2);
title('Apparent pseudo-diffusion coefficient sine');
xlabel('N');
ylabel('AD^{*}C');
legend('sin 1 osc', 'sin 2 osc', 'sin 4 osc', 'sin 6 osc', 'sin 10 osc', 'sin 16 osc', 'Location', 'southeast');

%% Second part: Simulated attenuation curves (ISMRM abstract)
global generic
load('generic.mat');

%simulation parameters:

%maximum b-value for attenuation curves:
bmax = 300;    % in s/mm²
bsteps = bmax; % number of data points

% physiologic parameters are given in intervals, e.g. v_cap = 1-2 mm/s. 
vsteps = 25; % vsteps = interpolation steps for physiologic velocities.  

% Duration of diffusion gradient profile
T = 70; % ms
T = T * ones(bmax+1, vsteps);

%physiologic parameters according to following sources:
% - Le Bihan et al. Radiology 1986;161:401-407
% - Le Bihan, Turner Magn Reson Med 1992;27:171-178
% - Maki et al. Magn Reson Med 1991;17:95-107
% - Klinke et al. "Physiologie. 6. Auflage ed. New york: Georg Thieme Verlag; 2010. 191 p.
% - Kunsch K, Kunsch S. Der Mensch in Zahlen - Eine Datensammlung in Tabellen mit über 20.000 Einzelwerten. 2. Auflage ed. Berlin: Spektrum Akademischer Verlag; 2000
% - Pawlik et al. Brain Res 1981;208:35-58

% mean velocity in capillaries
vcap_max = 2;  % mm/s
vcap_min = 1;  % mm/s

% mean velocity in arterioles
vatl_max = 30; % mm/s
vatl_min = 5;  % mm/s

% mean velocity in venules
vvnl_max = 10; % mm/s
vvnl_min = 5;  % mm/s

% length of vessel segments:
lcap = 0.075; % mm
latl = 2;     % mm
lvnl = 2;     % mm

% prepare matrices for calculation of IVIM-signal
lcap = lcap * ones(bmax+1, vsteps);
latl = latl * ones(bmax+1, vsteps);
lvnl = lvnl * ones(bmax+1, vsteps);
vcap = ones(bmax+1, 1) * (vcap_min + vcap_max * (0.5:vsteps) / vsteps);
vatl = ones(bmax+1, 1) * (vatl_min + vatl_max * (0.5:vsteps) / vsteps);
vvnl = ones(bmax+1, 1) * (vvnl_min + vvnl_max * (0.5:vsteps) / vsteps);
b = (0:bmax)' * ones(1, vsteps);

% plot attenuation curves for capillaries (assuming laminar velocity profile)
figure(1)
hold off
plot(0:bmax, sum(get_IVIM_laminar(b, T, lcap, vcap, 'bip'), 2) / vsteps);
hold on
plot(0:bmax, sum(get_IVIM_laminar(b, T, lcap, vcap, 'fc0'), 2) / vsteps, 'c');
plot(0:bmax, sum(get_IVIM_laminar(b, T, lcap, vcap, 'fc1'), 2) / vsteps, 'r');
plot(0:bmax, sum(get_IVIM_laminar(b, T, lcap, vcap, 'sin1'), 2) / vsteps, 'g');
plot(0:bmax, sum(get_IVIM_laminar(b, T, lcap, vcap, 'sin16'), 2) / vsteps, 'g');
plot(0:bmax, sum(get_IVIM_laminar(b, T, lcap, vcap, 'cos0_1'), 2) / vsteps, 'k');
plot(0:bmax, sum(get_IVIM_laminar(b, T, lcap, vcap, 'cos0_16'), 2) / vsteps, 'k');
title('IVIM Signal attenuation in capillaries (v=1-2mm/s, l=75µm)');
xlabel('b (s/mm²)');
ylabel('F(b,T,v,l)');
legend('bip', 'fc0', 'fc1', 'sin1', 'sin16', 'cos1', 'cos16');

% plot attenuation curves for arterioles (assuming laminar velocity profile)
figure(2)
hold off
plot(0:bmax, sum(get_IVIM_laminar(b, T, latl, vatl, 'bip'), 2) / vsteps);
hold on
plot(0:bmax, sum(get_IVIM_laminar(b, T, latl, vatl, 'fc0'), 2) / vsteps, 'c');
plot(0:bmax, sum(get_IVIM_laminar(b, T, latl, vatl, 'fc1'), 2) / vsteps, 'r');
plot(0:bmax, sum(get_IVIM_laminar(b, T, latl, vatl, 'sin1'), 2) / vsteps, 'g');
plot(0:bmax, sum(get_IVIM_laminar(b, T, latl, vatl, 'sin16'), 2) / vsteps, 'g');
plot(0:bmax, sum(get_IVIM_laminar(b, T, latl, vatl, 'cos0_1'), 2) / vsteps, 'k');
plot(0:bmax, sum(get_IVIM_laminar(b, T, latl, vatl, 'cos0_16'), 2) / vsteps, 'k');
title('IVIM Signal attenuation in arterioles (v=5-30mm/s, l=2mm)');
xlabel('b (s/mm²)');
ylabel('F(b,T,v,l)');
legend('bip', 'fc0', 'fc1', 'sin1', 'sin16', 'cos1', 'cos16');

% plot attenuation curves for venules (assuming laminar velocity profile)
figure(3)
hold off
plot(0:bmax, sum(get_IVIM_laminar(b, T, lvnl, vvnl, 'bip'), 2) / vsteps);
hold on
plot(0:bmax, sum(get_IVIM_laminar(b, T, lvnl, vvnl, 'fc0'), 2) / vsteps, 'c');
plot(0:bmax, sum(get_IVIM_laminar(b, T, lvnl, vvnl, 'fc1'), 2) / vsteps, 'r');
plot(0:bmax, sum(get_IVIM_laminar(b, T, lvnl, vvnl, 'sin1'), 2) / vsteps, 'g');
plot(0:bmax, sum(get_IVIM_laminar(b, T, lvnl, vvnl, 'sin16'), 2) / vsteps, 'g');
plot(0:bmax, sum(get_IVIM_laminar(b, T, lvnl, vvnl, 'cos0_1'), 2) / vsteps, 'k');
plot(0:bmax, sum(get_IVIM_laminar(b, T, lvnl, vvnl, 'cos0_16'), 2) / vsteps, 'k');
title('IVIM Signal attenuation in venules (v=5-10mm/s, l=2mm)');
xlabel('b (s/mm²)');
ylabel('F(b,T,v,l)');
legend('bip', 'fc0', 'fc1', 'sin1', 'sin16', 'cos1', 'cos16');

%% Third part: fit to in-vivo ROI data:

global generic
load('generic.mat'); % phase distributions
load('data.mat');    % roi-data (start with liver ...)

% select data (betta = gradient profile, 2 = fc1, 1 = bip)
FC40 = betta == 2 & T == 40;
FC70 = betta == 2 & T == 70;
FC100 = betta == 2 & T == 100;
BP40 = betta == 1 & T == 40;
BP70 = betta == 1 & T == 70;
BP100 = betta == 1 & T == 100;

% plot liver roi data
figure(1)
hold off
errorbar(b(FC40), liver(FC40), s_liver(FC40), 'ro');
hold on
errorbar(b(FC70), liver(FC70), s_liver(FC70), 'g^');
errorbar(b(FC100), liver(FC100), s_liver(FC100), 'ks');
errorbar(b(BP40), liver(BP40), s_liver(BP40), 'rh');
errorbar(b(BP70), liver(BP70), s_liver(BP70), 'gp');
errorbar(b(BP100), liver(BP100), s_liver(BP100), 'kd');
xlabel('b (s/mm²)');
ylabel('S(b) / S(0)');
title('Signal attenuation for liver');
legend('FC40', 'FC70', 'FC100', 'BP40', 'BP70', 'BP100');

% prefit (only high non-flowcompensated b-values) to get an initial 
% estimate for f, D.
biphighb = ((b >= 200) & (betta == 1));
options = optimset('Display', 'off');
[param(1:2), param(3)] = fminsearch(@(x) norm(liver(biphighb) - (1-x(2))*exp(-b(biphighb) * x(1) / 1000))^2, [0.2, 2], options);

% "real fit"
[param(4:7), param(8)] = fminsearch(@(x) norm(liver - generic_model_laminar(x, b, T, betta))^2, [param(1:2), 2, 4]);

% estimation of fit uncertainty
resi = param(8) / (numel(b) - 4);
jacky = jacobianest(@(x) (liver - generic_model_laminar(x, b, T, betta)), param(4:7));
param(9:12) = sqrt(diag(resi * inv(jacky'*jacky))); %#ok<MINV>

% liver fit - 
fprintf('Best fit liver: D=%f±%f, f=%f±%f, l=%f±%f, v=%f±%f\n', param(4), ...
    param(9), param(5), param(10), param(6), param(11), param(7), param(12));

% create attenuation curves based on fit-parameters (unoptimized, takes a
% while).
b = (0:550)';
betta = ones(size(b)) * 2;
T = ones(size(b)) * 40;
FC40_fit = generic_model_laminar(param(4:7), b, T, betta);
T = ones(size(b)) * 70;
FC70_fit = generic_model_laminar(param(4:7), b, T, betta);
T = ones(size(b)) * 100;
FC100_fit = generic_model_laminar(param(4:7), b, T, betta);
betta = ones(size(b));
MP100_fit = generic_model_laminar(param(4:7), b, T, betta);
T = ones(size(b)) * 70;
MP70_fit = generic_model_laminar(param(4:7), b, T, betta);
T = ones(size(b)) * 40;
MP40_fit = generic_model_laminar(param(4:7), b, T, betta);

plot(FC40_fit, 'r');
plot(FC70_fit, 'g');
plot(FC100_fit, 'k');
plot(MP40_fit, 'r');
plot(MP70_fit, 'g');
plot(MP100_fit, 'k');

% same agein for pancreas data ....
load('data.mat')

figure(2)
hold off
errorbar(b(FC40), pancreas(FC40), s_pancreas(FC40), 'ro');
hold on
errorbar(b(FC70), pancreas(FC70), s_pancreas(FC70), 'g^');
errorbar(b(FC100), pancreas(FC100), s_pancreas(FC100), 'ks');
errorbar(b(BP40), pancreas(BP40), s_pancreas(BP40), 'rh');
errorbar(b(BP70), pancreas(BP70), s_pancreas(BP70), 'gp');
errorbar(b(BP100), pancreas(BP100), s_pancreas(BP100), 'kd');
xlabel('b (s/mm²)');
ylabel('S(b) / S(0)');
title('Signal attenuation for pancreas');
legend('FC40', 'FC70', 'FC100', 'BP40', 'BP70', 'BP100');

biphighb = ((b >= 200) & (betta == 1));
options = optimset('Display', 'off');
[param(1:2), param(3)] = fminsearch(@(x) norm(pancreas(biphighb) - (1-x(2))*exp(-b(biphighb) * x(1) / 1000))^2, [0.2, 2], options);
[param(4:7), param(8)] = fminsearch(@(x) norm(pancreas - generic_model_laminar(x, b, T, betta))^2, [param(1:2), 2, 4]);

resi = param(8) / (numel(b) - 4);
jacky = jacobianest(@(x) (pancreas - generic_model_laminar(x, b, T, betta)), param(4:7));
param(9:12) = sqrt(diag(resi * inv(jacky'*jacky))); %#ok<MINV>

fprintf('Best fit pancreas: D=%f±%f, f=%f±%f, l=%f±%f, v=%f±%f\n', param(4), ...
    param(9), param(5), param(10), param(6), param(11), param(7), param(12));

[param(4:7), param(8)] = fminsearch(@(x) norm(pancreas - generic_model_laminar(x, b, T, betta))^2, [param(1:2), 2, 4]);

b = (0:550)';
betta = ones(size(b)) * 2;
T = ones(size(b)) * 40;
FC40_fit = generic_model_laminar(param(4:7), b, T, betta);
T = ones(size(b)) * 70;
FC70_fit = generic_model_laminar(param(4:7), b, T, betta);
T = ones(size(b)) * 100;
FC100_fit = generic_model_laminar(param(4:7), b, T, betta);
betta = ones(size(b));
MP100_fit = generic_model_laminar(param(4:7), b, T, betta);
T = ones(size(b)) * 70;
MP70_fit = generic_model_laminar(param(4:7), b, T, betta);
T = ones(size(b)) * 40;
MP40_fit = generic_model_laminar(param(4:7), b, T, betta);

figure(2)
plot(FC40_fit, 'r');
plot(FC70_fit, 'g');
plot(FC100_fit, 'k');
plot(MP40_fit, 'r');
plot(MP70_fit, 'g');
plot(MP100_fit, 'k');
