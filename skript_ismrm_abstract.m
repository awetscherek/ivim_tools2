

%% First part: Simulated attenuation curves (ISMRM abstract)
global generic
load('generic.mat');

bmax = 300;    % s/mm²
vsteps = 25;
bsteps = bmax;
T = 70 * ones(bmax+1, vsteps);    % ms

vcap_max = 2;  % mm/s
vcap_min = 1;  % mm/s

vatl_max = 30; % mm/s
vatl_min = 5;  % mm/s

vvnl_max = 10; % mm/s
vvnl_min = 5;  % mm/s

lcap = 0.075 * ones(bmax+1, vsteps); % mm
latl = 2 * ones(bmax+1, vsteps);     % mm
lvnl = 2 * ones(bmax+1, vsteps);     % mm

vcap = ones(bmax+1, 1) * (vcap_min + vcap_max * (0.5:vsteps) / vsteps);
vatl = ones(bmax+1, 1) * (vatl_min + vatl_max * (0.5:vsteps) / vsteps);
vvnl = ones(bmax+1, 1) * (vvnl_min + vvnl_max * (0.5:vsteps) / vsteps);

b = (0:bmax)' * ones(1, vsteps);

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

%% Second part: fit to in-vivo ROI data:

global generic
load('generic.mat');
load('data.mat');

FC40 = betta == 2 & T == 40;
FC70 = betta == 2 & T == 70;
FC100 = betta == 2 & T == 100;
BP40 = betta == 1 & T == 40;
BP70 = betta == 1 & T == 70;
BP100 = betta == 1 & T == 100;

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

biphighb = ((b >= 200) & (betta == 1));
options = optimset('Display', 'off');
[param(1:2), param(3)] = fminsearch(@(x) norm(liver(biphighb) - (1-x(2))*exp(-b(biphighb) * x(1) / 1000))^2, [0.2, 2], options);
[param(4:7), param(8)] = fminsearch(@(x) norm(liver - generic_model_laminar(x, b, T, betta))^2, [param(1:2), 2, 4]);

resi = param(8) / (numel(b) - 4);
jacky = jacobianest(@(x) (liver - generic_model_laminar(x, b, T, betta)), param(4:7));
param(9:12) = sqrt(diag(resi * inv(jacky'*jacky))); %#ok<MINV>

fprintf('Best fit liver: D=%f±%f, f=%f±%f, l=%f±%f, v=%f±%f\n', param(4), ...
    param(9), param(5), param(10), param(6), param(11), param(7), param(12));

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
