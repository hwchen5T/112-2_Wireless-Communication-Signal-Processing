clear all; clc;

N = 16; % number of antennas
theta0 = 0;
% ============ load "signal.mat" file =============
load('signal.mat') % We have sp belong to 1*10000 & x belong to 16*10000
% ====================== end ======================

% ======= compute steering vector of theta0 =======
a_t0 = steervec(N,theta0); % It is a 16*1 matrix
% ====================== end ======================

Pn = db2pow(-10);
Ns = 10000;

% compute the MVDR beamformer
wmvdr = mvdr(x,sp,a_t0);
% compute the MSINR beamformer
wmsinr = msinr(x,sp,a_t0);

% compute the MMSE beamformer
wmmse = mmse(x,sp,a_t0);

% plot
Thetas = -90 : 0.1 : 89.9;
a = steervec(N, Thetas);
bpmvdr = wmvdr' * a;
bpmsinr = wmsinr' * a;
bpmmse = wmmse' * a;

%% MVDR
figure;
semilogy(Thetas, abs(bpmvdr).^2, 'linewidth', 1.5)
axis([-90, 90, 1e-8, 1e2])
title('MVDR', 'fontsize', 14) 
xlabel('angle', 'fontsize', 11)
ylabel('P(\theta)', 'fontsize', 11)
%% MSINR
figure;
semilogy(Thetas, abs(bpmsinr).^2, 'linewidth', 1.5)
axis([-90, 90, 1e-8, 1e2])
title('MSINR', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('P(\theta)', 'fontsize', 11)
%% MMSE
figure;
semilogy(Thetas, abs(bpmmse).^2, 'linewidth', 1.5)
axis([-90, 90, 1e-8, 1e2])
title('MMSE', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('P(\theta)', 'fontsize', 11)
%% MVDR
figure;
plot(Thetas, angle(bpmvdr))
title('MVDR', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('\angle P(\theta)', 'fontsize', 11)
figure;
%% MSINR
plot(Thetas, angle(bpmsinr))
title('MSINR', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('\angle P(\theta)', 'fontsize', 11)
figure;
%% MMSE
plot(Thetas, angle(bpmmse))
title('MMSE', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('\angle P(\theta)', 'fontsize', 11)