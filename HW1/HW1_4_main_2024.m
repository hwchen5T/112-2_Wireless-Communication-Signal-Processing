clear;
clc;
close all;

%% Parameters of Jakes model

f_max = (50/3)*2*(10^9)/(3*10^8);
BW = 1000;
Ts = 1/BW;
M = 60;
N_FFT = 1024;

%% Jakes model

ch_Jakes = HW1_4_Jakes_2024(M, f_max, Ts);

% Draw the power spectrum density for Jakes model

%/************
%     Code 
ch_Jakes_fft = fft(ch_Jakes , N_FFT);
ch_Jakes_ffts = fftshift(ch_Jakes_fft);
ch_Jakes_psd = abs(ch_Jakes_ffts).^2 / N_FFT;
%f_idx = 2e9 - f_max : 2 * f_max/(N_FFT-1) :2e9 + f_max ;
f_idx=-BW/2:BW/(length(ch_Jakes_psd)-1):BW/2;
%*************/

H1 = figure(1);
plot(f_idx, ch_Jakes_psd);
%xticks([2e9-500,2e9,2e9+500]);
ylabel('PSD');
xlabel('freq. (Hz)');
grid;

% Evaluate the RMS frequency Doppler spread

%/************
%     Code 
f_avg = (f_idx * ch_Jakes_psd.')/sum(ch_Jakes_psd) ;

f_RMS = sqrt( ( ( (f_idx-f_avg) .^2 )*ch_Jakes_psd.' ) / sum(ch_Jakes_psd) );
%*************/


% Evaluate the PDF of the magitude fading coefficient for Jakes model

%/************
%     Code 
ch_mag = abs(ch_Jakes);
x_m = 0:0.01:7;
[pdf_env , x] = ksdensity(ch_mag,x_m);
pdf_theory = pdf('Rayleigh',x,sqrt(1/2)); 

%*************/

H2 = figure(2);
plot(x,pdf_env,'-',x,pdf_theory,'*');
legend('Simulated','Theoretic');
xlabel('mag.');
ylabel('PDF of mag.');
grid;

% Evaluate the PDF of the phase fading coefficient for Jakes model

%/************
%     Code 
ch_phase = angle(ch_Jakes);
x_p = -pi:0.05:pi;
[pdf_phase , x] = ksdensity(ch_phase,x_p);
%*************/

H3 = figure(3);
plot(x,pdf_phase,'-');
legend('Simulated');
xlabel('phase');
ylabel('PDF of phase');
xlim([-pi,pi]);
ylim([0,1]);
grid;
