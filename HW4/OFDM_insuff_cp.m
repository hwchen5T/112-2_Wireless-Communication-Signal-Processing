clear all;
close all;

SP.FFTsize = 128;           % size of the transmitter IFFT and the receiver FFT
SP.mod_size = 4;            % adopt QPSK modulation
SP.inputBlockSize = 128;    % input data block size
SP.SNR = [0:5:30];          % simulated SNR range in dB
SP.numRun = 10;             % number of simulation iterations
SP.CFO = 0;                 % CFO is zero
SP.channel = ones(1,10);    % initial channel
Packet_numRun = 1000;       % number of packets

% set CP size
CPsize = [2 4 8 16];
SER = zeros(length(SP.SNR), length(CPsize));

% Run the simulation for OFDM with different channel
for k = 1:length(CPsize)
    SP.CPsize = CPsize(k);
    for n = 1:Packet_numRun
        SER(:, k) = OFDM_system(SP) + SER(:, k);
    end
    SER(:, k) = SER(:, k) / Packet_numRun;
end

semilogy(SP.SNR, SER(:, 1), '-+', SP.SNR, SER(:, 2), '-x', SP.SNR, SER(:, 3), '-*', SP.SNR, SER(:, 4));
legend('CP = 2','CP = 4','CP = 8','CP = 16');
xlabel('SNR(dB)');
ylabel('SER');
title('QPSK OFDM with ZF equalizer');
