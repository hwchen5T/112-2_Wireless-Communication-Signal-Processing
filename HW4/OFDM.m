clear all;
close all;

SP.FFTsize = 128;           % size of the transmitter IFFT and the receiver FFT
SP.mod_size = 4;            % adopt QPSK modulation
SP.inputBlockSize = 128;    % input data block size
SP.CPsize = 16;             % CP length
SP.SNR = [0:5:30];          % simulated SNR range in dB
SP.numRun = 10;             % number of simulation iterations
SP.CFO = 0;                 % CFO is zero
SP.channel = ones(1,10);    % initial channel
Packet_numRun = 1000;       % number of packets

SER = 0;

% Run the simulation for OFDM with different channel
for n = 1:1:Packet_numRun
    SER = OFDM_system(SP) + SER;
end
SER = SER / Packet_numRun;

semilogy(SP.SNR, SER);
xlabel('SNR (dB)');
ylabel('SER');
title('QPSK OFDM with ZF equalizer');
