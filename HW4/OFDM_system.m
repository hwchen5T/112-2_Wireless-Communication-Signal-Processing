%===================================
% This is the simulator for OFDM.
%===================================
%
%  The function contains a completed SCFDMA transceiver, some parameters need to
%  be defined before running. These parameters are described for instance as the follows
%
%   SP.FFTsize = 128;           % size of the transmitter IFFT and the receiver FFT
%   SP.mod_size = 4;            % specified QPSK transmission
%   SP.inputBlockSize = 128;    % input data block size, i.e. number of sub-carrier within chunk
%   SP.CPsize = 16;             % CP length
%   SP.SNR = [0:5:30];          % simulated SNR range in dB
%   SP.numRun = 10;             % number of simulation iterations
%   SP.channel = pedAchannel;   % defined the channel response
%
%   The final SER will be returned.

function SER = OFDM_system(SP)

SER = zeros(length(SP.SNR),1);
numSymbols = SP.FFTsize;                % number of transmitted OFDM symbols
numCP = SP.CPsize;

% initial for M-QPSK symbol generation
BitPerSyml = log2(SP.mod_size);         % calculate bit number per symbol
SgnlPwr = sum(abs(qammod(0: SP.mod_size-1, SP.mod_size)).^2)/SP.mod_size;
WegVect = pow2([log2(SP.mod_size)-1:-1:0]);

% construct channel
% IID Rayleigh faing channel
tmp_fad = randn(2, length(SP.channel));
raychannel = (tmp_fad(1, :) + 1i*tmp_fad(2, :))/sqrt(2);

% normalize the channel
raychannel = raychannel/sqrt(sum(raychannel.^2));

SP.channel = raychannel;
H_channel = fft(SP.channel, SP.FFTsize);    % frequency domain version of the channel response


for n = 1:length(SP.SNR)
    
    CurrentSNR = SP.SNR(n);                 % show cuurent status (SNR)
    NPW = SgnlPwr*(SP.inputBlockSize/(SP.inputBlockSize+SP.CPsize))*10^(-CurrentSNR/10);             % show noise power status (SNR)
    
    % initialize the error count
    ErrCnt = 0;                             % symbol error count
    
    for k = 1:SP.numRun
        
        correct = 0;
        raw_bits = round(rand(BitPerSyml,SP.inputBlockSize));       % generate N QPSK symbol, each symbol contains n bits
        inputSymbols = qammod(WegVect * raw_bits, SP.mod_size);     % constellation mapping
        inputSamples_ofdm = ifft(inputSymbols, SP.FFTsize)*sqrt(numSymbols);         % convert the signal back to time domain
        
        % add CP
        theta = [zeros(numCP,(numSymbols-numCP)),eye(numCP);eye(numSymbols)];
        Tx_u = inputSamples_ofdm * (theta.');
        
        % propagate through multi-path channel. i.e take the convolution
        Rx_v = conv(Tx_u , SP.channel);
        Rx_v = Rx_v(1:numSymbols+numCP);
      
        
        % generate AWGN with appropriate noise power and add it into RxSmples
        tmp_noise = randn(1,numSymbols+numCP) + 1i*randn(1,numSymbols+numCP);
        noise = sqrt(NPW)*tmp_noise/sqrt(2);
        Rx_v = Rx_v + noise;
        
        % remove CP
        gamma = [zeros(numSymbols,numCP),eye(numSymbols)];
        NoCP = Rx_v * (gamma.');
        
        % add CFO
        CFO = exp(1i*2*pi*SP.CFO*(0:numSymbols-1)/numSymbols);
        output = NoCP.*CFO;
        
        % convert the received signal into frequency domain
        output = fft(output,numSymbols)/sqrt(numSymbols);
        
        % find the channel response for the subcarriers
        H = diag(H_channel);
        
        
        % perform zero forcing equalization in the frequency domain
        E_ZF = inv(H);
        X_zf_ofdm = output*(E_ZF.');
        
        
        % perform hard decision detection
        Det_Symbols = HardDet(X_zf_ofdm, SP.mod_size);
        
        
        % find and count errors
        ErrCnt = ErrCnt + sum(inputSymbols ~= Det_Symbols);
        
    end
    
    % calculate the symbol error rate (SER)
    SER(n) = ErrCnt/numSymbols/SP.numRun;
    
end


end
