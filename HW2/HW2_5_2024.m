clear all;
close all;

%% --- Parameter Setting ---

Trial = 500;  %% trials per SNR
Ns = 300;    %% samples per trial
Q = 1;
P = 2;
L = 7;
K = 9;
Kf = P*K;
Kb = 5;
Kb_2 = 7;
J = 3;
d = 8;

ref2_dB = [-20 -10 0 10 20 30 40 50 60 70];

%--------------------
SNRo_ZF = zeros(1,length(ref2_dB));
SNRo_MMSE = zeros(1,length(ref2_dB));
SNRo_DFE_5 = zeros(1,length(ref2_dB));
SNRo_DFE_7 = zeros(1,length(ref2_dB));

for n = 1:length(ref2_dB)
    
    SNR_dB = ref2_dB(n);    %% signal power is unit
    % Construct noise power with signal power is unit
    NPW = 10^(-SNR_dB/10);
    
    for t = 1:Trial
        
        %--- Fading Gain Generation ---
        fad(1,1:J) = randn(1,J)+1i*randn(1,J); % fading gain
        fad = (1/sqrt(2))*fad;
        %--- Raise Cosine filter ---
        g = rcosflt([0 0 1 0 0 0],1,10,'iir',0.25).';
        g = g(20:120);
        % oversample period P = 2, we sample raise-cosine function with 5 time slot.
        g = g(1:5:101);
        %--- Received signal generation ---
        H_up = [g*fad(1,1) 0 0 0 0]+[0 0 g*fad(1,2) 0 0]+[0 0 0 0 g*fad(1,3)];
        hp1 = H_up(7:2:21)';
        hp2 = H_up(8:2:22)';
        % Setup channel matrix with oversampling factor 2
        H1 = [hp1'; hp2'];
        
        % Construct concatenating channel matrix
        Hc_row = P*K;
        Hc_col = L+K;
        Hc = zeros(Hc_row,Hc_col);
        shift = 1;
        for i = 1:P:Hc_row
            Hc(i:i+1 ,shift:shift+7 ) = H1;
            shift = shift+1;
        end
        
        % Subroutine for generating transmit symbol vector: S and Sc
      
        % 4QAM source signal
       
        % Generate transmit symbol block signal
        
        % Generate signal(Block signal)
        S_iq = randi([0,1],2,300)*2-1;
        S = sqrt(1/2)*(S_iq(1,:)+1i*S_iq(2,:));
        Sc = zeros(L+K , Ns);
        for i=1:L+K
            Sc(i,:) = [zeros(1,i-1) , S(1:Ns-i+1)];
        end
        
        % Construct received signal model
        nc = normrnd(0,sqrt(NPW/2),size(Hc*Sc)) + 1i*normrnd(0,sqrt(NPW/2),size(Hc*Sc));
        Xc = Hc*Sc + nc;
        
        %--- Correlation Matrix Generation --- Part I.
        Rcin = (1/Ns)*(Xc-Hc*Sc)*(Xc-Hc*Sc)';
        Rcx = Hc*Hc' + Rcin;
        hcpr = Hc(:,d+1);
        rcxs = hcpr;
        
        %--- Weight Generation ---Part I.
        % Cobstruct solution of ZF LE & MMSE LE
        % ZF solution
        ed = [zeros(1,d) 1 zeros(1,(L+K)-d-1)]';
        qZF = Hc* inv(Hc' * Hc)*ed;
        
        % MMSE solution
        qMS = inv(Rcx)*rcxs;
        %--- Correlation matrix generation ---Part II.
        % Subroutine for generating FFF signal and interference-plus-noise
        Hcf = Hc(:,1:d+1);
        hcpr = Hc(:,d+1);
        Hcpa_1 = Hc(:,d+2:d+1+Kb);
        Hcpa_2 = Hc(:,d+2:d+1+Kb_2);
        Rcf = Hcf*Hcf' + Rcin;

        
        %--- Weight generation ---Part II.
        % Subroutine for generating FFF and FBF weight vector: qf and qb
        qf = inv(Rcf)*hcpr;
        qb_1 = Hcpa_1'*qf;
        qb_2 = Hcpa_2'*qf;
        
        
        %--- ZF, MMSE and DFE output SINR calculation ---
        %  Hint: use the formulation of page 33 and 39 with s(k-d) is known
        % Subroutine for calculating SINRo
        SNRo_ZF(n) = SNRo_ZF(n) + 1/(NPW*(qZF')*qZF);

        MMSE_S = sum(abs(qMS'*hcpr*Sc(d+1,:)).^2)/Ns;
        MMSE_N = sum(abs(qMS'*(Xc-hcpr*Sc(d+1,:))).^2)/Ns;
        SNRo_MMSE(n) = SNRo_MMSE(n) + MMSE_S/MMSE_N;

        DFE_5_S = sum(abs(qf'*hcpr*Sc(d+1,:)).^2)/Ns;
        DFE_5_N = sum(abs(qf'*(Xc-hcpr*Sc(d+1,:)-Hcpa_1*Sc(d+2:(d+2)+Kb-1,:))).^2)/Ns;
        SNRo_DFE_5(n) = SNRo_DFE_5(n) + DFE_5_S/DFE_5_N;

        DFE_7_S = sum(abs(qf'*hcpr*Sc(d+1,:)).^2)/Ns;
        DFE_7_N = sum(abs(qf'*(Xc-hcpr*Sc(d+1,:)-Hcpa_2*Sc(d+2:(d+2)+Kb_2-1,:))).^2)/Ns;
        SNRo_DFE_7(n) = SNRo_DFE_7(n) + DFE_7_S/DFE_7_N;


        
    end
    
end

SNRo_ZF = 10*log10(SNRo_ZF / Trial);
SNRo_MMSE = 10*log10(SNRo_MMSE / Trial);
SNRo_DFE_5 = 10*log10(SNRo_DFE_5 / Trial);
SNRo_DFE_7 = 10*log10(SNRo_DFE_7 / Trial);



%% plot should not be removed
plot(ref2_dB,SNRo_ZF,'bo--',ref2_dB,SNRo_MMSE,'r^--',ref2_dB,SNRo_DFE_5,'k*--',ref2_dB,SNRo_DFE_7,'md--')
xlabel('input SNR');
ylabel('output SNR');
legend('ZF','MMSE', 'DFE kb = 5','DFE kb = 7');
