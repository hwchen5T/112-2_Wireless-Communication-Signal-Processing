clear all;
close all;
clc;

%% Parameter setting
Ns = 2000;                       	% samples for estimation related parameters
Q = 3;                        		% number of active users
L = 4;                        		% number of paths (fingers)
N = 31;                           	% random code length
trial = 100;                      	% number of Monte Carlo runs
path_delay = [0 1 2 3 4];
J = length(path_delay);

SNR = -10:5:20;                     % signal-to-noise ratio

BER = zeros(2,length(SNR));

%% Simulation experiments

for kk = 1:length(SNR)
    SNR_list = [SNR(kk) SNR(kk)-6 SNR(kk)-12];
    for jj = 1:trial
        
        %-----------------symbol matrix----------------------
        Sa_detect = zeros (Q,Ns);
        Sb_detect = zeros (Q,Ns);
        S = 2*randi([0 1],Q,Ns)-1;
        SPW = 1;
        
        %-----------------noise matrix----------------------
        NPW = 10.^(-SNR_list/10);                                      % noise power
        noise = zeros(N+L,Ns,Q);
        noise(:,:,1) = sqrt(NPW(1)/2)*(randn(N+L,Ns)+1j*randn(N+L,Ns));         % noise
        noise(:,:,2) = sqrt(NPW(2)/2)*(randn(N+L,Ns)+1j*randn(N+L,Ns));
        noise(:,:,3) = sqrt(NPW(3)/2)*(randn(N+L,Ns)+1j*randn(N+L,Ns));
        for ii = 1:Ns
            
            %------------fading matrix----------------------
            for q = 1:Q
                fad(q,1:L+1) = randn(1,L+1);                          % fading gain
            end
            for qq = 1:Q                                              % fading normalization
                fad(qq,:) = fad(qq,:)./sqrt(sum(abs(fad(qq,:)).^2));
            end
            
            %------------------spreading code------------------------
            c = zeros(Q,N);
            c = randi([0 1],Q,N)*2-1;
            
            %------------------perfect channel estimation for H------
            C = zeros(N+L,J,Q);
            for qidx = 1:Q
                for jidx = 1:J
                      C(jidx:jidx+N-1,jidx,qidx) = c(qidx,:);
%                     C(:,jidx,qidx) = [zeros(1,jidx-1) c(qidx,:) zeros(1,N+L-(jidx-1)-N)]';
                end
            end
            H = zeros(N+L,Q);
            for qidx=1:Q
                H(:,qidx) = C(:,:,qidx)*fad(qidx,:)' ;
            end
            
            %------------------transmitted symbol--------------------
            x = H*S ; 
            y = x + noise; 
            
            %------------------detection SIC(a)----------------------
            y_descend = y;
            for qidx=1:Q
                z = H(:,qidx)'*y_descend(:,:,qidx);
                z((real(z)>0)) = 1;
                z((real(z)<0)) = -1;
                Sa_detect(qidx,:) = z;
                y_descend(:,:,qidx) = y_descend(:,:,qidx) - H(:,qidx)*Sa_detect(qidx,:);
            end
            
            %------------------detection SIC(b)----------------------
            y_ascend = y;
            for qidx=Q:-1:1
                z = H(:,qidx)'*y_ascend(:,:,qidx);
                z((real(z)>0)) = 1;
                z((real(z)<0)) = -1;
                Sb_detect(qidx,:) = z;
                y_ascend(:,:,qidx) = y_ascend(:,:,qidx) - H(:,qidx)*Sb_detect(qidx,:);
            end
        end
        
        %----------------------------BER calculation----------------------
        error = 0;
        for k = 1:Ns
            for n = 1:Q
                if (S(n,k) ~= Sa_detect(n,k))
                    error = error+1;
                end
            end
        end
        BER(1,kk) = error/(Q*Ns) + BER(1,kk);
        
        error = 0;
        for k = 1:Ns
            for n = 1:Q
                if (S(n,k) ~= Sb_detect(n,k))
                    error = error+1;
                end
            end
        end
        BER(2,kk) = error/(Q*Ns) + BER(2,kk);
        
    end
end

BER = BER/trial;

%% Simulation results

H = figure;

semilogy(SNR, BER(1,:), '-bo', SNR, BER(2,:), '-rx');
grid on;
legend('SIC_a', 'SIC_b');
xlabel('SNR (dB)');
ylabel('BER');


