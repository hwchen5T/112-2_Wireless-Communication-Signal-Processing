clear all;
close all;
clc;

%% Parameter setting
Ns = 2000;                       	% samples for estimation related parameters
Q = 3;                        		% number of active users
L = 4;                        		% number of paths (fingers)
N = 31;                           	% random code length
trial = 100;                      	% number of Monte Carlo runs

SNR = -10:5:20;                     % signal-to-noise ratio

BER = zeros(2,length(SNR));

%% Simulation experiments

for kk = 1:length(SNR)
    
    for jj = 1:trial
        
        %-----------------symbol matrix----------------------
        Sa_detect = zeros (Q,Ns);
        Sb_detect = zeros (Q,Ns);
        S = 2*randi([0 1],Q,Ns)-1;
        
        %-----------------noise matrix----------------------
            NPW=10^(-SNR(kk)/10);          % noise power
             noise=(NPW/2)^0.5 * (randn(Q,N+L) +1i*randn(Q,N+L));      % noise
        
        for ii = 1:Ns
            
            %------------fading matrix----------------------
            for q = 1:Q
                fad(q,1:L+1) = randn(1,L+1);                          % fading gain
            end
            for qq = 1:Q                                              % fading normalization
                fad(qq,:) = fad(qq,:)./sqrt(sum(abs(fad(qq,:)).^2));
            end
            
            %------------------spreading code------------------------
            spread_code = 2*round(rand(Q,N)) - 1;
            
            %------------------perfect channel estimation for H------
             H=[];
            for q=1:Q
                C=zeros(N+L,L+1);
                for index=1:L+1
                    for n =1:N
                        C( index+n-1 ,index) =  spread_code( q,n );
                    end
                end
                temp=C*transpose(fad(q,:));
                H=[H,temp];
            end      
            
            %------------------transmitted symbol--------------------
                SNR1=10^(SNR(kk)/10);
                SNR2=10^((SNR(kk)-6)/10);
                SNR3=10^((SNR(kk)-12)/10);
                SPW1 = sqrt(SNR1*NPW);
                SPW2 = sqrt(SNR2*NPW);
                SPW3 = sqrt(SNR3*NPW);
                SPW = [SPW1,SPW2,SPW3];
                A=[SPW1 0 0; 0 SPW2 0; 0 0 SPW3];
                x=H*A*S(:,ii)+noise';
            %------------------detection SIC(a)----------------------
            for q=1:Q
                d_r=H(:,q);
                z=d_r'*x;
                if real(z)>=0
                     Sa_detect(q,ii)=1;
                else
                     Sa_detect(q,ii)=-1;
                end
                x=H*A*S(:,ii)+noise';
                x = x - H*A* Sa_detect(:,ii) ;
            end
             x=H*A*S(:,ii)+noise';
            %------------------detection SIC(b)----------------------
               for q=1:Q
                d_r=H(:,Q-q+1);
                z=d_r'*x;
                if real(z)>=0
                     Sb_detect(Q-q+1,ii)=1;
                else
                     Sb_detect(Q-q+1,ii)=-1;
                end
                x=H*A*S(:,ii)+noise';
                x = x - H*A* Sb_detect(:,ii) ;
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

H_f = figure;

semilogy(SNR, BER(1,:), '-bo', SNR, BER(2,:), '-rx');
grid on;
legend('SIC_a', 'SIC_b');
xlabel('SNR (dB)');
ylabel('BER');


