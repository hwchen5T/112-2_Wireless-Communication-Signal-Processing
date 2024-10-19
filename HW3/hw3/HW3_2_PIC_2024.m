clear all;
close all;
clc;

%% Parameter setting
Ns = 2000;                       	% samples for estimation related parameters
Q = 3;                        		% number of active users
L = 4;                        		% number of paths (fingers)
N = 31;                           	% random code length
trial = 100;                      	% number of Monte Carlo runs
path = [0 1 2 3 4];
J = length(path);                   % J=5

SNR = -10:5:20;                     % signal-to-noise ratio

BER = zeros(2,length(SNR));

%% Simulation experiments

for kk = 1:length(SNR)
    
    for jj = 1:trial
        
        %-----------------symbol matrix----------------------
        S_detect = zeros(Q,Ns);
        Sa_detect = zeros (Q,Ns);
        Sb_detect = zeros (Q,Ns);
        S = 2*randi([0 1],Q,Ns)-1;
        
        %-----------------noise matrix----------------------
        NPW = 1;                                                            % noise power
        noise = sqrt(NPW/2)*(randn(N+L,Ns)+1i*randn(N+L,Ns));               % noise
        
        for ii = 1:Ns
            
            %------------fading matrix----------------------
            for q = 1:Q
                fad(q,1:L+1) = randn(1,L+1);                          % fading gain
            end
            for qq = 1:Q                                              % fading normalization
                fad(qq,:) = fad(qq,:)./sqrt(sum(abs(fad(qq,:)).^2));
            end
            
            %------------------spreading code------------------------
            c = randi([0,1],Q,N)*2-1;
            
            %------------------perfect channel estimation for H------
            C_q = zeros(N+L,J,Q);
            for qth = 1:Q
                for jth = 1:J
                    C_q(:,jth,qth) = [zeros(1,jth-1),c(qth,:),zeros(1,N+L-N-(jth-1))]';
                end
            end
            
            H = zeros(N+L,Q);

            for qth = 1:Q
                H(:,qth) = C_q(:,:,qth)*fad(qth,:)';
            end
            %------------------transmitted symbol--------------------
            SNR1 = SNR(kk);
            SNR2 = SNR(kk)-6;
            SNR3 = SNR(kk)-12;
            SPW1 = 10^(SNR1/10);
            SPW2 = 10^(SNR2/10);
            SPW3 = 10^(SNR3/10);

            A = [sqrt(SPW1) 0 0 ; 0 sqrt(SPW2) 0 ; 0 0 sqrt(SPW3)];

            x = H*A*S(:,ii) + noise(:,ii); 
            
            
            %------------------detection PIC----------------------
            x_PIC = x ;
            D_rake = H;
            z = D_rake'*x_PIC;
            for qth = 1:Q
                if real(z(qth))>=0
                    S_detect(qth,ii)=1;
                else
                    S_detect(qth,ii)=-1;
                end
            end
            y1 = x_PIC - H(:,2)*A(2,2)*S_detect(2,ii)-H(:,3)*A(3,3)*S_detect(3,ii);
            y2 = x_PIC - H(:,1)*A(1,1)*S_detect(1,ii)-H(:,3)*A(3,3)*S_detect(3,ii);
            y3 = x_PIC - H(:,2)*A(2,2)*S_detect(2,ii)-H(:,1)*A(1,1)*S_detect(1,ii);
           
            z_re1 = D_rake(:,1)'*y1;
            z_re2 = D_rake(:,2)'*y2;
            z_re3 = D_rake(:,3)'*y3;
            z_re = [z_re1,z_re2,z_re3];
            for qth = 1:Q
                if real(z_re(qth))>=0
                    Sb_detect(qth,ii)=1;
                else
                    Sb_detect(qth,ii)=-1;
                end
            end
            %test
            y1t = x_PIC - H(:,2)*A(2,2)*Sb_detect(2,ii)-H(:,3)*A(3,3)*Sb_detect(3,ii);
            y2t = x_PIC - H(:,1)*A(1,1)*Sb_detect(1,ii)-H(:,3)*A(3,3)*Sb_detect(3,ii);
            y3t = x_PIC - H(:,2)*A(2,2)*Sb_detect(2,ii)-H(:,1)*A(1,1)*Sb_detect(1,ii);
            z_re1t = D_rake(:,1)'*y1t;
            z_re2t = D_rake(:,2)'*y2t;
            z_re3t = D_rake(:,3)'*y3t;
            z_ret = [z_re1t,z_re2t,z_re3t];
            for qth = 1:Q
                if real(z_ret(qth))>=0
                    Sa_detect(qth,ii)=1;
                else
                    Sa_detect(qth,ii)=-1;
                end
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

H_a = figure;

semilogy(SNR, BER(1,:), '-bo', SNR, BER(2,:), '-rx');
grid on;
legend('PIC_a', 'PIC_b');
xlabel('SNR (dB)');
ylabel('BER');


