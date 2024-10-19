clear all;
close all;
clc;

%% Parameter setting

SIR = 0;                      		% signal-to-CCI ratio
IPW = 10^(-SIR/10);                   % interference power
L = 4;                        		% number of paths (fingers)
N = 31;                               % random code length
trial = 5000;                         % number of Monte Carlo runs
path = [0 1 2 3 4];                   % path delays assume uniform[0,4Tc]
J = length(path);                     % Number of paths

SNR = 0;
NPW = 10^(-SNR/10);

QQ = 1:10;                            % number of active users

%% Simulation experiments

rMS = zeros(trial,length(QQ));
rZF = zeros(trial,length(QQ));
rRAKE = zeros(trial,length(QQ));

for ii = 1:length(QQ)
    for jj = 1:trial
        
        Q = QQ(ii);
        
        %-------------------------fading matrix
        for q = 1:Q
            fad(q,1:L+1) = randn(1,L+1);                          % fading gain
        end
        for qq = 1:Q                                            %fading  normalization
            fad(qq,:) = fad(qq,:)./sqrt(sum(abs(fad(qq,:)).^2));
        end
        
        %----------------------------spreading code
        c = randi([0,1],Q,N)*2-1;
        
        
        %----------------------------perfect channel estimation for H and h1
        C_q = zeros(N+L,J,Q);

        for qth = 1:Q
            for jth = 1:J
                C_q(:,jth,qth) = [zeros(1,jth-1) , c(qth,:) , zeros(1,(N+L)-N-(jth-1)) ]' ;
            end
        end

        H = zeros(N+L,Q);

        for qth = 1:Q
            H(:,qth) = C_q(:,:,qth)*fad(qth,:)' ;
        end
        
        %----------------------------correlation matrix generation

        Rin = IPW*H(:,2:Q)*(H(:,2:Q)') + NPW*eye(N+L) ;
        Rxx = IPW*H(:,:)*(H(:,:)') + NPW*eye(N+L) ;
        rxs = IPW*H(:,1);
        
        
        %----------------------------SINR for RAKE
        % RAKE receiver
        h1 = H(:,1);
        d_RAKE = h1;
        % RAKE receiver output SINR
        rRAKE(jj,ii) = rRAKE(jj,ii) + abs(IPW*d_RAKE'*h1)^2 / (d_RAKE'*Rin*d_RAKE);
        
        %----------------------------SINR for ZF
        % ZF receiver
        e0 = [1 , zeros(1,Q-1)]';
        d_ZF = H*inv(H'*H)*e0;
        % ZF receiver output SINR
        rZF(jj,ii) = rZF(jj,ii) + abs(IPW*d_ZF'*h1)^2 / (d_ZF'*Rin*d_ZF) ;
        
        %----------------------------SINR for MMSE
        % MMSE detector
        d_MS = inv(Rxx)*rxs;
        % MMSE detector output SINR
        rMS(jj,ii) = rMS(jj,ii) + abs(IPW*d_MS'*h1)^2 / (d_MS'*Rin*d_MS);
        
    end
end

rRAKE_AV = sum(rRAKE)/trial;
rZF_AV = sum(rZF)/trial;
rMS_AV = sum(rMS)/trial;

%% Simulation results

H1 = figure;

plot(QQ,10*log10(rMS_AV),'-b+',QQ,10*log10(rZF_AV),'-r+',QQ,10*log10(rRAKE_AV),'-g+');
h = legend('MMSE','ZF','RAKE','location','SouthWest');
set(h,'Fontsize',15);
xlabel('Q');
ylabel('SINRo (dB)');
