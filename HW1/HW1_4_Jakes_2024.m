% hn = Jakes(M,f_max,Ts)
% It's a subroutine which is about to generate
% the Doppler spread fading channel based on Jakes model
%
% Input:
%       M: a parameter ralated to no. of paths
%          (no. of path N = 4*M + 2)
%       f_max: maximum Doppler shift (Hz)
%       Ts: sampling rate (sec)
% Output:
%       hn = normalized channel in the interval, [0,10^3*Ts]

function hn = HW1_4_Jakes_2024(M,f_max,Ts)

N = 4*M + 2;    % no. of paths
alpha = 0;      % a parameter ralated to the inner product properties of the output channel
n = 1:1:M;
t = 0:Ts:10^3*Ts;   % time interval
beta = n*pi/M;

% construct h using the parameters above

h_i = 0;
h_q = 0;
f_n= f_max * cos(2*pi*n/N) ;

for j=1:M

    h_i_temp = cos(2*pi*f_n(j)*t)*2*cos(beta(j));
    h_q_temp = cos(2*pi*f_n(j)*t)*2*sin(beta(j));
    h_i = h_i + h_i_temp;
    h_q = h_q + h_q_temp;
end

h_i = h_i + 2*(1/sqrt(2))*cos(2*pi*f_max*t)*cos(alpha);
h_q = h_q + 2*(1/sqrt(2))*cos(2*pi*f_max*t)*sin(alpha);

h = h_i + 1i*h_q;


% Output of the function

power = sum((abs(h)).^2)/length(h);
hn = h/sqrt(power); % normalized channel

end
