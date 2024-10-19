function [ symbol_MMSE ] = MMSE(Da_Str,y,H,NPW)

%% chapter5 p. 119 MMSE

symbol_MMSE = zeros(Da_Str,1);
W = (inv(H'*H + Da_Str*NPW*eye(Da_Str) ) * H')' ;
symbol_MMSE = qamdemod(W'*y,4);
symbol_MMSE = qammod(symbol_MMSE,4)/sqrt(2);

end
