function [ symbol_MMSE_OSIC ] = MMSE_OSIC(Da_Str,y,H,NPW)

F = H'*H;
in_match = H'*y;
symbol_MMSE_OSIC = zeros(Da_Str,1);
dec_table = 1:Da_Str;

for i = 1:Da_Str

    G = H *inv( F + (NPW*Da_Str)*eye(size(F)) );
    [~ , index] = sort(diag(G'*G));
    w = G(:,index(1));
    z = w' * y ;
    x_est = qamdemod(z,4);
    x_est = qammod(x_est,4)/sqrt(2);
    symbol_MMSE_OSIC(dec_table(index(1))) = x_est;
    
    h = H(:,index(1));
    y = y - h * x_est;

    dec_table(index(1)) = [];
    H(:,index(1)) = [];
    in_match = H'*y;
    F = H'*H;

end


end
