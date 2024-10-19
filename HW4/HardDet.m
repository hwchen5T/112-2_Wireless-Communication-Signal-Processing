function y = HardDet(x,mod)  %% Design your own desion maker
    
    de_r = real(x);
    de_i = imag(x);
    de_r((de_r > 0)) = 1;
    de_r((de_r < 0)) = -1;
    de_i((de_i > 0)) = 1;
    de_i((de_i < 0)) = -1;
    y = de_r + 1j*de_i;

end
