function w = mmse(x,sp,a_t0)

Rxx = x*x';
rxs = x*sp';

w = inv(Rxx)*rxs;

end