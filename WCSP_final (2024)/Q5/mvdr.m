function w = mvdr(x,sp,a_t0)


Rxx = x*x';
w = inv(Rxx)*a_t0/(a_t0'*inv(Rxx)*a_t0);

end