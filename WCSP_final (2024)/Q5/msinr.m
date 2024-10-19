function w = msinr(x,sp,a_t0)


Rin = (x - a_t0*sp)*(x - a_t0*sp)';

w = inv(Rin)*a_t0 / (a_t0'*inv(Rin)*a_t0);

end