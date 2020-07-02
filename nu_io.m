function nu = nu_io(n_n,sig,v)
% n_n in cm^-3, sig in cm^2, v in cm/s
% nu in 1/s

nu = n_n.*sig.*v;