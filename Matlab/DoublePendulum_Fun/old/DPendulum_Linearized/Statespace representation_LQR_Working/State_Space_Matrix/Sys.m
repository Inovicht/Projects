function sys = Sys(params) 
[sys.A,sys.K,sys.D,sys.M] = Amat(params);             %Generation of state matrix
[sys.Cc,sys.Ca] = Cmat(sys.M,sys.K,sys.D,params.d,params.ddot,params.dddot);    %Generation of output matrix
[sys.B,sys.B2] = Bmat(sys.M,params.q);                                          %Generation of input matrix
if (isempty(params.dddot))                                                      %Generation of Transmission matrix
    sys.H = 0;
else
    sys.H = Hmat(sys.Ca,sys.M,sys.B2);
end