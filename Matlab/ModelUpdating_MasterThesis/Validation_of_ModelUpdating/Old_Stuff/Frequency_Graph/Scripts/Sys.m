function sys = Sys(params,theta)
[params.K_beam,params.Kg_beam,sys.M,params.freedof] =beam2D(params,theta);
[sys.A,sys.K,sys.D] = Amat(params,sys.M,theta);
[sys.Cc,sys.Ca] = Cmat(sys.M,sys.K,sys.D,params.d,params.ddot,params.dddot);
[sys.B,sys.B2] = Bmat(sys.M,params.q);
if (isempty(params.dddot))
    sys.H = 0;
else
    sys.H = Hmat(sys.Ca,sys.M,sys.B2);
end
sys.lambda = eig(sys.A);
