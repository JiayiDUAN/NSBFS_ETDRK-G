function mu_rho=mu_rho(phi,rho)
global eta2 theta
mu_rho=-eta2*Laplace2d_PNII(rho)+hd(rho)-theta*g2(phi);
end