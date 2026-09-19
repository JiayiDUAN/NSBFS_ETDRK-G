function mu_phi=mu_phi(phi,rho)
global epsilon eta1 theta
mu_phi=eta1*Laplace2d_PNII(Laplace2d_PNII(phi))-epsilon*Laplace2d_PNII(phi)...
                 +fd(phi)/epsilon+2*theta*g1(phi,rho);
end