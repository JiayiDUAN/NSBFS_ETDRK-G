function F2= F2(u,v,phi,rho)
global theta B M2
F2=M2*Laplace2d_PNII(hd(rho)-B*rho-theta*g2(phi))...
       -gradx(u.*averx_m(rho))-gradyDtoN(v.*averyNtoD(rho));  
end
