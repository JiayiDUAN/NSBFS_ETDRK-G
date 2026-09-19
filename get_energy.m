%get_energy.m - compute the energy functional
function G = get_energy(u,v,phi,rho)
global  hx hy epsilon theta eta1 eta2 We
f1=(eta1/2)*Laplace2d_PNII(phi).^2;
f2=-(epsilon/2)*phi.*Laplace2d_PNII(phi);
f3=f(phi)/epsilon;
f4=-(eta2/2)*rho.*Laplace2d_PNII(rho);
f5=hrho(rho);
f6=-theta*rho.*g2(phi);
f7=(We/2)*(averx(u.^2)+averyDtoN(v.^2));
G=hx*hy*sum(sum(f1+f2+f3+f4+f5+f6+f7));
end