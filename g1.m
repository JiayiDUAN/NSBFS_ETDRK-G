function g1=g1(phi,rho)
g1= gradx_m(averx(rho).*gradx(phi))...        
         +gradyDtoN(averyNtoD(rho).*gradyNtoD(phi));
end
