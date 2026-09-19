function hd=hd(rho)
global beta beta1
hd=beta*(log(rho)-log(1-rho))+(beta1/2)*(1-2*rho);
end

