function hrho=hrho(rho)
global beta beta1
hrho=beta*(rho.*log(rho)+(1-rho).*log(1-rho))+(beta1/2)*rho.*(1-rho);
end

