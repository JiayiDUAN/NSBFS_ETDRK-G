function varphi32 = varphi32(x)
varphi32 = (exp(x)-1-x)./x.^2 - (1/2)*(exp(x)-1)./x;
end