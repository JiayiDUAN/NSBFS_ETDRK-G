function LaplaceA=Laplace2d_PNII(A)
global hx hy
LaplaceA=(A(:,[2:end 1])+A(:,[end 1:end-1])-2*A)/(hx^2)...
                   +(A([2:end end],:)+A([1 1:end-1],:)-2*A)/(hy^2);
end