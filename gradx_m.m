function gradxA=gradx_m(A)
global hx
gradxA=(A-A(:,[end 1:end-1]))/hx;
%LaplaceA=(A([2:end 1],:)+A([end 1:end-1],:)+A(:,[2:end 1])+A(:,[end 1:end-1])-4*A)/h^2;
end