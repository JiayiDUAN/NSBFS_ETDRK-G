function gradxA=gradx(A)
global hx
gradxA=(A(:,[2:end 1])-A)/hx;
end