function gradyA=gradyNtoD(A)
global hy
gradyA=(A([1:end end],:)-A([1 1:end],:))/hy;
end