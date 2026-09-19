function gradyA=gradyNtoDin(A)
global hy
gradyA=(A([2:end],:)-A([1:end-1],:))/hy;
end