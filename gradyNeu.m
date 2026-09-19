function gradyA=gradyNeu(A)
global hy
gradyA=(A([2:end end],:)-A)/hy;
end