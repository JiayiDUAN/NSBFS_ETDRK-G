function gradyA=grady_mNeu(A)
global hy
gradyA=(A-A([1 1:end-1],:))/hy;
end