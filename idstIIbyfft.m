function v=idstIIbyfft(vhat,Omebar,Nx,Ny)
ymy=flip(vhat(1:end-1,:))-1i*vhat(1:end-1,:);
ff=[vhat(end,:); Omebar.*ymy];
uu=real(fft(conj(ff))/Ny);
v=zeros(Ny,Nx);
v([1:2:Ny 2:2:Ny],:)=[uu(1:Ny/2,:); -flip(uu(Ny/2+1:end,:),1)]; 
end