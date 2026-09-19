function u=idstIbyfft(uhat,Nx,Ny)
vhat=[zeros(1,Nx); uhat; zeros(1,Nx); -flip(uhat,1)];
v=fft(vhat)/(-2i);
u=real(v(2:Ny,:))*(2/Ny);
end

