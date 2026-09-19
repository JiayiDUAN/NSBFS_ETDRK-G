function uhat=dstIbyfft(u,Nx,Ny)
v=[zeros(1,Nx); u; zeros(1,Nx); -flip(u,1)];
vhat=fft(v)/(-2i);
uhat=real(vhat(2:Ny,:));
end
