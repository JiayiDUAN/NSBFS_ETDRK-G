function vhat=dstIIbyfft(v,Ome,Ny)
ff=fft([v(1:2:Ny,:); -flip(v(2:2:Ny,:),1)]);
vhat=[-imag(Ome.*ff(2:end,:)) ; ff(1,:)];
end