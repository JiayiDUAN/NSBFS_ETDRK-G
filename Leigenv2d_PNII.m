function lambda = Leigenv2d_PNII(Ny,Nx)
global hx hy 
kx = [0:Nx/2 -Nx/2+1:-1]; 
ky = [0:Ny-1]';  
[Kx,Ky] = meshgrid(kx,ky);
lambda = - (4/hx^2)*(sin(Kx*pi/Nx)).^2 - (4/hy^2)*(sin(Ky*pi/(2*Ny))).^2;
end