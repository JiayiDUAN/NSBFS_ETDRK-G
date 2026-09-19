function draw_figure(N,T,dt,dtn)
% clear, clc, close all
% N=128; T=0.9;  dt=0.001; dtn=0.1;
tvec=0:dt:T; tn=0:dtn:T; n=size(tn,2)
Lx = 2*pi; Ly=pi; Nx=2*N; Ny=N;
hx = Lx/Nx; hy = Ly/Ny; 
xnode = Lx/Nx*[0:Nx-1]'+hx/2;
ynode = Ly/Ny*[0:Ny-1]'+hy/2; 
[X,Y] = meshgrid(xnode,ynode);
load('energy.mat')
load('max_min.mat')
max(phimax)
min(phimin)
max(rhomax)
min(rhomin)
figure(1), plot(tvec(1:end),phimax(1:end),'k','Linewidth',1), title('Maximum and minimum of phi'),xlabel('t'), ylabel('phi')
hold on, plot(tvec(1:end),phimin(1:end),'b','Linewidth',1)
figure(2), plot(tvec(1:end),rhomax(1:end),'Linewidth',1), title('Maximum and minimum of rho'),xlabel('t'), ylabel('rho')
hold on, plot(tvec(1:end),rhomin(1:end),'b','Linewidth',1)
figure(3);
plot(tvec(1:end), Energy(1:end),'k','Linewidth',1), title('Energy'),xlabel('t'), ylabel('G');

load('phase_seperation.mat')
for i=n:-1:1
    figure,h=mesh(X,Y,phitn(:,(i-1)*Nx+1:i*Nx)),axis([0 Lx 0 Ly]), axis equal, colorbar, view([0,-90])
    set(h, 'FaceColor', 'interp');
end

diff_Energy=Energy(2:end)-Energy(1:end-1);
energy_increase_step=find(diff_Energy>0)+1

