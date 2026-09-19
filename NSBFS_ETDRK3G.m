%clear, clc, close all
function  [time,u,v,phi,rho]=NSBFS_ETDRK3G(N, T, dt, dtn)
global Lx Ly hx hy A1 A2 B eta1 epsilon theta eta2 beta beta1 nu M1 M2 We 
Lx = 2*pi; Ly=pi; Nx=2*N; Ny=N;
hx = Lx/Nx; hy = Ly/Ny;  
M1=1; M2=1;
eta1=0.0002; epsilon = 0.05; theta=0.002;
eta2=0.0004; beta=0.05; beta1=0.15;
A1=4; B=6; A2=0.05;     %clear, clc, close all; NSBFS_ETDRK3(128, 0.9, 0.001, 0.1)
nu=0.01; We=1;  g1value=3.5; 
xnode = Lx/Nx*[0:Nx-1]'+hx/2;
ynode = Ly/Ny*[0:Ny-1]'+hy/2; 
[X,Y] = meshgrid(xnode,ynode); 
uxnode = Lx/Nx*[0:Nx-1]';
uynode = Ly/Ny*[0:Ny-1]'+hy/2; 
[Xu,Yu] = meshgrid(uxnode,uynode); 
vxnode = Lx/Nx*[0:Nx-1]'+hx/2;
vynode = Ly/Ny*[0:Ny]'; 
[Xv,Yv] = meshgrid(vxnode,vynode); 
u0=zeros(Ny,Nx); v0=zeros(Ny+1,Nx); p0=zeros(Ny,Nx);
Ov=zeros(1,Nx); Oin=zeros(Ny-2,Nx);
g1=g1value*ones(1,Nx); gup=g1; glow=-g1;
phi0=tanh((0.2*pi-sqrt((X-pi-0.19*pi).^2+(Y-pi/2-0.22*pi).^2))/(epsilon*2^0.5))...
          +tanh((0.2*pi-sqrt((X-pi+0.19*pi).^2+(Y-pi/2+0.22*pi).^2))/(epsilon*2^0.5))+1;
rho0=0.005+0.2*(1-abs(phi0));
% rho0=0.02+0.8*(1-abs(phi0));
ve=g1value; 
u0=ve*(Yu-Ly/2)/(Ly/2);
u0=-u0;

phi0_mean=sum(sum(phi0))/Nx/Ny
rho0_mean=sum(sum(rho0))/Nx/Ny

lambdaPNII = Leigenv2d_PNII(Ny,Nx);
lambdaPDI = Leigenv2d_PDI(Ny,Nx);
lambdaPDII = Leigenv2d_PDII(Ny,Nx);
L1hat= M1*(eta1*lambdaPNII.^3-(epsilon+A2)*lambdaPNII.^2 +(A1/epsilon)*lambdaPNII); 
L2hat= M2*(-eta2*lambdaPNII.^2+B*lambdaPNII); 
L3uhat=nu*lambdaPDII;
L3vhat=nu*lambdaPDI;
varphi0dtL1= varphi0(dt*L1hat);
varphi0dtL2= varphi0(dt*L2hat);
varphi0dtL3u= varphi0(dt*L3uhat);
varphi0dtL3v= varphi0(dt*L3vhat);
varphi023dtL1= varphi0((2/3)*dt*L1hat);
varphi023dtL2= varphi0((2/3)*dt*L2hat);
varphi023dtL3u= varphi0((2/3)*dt*L3uhat);
varphi023dtL3v= varphi0((2/3)*dt*L3vhat);
varphi1dtL1= varphi1(dt*L1hat);
varphi1dtL2= varphi1(dt*L2hat);
varphi1dtL3u= varphi1(dt*L3uhat);
varphi1dtL3v= varphi1(dt*L3vhat);
varphi2123dtL1= varphi21((2/3)*dt*L1hat);
varphi2123dtL2= varphi21((2/3)*dt*L2hat);
varphi2123dtL3u= varphi21((2/3)*dt*L3uhat);
varphi2123dtL3v= varphi21((2/3)*dt*L3vhat);
varphi2223dtL1= varphi22((2/3)*dt*L1hat);
varphi2223dtL2= varphi22((2/3)*dt*L2hat);
varphi2223dtL3u= varphi22((2/3)*dt*L3uhat);
varphi2223dtL3v= varphi22((2/3)*dt*L3vhat);
varphi31dtL1= varphi31(dt*L1hat);
varphi31dtL2= varphi31(dt*L2hat);
varphi31dtL3u= varphi31(dt*L3uhat);
varphi31dtL3v= varphi31(dt*L3vhat);
varphi32dtL1= varphi32(dt*L1hat);
varphi32dtL2= varphi32(dt*L2hat);
varphi32dtL3u= varphi32(dt*L3uhat);
varphi32dtL3v= varphi32(dt*L3vhat);
varphi33dtL1= varphi33(dt*L1hat);
varphi33dtL2= varphi33(dt*L2hat);
varphi33dtL3u= varphi33(dt*L3uhat);
varphi33dtL3v= varphi33(dt*L3vhat);
varphi1dtL1(1,1)=1; 
varphi1dtL2(1,1)=1; 
varphi2123dtL1(1,1)=4/9;
varphi2123dtL2(1,1)=4/9;
varphi2223dtL1(1,1)=2/9;
varphi2223dtL2(1,1)=2/9;
varphi31dtL1(1,1)=1/4;
varphi31dtL2(1,1)=1/4;
varphi32dtL1(1,1)=0;
varphi32dtL2(1,1)=0;
varphi33dtL1(1,1)=3/4;
varphi33dtL2(1,1)=3/4;

k=[1:Ny-1]';
omega=exp(-pi*1i*k/(2*Ny));
Ome=repmat(omega,1,Nx);
Omebar=conj(Ome);

q0=zeros(Ny,Nx); 
t = 0; u=u0; v=v0;  %p=p0;
m_u=u0; m_v=v0;  
phi = phi0; rho=rho0;
q=q0; dxqonDiub=gradx_m(averyNtoD(q)); 

tvec = t;
Energy=get_energy(u,v, phi,rho);
phimax = max(max(phi));
rhomax = max(max(rho));
phimin = min(min(phi));
rhomin = min(min(rho));

tn=0:dtn:T; %the time saving snapshots and velocity
n=2;
phitn=phi0;
rhotn=rho0;
utn=u0;
vtn=v0;

lambdaa=lambdaPNII;
lambdaa(1,1)=1e-16;
nstepend=round(T/dt)-4;
for nstep = 1: nstepend  % from t=0 to t=T-4*dt
    t = t + dt
    m_vin=m_v(2:end-1,:);
    m_uhat=fft(dstIIbyfft(m_u,Ome,Ny)')';
    m_vinhat=fft(dstIbyfft(m_vin,Nx,Ny)')';
    gqupO=gup+dxqonDiub(1,:);
    gqlowO=glow+dxqonDiub(end,:);
    gqO=[gqupO; Oin; gqlowO];
    F3uhatO = fft(dstIIbyfft(nu*2*gqO/hy/hy+F3u(gup,glow,u,v,phi,rho),Ome,Ny)')';
    F3vhatO = fft(dstIbyfft(F3v(u,v,phi,rho),Nx,Ny)')';
    m_uhatA = varphi0dtL3u .*m_uhat + dt*varphi1dtL3u .* F3uhatO ;  
    m_vinhatA = varphi0dtL3v .*m_vinhat + dt*varphi1dtL3v .*F3vhatO ; 
    m_uA = idstIIbyfft(real(ifft(m_uhatA')'),Omebar,Nx,Ny); 
    m_vinA = idstIbyfft(real(ifft(m_vinhatA')'),Nx,Ny); 
    m_vA=[Ov; m_vinA ;Ov];
    bA=gradx(m_uA)+gradyDtoN(m_vA);
    bAhat=fft(dct(bA,'Type',2)')'; bAhat(1,1)=0;
    qA=idct(real(ifft((bAhat./lambdaa)')'),'Type',2);  
    dxqonDiubA=gradx_m(averyNtoD(qA)); 
    uA=m_uA-gradx_m(qA); 
    vA=m_vA-gradyNtoD(qA); vA(1,:)=Ov; vA(end,:)=Ov;

    phihat=fft(dct(phi,'Type',2)')';
    rhohat=fft(dct(rho,'Type',2)')';
    F1hatO = fft(dct(F1(u,v,phi,rho),'Type',2)')';
    F2hatO = fft(dct(F2(u,v,phi,rho),'Type',2)')';
    phihatA = varphi0dtL1.*phihat + dt*varphi1dtL1.*F1hatO ;  
    rhohatA = varphi0dtL2.*rhohat + dt*varphi1dtL2.*F2hatO ; 
    phiA = idct(real(ifft(phihatA')'),'Type',2); 
    rhoA = idct(real(ifft(rhohatA')'),'Type',2); 
    
    gqupA=gup+dxqonDiubA(1,:);
    gqlowA=glow+dxqonDiubA(end,:);
    gqA=[gqupA; Oin; gqlowA];
    F3uhatA = fft(dstIIbyfft(nu*2*gqA/hy/hy+F3u(gup,glow,uA,vA,phiA,rhoA),Ome,Ny)')';
    F3vhatA = fft(dstIbyfft(F3v(uA,vA,phiA,rhoA),Nx,Ny)')';
    m_uhatB = varphi023dtL3u.*m_uhat + dt*(varphi2123dtL3u.*F3uhatO + varphi2223dtL3u.*F3uhatA) ;  
    m_vinhatB = varphi023dtL3v.*m_vinhat + dt*(varphi2123dtL3v.*F3vhatO + varphi2223dtL3v.*F3vhatA) ;  
    m_uB = idstIIbyfft(real(ifft(m_uhatB')'),Omebar,Nx,Ny); 
    m_vinB = idstIbyfft(real(ifft(m_vinhatB')'),Nx,Ny); 
    m_vB=[Ov; m_vinB ;Ov];
    bB=gradx(m_uB)+gradyDtoN(m_vB);
    bBhat=fft(dct(bB,'Type',2)')'; bBhat(1,1)=0;
    qB=idct(real(ifft((bBhat./lambdaa)')'),'Type',2);  
    dxqonDiubB=gradx_m(averyNtoD(qB)); 
    uB=m_uB-gradx_m(qB); 
    vB=m_vB-gradyNtoD(qB); vB(1,:)=Ov; vB(end,:)=Ov;
    

    F1hatA = fft(dct(F1(uA,vA,phiA,rhoA),'Type',2)')';
    F2hatA = fft(dct(F2(uA,vA,phiA,rhoA),'Type',2)')';
    phihatB = varphi023dtL1.*phihat + dt*(varphi2123dtL1.*F1hatO + varphi2223dtL1.*F1hatA);  
    rhohatB = varphi023dtL2.*rhohat + dt*(varphi2123dtL2.*F2hatO + varphi2223dtL2.*F2hatA);  
    phiB = idct(real(ifft(phihatB')'),'Type',2); 
    rhoB = idct(real(ifft(rhohatB')'),'Type',2); 

    gqupB=gup+dxqonDiubB(1,:);
    gqlowB=glow+dxqonDiubB(end,:);
    gqB=[gqupB; Oin; gqlowB];
    F3uhatB = fft(dstIIbyfft(nu*2*gqB/hy/hy+F3u(gup,glow,uB,vB,phiB,rhoB),Ome,Ny)')';
    F3vhatB = fft(dstIbyfft(F3v(uB,vB,phiB,rhoB),Nx,Ny)')';
    m_uhat3 = varphi0dtL3u.*m_uhat + dt*(varphi31dtL3u .* F3uhatO + varphi32dtL3u .*F3uhatA + varphi33dtL3u .*F3uhatB);  
    m_vinhat3 = varphi0dtL3v.*m_vinhat + dt*(varphi31dtL3v .*F3vhatO + varphi32dtL3v .*F3vhatA + varphi33dtL3v .*F3vhatB) ; 
    m_u3 = idstIIbyfft(real(ifft(m_uhat3')'),Omebar,Nx,Ny); 
    m_vin3 = idstIbyfft(real(ifft(m_vinhat3')'),Nx,Ny); 
    m_v3=[Ov; m_vin3 ;Ov];
    b3=gradx(m_u3)+gradyDtoN(m_v3);
    b3hat=fft(dct(b3,'Type',2)')'; b3hat(1,1)=0;
    q3=idct(real(ifft((b3hat./lambdaa)')'),'Type',2); 
    dxqonDiub3=gradx_m(averyNtoD(q3)); 
    u3=m_u3-gradx_m(q3); 
    v3=m_v3-gradyNtoD(q3); v3(1,:)=Ov; v3(end,:)=Ov;

    F1hatB = fft(dct(F1(uB,vB,phiB,rhoB),'Type',2)')';
    F2hatB = fft(dct(F2(uB,vB,phiB,rhoB),'Type',2)')';
    phihat3 = varphi0dtL1.*phihat + dt*(varphi31dtL1.*F1hatO + varphi32dtL1.*F1hatA + varphi33dtL1.*F1hatB);   
    rhohat3 = varphi0dtL2.*rhohat + dt*(varphi31dtL2.*F2hatO + varphi32dtL2.*F2hatA + varphi33dtL2.*F2hatB); 
    phi3 = idct(real(ifft(phihat3')'),'Type',2); 
    rho3 = idct(real(ifft(rhohat3')'),'Type',2); 

    phi=phi3; rho=rho3;
    m_u=m_u3; m_v=m_v3; 
    u=u3; v=v3;
    dxqonDiub=dxqonDiub3;
    q=q3;

    if abs(t-tn(n))<0.1*dt %comparing float point numbers dircetly yields error
        phitn=[phitn phi];
        rhotn=[rhotn rho];
        utn=[utn u];
        vtn=[vtn v];
        if n<size(tn,2)
            n=n+1;
        end
    end

    tvec = [tvec t];
    maxmin=[max(max(phi)) max(max(rho));min(min(phi)) min(min(rho))]
    G = get_energy(u,v,phi,rho)
    Energy=[Energy;G];
    phimax = [phimax max(max(phi))];
    rhomax = [rhomax max(max(rho))];
    phimin = [phimin min(min(phi))];
    rhomin = [rhomin min(min(rho))];
end

qq=[];
for s = 1 : 4
    t = t + dt %from t=T-3dt to t=T
    m_vin=m_v(2:end-1,:);
    m_uhat=fft(dstIIbyfft(m_u,Ome,Ny)')';
    m_vinhat=fft(dstIbyfft(m_vin,Nx,Ny)')';
    gqupO=gup+dxqonDiub(1,:);
    gqlowO=glow+dxqonDiub(end,:);
    gqO=[gqupO; Oin; gqlowO];
    F3uhatO = fft(dstIIbyfft(nu*2*gqO/hy/hy+F3u(gup,glow,u,v,phi,rho),Ome,Ny)')';
    F3vhatO = fft(dstIbyfft(F3v(u,v,phi,rho),Nx,Ny)')';
    m_uhatA = varphi0dtL3u .*m_uhat + dt*varphi1dtL3u .* F3uhatO ;  
    m_vinhatA = varphi0dtL3v .*m_vinhat + dt*varphi1dtL3v .*F3vhatO ; 
    m_uA = idstIIbyfft(real(ifft(m_uhatA')'),Omebar,Nx,Ny); 
    m_vinA = idstIbyfft(real(ifft(m_vinhatA')'),Nx,Ny); 
    m_vA=[Ov; m_vinA ;Ov];
    bA=gradx(m_uA)+gradyDtoN(m_vA);
    bAhat=fft(dct(bA,'Type',2)')'; bAhat(1,1)=0;
    qA=idct(real(ifft((bAhat./lambdaa)')'),'Type',2);  
    dxqonDiubA=gradx_m(averyNtoD(qA)); 
    uA=m_uA-gradx_m(qA); 
    vA=m_vA-gradyNtoD(qA); vA(1,:)=Ov; vA(end,:)=Ov;

    phihat=fft(dct(phi,'Type',2)')';
    rhohat=fft(dct(rho,'Type',2)')';
    F1hatO = fft(dct(F1(u,v,phi,rho),'Type',2)')';
    F2hatO = fft(dct(F2(u,v,phi,rho),'Type',2)')';
    phihatA = varphi0dtL1.*phihat + dt*varphi1dtL1.*F1hatO ;  
    rhohatA = varphi0dtL2.*rhohat + dt*varphi1dtL2.*F2hatO ; 
    phiA = idct(real(ifft(phihatA')'),'Type',2); 
    rhoA = idct(real(ifft(rhohatA')'),'Type',2); 
    
    gqupA=gup+dxqonDiubA(1,:);
    gqlowA=glow+dxqonDiubA(end,:);
    gqA=[gqupA; Oin; gqlowA];
    F3uhatA = fft(dstIIbyfft(nu*2*gqA/hy/hy+F3u(gup,glow,uA,vA,phiA,rhoA),Ome,Ny)')';
    F3vhatA = fft(dstIbyfft(F3v(uA,vA,phiA,rhoA),Nx,Ny)')';
    m_uhatB = varphi023dtL3u.*m_uhat + dt*(varphi2123dtL3u.*F3uhatO + varphi2223dtL3u.*F3uhatA) ;  
    m_vinhatB = varphi023dtL3v.*m_vinhat + dt*(varphi2123dtL3v.*F3vhatO + varphi2223dtL3v.*F3vhatA) ;  
    m_uB = idstIIbyfft(real(ifft(m_uhatB')'),Omebar,Nx,Ny); 
    m_vinB = idstIbyfft(real(ifft(m_vinhatB')'),Nx,Ny); 
    m_vB=[Ov; m_vinB ;Ov];
    bB=gradx(m_uB)+gradyDtoN(m_vB);
    bBhat=fft(dct(bB,'Type',2)')'; bBhat(1,1)=0;
    qB=idct(real(ifft((bBhat./lambdaa)')'),'Type',2);  
    dxqonDiubB=gradx_m(averyNtoD(qB)); 
    uB=m_uB-gradx_m(qB); 
    vB=m_vB-gradyNtoD(qB); vB(1,:)=Ov; vB(end,:)=Ov;
    

    F1hatA = fft(dct(F1(uA,vA,phiA,rhoA),'Type',2)')';
    F2hatA = fft(dct(F2(uA,vA,phiA,rhoA),'Type',2)')';
    phihatB = varphi023dtL1.*phihat + dt*(varphi2123dtL1.*F1hatO + varphi2223dtL1.*F1hatA);  
    rhohatB = varphi023dtL2.*rhohat + dt*(varphi2123dtL2.*F2hatO + varphi2223dtL2.*F2hatA);  
    phiB = idct(real(ifft(phihatB')'),'Type',2); 
    rhoB = idct(real(ifft(rhohatB')'),'Type',2); 

    gqupB=gup+dxqonDiubB(1,:);
    gqlowB=glow+dxqonDiubB(end,:);
    gqB=[gqupB; Oin; gqlowB];
    F3uhatB = fft(dstIIbyfft(nu*2*gqB/hy/hy+F3u(gup,glow,uB,vB,phiB,rhoB),Ome,Ny)')';
    F3vhatB = fft(dstIbyfft(F3v(uB,vB,phiB,rhoB),Nx,Ny)')';
    m_uhat3 = varphi0dtL3u.*m_uhat + dt*(varphi31dtL3u .* F3uhatO + varphi32dtL3u .*F3uhatA + varphi33dtL3u .*F3uhatB);  
    m_vinhat3 = varphi0dtL3v.*m_vinhat + dt*(varphi31dtL3v .*F3vhatO + varphi32dtL3v .*F3vhatA + varphi33dtL3v .*F3vhatB) ; 
    m_u3 = idstIIbyfft(real(ifft(m_uhat3')'),Omebar,Nx,Ny); 
    m_vin3 = idstIbyfft(real(ifft(m_vinhat3')'),Nx,Ny); 
    m_v3=[Ov; m_vin3 ;Ov];
    b3=gradx(m_u3)+gradyDtoN(m_v3);
    b3hat=fft(dct(b3,'Type',2)')'; b3hat(1,1)=0;
    q3=idct(real(ifft((b3hat./lambdaa)')'),'Type',2);  
    dxqonDiub3=gradx_m(averyNtoD(q3)); 
    u3=m_u3-gradx_m(q3); 
    v3=m_v3-gradyNtoD(q3); v3(1,:)=Ov; v3(end,:)=Ov;

    F1hatB = fft(dct(F1(uB,vB,phiB,rhoB),'Type',2)')';
    F2hatB = fft(dct(F2(uB,vB,phiB,rhoB),'Type',2)')';
    phihat3 = varphi0dtL1.*phihat + dt*(varphi31dtL1.*F1hatO + varphi32dtL1.*F1hatA + varphi33dtL1.*F1hatB);   
    rhohat3 = varphi0dtL2.*rhohat + dt*(varphi31dtL2.*F2hatO + varphi32dtL2.*F2hatA + varphi33dtL2.*F2hatB); 
    phi3 = idct(real(ifft(phihat3')'),'Type',2); 
    rho3 = idct(real(ifft(rhohat3')'),'Type',2); 

    phi=phi3; rho=rho3;
    m_u=m_u3; m_v=m_v3; 
    u=u3; v=v3;
    dxqonDiub=dxqonDiub3;
    q=q3; qq=[qq q];
    
    if abs(t-tn(n))<0.1*dt %comparing float point numbers dircetly yields error
        phitn=[phitn phi];
        rhotn=[rhotn rho];
        utn=[utn u];
        vtn=[vtn v];
        if n<size(tn,2)
            n=n+1;
        end
    end
    
    tvec = [tvec t];
    G = get_energy(u,v,phi,rho)
    Energy=[Energy;G];
    phimax = [phimax max(max(phi))];
    rhomax = [rhomax max(max(rho))];
    phimin = [phimin min(min(phi))];
    rhomin = [rhomin min(min(rho))];
end
save('qq.mat','qq')
qq1=qq(:,1:Nx); qq2=qq(:,Nx+1:2*Nx); qq3=qq(:,2*Nx+1:3*Nx); qq4=qq(:,3*Nx+1:4*Nx); 
p=((11/6)*qq4-3*qq3+1.5*qq2-(1/3)*qq1)/dt-nu*Laplace2d_PNII(qq4);
save('pressure_T.mat','p')

Energy
save('energy.mat','Energy')
save('max_min.mat','phimax','rhomax','phimin','rhomin')
save('phase_seperation.mat','phitn','rhotn')
save('velocity.mat','utn','vtn')

