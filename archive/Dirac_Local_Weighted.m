clear all
N=60;
NL0=sqrt(N);
Omega1=9;
omegaedge=4*randn(N,N);
for i2=1:NL0^2, 
    x1=floor((i2-1)/NL0);
    y1=mod((i2-1),NL0);
    x(i2)=x1;
    y(i2)=y1;
    if(y1+1<NL0)
    i=i2;
    j=x1*NL0+y1+2;
    a(i,j)=1;
    a(j,i)=1;    
    if(x1==0),
        omegaedge(i,j)=Omega1;
        omegaedge(j,i)=-Omega1;
    end
    if (x1==NL0-1),
        omegaedge(i,j)=-Omega1;
        omegaedge(j,i)=Omega1;
    end
    end
    if(x1+1<NL0)
    j=(x1+1)*NL0+y1+1;
    i=i2;
    a(i,j)=1;
    a(j,i)=1;
    if((y1==NL0-1))
        omegaedge(i,j)=Omega1;
        omegaedge(j,i)=-Omega1;
    end
    if (y1==0),
            omegaedge(i,j)=-Omega1;
        omegaedge(j,i)=Omega1;
    end
    end
end
k=sum(a,2);
Omega0=9;
tau0=1;
z=0;
w = Omega0 + randn(N,1)/tau0; %Gaussian node frequencies
%w = w - mean(w);%enforce mean 0
epsilon=0.01;
[I,J,V]=find(triu(a,1));
L=numel(V);
for n=1:L,
    B(I(n),n)=-1;
    B(J(n),n)=1;
end
w=rand(L,1);
W=diag(w);
strength=abs(B)*W;
D=sparse(N+L,N+L);
D([1:N],[N+1:N+L])=diag(diag(strength.^(-1)))*B*W/2;
D([N+1:N+L],[1:N])=transpose(B);
gamma=sparse(N+L,N+L);
gamma([1:N],[1:N])=eye(N);
gamma([N+1:N+L],[N+1:N+L])=-eye(L);

D2=z*gamma*D*D;
%D2=zeros(N+L,N+L);

for n=1:numel(V),
    wedge(n)=omegaedge(I(n),J(n));
end

Omega([1:N])=W;
Omega([N+1:N+L])=wedge;
%Omega([N+1:N+L])=10+20*randn(L,1);
Omega=Omega';
Psi0=2*pi*rand(N+L,1);


func = @Dirac_metric; %sets the model between ddt_I and ddt_II


Psi=Psi0;


Tmax = 10; % max time to integrate
dt = 0.01; %time step
sigma=12;
dsigma=0.03;
vidfile = VideoWriter('squaremovie_z_0_sigma12_omega_4_9.mp4','MPEG-4');
vidfile.FrameRate = 20.96;
open(vidfile);
nframe=1;

for t = 0:dt:Tmax %time loop

    %RK4 method
    [k1{1},w1{1}] = func(Psi, w,sigma,B,Omega); %Eqs (32) in notes 
    [k2{1},w2{1}] = func(Psi+0.5*dt*k1{1},w+0.5*dt*w1{1},sigma,B,Omega);
    [k3{1},w3{1}] = func(Psi+0.5*dt*k2{1},w+0.5*dt*w2{1},sigma,B,Omega);
    [k4{1},w{4}] = func(Psi+dt*k3{1},w+dt*w3{1}, sigma,B,Omega);
    
    %Final update
    Psi = Psi + (dt/6)*(k1{1}+2*k2{1}+2*k3{1}+k4{1});
    w=w+(dt/6)*(w1{1}+2*w2{1}+2*w3{1}+w4{1});

    theta=Psi([1:N]);
    psi=(B*Psi([N+1:N+L]))./k;
    
    phi=Psi([N+1:N+L]);
    
    %Calculates averages over the last fifth of the time series   
    

if(t>7)
figure(1)
cc=colormap(hsv(100));
mc=-1;
Mac=1;
for n=1:L,
    c2=(floor((real(cos(phi(n))-mc)*(90)/(Mac-mc))))+1;
    plot([x(I(n)),x(J(n))],[y(I(n)),y(J(n))],'Color',cc(c2,:),'LineWidth',3)
    hold on
end
for i=1:N,
    c2=(floor((real(cos(theta(i))-mc)*(90)/(Mac-mc))))+1;
    plot(x(i),y(i),'o','Color',cc(c2,:),'MarkerFaceColor',cc(c2,:))
    hold on
end
axis off
    
F(nframe) = getframe(gcf); 

    writeVideo(vidfile,F(nframe));
    nframe=nframe+1;
    hold off
end

end

close(vidfile)
