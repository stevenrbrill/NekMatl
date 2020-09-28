clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navier-Stokes Solver Demo (2D Channel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact; format longe; clear all
Re = 1; Pr=0.8; Pe=Re*Pr; 

%N=16; E=5; N1=N+1; nL=N1*N1*E;  % 16th order
N=10;E=5; N1=N+1; nL=N1*N1*E; % 10th order 

Q=makeq(E,N);
R=maker(Q,E,N);
[X,Y]=make_geom_channel(E,N);      % Geometry in local form
[G,J,ML,RX]=make_coef(X,Y);
[Ah,Bh,Ch,Dh,z,w]=semhat(N);

[n,nb]=size(R); 
nL=N1*N1*E;


Ab=spalloc(nb,nb,nb*N);    %  Generate Abar
for j=1:nb;
    x=zeros(nb,1); 
    x(j)=1; 
    Ab(:,j)=abar(x,Dh,G,Q); 
end;
A=R*Ab*R';
Bb=reshape(ML,nL,1); 
Bb=diag(Bb); 
Bb=sparse(Bb); 
Bb=Q'*Bb*Q; 
Ma=R*Bb*R';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CFL=0.3; dxmin=pi*(z(N1)-z(N))/(2*E); % Get min dx for CFL constraint
Tfinal=150; dt=CFL*dxmin; nstep=ceil(Tfinal/dt); dt=Tfinal/nstep; nstep


u=ones(size(ML)); u1=u; u2=u; u3=u; fx3=u; fx2=u; fx1=u;
v=0*ML; v1=v; v2=v; v3=v; fy3=v; fy2=v; fy1=v;

F = ones(size(ML));
Fb=reshape(ML.*F,nL,1); 
Fb=Bb\(Q'*Fb);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for step=1:nstep; time=step*dt;

    if step==1; b0=1.0;    b= [ -1 0 0 ]';       a=[ 1  0 0 ]'; end;
    if step==2; b0=1.5;    b=([ -4 1 0 ]')./2;   a=[ 2 -1 0 ]'; end;
    if step==3; b0=11./6.; b=([ -18 9 -2 ]')./6; a=[ 3 -3 1 ]'; end;
    if step<=3; H=(Ma+ A*dt/(b0*Re)); [LH,UH]=lu(H); b0i=1./b0; end; % Viscous op
    if step<=3; Hbar=(Bb+ Ab*dt/(b0*Re)); end; % Viscous op

%   Nonlinear step - unassembled, not multiplied by mass matrix

    fx1 = -convl(u,RX,Dh,u,v) + F; % du = Cu  (without mass matrix)
    fy1 = -convl(v,RX,Dh,u,v); % dv = Cv

    rx  = a(1)*fx1+a(2)*fx2+a(3)*fx3; % kth-order extrapolation
    ry  = a(1)*fy1+a(2)*fy2+a(3)*fy3;

    fx3=fx2; fx2=fx1; 
    fy3=fy2; fy2=fy1; 

    rx  = dt*rx - (b(1)*u+b(2)*u2+b(3)*u3); u3=u2; u2=u; % Add BDF terms
    ry  = dt*ry - (b(1)*v+b(2)*v2+b(3)*v3); v3=v2; v2=v; %     and

    ut  = b0i*rx; 
    vt = b0i*ry; 

%   uL=ut; vL=vt;
    [uL,vL,pr]=pressure_project(ut,vt,Ab,Q,ML,RX,Dh); % Div-free velocity
    pr = (b0/dt)*pr;

    %   Set RHS.                 %Viscous update. %  Convert to local form.
    u=R*(Q'*reshape(ML.*uL,nL,1));
    u=UH\(LH\u);u=Q*(R'*u);
    u=reshape(u,N1,N1,E);
    
    v=R*(Q'*reshape(ML.*vL,nL,1));
    v=UH\(LH\v);v=Q*(R'*v);
    v=reshape(v,N1,N1,E);

    if mod(step,20)==0;
        figure(1)
        plotit(u,X,Y); 
        figure(2)
        plot(u(1,:,1),Y(1,:,1))
        xlim([0,1])
        ylim([-1,1])
        xlabel('u')
        ylabel('y')
        pause(.1); 
        [step time glmax(u)] 
    end;

end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
