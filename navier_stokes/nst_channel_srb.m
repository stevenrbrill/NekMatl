clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navier-Stokes Solver Demo (2D Channel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact; format short; clear all
Re = 200; Pr=0.8; Pe=Re*Pr; 

%N=16; E=5; N1=N+1; nL=N1*N1*E;  % 16th order
N=11; % polynomial order  
Ex=5; % Number of elements in x
Ey=5; % Number of elements in y
CFL=0.1;
u_ic = Re/2;
pert = 0.1;
f_ic = @(x,y) u_ic*(1-y.^2);


E=Ex*Ey; % Total number of elements
N1=N+1; 

% Enrichment information
psi = @(x,y) y.^2 + 0.*x;
gpsi = {@(x,y) 0.*y + 0.*x, @(x,y) 2.*y + 0.*x};
hpsi = {@(x,y) 0.*y + 0.*x, @(x,y) 2 + 0.*y + 0.*x};

disp("Generating Matrices")
Q=makeq(Ex,Ey,N); % Global continuity
R=maker(Q,Ex,N); % Restriction matrix, applies Dirichlet conditions
[X,Y]=make_geom_channel(Ex,Ey,N);      % Geometry in local form
[G,J,ML,RX]=make_coef(X,Y);
[Ah,Bh,Ch,Dh,z,w]=semhat(N);

[n,nb]=size(R); 
nL=N1*N1*E;

Ab=spalloc(nb,nb,nb*N);    %  Generate Abar
for j=1:nb;
    x=zeros(nb,1); 
    x(j)=1; 
    Ab(:,j)=abar(x,Dh,G,Q);  % assembled Neumann Operator
end;
A=R*Ab*R'; % Full stiffness matrix
Bb=reshape(ML,nL,1); % Form Mass matrix
Bb=diag(Bb); 
Bb=sparse(Bb); 
Bb=Q'*Bb*Q; % Assembling mass matrix
Ma=R*Bb*R'; % Full mass matrix

% Assemble enrichment matrices
disp("Computing enrichment")
[Mp,Sp,T1,T2,z_en,w_en] = enrich_mats(X,Y,E,N,psi,gpsi,hpsi);
nb = (N+1)^2; 
for k=1:2
    Mp_all{k} = zeros(nb*E,nb*E);
    Sp_all{k} = zeros(nb*E,nb*E);
    T1_all{k} = zeros(nb*E,1);
    T2_all{k} = zeros(nb*E,1);
    for i = 1:E
        Mp_all{k}((i-1)*nb+1:i*nb,(i-1)*nb+1:i*nb) = Mp{k}(:,:,i);
        Sp_all{k}((i-1)*nb+1:i*nb,(i-1)*nb+1:i*nb) = Sp{k}(:,:,i);
        T1_all{k}((i-1)*nb+1:i*nb) = T1{k}(:,i);
        T2_all{k}((i-1)*nb+1:i*nb) = T2{k}(:,i);
    end
    % TODO: Remove small entries
    Mp_all{k} = sparse(Mp_all{k});
    Mp_all{k} = Q'*Mp_all{k}*Q;
    Sp_all{k} = sparse(Sp_all{k});
    Sp_all{k} = Q'*Sp_all{k}*Q;
    
    T1_all{k} = sparse(T1_all{k});
    T1_all{k} = Q'*T1_all{k};
    T2_all{k} = sparse(T2_all{k});
    T2_all{k} = Q'*T2_all{k};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dxmin=pi*(z(N1)-z(N))/(2*Ex); % Get min dx for CFL constraint
Tfinal=150; 
dt=CFL*dxmin/u_ic; 
nstep=ceil(Tfinal/dt); 
dt=Tfinal/nstep; 
% Print information
nstep
dt
yp1 = (1+Y(1,2))*Re

disp("Computing inverse")
Ai=pinv(full(Ab)); 

disp("Setting Initial Conditions")
% Initial conditions
u=u_ic*ones(size(ML)); 
v=0*ML;

for e = 1:E
    for i = 1:N+1
        for j = 1:N+1
            u(i,j,e) = f_ic(X(i,j,e),Y(i,j,e))+pert*rand()*u_ic;
            v(i,j,e) = v(i,j,e)+pert*rand()*u_ic;
        end
    end
end


u1=u; u2=u; u3=u; fx3=u; fx2=u; fx1=u;
v1=v; v2=v; v3=v; fy3=v; fy2=v; fy1=v;

F = ones(size(ML));
Fb=reshape(ML.*F,nL,1); 
Fb=Bb\(Q'*Fb);  


disp("Timestepping")
plot1 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for step=1:nstep; time=step*dt;

    if step==1; b0=1.0;    b= [ -1 0 0 ]';       a=[ 1  0 0 ]'; end
    if step==2; b0=1.5;    b=([ -4 1 0 ]')./2;   a=[ 2 -1 0 ]'; end
    if step==3; b0=11./6.; b=([ -18 9 -2 ]')./6; a=[ 3 -3 1 ]'; end
    if step<=3 
        H=(Ma+ A*dt/(b0*Re)); 
        [LH,UH]=lu(H); 
        b0i=1./b0; 
    end % Viscous op
    if step<=3
        Hbar=(Bb+ Ab*dt/(b0*Re)); 
    end % Viscous op

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
    [uL,vL,pr]=pressure_project(ut,vt,Ai,Q,ML,RX,Dh); % Div-free velocity
    pr = (b0/dt)*pr;

    %   Set RHS.                 %Viscous update. %  Convert to local form.
    u=R*(Q'*reshape(ML.*uL,nL,1));
    u=UH\(LH\u);u=Q*(R'*u);
    u=reshape(u,N1,N1,E);
    
    v=R*(Q'*reshape(ML.*vL,nL,1));
    v=UH\(LH\v);v=Q*(R'*v);
    v=reshape(v,N1,N1,E);

    if mod(step,1000)==0
        if plot1
            figure(1);
        else
            figure(3);
        end
        plotit(u,X,Y); 
        if plot1
            figure(2);
            plot1 = 0;
        else
            figure(4);
            plot1 = 1;
        end
        for i = 1:Ey
            plot(u(1,:,1+Ex*(i-1)),Y(1,:,1+Ex*(i-1)),'b-o')
            hold on
        end
        hold off
%         xlim([0,1])
        ylim([-1,1])
        xlabel('u')
        ylabel('y')
%         pause(.1); 
        [time glmax(u)] 
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
