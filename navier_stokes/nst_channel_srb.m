% clc
clear all
close all


set(0,'DefaultAxesFontSize',14)
set(0,'defaultLineLineWidth',2.5)
% set(0,'DefaultFigureWindowStyle','docked')


linestyles = {'k-','b-','r-','g-','c-','m-','k--','b--','r--','g--','c--','m--','k-.','b-.','r-.','g-.','c-.','m-.','k:','b:','r:','g:','c:','m:'};
pointstyles = {'ko','bo','ro','go','co','mo','k^','b^','r^','g^','c^','m^','ks','bs','rs','gs','cs','ms','k*','b*','r*','g*','c*','m*','k+','b+','r+','g+','c+','m+'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navier-Stokes Solver Demo (2D Channel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact;
format short; 
Re = 1; Pr=0.8; Pe=Re*Pr; 

%N=16; E=5; N1=N+1; nL=N1*N1*E;  % 16th order
N=5; % polynomial order  
Ex=1; % Number of elements in x
Ey=3; % Number of elements in y
CFL=0.1;
u_ic = Re;
pert = 0.0;
f_ic = @(x,y) u_ic*(1-y.^2)/2;

%% Enrichment information
en_on = 1;
N_en_y = 1; 
psi = {@(x,y) (0.5*(1 - y.^2) + 0.*x), @(x,y) 0.*y + 0.*x};
gpsi = {@(x,y) 0.*y + 0.*x, @(x,y) (-1.*y + 0.*x), ...
        @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
hpsi = {@(x,y) 0.*y + 0.*x, @(x,y) (-1 - 0.*y + 0.*x), ...
        @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
    
% No enrichment on one half
% psi = {@(x,y) (y<=0).*(0.5*(1 - y.^2)) + (y>0).*0.5 + 0.*x, @(x,y) 0.*y + 0.*x};
% gpsi = {@(x,y) 0.*y + 0.*x, @(x,y) (y<=0).*(-1.*y)+ (y>0).*0 + 0.*x, ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
% hpsi = {@(x,y) 0.*y + 0.*x, @(x,y) (y<=0).*(-1.)+ (y>0).*0 + 0.*y + 0.*x, ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};

% Blended enrichment
% En - Blend - normal .. - Blend - En
% N_en_form = 1;
% Ly = 2/Ey;
% ye = 1-Ly*N_en_form; % end of enrichment
% yb = 1-Ly*(N_en_form+1); % end of blending
% 
% blend = @(x,y) 0.*x + ...
%     (abs(y) >= (ye-eps)).*1 + ...
%     ((-(ye-eps) < y).*(y <= -(yb-eps))).*(y-(-yb))/((-ye)-(-yb)) + ...
%     + y.*0 + ...
%     (((ye-eps) > y).*(y >= (yb-eps))).*(y-(yb))/((ye)-(yb));
% dblend = @(x,y) 0.*x + ...
%     (abs(y) >= (ye-eps)).*0 + ...
%     ((-(ye-eps) < y).*(y <= -(yb-eps))).*1./((-ye)-(-yb)) + ...
%     + y.*0 + ...
%     (((ye-eps) > y).*(y >= (yb-eps))).*1./((ye)-(yb));
% hblend = @(x,y) 0.*x + ...
%     (abs(y) >= (ye-eps)).*0 + ...
%     ((-(ye-eps) < y).*(y <= -(yb-eps))).*0 + ...
%     + y.*0 + ...
%     (((ye-eps) > y).*(y >= (yb-eps))).*0;
% 
% % c1 = 1./(-ye^2-yb^2+2*ye*yb);
% % c2 = -2*ye*c1;
% % c3 = -(yb^2-2*yb*ye)*c1;
% % c1b = 1./(-(-ye)^2-(-yb)^2+2*(-ye)*(-yb));
% % c2b = -2*(-ye)*c1;
% % c3b = -((-yb)^2-2*(-yb)*(-ye))*c1;
% % blend = @(x,y) 0.*x + ...
% %     (abs(y) >= ye-eps).*1 + ...
% %     ((-(ye-eps) < y).*(y <= -(yb-eps))).*(c1b*y.^2+c2b*y+c3b) + ...
% %     + y.*0 + ...
% %     (((ye-eps) > y).*(y >= (yb-eps))).*(c1*y.^2+c2*y+c3);
% % 
% % dblend = @(x,y) 0.*x + ...
% %     (abs(y) >= (ye-eps)).*0 + ...
% %     ((-(ye-eps) < y).*(y <= -(yb-eps))).*(2*c1b*y+c2b) + ...
% %     + y.*0 + ...
% %     (((ye-eps) > y).*(y >= (yb-eps))).*(2*c1*y+c2);
% % 
% % hblend = @(x,y) 0.*x + ...
% %     (abs(y) >= (ye-eps)).*0 + ...
% %     ((-(ye-eps) < y).*(y <= -(yb-eps))).*2*c1b + ...
% %     + y.*0 + ...
% %     (((ye-eps) > y).*(y >= (yb-eps))).*2*c1;
% 
% en = @(x,y)(0.5*(1 - y.^2)) + 0.*x;
% den = @(x,y) (-1.*y + 0.*x);
% hen = @(x,y) (-1 - 0.*y + 0.*x);
% 
% psi = {@(x,y) blend(x,y).*en(x,y), @(x,y) 0.*y + 0.*x};
% gpsi = {@(x,y) 0.*y + 0.*x, ...
%         @(x,y) blend(x,y).*den(x,y) + dblend(x,y).*en(x,y), ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
% hpsi = {@(x,y) 0.*y + 0.*x, ...
%         @(x,y) blend(x,y).*hen(x,y) + 2*dblend(x,y).*den(x,y) + hblend(x,y).*en(x,y), ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};

% Flat enrichment in middle
% Middle is higher than normal enrichment
% en_loc = 2/Ey+0.0001;
% psi = {@(x,y) ((1-abs(y))<=en_loc).*(0.5*(1 - y.^2) + 0.*x)+((1-abs(y))>en_loc).*0.5*(1-(1-2/Ey)^2), @(x,y) 0.*y + 0.*x};
% gpsi = {@(x,y) 0.*y + 0.*x, @(x,y) ((1-abs(y))<=en_loc).*(-1.*y + 0.*x), ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
% hpsi = {@(x,y) 0.*y + 0.*x, @(x,y) ((1-abs(y))<=en_loc).*(-1 - 0.*y + 0.*x), ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};

% % 0 enrichment in middle  
% Middle is slightly below normal enrichment
% en_loc = 2/Ey+0.000000001;
% psi = {@(x,y) ((1-abs(y))<=en_loc).*(0.5*(1 - y.^2) + 0.*x), @(x,y) 0.*y + 0.*x};
% gpsi = {@(x,y) 0.*y + 0.*x, @(x,y) ((1-abs(y))<=en_loc).*(-1.*y + 0.*x), ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
% hpsi = {@(x,y) 0.*y + 0.*x, @(x,y) ((1-abs(y))<=en_loc).*(-1 - 0.*y + 0.*x), ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};

% % 0 enrichment in middle not full element
% Middle is slightly below normal enrichment
% [X,Y]=make_geom_channel(Ex,Ey,N);      % Geometry in local form
% en_loc = 1-abs(Y(1,end-1,1))+0.000000001;
% psi = {@(x,y) ((1-abs(y))<=en_loc).*(0.5*(1 - y.^2) + 0.*x), @(x,y) 0.*y + 0.*x};
% gpsi = {@(x,y) 0.*y + 0.*x, @(x,y) ((1-abs(y))<=en_loc).*(-1.*y + 0.*x), ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
% hpsi = {@(x,y) 0.*y + 0.*x, @(x,y) ((1-abs(y))<=en_loc).*(-1 - 0.*y + 0.*x), ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};

% 0 enrichment
% psi = {@(x,y) 0+0.*y + 0.*x, @(x,y) 0+0.*y + 0.*x;};
% gpsi = {@(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
% hpsi = {@(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};

%% Plot psi
figure
ys_plot = linspace(-1,1,1000);
plot(psi{1}(0,ys_plot),ys_plot)
xlabel('\psi')
ylabel('y')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
    yticks(linspace(-1,1,Ey+1));

psi_p = psi{1}(0,N_en_y*2/Ey-1);

%% Begin Solve
E=Ex*Ey; % Total number of elements
N1=N+1;
disp("Generating Matrices")
% if en_on
%     Q=makeq_en(Ex,Ey,N,N_en_y);
% else
    Q=makeq(Ex,Ey,N); % Global continuity
% end
R=maker(Q,Ex,N); % Restriction matrix, applies Dirichlet conditions
[X,Y]=make_geom_channel(Ex,Ey,N);      % Geometry in local form
en_b_nodes = get_en_bound_nodes(Ex,Ey,N,N_en_y);

Ys = zeros(Ey*length(Y(1,:,1)),1);
for i = 1:Ey
    Ys((i-1)*(N+1)+1:(i)*(N+1)) = Y(1,:,(i-1)*Ex+1);
end
        
[G,J,ML,RX]=make_coef(X,Y);
[Ah,Bh,Ch,Dh,z,w]=semhat(N);

[n,nb]=size(R); 
nL=N1*N1*E;

Ab=spalloc(nL,nL,nL*N);    %  Generate Abar
for j=1:nL
    x=zeros(nL,1); 
    x(j)=1; 
    Ab(:,j)=abar_c(x,Dh,G,Q);  % assembled Neumann Operator
end
A_c=R*Q'*apply_en_cont(Ab,en_b_nodes,psi_p);
S_check=full(Ab);
Ab = Q'*Ab*Q;
A=R*Ab*R'; % Full stiffness matrix
Bb=reshape(ML,nL,1); % Form Mass matrix
Bb=diag(Bb); 
Bb=sparse(Bb); 
Ma=R*Q'*Bb*Q*R'; % Assembling mass matrix  % Full mass matrix
M_c = R*Q'*apply_en_cont(Bb,en_b_nodes,psi_p);


%%
M_check=full(Bb);

Ma_uv_check = zeros(2*nL);
Ma_uv_check(1:nL,1:nL) = M_check;
Ma_uv_check(nL+1:2*nL,nL+1:2*nL) = M_check;
A_uv_check = zeros(2*nL);
A_uv_check(1:nL,1:nL) = S_check;
A_uv_check(nL+1:2*nL,nL+1:2*nL) = S_check;

M_q = full(Ma);
S_q = full(A);


%%


% Assemble overall system matrices
nn = length(Ma); % number of nodes in the domain
Ma_uv = zeros(2*nn);
Ma_uv(1:nn,1:nn) = Ma;
Ma_uv(nn+1:2*nn,nn+1:2*nn) = Ma;
Ma_uv = sparse(Ma_uv);
A_uv = zeros(2*nn);
A_uv(1:nn,1:nn) = A;
A_uv(nn+1:2*nn,nn+1:2*nn) = A;
A_uv = sparse(A_uv);
ML_uv = zeros(2*(N+1),2*(N+1),E);
ML_uv(1:N+1,1:N+1,1:E) = ML;
ML_uv(N+2:2*(N+1),N+2:2*(N+1),1:E) = ML;


%% Assemble enrichment matrices
psi_xy{1} = psi{1}(X,Y);
psi_xy{2} = psi{2}(X,Y);
if en_on
    disp("Computing enrichment")
    [Mp,Sp,T1,T2,z_en,w_en] = enrich_mats(X,Y,E,N,psi,gpsi,hpsi);
    nb = N1*N1;
        
    for k=1:2
        Mp_all{k} = zeros(nb*E,nb*E);
        Sp_all{k} = zeros(nb*E,nb*E);
%         T1_all{k} = zeros(nb*E,1);
%         T2_all{k} = zeros(nb*E,1);
        T1_all{k} = zeros(N+1,N+1,E);
        T2_all{k} = zeros(N+1,N+1,E);
        
        T1_rs{k} = reshape(T1{k},[N+1,N+1,E]);
        T2_rs{k} = reshape(T2{k},[N+1,N+1,E]);
        
        % TODO: Assemble differently for different elements
        for iy = 1:Ey
            for ix = 1:Ex
                i = (iy-1)*Ex+ix;
                if ((iy <= N_en_y) || (iy > Ey-N_en_y))
                    Mp_all{k}((i-1)*nb+1:i*nb,(i-1)*nb+1:i*nb) = Mp{k}(:,:,i);
                    Sp_all{k}((i-1)*nb+1:i*nb,(i-1)*nb+1:i*nb) = Sp{k}(:,:,i);
%                     T1_all{k}((i-1)*nb+1:i*nb) = T1{k}(:,i);
%                     T2_all{k}((i-1)*nb+1:i*nb) = T2{k}(:,i);
%                     T1{k}(:,i) = T1{k}(:,i);
%                     T2{k}(:,i) = T2{k}(:,i);
                    T1_all{k}(:,:,i) = T1_rs{k}(:,:,i); 
                    T2_all{k}(:,:,i) = T2_rs{k}(:,:,i);
                end
            end
        end
        %%
        if k == 1
            Mp_check_1 = full(Mp_all{k});
            Sp_check_1 = full(Sp_all{k});
            T1_check_1 = full(T1_all{k});
            T2_check_1 = full(T2_all{k});
            Mp_check_2 = zeros(size(Mp_all{k}));
            Sp_check_2 = zeros(size(Mp_all{k}));
            T1_check_2 = zeros(size(Mp_all{k}));
            T2_check_2 = zeros(size(Mp_all{k}));
        elseif k == 2
            Mp_check_2 = full(Mp_all{k});
            Sp_check_2 = full(Sp_all{k});
            T1_check_2 = full(T1_all{k});
            T2_check_2 = full(T2_all{k});
        end
        %%
        Mp_all{k} = sparse(Mp_all{k});        
        Mp_all_c{k} = R*Q'*apply_en_cont(Mp_all{k},en_b_nodes,psi_p);
        Mp_all{k} = R*Q'*Mp_all{k}*Q*R';
        Sp_all{k} = sparse(Sp_all{k});
        Sp_all_c{k} = R*Q'*apply_en_cont(Sp_all{k},en_b_nodes,psi_p);
        Sp_all{k} = R*Q'*Sp_all{k}*Q*R';
        
%         T1_all{k} = sparse(T1_all{k});
%         T2_all{k} = sparse(T2_all{k});
%         T2_all{k} = Q'*T2_all{k};
   
    end
    
    Mp_uv = zeros(2*nn);
    Sp_uv = zeros(2*nn);
    Mp_uv(1:nn,1:nn) = Mp_all{1};
    Mp_uv(1:nn,nn+1:2*nn) = Mp_all{2};
    Sp_uv(1:nn,1:nn) = Sp_all{1};
    Sp_uv(nn+1:2*nn,nn+1:2*nn) = Sp_all{1};
    
    %%
    Mp_uv_check = zeros(2*nL);
    Mp_uv_check(1:nL,1:nL) = Mp_check_1;
    Mp_uv_check(1:nL,nL+1:2*nL) = Mp_check_2;
    Sp_uv_check(1:nL,1:nL) = Sp_check_1;
    Sp_uv_check(nL+1:2*nL,nL+1:2*nL) = Sp_check_1;
    Mp_q_1 = full(Mp_all{1});
    Sp_q_1 = full(Sp_all{1});
    
    Mp_q_2 = full(Mp_all{2});
    Sp_q_2 = full(Sp_all{2});
    %%

    Mp_uv = sparse(Mp_uv);
    Sp_uv = sparse(Sp_uv);   
end

%%

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
            if en_on
                if ((e <= Ex*N_en_y) || (e > (Ey-N_en_y)*Ex))
                    u(i,j,e) = u(i,j,e) - psi{1}(X(i,j,e),Y(i,j,e));
                    v(i,j,e) = v(i,j,e) - psi{2}(X(i,j,e),Y(i,j,e));
                end
            end
        end
    end
end

% u = apply_en_cont_soln(u,en_b_nodes,psi_p);
u1=u; u2=u; u3=u; fx3=u; fx2=u; fx1=u;
v1=v; v2=v; v3=v; fy3=v; fy2=v; fy1=v;

F = ones(size(ML));
% Fb=reshape(ML.*F,nL,1); %% Output
% 
% Fb=Bb\(Q'*Fb);  

uv = [u;v];

%%
uv_ic_check = uv;
%%

%% Setup BC if nonzero
% u_bc = ones(size(u)).*((Y==1)+(Y==-1));
% u_bc = reshape(ML.*u_bc,nL,1);
% u_bc = Bb\(Q'*u_bc);

%%
disp("Timestepping")
plot1 = 1;
time = 0;
plot1 = post_channel(N,Ex,Ey,w,X,Y,Ys,en_on,time,u,psi_xy,N_en_y,plot1);
%%
% psi_len = length(M_c);
% psi_c = zeros(psi_len,1);
% psi_c(1:psi_len/2) = -1*ones(psi_len/2,1)*psi_p;
% psi_c(psi_len/2+1:psi_len) = 1*ones(psi_len/2,1)*psi_p;
for step=1:nstep 
    time=step*dt;
    if step==1; b0=1.0;    b= [ -1 0 0 ]';       a=[ 1  0 0 ]'; end
    if step==2; b0=1.5;    b=([ -4 1 0 ]')./2;   a=[ 2 -1 0 ]'; end
    if step==3; b0=11./6.; b=([ -18 9 -2 ]')./6; a=[ 3 -3 1 ]'; end
    if step<=3
        if en_on
            H_x=(Ma + A*dt/(b0*Re) + dt/b0*(Mp_all{1} + Sp_all{1}));
            H_y=(Ma + A*dt/(b0*Re));
            [LH_x,UH_x]=lu(H_x);
            [LH_y,UH_y]=lu(H_y);
            terms_x = 1/Re*(T1_all{1})+T2_all{1};
            terms_y = 1/Re*(T1_all{2})+T2_all{2};
            
            H_uv = (Ma_uv + A_uv*dt/(b0*Re) + dt/b0*(Mp_uv + Sp_uv));
            H_c = (M_c + A_c*dt/(b0*Re) + dt/b0*(Mp_all_c{1}+Sp_all_c{1}));
            H_check = (Ma_uv_check + A_uv_check*dt/(b0*Re) + dt/b0*(Mp_uv_check + Sp_uv_check));
            H_q_check = R*Q'*(apply_en_cont(H_check(1:nL,1:nL),en_b_nodes,psi_p)+apply_en_cont(H_check(1:nL,nL+1:2*nL),en_b_nodes,psi_p));
            H_q = full(H_uv);
            rhs_c = (H_c);
            [LH_uv,UH_uv]=lu(H_uv);
        else
            H=(Ma + A*dt/(b0*Re));
            [LH_x,UH_x]=lu(H);
            LH_y = LH_x;
            UH_y = UH_x;

            terms_x = zeros(N+1,N+1,E);
            terms_y = zeros(N+1,N+1,E);
            
            H_uv = (Ma_uv + A_uv*dt/(b0*Re));
            rhs_c = zeros(size(Ma(:,1)));
            [LH_uv,UH_uv]=lu(H_uv);
%             Hbar=(Bb+ Ab*dt/(b0*Re));
        end
        
        b0i=1./b0;
    end % Viscous op
 
    %% uv version
    
%   Nonlinear step - unassembled, not multiplied by mass matrix

    fx1 = -convl(u,RX,Dh,u,v) + F + terms_x; % du = Cu  
    fy1 = -convl(v,RX,Dh,u,v) + terms_y; % dv = Cv
    
    %%
    convl_x_check = -convl(u,RX,Dh,u,v);
    convl_y_check = -convl(v,RX,Dh,u,v);
    %%

    rx  = a(1)*fx1+a(2)*fx2+a(3)*fx3; % kth-order extrapolation
    ry  = a(1)*fy1+a(2)*fy2+a(3)*fy3;

    fx3=fx2; fx2=fx1; 
    fy3=fy2; fy2=fy1; 

    rx  = dt*rx - (b(1)*u+b(2)*u2+b(3)*u3); u3=u2; u2=u; % Add BDF terms
    ry  = dt*ry - (b(1)*v+b(2)*v2+b(3)*v3); v3=v2; v2=v; %     and

    ut  = b0i*rx; 
    vt  = b0i*ry; 

%   uL=ut; vL=vt;
    [uL,vL,pr]=pressure_project(ut,vt,Ai,Q,ML,RX,Dh); % Div-free velocity
    pr = (b0/dt)*pr;
    
    %%
    uL_check = uL;
    vL_check = vL;
    %%

    %   Set RHS.                 %Viscous update. %  Convert to local form.
%     u=R*(Q'*reshape(ML.*uL,nL,1)-Hbar*u_bc);
    u_rhs=R*(Q'*reshape(ML.*uL,nL,1));
    v_rhs=R*(Q'*reshape(ML.*vL,nL,1));
    
        %%
    u_rhs_check = u_rhs;
    v_rhs_check = v_rhs;
    %%
    
    uv = [u_rhs+rhs_c;v_rhs];
    uv_check = uv;
%     uv = [u;v];
    uv=UH_uv\(LH_uv\uv);
    u = uv(1:nn);
    v = uv(nn+1:2*nn);
%     u=Q*(R'*u+u_bc);
    u=Q*(R'*u);
    if en_on
        u = apply_en_cont_soln(Ey,N_en_y,en_b_nodes,u,psi_p);
    end
    u=reshape(u,N1,N1,E);
    v=Q*(R'*v);
    v=reshape(v,N1,N1,E);
    
    
    
%% Output
    if mod(step,100)==0
        plot1 = post_channel(N,Ex,Ey,w,X,Y,Ys,en_on,time,u,psi_xy,N_en_y,plot1);
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
