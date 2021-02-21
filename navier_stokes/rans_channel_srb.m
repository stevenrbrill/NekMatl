% clc
clear all
close all


set(0,'DefaultAxesFontSize',14)
set(0,'defaultLineLineWidth',2.5)
% set(0,'DefaultFigureWindowStyle','docked')


linestyles = {'k-','b-','r-','g-','c-','m-','k--','b--','r--','g--','c--','m--','k-.','b-.','r-.','g-.','c-.','m-.','k:','b:','r:','g:','c:','m:'};
pointstyles = {'ko','bo','ro','go','co','mo','k^','b^','r^','g^','c^','m^','ks','bs','rs','gs','cs','ms','k*','b*','r*','g*','c*','m*','k+','b+','r+','g+','c+','m+'};
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navier-Stokes Solver Demo (2D Channel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact;
format short; 
mu = 1/100000; %1/395; 1/6874;
rho = 1;
Re = rho/mu; 
dpdx = 1;
k_bc_val = 0;
omg_bc_val = 0;

N=6; % polynomial order  
Ex=1; % Number of elements in x
Ey=8; % Number of elements in y
Tfinal=300; 
CFL=500;
head = '';

rans_on = 1;
exp_mesh = 1;

dir_name = [head,'re',num2str(ceil(Re)),'_p',num2str(N),'_e',num2str(Ey)];
soln_dir = dir_name;
plot_soln = 1;
save_soln = 0;
plot_int = 5000;
save_soln_int = 5000;
restart = 0;
rst_step = 300000;

u_ic = Re;
pert = 0.0;
f_ic = @(x,y) 3/2*(1-y.^2);

%% Enrichment information
en_on = 0;
N_en_y = 1; 
delay_en = 0;
en_start_time = 100;



en_mag = 3;
psi = {@(x,y) en_mag*(0.5*(1 - y.^2) + 0.*x), @(x,y) 0.*y + 0.*x};
gpsi = {@(x,y) 0.*y + 0.*x, @(x,y) en_mag*(-1.*y + 0.*x), ...
        @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
hpsi = {@(x,y) 0.*y + 0.*x, @(x,y) en_mag*(-1 - 0.*y + 0.*x), ...
        @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
    
% psi = {@(x,y) (0.5*(1 - y.^4) + 0.*x), @(x,y) 0.*y + 0.*x};
% gpsi = {@(x,y) 0.*y + 0.*x, @(x,y) (-2.*y.^3 + 0.*x), ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
% hpsi = {@(x,y) 0.*y + 0.*x, @(x,y) (-6.*y.^2 + 0.*x), ...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
%    

    
% Law of the wall
u_tau = sqrt(0.316./(Re.^0.25)/8);
Re_t = u_tau/mu; %Re;
% u_tau = 1;
nu = 1/Re_t;
kap = 0.41;
beta = 5.2;
dypdy = u_tau/nu;
ypb = 11.062299784340414;
yp = @(y) (1-abs(y))*Re_t;
u_tau = u_tau;
lotw = @(yp) ((yp <= ypb).*yp + (yp > ypb).*(1./kap.*log(yp+eps)+beta));
psi = {@(x,y) u_tau*((yp(y) <= ypb).*yp(y) + (yp(y) > ypb).*(1./kap.*log(yp(y)+eps)+beta) + 0.*x), @(x,y) 0.*y + 0.*x};
gpsi = {@(x,y) 0.*y + 0.*x,  @(x,y) -1*sign(y).*u_tau.*(((yp(y) <= ypb).*1 + (yp(y) > ypb).*1./(kap*(yp(y)+eps)))*dypdy + 0.*x),...
        @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
hpsi = {@(x,y) 0.*y + 0.*x,  @(x,y) u_tau*(((yp(y) <= ypb).*0 + (yp(y) > ypb).*-1./(kap*(yp(y)+eps).^2))*dypdy*dypdy + 0.*x),...
        @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
    

% 1/7th
% mag = 1;
% yd = @(y) 1-abs(y)+eps;
% psi = {@(x,y) mag*((yd(y)).^(1/7).*(yd(y)>eps*10)+ 0.*x), @(x,y) 0.*y + 0.*x};
% gpsi = {@(x,y) 0.*y + 0.*x,  @(x,y) mag*(-1*sign(y).*1./7.*(yd(y)).^(-6/7).*(yd(y)>eps*10) + 0.*x),...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
% hpsi = {@(x,y) 0.*y + 0.*x,  @(x,y) mag*(-6/49*(yd(y)).^(-13/7).*(yd(y)>eps*10) + 0.*x),...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
    

%% Plot psi
figure(5)
ys_plot = linspace(-1,1,1000);
plot(psi{1}(0,ys_plot),ys_plot)
xlabel('\psi')
ylabel('y')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
yticks(linspace(-1,1,Ey+1));

%% Begin Solve
E=Ex*Ey; % Total number of elements
N1=N+1;
disp("Generating Matrices")
Q=makeq(Ex,Ey,N); % Global continuity
R=maker(Q,Ex,N); % Restriction matrix, applies Dirichlet conditions
[X,Y]=make_geom_channel(Ex,Ey,N);      % Geometry in local form
if exp_mesh
    [X,Y]=make_geom_channel_exp(Ex,Ey,N);      % Geometry in local form
end
en_b_nodes = get_en_bound_nodes(Ex,Ey,N,N_en_y);
dom_vol = (max(max(max(X)))-min(min(min(X))))*(max(max(max(Y)))-min(min(min(Y))));

Ys = zeros(Ey*length(Y(1,:,1)),1);
for i = 1:Ey
    Ys((i-1)*(N+1)+1:(i)*(N+1)) = Y(1,:,(i-1)*Ex+1);
end

psi_p = psi{1}(0,Ys(N_en_y*(N+1)));
dypsi_p1 = gpsi{2}(0,Ys(N_en_y*(N+1)));
dypsi_p2 = gpsi{2}(0,Ys(end-N_en_y*(N+1)));
        
[G,J,ML,RX]=make_coef(X,Y);
[Ah,Bh,Ch,Dh,z,w]=semhat(N);
[Gs, Gs_x, Gs_y] =  face_mats(X,Y,E,N);
Gs_c = R*Q'*apply_en_cont(Gs,en_b_nodes,psi_p);
Gsb = R*Q'*Gs*Q*R';


[n,nb]=size(R); 
nL=N1*N1*E;

Ab=spalloc(nL,nL,nL*N);    %  Generate Abar
for j=1:nL
    x=zeros(nL,1); 
    x(j)=1; 
    Ab(:,j)=abar_c(x,Dh,G,Q);  % assembled Neumann Operator
end
Ab_save = Ab;
A_c=R*Q'*apply_en_cont(Ab,en_b_nodes,psi_p);
S_check=full(Ab);
Ab = Q'*Ab*Q;
A=R*Ab*R'; % Full stiffness matrix
Bb=reshape(ML,nL,1); % Form Mass matrix
Bb=diag(Bb); 
Bb=sparse(Bb); 
Ma=R*Q'*Bb*Q*R'; % Assembling mass matrix  % Full mass matrix
M_c = R*Q'*apply_en_cont(Bb,en_b_nodes,psi_p);


%% Form matrices for rans

[L_x,L_y,J_x,J_y] = dir_jac(E,X,Y);
w2d = w*w';
w1d = reshape(w2d,[N1*N1,1])';
[dphi_dxi, dphi_deta, dphi_dx, dphi_dy, dpdx_dpdx, dpdy_dpdy, dpdx_dpdy, ...
        dpdx_dpdx_flat, dpdy_dpdy_flat, dpdx_dpdy_flat] ...
    = get_phi_grads(N1,Dh,E,J_x,J_y);
dphi_dx_flat = reshape(dphi_dx,N1*N1,N1*N1,E);
dphi_dy_flat = reshape(dphi_dy,N1*N1,N1*N1,E);

w2d_e = zeros(N+1,N+1,E);
for ie = 1:E
    w2d_e(:,:,ie) = w2d*L_x(ie)/2*L_y(ie)/2;
end

A_x = form_Ax(N1,E,w2d,J_x,J,dpdx_dpdx,ones(size(J)));
A_y = form_Ay(N1,E,w2d,J_y,J,dpdy_dpdy,ones(size(J)));
A_xy = form_Axy(N1,E,w2d,J_x,J_y,J,dpdx_dpdy,ones(size(J)));


SS = zeros(N1,N1,E,2,2);
OS = zeros(N1,N1,E,2,2);
dudx = zeros(N1,N1,E);
dudy = zeros(N1,N1,E);
dvdx = zeros(N1,N1,E);
dvdy = zeros(N1,N1,E);
dkdx = zeros(N1,N1,E);
dkdy = zeros(N1,N1,E);
domgdx = zeros(N1,N1,E);
domgdy = zeros(N1,N1,E);

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
G_uv = zeros(2*nn);
G_uv(1:nn,1:nn) = Gsb;
G_uv(nn+1:2*nn,nn+1:2*nn) = Gsb;
G_uv = sparse(G_uv);



%% Assemble enrichment matrices
psi_xy = zeros(N+1,N+1,E);
if en_on
    [psi_xy,psi_xy_act,gpsi_xy_act,Sp_all,Sp_all_Q,Jac_e_flat,gpsi_e_flat,Mp_uv,Sp_uv,Mp_all_c,Sp_all_c,Mp_full,Sp_full] = ...
        assemble_enrichment(X,Y,Ex,Ey,E,N,N1,nn,nL,J,Q,R,N_en_y,psi,gpsi,hpsi,en_b_nodes,psi_p);
end

%%

dxmin=pi*(z(N1)-z(N))/(2*Ex); % Get min dx for CFL constraint
dt=CFL*dxmin/u_ic; 

% dt = 1e-3;
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
k=k_bc_val*ones(size(ML));
omg=omg_bc_val*ones(size(ML));

for e = 1:E
    for i = 1:N+1
        for j = 1:N+1
            u(i,j,e) = f_ic(X(i,j,e),Y(i,j,e))+pert*rand()*u_ic;
            v(i,j,e) = v(i,j,e)+pert*rand()*u_ic;
        end
    end
end
if en_on
    [u, v] = project_to_enrich(X,Y,E,Ex,Ey,N,N_en_y,u,v,psi);
end

darcy = 0.316./(Re.^0.25);
u_tau = sqrt(darcy/8);
Yp = max((1-abs(Y))*u_tau*Re,1e-3)+eps;
sigma = 0.6;
fact = exp((-(log10(Yp)-1).^2)./(2*sigma.^2));
k = k_bc_val + 4.5*u_tau*u_tau*fact;

eps_s = 3;
omg_bc_val = 0; 
omg = omg_bc_val + 0.5*Re*u_tau*u_tau*fact;

% u = apply_en_cont_soln(u,en_b_nodes,psi_p);
u1_0=u; u2_0=u; u3_0=u; fx3_0=u; fx2_0=u; fx1_0=u; u_0 = u;
v1_0=v; v2_0=v; v3_0=v; fy3_0=v; fy2_0=v; fy1_0=v; v_0 = v;
u1=u; u2=u; u3=u; fx3=u; fx2=u; fx1=u;
v1=v; v2=v; v3=v; fy3=v; fy2=v; fy1=v;
k1=k; k2=k; k3=k; fk3=k; fk2=k; fk1=k;
omg1=omg; omg2=omg; omg3=omg; fomg3=omg; fomg2=omg; fomg1=omg;

F = dpdx*ones(size(ML));
% Fb=reshape(ML.*F,nL,1); %% Output
% 
% Fb=Bb\(Q'*Fb);  

uv = [u;v];

%%
uv_ic_check = uv;

%% Setup BC if nonzero
% u_bc = ones(size(u)).*((Y==1)+(Y==-1));
% u_bc = reshape(ML.*u_bc,nL,1);
% u_bc = Bb\(Q'*u_bc);

k_bc = k_bc_val.*ones(size(k));
k_bc = reshape(ML.*k_bc,nL,1);
k_bc = (Q'*Bb*Q)\(Q'*k_bc);

omg_bc = omg_bc_val.*ones(size(omg));
omg_bc = reshape(ML.*omg_bc,nL,1);
omg_bc = (Q'*Bb*Q)\(Q'*omg_bc);

%%
disp("Timestepping")
step = 1;
if (restart == 1)
   fname = strcat(soln_dir,"/soln_Re_",num2str(Re),"_P_",num2str(N),"_",num2str(Ex),"x",num2str(Ey),"_step_",num2str(rst_step),".mat");
   load(fname);
   step = rst_step;
end

plot1 = 1;
time = 0;
[plot1,u_mean] = post_channel(N,Ex,Ey,w,X,Y,Ys,en_on,step,time,u,psi_xy,N_en_y,plot1);
if save_soln
    mkdir(soln_dir);
end

%%
while step <= nstep
    time=step*dt;
    
    if (time > en_start_time) && (~en_on) && (delay_en)
        en_on = 1;
        [psi_xy,psi_xy_act,gpsi_xy_act,Sp_all,Sp_all_Q,Jac_e_flat,gpsi_e_flat,Mp_uv,Sp_uv,Mp_all_c,Sp_all_c,Mp_full,Sp_full] = ...
            assemble_enrichment(X,Y,Ex,Ey,E,N,N1,nn,nL,J,Q,R,N_en_y,psi,gpsi,hpsi,en_b_nodes,psi_p);
        [u, v] = project_to_enrich(X,Y,E,Ex,Ey,N,N_en_y,u,v,psi);
    end
    
    % Form combined u for RANS terms
    u_comb = u;
    if en_on
        u_comb = u+psi_xy_act{1};
    end
    
    %% Form S
    u_flat = reshape(u,N1*N1,E);
    v_flat = reshape(v,N1*N1,E);
    k_flat = reshape(k,N1*N1,E);
    omg_flat = reshape(omg,N1*N1,E);

    
    for ie = 1:E
        dudx(:,:,ie) = reshape(dphi_dx_flat(:,:,ie)*u_flat(:,ie),N1,N1);
        dvdx(:,:,ie) = reshape(dphi_dx_flat(:,:,ie)*v_flat(:,ie),N1,N1);
        dudy(:,:,ie) = reshape(dphi_dy_flat(:,:,ie)*u_flat(:,ie),N1,N1);
        dvdy(:,:,ie) = reshape(dphi_dy_flat(:,:,ie)*v_flat(:,ie),N1,N1);
        dkdx(:,:,ie) = reshape(dphi_dx_flat(:,:,ie)*k_flat(:,ie),N1,N1);
        dkdy(:,:,ie) = reshape(dphi_dy_flat(:,:,ie)*k_flat(:,ie),N1,N1);
        domgdx(:,:,ie) = reshape(dphi_dx_flat(:,:,ie)*omg_flat(:,ie),N1,N1);
        domgdy(:,:,ie) = reshape(dphi_dy_flat(:,:,ie)*omg_flat(:,ie),N1,N1);
    end
    
    if en_on
        dudx = dudx + gpsi_xy_act{1};
        dudy = dudy + gpsi_xy_act{2};
        dvdx = dvdx + gpsi_xy_act{3}; % Doesn't matter given assumptions
        dvdy = dvdy + gpsi_xy_act{4}; % Doesn't matter given assumptions
    end
    
    % Modify enrichment
%     if (mod(step,2000)==0) && (time > 50)
%         u_tau_calc = sqrt(mu*dudy(:,1,1)/rho);
%         u_tau = u_tau_calc(1);
%         psi = {@(x,y) u_tau*((yp(y) <= ypb).*yp(y) + (yp(y) > ypb).*(1./kap.*log(yp(y)+eps)+beta) + 0.*x), @(x,y) 0.*y + 0.*x};
%         gpsi = {@(x,y) 0.*y + 0.*x,  @(x,y) -1*sign(y).*u_tau.*(((yp(y) <= ypb).*1 + (yp(y) > ypb).*1./(kap*(yp(y)+eps)))*dypdy + 0.*x),...
%             @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
%         hpsi = {@(x,y) 0.*y + 0.*x,  @(x,y) u_tau*(((yp(y) <= ypb).*0 + (yp(y) > ypb).*-1./(kap*(yp(y)+eps).^2))*dypdy*dypdy + 0.*x),...
%             @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
%         
%         psi_p = psi{1}(0,Ys(N_en_y*(N+1)));
%         dypsi_p1 = gpsi{2}(0,Ys(N_en_y*(N+1)));
%         dypsi_p2 = gpsi{2}(0,Ys(end-N_en_y*(N+1)));
%         
%         A_c=R*Q'*apply_en_cont(Ab_save,en_b_nodes,psi_p);
%         M_c = R*Q'*apply_en_cont(Bb,en_b_nodes,psi_p);
%         
%         [psi_xy,psi_xy_act,gpsi_xy_act,Sp_all,Sp_all_Q,Jac_e_flat,gpsi_e_flat,Mp_uv,Sp_uv,Mp_all_c,Sp_all_c,Mp_full,Sp_full] = ...
%             assemble_enrichment(X,Y,Ex,Ey,E,N,N1,nn,nL,J,Q,R,N_en_y,psi,gpsi,hpsi,en_b_nodes,psi_p);
%         [u, v] = project_to_enrich(X,Y,E,Ex,Ey,N,N_en_y,u_comb,v,psi);
%         
%         en_on = 2;
%         
%         for ie = 1:E
%             dudx(:,:,ie) = reshape(dphi_dx_flat(:,:,ie)*u_flat(:,ie),N1,N1);
%             dvdx(:,:,ie) = reshape(dphi_dx_flat(:,:,ie)*v_flat(:,ie),N1,N1);
%             dudy(:,:,ie) = reshape(dphi_dy_flat(:,:,ie)*u_flat(:,ie),N1,N1);
%             dvdy(:,:,ie) = reshape(dphi_dy_flat(:,:,ie)*v_flat(:,ie),N1,N1);
%             dkdx(:,:,ie) = reshape(dphi_dx_flat(:,:,ie)*k_flat(:,ie),N1,N1);
%             dkdy(:,:,ie) = reshape(dphi_dy_flat(:,:,ie)*k_flat(:,ie),N1,N1);
%             domgdx(:,:,ie) = reshape(dphi_dx_flat(:,:,ie)*omg_flat(:,ie),N1,N1);
%             domgdy(:,:,ie) = reshape(dphi_dy_flat(:,:,ie)*omg_flat(:,ie),N1,N1);
%         end
%         
%         if en_on
%             dudx = dudx + gpsi_xy_act{1};
%             dudy = dudy + gpsi_xy_act{2};
%             dvdx = dvdx + gpsi_xy_act{3}; % Doesn't matter given assumptions
%             dvdy = dvdy + gpsi_xy_act{4}; % Doesn't matter given assumptions
%         end
%     end
    
    
    SS(:,:,:,1,1) = 1/2*(dudx+dudx);
    SS(:,:,:,1,2) = 1/2*(dudy+dvdx);
    SS(:,:,:,2,1) = 1/2*(dvdx+dudy);
    SS(:,:,:,2,2) = 1/2*(dvdy+dvdy);
    OS(:,:,:,1,1) = 1/2*(dudx-dudx);
    OS(:,:,:,1,2) = 1/2*(dudy-dvdx);
    OS(:,:,:,2,1) = 1/2*(dvdx-dudy);
    OS(:,:,:,2,2) = 1/2*(dvdy-dvdy);
    

    %% Get rans values
    [mu_t,gam_k,gam_omg,G_k,G_omg,Y_k,Y_omg,S_k,S_omg,R1,R2,R3,R4,omg_w] ...
        = get_rans_coeffs(rho,mu,k,omg,SS,OS,dkdx,dkdy,domgdx,domgdy,Y,u_comb,v);
    if ~rans_on
        mu_t = zeros(size(mu_t));
    end
    Re_t = rho./(mu_t+eps);
    Re_comb = rho./(mu*ones(size(mu_t))+mu_t);
    Re_k = rho./gam_k;
    Re_omg = rho./gam_omg;
    %% Timestepping
    
    if step==1; b0=1.0;    b= [ -1 0 0 ]';       a=[ 1  0 0 ]'; end
    if step==2; b0=1.5;    b=([ -4 1 0 ]')./2;   a=[ 2 -1 0 ]'; end
    if step>=3; b0=11./6.; b=([ -18 9 -2 ]')./6; a=[ 3 -3 1 ]'; end

    if ~rans_on
        % Compute A
        [A_x_full,A_y_full,A_xy_full] = form_Ax_Ay_Axy(N1,E,w1d,J,dpdx_dpdx_flat,dpdy_dpdy_flat,dpdx_dpdy_flat,Re_comb);
        A_x = R*Q'*A_x_full*Q*R';
        A_y = R*Q'*A_y_full*Q*R';
        A_xy = R*Q'*A_xy_full*Q*R';
        A_uv = zeros(2*nn);
        A_uv(1:nn,1:nn) = 2*A_x+A_y+A_xy;
        A_uv(nn+1:2*nn,nn+1:2*nn) = A_x+2*A_y+A_xy';
        A_uv = sparse(A_uv);
        
        T1_rhs = 0;
        
        if en_on
            A_full = 2*A_x_full+A_y_full+A_xy_full;
            A_c=R*Q'*apply_en_cont(A_full,en_b_nodes,psi_p);
            [T1_new] = form_T1_psi(E,N,w1d,Jac_e_flat,dphi_dy_flat,gpsi_e_flat,Re_comb,N_en_y);
            T1_rhs = (dt/b0)*R*Q'*reshape(T1_new,nL,1);
            
            H_uv = (Ma_uv + (A_uv)*dt/(b0) + dt/b0*(Mp_uv + Sp_uv));
            H_c = (M_c + (A_c)*dt/(b0) + dt/b0*(Mp_all_c{1}+Sp_all_c{1}));
            if en_on == 2
                H_uv = (Ma_uv + (A_uv)*dt/(b0)); % + dt/b0*(Mp_uv + Sp_uv));
                H_c = (M_c + (A_c)*dt/(b0)); % + dt/b0*(Mp_all_c{1}+Sp_all_c{1}));
            end
            %             H_q = full(H_uv);
            rhs_c = (H_c);
            [LH_uv,UH_uv]=lu(H_uv);
        end
    end
    
%     if step>=3
        if en_on          
            if rans_on
                [A_x_k,A_y_k] = form_Ax_Ay(N1,E,w1d,J,dpdx_dpdx_flat,dpdy_dpdy_flat,Re_k);
                A_k = R*Q'*(A_x_k + A_y_k)*Q*R';
                Ab_k = Q'*(A_x_k + A_y_k)*Q;
                H_k=(Ma+A_k*dt/(b0)+dt/b0*Sp_all{1});
                H_k_bar = (Q'*Bb*Q+ Ab_k*dt/(b0)+dt/b0*Sp_all_Q{1});
                if en_on == 2
                    H_k=(Ma+A_k*dt/(b0));
                    H_k_bar = (Q'*Bb*Q+ Ab_k*dt/(b0));
                end
        
                [LH_k,UH_k]=lu(H_k);
        
                % Assuming Re_omg=Re_k
                H_omg = H_k;
                H_omg_bar = H_k_bar;
                LH_omg = LH_k;
                UH_omg = UH_k;
            end
        else
            terms_x = zeros(N+1,N+1,E);
            terms_y = zeros(N+1,N+1,E);
            T1_rhs = 0;
            
            H_uv = (Ma_uv + (A_uv)*dt/(b0));
            rhs_c = zeros(size(Ma(:,1)));
            [LH_uv,UH_uv]=lu(H_uv);
            
            if rans_on
                [A_x_k,A_y_k] = form_Ax_Ay(N1,E,w1d,J,dpdx_dpdx_flat,dpdy_dpdy_flat,Re_k);
                A_k = R*Q'*(A_x_k + A_y_k)*Q*R';
                Ab_k = Q'*(A_x_k + A_y_k)*Q;
                H_k=(Ma+A_k*dt/(b0));
                H_k_bar = (Q'*Bb*Q+ Ab_k*dt/(b0));
        
                [LH_k,UH_k]=lu(H_k);
        
                % Assuming Re_omg=Re_k
                H_omg = H_k;
                H_omg_bar = H_k_bar;
                LH_omg = LH_k;
                UH_omg = UH_k;
            end
        end
        
        b0i=1./b0;
%     end % Viscous op

    
%   Nonlinear step - unassembled, not multiplied by mass matrix   
    %% Setup Convection terms and forcing terms
    fx1 = -convl(u,RX,Dh,u,v); % + F; % du = Cu  
    fy1 = -convl(v,RX,Dh,u,v); % dv = Cv
    fk1 = -convl(k,RX,Dh,u,v) + G_k - Y_k + S_k; % dk = Ck
    fomg1 = -convl(omg,RX,Dh,u,v) + G_omg - Y_omg + S_omg + R1 + R2 + R3 + R4; % domg = Comg
    fx1_0 = F;
    fy1_0 = zeros(size(F));
    
    if en_on == 2
        en_u = dt/b0*((Sp_full{1}+Mp_full{1})*reshape(u,nL,1)+Mp_full{2}*reshape(v,nL,1));
        en_u = reshape(en_u,N1,N1,E);
        en_u_0 = dt/b0*((Sp_full{1}+Mp_full{1})*reshape(u_0,nL,1)+Mp_full{2}*reshape(v_0,nL,1));
        en_u_0 = reshape(en_u_0,N1,N1,E);
        en_v = dt/b0*Mp_full{2}*reshape(v,nL,1);
        en_v = reshape(en_v,N1,N1,E);
        
        en_k = dt/b0*Sp_full{1}*reshape(k,nL,1);
        en_k = reshape(en_k,N1,N1,E); 
        en_omg = dt/b0*Sp_full{1}*reshape(omg,nL,1);
        en_omg = reshape(en_omg,N1,N1,E); 
        
        
        % SRB
        fx1 = fx1 - en_u;
        fx1_0 = fx1_0;% - en_u_0; 
        fy1 = fy1 - en_v;
        fk1 = fk1 - en_k;
        fomg1 = fomg1 - en_omg;
    end
    
    %% Setup residual terms in bdf

    rx   = a(1)*fx1   + a(2)*fx2   + a(3)*fx3; % kth-order extrapolation
    ry   = a(1)*fy1   + a(2)*fy2   + a(3)*fy3;
    rk   = a(1)*fk1   + a(2)*fk2   + a(3)*fk3;
    romg = a(1)*fomg1 + a(2)*fomg2 + a(3)*fomg3;
    rx_0 = a(1)*fx1_0   + a(2)*fx2_0   + a(3)*fx3_0;
    ry_0 = a(1)*fy1_0   + a(2)*fy2_0   + a(3)*fy3_0;

    fx3=fx2; fx2=fx1; 
    fy3=fy2; fy2=fy1; 
    fk3=fk2; fk2=fk1;
    fomg3=fomg2; fomg2=fomg1;
    fx3_0=fx2_0; fx2_0=fx1_0; 
    fy3_0=fy2_0; fy2_0=fy1_0; 

    %% Add residual terms to bdf
    rx  = dt*rx - (b(1)*u+b(2)*u2+b(3)*u3); u3=u2; u2=u; % Add BDF terms
    ry  = dt*ry - (b(1)*v+b(2)*v2+b(3)*v3); v3=v2; v2=v; %     and
    rk  = dt*rk - (b(1)*k+b(2)*k2+b(3)*k3); k3=k2; k2=k;
    romg = dt*romg - (b(1)*omg+b(2)*omg2+b(3)*omg3); omg3=omg2; omg2=omg;
    rx_0  = dt*rx_0 - (b(1)*u_0+b(2)*u2_0+b(3)*u3_0); u3_0=u2_0; u2_0=u_0; % Add BDF terms
    ry_0  = dt*ry_0 - (b(1)*v_0+b(2)*v2_0+b(3)*v3_0); v3_0=v2_0; v2_0=v_0; % Add BDF terms
    
    ut  = b0i*rx; 
    vt  = b0i*ry; 
    k   = b0i*rk;
    omg = b0i*romg;
    ut_0  = b0i*rx_0; 
    vt_0  = b0i*ry_0; 
 
    [uL,vL,pr]=pressure_project(ut,vt,Ai,Q,ML,RX,Dh); % Div-free velocity
    [uL_0,vL_0,pr_0]=pressure_project(ut_0,vt_0,Ai,Q,ML,RX,Dh); % Div-free velocity
    pr = (b0/dt)*pr;
    pr_0 = (b0/dt)*pr_0;
    

    %%

    %   Set RHS.                 %Viscous update. %  Convert to local form.
%     u=R*(Q'*reshape(ML.*uL,nL,1)-Hbar*u_bc);
    u_rhs=R*(Q'*reshape(ML.*uL,nL,1));
    v_rhs=R*(Q'*reshape(ML.*vL,nL,1));
    u_rhs_0=R*(Q'*reshape(ML.*uL_0,nL,1));
    v_rhs_0=R*(Q'*reshape(ML.*vL_0,nL,1));
    
    if rans_on
        k_rhs=R*(Q'*reshape(ML.*k,nL,1));%-H_k_bar*k_bc);
        omg_rhs=R*(Q'*reshape(ML.*omg,nL,1));%-H_omg_bar*omg_bc);
    end
    
    %%
    
    % Solve k and omg first
    if rans_on
        k = UH_k\(LH_k\k_rhs);
        omg = UH_omg\(LH_omg\omg_rhs);
        k=Q*(R'*k);%+k_bc);
        omg=Q*(R'*omg);%+omg_bc);
        k=reshape(k,N1,N1,E);
        omg=reshape(omg,N1,N1,E);
        mu_t = rho.*k./(omg_w+omg);
        
        % recompute matrices for uv 
        Re_t = rho./(mu_t+eps);
        Re_comb = rho./(mu*ones(size(mu_t))+mu_t*rans_on);
        
        [A_x_full,A_y_full,A_xy_full] = form_Ax_Ay_Axy(N1,E,w1d,J,dpdx_dpdx_flat,dpdy_dpdy_flat,dpdx_dpdy_flat,Re_comb);
        A_x = R*Q'*A_x_full*Q*R';
        A_y = R*Q'*A_y_full*Q*R';
        A_xy = R*Q'*A_xy_full*Q*R';
        A_uv = zeros(2*nn);
        A_uv(1:nn,1:nn) = 2*A_x+A_y+A_xy;
        A_uv(nn+1:2*nn,nn+1:2*nn) = A_x+2*A_y+A_xy';
        A_uv = sparse(A_uv);
                
        if ~en_on
            H_uv = (Ma_uv + (A_uv)*dt/(b0));
            rhs_c = zeros(size(Ma(:,1)));
            [LH_uv,UH_uv]=lu(H_uv);
        end
        
        if en_on           
            A_full = 2*A_x_full+A_y_full+A_xy_full;
            A_c=R*Q'*apply_en_cont(A_full,en_b_nodes,psi_p);
            [T1_new] = form_T1_psi(E,N,w1d,Jac_e_flat,dphi_dy_flat,gpsi_e_flat,Re_comb,N_en_y);
            T1_rhs = (dt/b0)*R*Q'*reshape(T1_new,nL,1);
            
            if en_on == 1
                H_uv = (Ma_uv + (A_uv)*dt/(b0) + dt/b0*(Mp_uv + Sp_uv));
                H_c = (M_c + (A_c)*dt/(b0) + dt/b0*(Mp_all_c{1}+Sp_all_c{1}));
            end
            if en_on == 2
                H_uv = (Ma_uv + (A_uv)*dt/(b0));
                H_c = (M_c + (A_c)*dt/(b0)); 
            end
            rhs_c = (H_c);
            [LH_uv,UH_uv]=lu(H_uv);
        end
    end
    
    % SRB
    u_rhs = u_rhs + T1_rhs;
    u_rhs_0 = u_rhs_0;% + T1_rhs;
    
    %% Solve for u_0, v_0, p_0
    % SRB
    uv_0 = [u_rhs_0;v_rhs_0];
    uv_0=UH_uv\(LH_uv\uv_0);
    
    u_0 = uv_0(1:nn);
    v_0 = uv_0(nn+1:2*nn);
    u_0=Q*(R'*u_0);
    %SRB
%     if en_on
%         u_0 = apply_en_cont_soln(Ey,N_en_y,en_b_nodes,u_0,psi_p);
%     end
    u_0=reshape(u_0,N1,N1,E);
    v_0=Q*(R'*v_0);
    v_0=reshape(v_0,N1,N1,E);
    
    
    %% Solve for u', v', p'
    %SRB
    uv = [u_rhs+rhs_c;v_rhs];
    uv=UH_uv\(LH_uv\uv); 
   
    u = uv(1:nn);
    v = uv(nn+1:2*nn);
    u=Q*(R'*u);
    %SRB
    if en_on
        u = apply_en_cont_soln(Ey,N_en_y,en_b_nodes,u,psi_p);
    end
    u=reshape(u,N1,N1,E);
    v=Q*(R'*v);
    v=reshape(v,N1,N1,E);
    
    %% Add in forcing term
    u_comb = u;
    u_0_comb = u_0;
    if en_on
        u_comb = u+psi_xy_act{1};
%         u_0_comb = u_0+psi_xy_act{1};
    end
    % Integrate domain
    avg_u = sum(u_comb.*w2d_e,'All')/dom_vol;
    avg_u_0 = sum(u_0_comb.*w2d_e,'All')/dom_vol;
    alpha_0 = (1-avg_u)/avg_u_0;
    
    % SRB 
    u = u + alpha_0*u_0;
%     u_comb = u_comb + alpha_0*u_comb;
%     u =  u_comb;
%     if en_on
%         u = u_comb - psi_xy_act{1};
%     end
    v = v + alpha_0*v_0;
    pr = pr + alpha_0*pr_0;
    
    
%% Output
    if mod(step,plot_int)==0 && plot_soln
        [plot1,u_mean] = post_channel(N,Ex,Ey,w,X,Y,Ys,en_on,step,time,u,psi_xy,N_en_y,plot1);
    end
    if mod(step,save_soln_int)==0 && save_soln
        disp(step)
        fname = strcat(soln_dir,"/soln_Re_",num2str(Re),"_P_",num2str(N),"_",num2str(Ex),"x",num2str(Ey),"_step_",num2str(step),".mat");
        save(fname); 
    end
    if isnan(u(1,2))
        step
        time
        return
    end
    
    step = step + 1;
end
if  plot_soln
    [plot1,u_mean] = post_channel(N,Ex,Ey,w,X,Y,Ys,en_on,step,time,u,psi_xy,N_en_y,plot1);
    get_Up(Ys,u_mean,u_tau,mu,lotw,Ey,N,psi,en_on);
end
if  save_soln
    fname = strcat(soln_dir,"/soln_Re_",num2str(Re),"_P_",num2str(N),"_",num2str(Ex),"x",num2str(Ey),"_step_",num2str(step),".mat");
    save(fname);
end