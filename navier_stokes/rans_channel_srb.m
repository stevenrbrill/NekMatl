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
mu = 1/395;
rho = 1;
Re = rho/mu; 
dpdx = 1;
k_bc_val = 0;
omg_bc_val = 1*k_bc_val;

%N=16; E=5; N1=N+1; nL=N1*N1*E;  % 16th order
N=5; % polynomial order  
Ex=1; % Number of elements in x
Ey=7; % Number of elements in y
CFL=0.1;
u_ic = Re;
pert = 0.0;
f_ic = @(x,y) 3/2*(1-y.^2);

rans_on = 1;


soln_dir = "test";
save_soln = 1;
plot_int = 100;
save_soln_int = 500;

%% Enrichment information
en_on = 0;
N_en_y = 1; 
en_mag = 1;
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
% Re_t = Re;
% u_tau = 1;
% nu = 1/Re_t;
% kap = 0.41;
% beta = 5.2;
% dypdy = u_tau/nu;
% ypb = 11.062299784340414;
% yp = @(y) (1-abs(y))*Re_t;
% psi = {@(x,y) (yp(y) <= ypb).*yp(y) + (yp(y) > ypb).*(1./kap.*log(yp(y)+eps)+beta) + 0.*x, @(x,y) 0.*y + 0.*x};
% gpsi = {@(x,y) 0.*y + 0.*x,  @(x,y) ((yp(y) <= ypb).*1 + (yp(y) > ypb).*1/(kap*(yp(y)+eps)))*dypdy + 0.*x,...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
% hpsi = {@(x,y) 0.*y + 0.*x,  @(x,y) ((yp(y) <= ypb).*0 + (yp(y) > ypb).*-1./(kap*(yp(y)+eps).^2))*dypdy*dypdy + 0.*x,...
%         @(x,y) 0.*y + 0.*x, @(x,y) 0.*y + 0.*x};
    

%% Plot psi
figure(5)
ys_plot = linspace(-1,1,1000);
plot(psi{1}(0,ys_plot),ys_plot)
xlabel('\psi')
ylabel('y')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
yticks(linspace(-1,1,Ey+1));

psi_p = psi{1}(0,N_en_y*2/Ey-1);
dypsi_p1 = gpsi{2}(0,N_en_y*2/Ey-1);
dypsi_p2 = gpsi{2}(0,1-N_en_y*2/Ey);

%% Begin Solve
E=Ex*Ey; % Total number of elements
N1=N+1;
disp("Generating Matrices")
Q=makeq(Ex,Ey,N); % Global continuity
R=maker(Q,Ex,N); % Restriction matrix, applies Dirichlet conditions
[X,Y]=make_geom_channel(Ex,Ey,N);      % Geometry in local form
en_b_nodes = get_en_bound_nodes(Ex,Ey,N,N_en_y);

Ys = zeros(Ey*length(Y(1,:,1)),1);
for i = 1:Ey
    Ys((i-1)*(N+1)+1:(i)*(N+1)) = Y(1,:,(i-1)*Ex+1);
end
        
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
[dphi_dxi, dphi_deta, dphi_dx, dphi_dy, dpdx_dpdx, dpdy_dpdy, dpdx_dpdy] ...
    = get_phi_grads(N1,Dh,E,J_x,J_y);
dphi_dx_flat = reshape(dphi_dx,N1*N1,N1*N1,E);
dphi_dy_flat = reshape(dphi_dy,N1*N1,N1*N1,E);

myA = zeros(N1*N1*E);
for ie=1:E
    for i=1:N1*N1
        for j=1:N1*N1
            myA((ie-1)*N1*N1+i,(ie-1)*N1*N1+j) = ... 
                sum(J(:,:,ie).*w2d.*dpdx_dpdx(:,:,i,j,ie) ...
               + J(:,:,ie).*w2d.*dpdy_dpdy(:,:,i,j,ie),'All');           
        end
    end    
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
psi_xy{1} = psi{1}(X,Y);
psi_xy{2} = psi{2}(X,Y);
if en_on
    disp("Computing enrichment")
    [Mp,Sp,T1,T2,T1_alt,T1_alt2,Mp_alt,Sp_alt,z_en,w_en] = enrich_mats(X,Y,E,N,psi,gpsi,hpsi,J);
    nb = N1*N1;
        
    for is=1:2
        Mp_all{is} = zeros(nb*E,nb*E);
        Sp_all{is} = zeros(nb*E,nb*E);
        T1_all{is} = zeros(N+1,N+1,E);
        T2_all{is} = zeros(N+1,N+1,E);
        T1_alt_all{is} = zeros(N+1,N+1,E);
        T1_alt2_all{is} = zeros(N+1,N+1,E);
        Mp_alt_all{is} = zeros(nb*E,nb*E);
        Sp_alt_all{is} = zeros(nb*E,nb*E);
        
        T1_rs{is} = reshape(T1{is},[N+1,N+1,E]);
        T2_rs{is} = reshape(T2{is},[N+1,N+1,E]);
        T1_alt_rs{is} = reshape(T1_alt{is},[N+1,N+1,E]);
        T1_alt2_rs{is} = reshape(T1_alt2{is},[N+1,N+1,E]);
        
        % TODO: Assemble differently for different elements
        for iy = 1:Ey
            for ix = 1:Ex
                i = (iy-1)*Ex+ix;
                if ((iy <= N_en_y) || (iy > Ey-N_en_y))
                    Mp_all{is}((i-1)*nb+1:i*nb,(i-1)*nb+1:i*nb) = Mp{is}(:,:,i);
                    Sp_all{is}((i-1)*nb+1:i*nb,(i-1)*nb+1:i*nb) = Sp{is}(:,:,i);
                    T1_all{is}(:,:,i) = T1_rs{is}(:,:,i); 
                    T2_all{is}(:,:,i) = T2_rs{is}(:,:,i);
                    T1_alt_all{is}(:,:,i) = T1_alt_rs{is}(:,:,i); 
                    T1_alt2_all{is}(:,:,i) = T1_alt2_rs{is}(:,:,i); 
                end                
            end
        end
        %%
        if is == 1
            Mp_check_1 = full(Mp_all{is});
            Sp_check_1 = full(Sp_all{is});
            T1_check_1 = full(T1_all{is});
            T2_check_1 = full(T2_all{is});
            Mp_check_2 = zeros(size(Mp_all{is}));
            Sp_check_2 = zeros(size(Mp_all{is}));
            T1_check_2 = zeros(size(Mp_all{is}));
            T2_check_2 = zeros(size(Mp_all{is}));
        elseif is == 2
            Mp_check_2 = full(Mp_all{is});
            Sp_check_2 = full(Sp_all{is});
            T1_check_2 = full(T1_all{is});
            T2_check_2 = full(T2_all{is});
        end
        %%
        Mp_all{is} = sparse(Mp_all{is});        
        Mp_all_c{is} = R*Q'*apply_en_cont(Mp_all{is},en_b_nodes,psi_p);
        Mp_all{is} = R*Q'*Mp_all{is}*Q*R';
        Sp_all{is} = sparse(Sp_all{is});
        Sp_all_c{is} = R*Q'*apply_en_cont(Sp_all{is},en_b_nodes,psi_p);
        Sp_all{is} = R*Q'*Sp_all{is}*Q*R';
        
        Mp_alt_all{is} = sparse(Mp_alt_all{is});        
        Mp_alt_all_c{is} = R*Q'*apply_en_cont(Mp_alt_all{is},en_b_nodes,psi_p);
        Mp_alt_all{is} = R*Q'*Mp_alt_all{is}*Q*R';   
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
Tfinal=50; 
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
k=k_bc_val*ones(size(ML));
omg=omg_bc_val*ones(size(ML));

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

% u_tau = 1;
u_tau = sqrt(0.316./(Re.^0.25)/8);
Yp = max((1-abs(Y))*u_tau*Re,1e-3);
sigma = 0.6;
fact = exp((-(log10(Yp)-1).^2)./(2*sigma.^2));
k = k_bc_val + 4.5*u_tau*u_tau*fact;

eps_s = 3;
omg_bc_val = 40000*u_tau.^2/(mu*eps_s.^2);
omg = omg_bc_val + 0.5*Re*u_tau*u_tau*fact;

% u = apply_en_cont_soln(u,en_b_nodes,psi_p);
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
plot1 = 1;
time = 0;
plot1 = post_channel(N,Ex,Ey,w,X,Y,Ys,en_on,time,u,psi_xy,N_en_y,plot1);
if save_soln
    mkdir(soln_dir);
end
%%
% psi_len = length(M_c);
% psi_c = zeros(psi_len,1);
% psi_c(1:psi_len/2) = -1*ones(psi_len/2,1)*psi_p;
% psi_c(psi_len/2+1:psi_len) = 1*ones(psi_len/2,1)*psi_p;
for step=1:nstep
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
    SS(:,:,:,1,1) = 1/2*(dudx+dudx);
    SS(:,:,:,1,2) = 1/2*(dudy+dvdx);
    SS(:,:,:,2,1) = 1/2*(dvdx+dudy);
    SS(:,:,:,2,2) = 1/2*(dvdy+dvdy);
    OS(:,:,:,1,1) = 1/2*(dudx-dudx);
    OS(:,:,:,1,2) = 1/2*(dudy-dvdx);
    OS(:,:,:,2,1) = 1/2*(dvdx-dudy);
    OS(:,:,:,2,2) = 1/2*(dvdy-dvdy);
    

    %% Get rans values
    [mu_t,gam_k,gam_omg,G_k,G_omg,Y_k,Y_omg,S_k,S_omg] ...
        = get_rans_coeffs(rho,mu,k,omg,SS,OS,dkdx,dkdy,domgdx,domgdy);
    Re_t = rho./(mu_t+eps);
    Re_comb = rho./(mu*ones(size(mu_t))+mu_t*rans_on);
    Re_k = rho./gam_k;
    Re_omg = rho./gam_omg;
    %% Timestepping
    time=step*dt;
    if step==1; b0=1.0;    b= [ -1 0 0 ]';       a=[ 1  0 0 ]'; end
    if step==2; b0=1.5;    b=([ -4 1 0 ]')./2;   a=[ 2 -1 0 ]'; end
    if step==3; b0=11./6.; b=([ -18 9 -2 ]')./6; a=[ 3 -3 1 ]'; end
%     if step>=3
        if en_on
            H_x=(Ma + A*dt/(b0*Re) + dt/b0*(Mp_all{1} + Sp_all{1}));
            H_y=(Ma + A*dt/(b0*Re));
            [LH_x,UH_x]=lu(H_x);
            [LH_y,UH_y]=lu(H_y);
            terms_x = 1/Re*(T1_all{1})+T2_all{1};
            terms_y = 1/Re*(T1_all{2})+T2_all{2};
            
            T1_rhs = (dt/b0)*1/Re*R*Q'*reshape(T1_alt_all{1},nL,1);
            
            H_uv = (Ma_uv + (A_uv)*dt/(b0*Re) + dt/b0*(Mp_uv + Sp_uv));
            H_c = (M_c + (A_c)*dt/(b0*Re) + dt/b0*(Mp_all_c{1}+Sp_all_c{1}));
            H_q = full(H_uv);
            rhs_c = (H_c);
            [LH_uv,UH_uv]=lu(H_uv);
        else
%             H=(Ma + A*dt/(b0*Re));
%             [LH_x,UH_x]=lu(H);
%             LH_y = LH_x;
%             UH_y = UH_x;

            terms_x = zeros(N+1,N+1,E);
            terms_y = zeros(N+1,N+1,E);
            T1_rhs = 0;

            [A_x_full,A_y_full,A_xy_full] = form_Ax_Ay_Axy(N1,E,w2d,J_x,J_y,J,dpdx_dpdx,dpdy_dpdy,dpdx_dpdy,Re_comb);
            A_x = R*Q'*A_x_full*Q*R';
            A_y = R*Q'*A_y_full*Q*R';
            A_xy = R*Q'*A_xy_full*Q*R';
            A_uv = zeros(2*nn);
            A_uv(1:nn,1:nn) = 2*A_x+A_y+A_xy;
            A_uv(nn+1:2*nn,nn+1:2*nn) = A_x+2*A_y+A_xy';
            A_uv = sparse(A_uv);
            
            H_uv = (Ma_uv + (A_uv)*dt/(b0));
            H_check = (Ma_uv_check + A_uv_check*dt/(b0*Re)); 
            rhs_c = zeros(size(Ma(:,1)));
            [LH_uv,UH_uv]=lu(H_uv);
            
            if rans_on
                [A_x_k,A_y_k] = form_Ax_Ay(N1,E,w2d,J_x,J_y,J,dpdx_dpdx,dpdy_dpdy,Re_k);
                A_k = R*Q'*(A_x_k + A_y_k)*Q*R';
                Ab_k = Q'*(A_x_k + A_y_k)*Q;
                H_k=(Ma+A_k*dt/(b0));
                H_k_bar = (Q'*Bb*Q+ Ab_k*dt/(b0));
                
                [A_x_omg,A_y_omg] = form_Ax_Ay(N1,E,w2d,J_x,J_y,J,dpdx_dpdx,dpdy_dpdy,Re_omg);
                A_omg = R*Q'*(A_x_omg + A_y_omg)*Q*R';
                Ab_omg = Q'*(A_x_omg + A_y_omg)*Q;
                H_omg=(Ma+A_omg*dt/(b0));
                H_omg_bar = (Q'*Bb*Q+ Ab_omg*dt/(b0));
                
                [LH_k,UH_k]=lu(H_k);
                [LH_omg,UH_omg]=lu(H_omg);
            end
%             Hbar=(Bb+ Ab*dt/(b0*Re));
        end
        
        b0i=1./b0;
%     end % Viscous op
    

    
 
    %% uv version
    
%   Nonlinear step - unassembled, not multiplied by mass matrix

    fx1 = -convl(u,RX,Dh,u,v) + F; % du = Cu  
    fy1 = -convl(v,RX,Dh,u,v); % dv = Cv
    fk1 = -convl(k,RX,Dh,u,v) + G_k - Y_k + S_k; % dk = Ck
    fomg1 = -convl(omg,RX,Dh,u,v) + G_omg - Y_omg + S_omg; % domg = Comg
    
    %%

    rx   = a(1)*fx1   + a(2)*fx2   + a(3)*fx3; % kth-order extrapolation
    ry   = a(1)*fy1   + a(2)*fy2   + a(3)*fy3;
    rk   = a(1)*fk1   + a(2)*fk2   + a(3)*fk3;
    romg = a(1)*fomg1 + a(2)*fomg2 + a(3)*fomg3;

    fx3=fx2; fx2=fx1; 
    fy3=fy2; fy2=fy1; 
    fk3=fk2; fk2=fk1;
    fomg3=fomg2; fomg2=fomg1;

    rx  = dt*rx - (b(1)*u+b(2)*u2+b(3)*u3); u3=u2; u2=u; % Add BDF terms
    ry  = dt*ry - (b(1)*v+b(2)*v2+b(3)*v3); v3=v2; v2=v; %     and
    rk  = dt*rk - (b(1)*k+b(2)*k2+b(3)*k3); k3=k2; k2=k;
    romg = dt*romg - (b(1)*omg+b(2)*omg2+b(3)*omg3); omg3=omg2; omg2=omg;
    
    ut  = b0i*rx; 
    vt  = b0i*ry; 
    k   = b0i*rk;
    omg = b0i*romg;

%   uL=ut; vL=vt;
    [uL,vL,pr]=pressure_project(ut,vt,Ai,Q,ML,RX,Dh); % Div-free velocity
    pr = (b0/dt)*pr;
    

    %%

    %   Set RHS.                 %Viscous update. %  Convert to local form.
%     u=R*(Q'*reshape(ML.*uL,nL,1)-Hbar*u_bc);
    u_rhs=R*(Q'*reshape(ML.*uL,nL,1));
    v_rhs=R*(Q'*reshape(ML.*vL,nL,1));
    
    if rans_on
        k_rhs=R*(Q'*reshape(ML.*k,nL,1)-H_k_bar*k_bc);
        omg_rhs=R*(Q'*reshape(ML.*omg,nL,1)-H_omg_bar*omg_bc);
    end
    
    u_rhs = u_rhs + T1_rhs;
    
    %%
    
    uv = [u_rhs+rhs_c;v_rhs];
    uv=UH_uv\(LH_uv\uv);
 
    if rans_on
        k = UH_k\(LH_k\k_rhs);
        omg = UH_omg\(LH_omg\omg_rhs);
    end
    
    u = uv(1:nn);
    v = uv(nn+1:2*nn);
    u=Q*(R'*u);
    if en_on
        u = apply_en_cont_soln(Ey,N_en_y,en_b_nodes,u,psi_p);
    end
    u=reshape(u,N1,N1,E);
    v=Q*(R'*v);
    v=reshape(v,N1,N1,E);
    
    
    if rans_on
        k=Q*(R'*k+k_bc);
        k=reshape(k,N1,N1,E);
        omg=Q*(R'*omg+omg_bc);
        omg=reshape(omg,N1,N1,E);
    end
    
    
%% Output
    if mod(step,plot_int)==0
        plot1 = post_channel(N,Ex,Ey,w,X,Y,Ys,en_on,time,u,psi_xy,N_en_y,plot1);
        if mod(step,save_soln_int)==0 && save_soln
            fname = strcat(soln_dir,"/soln_Re_",num2str(Re),"_P_",num2str(N),"_",num2str(Ex),"x",num2str(Ey),"_step_",num2str(step),".mat");
            save(fname,"u","v","pr","psi","gpsi","hpsi","N","Ex","Ey","en_on","time","psi_xy","N_en_y","w","X","Y","Ys","Re","pert","k","omg","mu_t","k_bc","omg_bc","rho");
        end
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
