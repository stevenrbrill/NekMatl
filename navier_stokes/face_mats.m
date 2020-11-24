function[G_total,Gs_x,Gs_y] =  face_mats(Xe,Ye,E,N)
%
%                                                 ^
%     Compute the single element 1D SEM Stiffness Mass, and Convection
%     matrices, as well as the points and weights for a polynomial
%     of degree N
%

[Ah,Bh,Ch,Dh,zb,wb] = semhat(N);
N_quad = N;
[z,w] = zwgll(N_quad);
nq = length(z);
basis = get_nodal_basis_coeffs(zb);
phi = zeros(nq,N+1);
phi_p = zeros(nq,N+1);
nb = (N+1)*(N+1);
phi2d = zeros(nq,nq,nb);
dxphi2d = zeros(nq,nq,nb);
dyphi2d = zeros(nq,nq,nb);

% Compute basis at quad points
for i = 1:N+1
    phi(:,i) = polyval(basis(i,:),z);
    phi_p(:,i) = polyval(polyder(basis(i,:)),z);
end

count = 1;
for j = 1:N+1
    for i = 1:N+1
        phi2d(:,:,count) = phi(:,i)*phi(:,j)';
        dxphi2d(:,:,count) = phi_p(:,i)*phi(:,j)';
        dyphi2d(:,:,count) = phi(:,i)*phi_p(:,j)';
        count = count + 1;
    end
end

G{1} = zeros(nb,nb,E);
G{2} = zeros(nb,nb,E);
G{3} = zeros(nb,nb,E);
G{4} = zeros(nb,nb,E);
G_total = zeros(nb*E);
Gs_x = zeros(nb*E);
Gs_y = zeros(nb*E);


for e = 1:E
    X_min = min(min(Xe(:,:,e)));
    X_max = max(max(Xe(:,:,e)));
    Y_min = min(min(Ye(:,:,e)));
    Y_max = max(max(Ye(:,:,e)));
    L_x = X_max-X_min;
    L_y = Y_max-Y_min;
    w2d = w*w';
    
    z_x = L_x/2*(z--1)+X_min;
    z_y = L_y/2*(z--1)+Y_min;
    
    zb_x = L_x/2*(zb--1)+X_min;
    zb_y = L_y/2*(zb--1)+Y_min;
   
    
    Jac = L_x/2*L_y/2*ones(size(w2d));
    Jac_x = L_x/2*ones(size(w));
    Jac_y = L_y/2*ones(size(w'));
    ijac_y = 1/Jac_y(1);
    ijac_x = 1/Jac_x(1);
    
%     for k = 1:4
        for i = 1:nb
            for j = 1:nb
                G{1}(j,i,e) = G{1}(j,i,e) + sum(-Jac_x.*w.*phi2d(:,1,j).*dyphi2d(:,1,i)*ijac_y,'All');
                G{2}(j,i,e) = G{2}(j,i,e) + sum(Jac_y.*w'.*phi2d(end,:,j).*dxphi2d(end,:,i)*ijac_x,'All');
                G{3}(j,i,e) = G{3}(j,i,e) + sum(Jac_x.*w.*phi2d(:,end,j).*dyphi2d(:,end,i)*ijac_y,'All');
                G{4}(j,i,e) = G{4}(j,i,e) + sum(-Jac_y.*w'.*phi2d(1,:,j).*dxphi2d(1,:,i)*ijac_x,'All');
            end
        end
        G_total(1+nb*(e-1):nb*(e),1+nb*(e-1):nb*(e)) = G{1}(:,:,e) + G{2}(:,:,e) + G{3}(:,:,e) + G{4}(:,:,e);
        Gs_y(1+nb*(e-1):nb*(e),1+nb*(e-1):nb*(e)) = G{1}(:,:,e) + G{3}(:,:,e);
        Gs_x(1+nb*(e-1):nb*(e),1+nb*(e-1):nb*(e)) = G{2}(:,:,e) + G{4}(:,:,e);
%     end
end

