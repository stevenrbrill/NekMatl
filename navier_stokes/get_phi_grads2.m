function [dphi_dxi, dphi_deta, dphi_dx, dphi_dy, dpdx_dpdx, dpdy_dpdy, dpdx_dpdy, ...
        dpdx_dpdx_flat, dpdy_dpdy_flat, dpdx_dpdy_flat,phi_2d_flat] ...
    = get_phi_grads2(N1,E,J_x,J_y,N_quad)

nq = N_quad+1;
N = N1-1;
[z_basis,w_basis] = zwgll(N);
[z,w] = zwgll(N_quad); % Quad weights and points
basis = get_nodal_basis_coeffs(z_basis);
% Compute basis at quad points
phi = zeros(nq,N+1);
phi_p = zeros(nq,N+1);
for i = 1:N1
    phi(:,i) = polyval(basis(i,:),z);
    phi_p(:,i) = polyval(polyder(basis(i,:)),z);
end

N1_quad = N_quad+1;

phi_2d = zeros(N1_quad,N1_quad,N1*N1);
phi_2d_flat = zeros(N1_quad*N1_quad,N1*N1);
ind=1;
for i=1:N1
    for j=1:N1       
        phi_2d(:,:,ind) = phi(:,j)*phi(:,i)';
        phi_2d_flat(:,ind) = reshape(phi_2d(:,:,ind),[N1_quad*N1_quad,1]);
        ind = ind +1;
    end
end


dphi_dxi = zeros(N1_quad,N1_quad,N1*N1);
dphi_deta = zeros(N1_quad,N1_quad,N1*N1);
ind=1;
for i=1:N1
    for j=1:N1       
        dphi_dxi(:,:,ind) = phi_p(:,j)*phi(:,i)';
        dphi_deta(:,:,ind) = phi(:,j)*phi_p(:,i)';
        ind = ind +1;
    end
end

dphi_dx = zeros(N1_quad,N1_quad,N1*N1,E);
dphi_dy = zeros(N1_quad,N1_quad,N1*N1,E);
for ie=1:E
    ind=1;
    for i=1:N1
        for j=1:N1        
            dphi_dx(:,:,ind,ie) = dphi_dxi(:,:,ind)*J_x(ie);
            dphi_dy(:,:,ind,ie) = dphi_deta(:,:,ind)*J_y(ie);
            ind = ind +1;
        end
    end
end

dpdx_dpdx = zeros(N1_quad,N1_quad,N1*N1,N1*N1,E);
dpdy_dpdy = zeros(N1_quad,N1_quad,N1*N1,N1*N1,E);
dpdx_dpdy = zeros(N1_quad,N1_quad,N1*N1,N1*N1,E);
dpdx_dpdx_flat = zeros(N1_quad*N1_quad,N1*N1,E);
dpdy_dpdy_flat = zeros(N1_quad*N1_quad,N1*N1,E);
dpdx_dpdy_flat = zeros(N1_quad*N1_quad,N1*N1,E);
for ie=1:E
for i=1:N1*N1
    for j=1:N1*N1
        dpdx_dpdx(:,:,i,j,ie) = dphi_dx(:,:,i,ie).*dphi_dx(:,:,j,ie);
        dpdy_dpdy(:,:,i,j,ie) = dphi_dy(:,:,i,ie).*dphi_dy(:,:,j,ie);
        dpdx_dpdy(:,:,i,j,ie) = dphi_dx(:,:,i,ie).*dphi_dy(:,:,j,ie);
        
        dpdx_dpdx_flat(:,(i-1)*N1*N1+j,ie) = reshape(dpdx_dpdx(:,:,i,j,ie),[N1_quad*N1_quad,1]);
        dpdy_dpdy_flat(:,(i-1)*N1*N1+j,ie) = reshape(dpdy_dpdy(:,:,i,j,ie),[N1_quad*N1_quad,1]);
        dpdx_dpdy_flat(:,(i-1)*N1*N1+j,ie) = reshape(dpdx_dpdy(:,:,i,j,ie),[N1_quad*N1_quad,1]);
    end
end   
end