function[Mp,Sp,T1,T2,T1_alt,T1_alt2,Mp_alt,Sp_alt,z,w] =  enrich_mats(Xe,Ye,E,N,fpsi,fgpsi,fhpsi,jac)
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

Mp{1} = zeros(nb,nb,E);
Mp{2} = zeros(nb,nb,E);
Sp{1} = zeros(nb,nb,E);
Sp{2} = zeros(nb,nb,E);
T1{1} = zeros(nb,E);
T1{2} = zeros(nb,E);
T2{1} = zeros(nb,E);
T2{2} = zeros(nb,E);
T1_alt{1} = zeros(nb,E);
T1_alt{2} = zeros(nb,E);
T1_alt2{1} = zeros(nb,E);
T1_alt2{2} = zeros(nb,E);
Mp_alt{1} = zeros(nb,nb,E);
Mp_alt{2} = zeros(nb,nb,E);
Sp_alt{1} = zeros(nb,nb,E);
Sp_alt{2} = zeros(nb,nb,E);

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
    
    psi{1} = fpsi{1}(z_x,z_y');
    psi{2} = fpsi{2}(z_x,z_y');
    gpsi{1} = fgpsi{1}(z_x,z_y');
    gpsi{2} = fgpsi{2}(z_x,z_y');
    gpsi{3} = fgpsi{3}(z_x,z_y');
    gpsi{4} = fgpsi{4}(z_x,z_y');
    hpsi{1} = fhpsi{1}(z_x,z_y');
    hpsi{2} = fhpsi{2}(z_x,z_y');
    hpsi{3} = fhpsi{3}(z_x,z_y');
    hpsi{4} = fhpsi{4}(z_x,z_y');
    
    zb_x = L_x/2*(zb--1)+X_min;
    zb_y = L_y/2*(zb--1)+Y_min;
    
    psib{1} = fpsi{1}(zb_x,zb_y');
    psib{2} = fpsi{2}(zb_x,zb_y');
    gpsib{1} = fgpsi{1}(zb_x,zb_y');
    gpsib{2} = fgpsi{2}(zb_x,zb_y');
    gpsib{3} = fgpsi{3}(zb_x,zb_y');
    gpsib{4} = fgpsi{4}(zb_x,zb_y');
    hpsib{1} = fhpsi{1}(zb_x,zb_y');
    hpsib{2} = fhpsi{2}(zb_x,zb_y');
    hpsib{3} = fhpsi{3}(zb_x,zb_y');
    hpsib{4} = fhpsi{4}(zb_x,zb_y');
    
    Jac = L_x/2*L_y/2*ones(size(w2d));
    ijac = 1./Jac;
    
    mask_y = zeros(size(z_y));
    mask_y(1) = 1;
    mask_y(end) = 1;
    mask = ones(size(z_x))*mask_y';
    
    for k = 1:4
        for i = 1:nb
            for j = 1:nb
                Mp{k}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*phi2d(:,:,j).*gpsi{k},'All');
            end
        end
    end
    
    for k = 1:4
        for i = 1:nb
            for j = 1:nb
                Mp_alt{k}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*phi2d(:,:,j).*gpsi{k}.*mask,'All');
            end
        end
    end

    for i = 1:nb
        for j = 1:nb
            Sp{1}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*dxphi2d(:,:,j).*2/L_x.*psi{1},'All');
            Sp{2}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*dyphi2d(:,:,j).*2/L_y.*psi{1},'All');
            Sp{3}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*dxphi2d(:,:,j).*2/L_x.*psi{2},'All');
            Sp{4}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*dyphi2d(:,:,j).*2/L_y.*psi{2},'All');
            
            Sp_alt{1}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*dxphi2d(:,:,j).*2/L_x.*psi{1}.*mask,'All');
            Sp_alt{2}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*dyphi2d(:,:,j).*2/L_y.*psi{1}.*mask,'All');
            Sp_alt{3}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*dxphi2d(:,:,j).*2/L_x.*psi{2}.*mask,'All');
            Sp_alt{4}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*dyphi2d(:,:,j).*2/L_y.*psi{2}.*mask,'All');
        end
    end
    
    for k = 1:2
        for j = 1:nb
            %sum(Jac.*w2d.*phi2d(:,:,j).*hpsi{k},'All');
            T1{k}(j,e) = hpsib{2*(k-1)+1}(j)+hpsib{2*(k-1)+2}(j);
            T2{k}(j,e) = psib{1}(j)*gpsib{2*(k-1)+1}(j)+psib{2}(j)*gpsib{2*(k-1)+2}(j);
            T1_alt{k}(j,e) = -sum(Jac.*w2d.*dyphi2d(:,:,j).*2/L_y.*gpsi{2*(k-1)+2},'All') + sum(-L_x/2.*w.*phi2d(:,1,j).*gpsi{2*(k-1)+2}(:,1),'All')  + sum(L_x/2.*w.*phi2d(:,end,j).*gpsi{2*(k-1)+2}(:,end),'All');
            T1_alt2{k}(j,e) = sum(Jac.*w2d.*phi2d(:,:,j).*hpsi{2*(k-1)+2},'All');
        end
    end
    
end