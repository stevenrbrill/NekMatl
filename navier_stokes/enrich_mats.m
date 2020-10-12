function[Mp,Sp,T1,T2,z,w] =  enrich_mats(Xe,Ye,E,N,fpsi,fgpsi,fhpsi)
%
%                                                 ^
%     Compute the single element 1D SEM Stiffness Mass, and Convection
%     matrices, as well as the points and weights for a polynomial
%     of degree N
%

[Ah,Bh,Ch,Dh,zb,wb] = semhat(N);
[z,w] = zwgll(4*N);
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
    
    psi = fpsi(z_x,z_y');
    gpsi{1} = fgpsi{1}(z_x,z_y');
    gpsi{2} = fgpsi{2}(z_x,z_y');
    hpsi{1} = fhpsi{1}(z_x,z_y');
    hpsi{2} = fhpsi{2}(z_x,z_y');
    
    
%     X=Xe(:,:,e); Y=Ye(:,:,e); 
%     Xr  = Dh*X; Xs = X*Dh';
%     Yr  = Dh*Y; Ys = Y*Dh';
    Jac = L_x*L_y*ones(size(w2d));
%     Jac = Xr.*Ys - Xs.*Yr;
%     rx  =  Ys./Jac; ry  = -Xs./Jac;
%     sx  = -Yr./Jac; sy  =  Xr./Jac;
    
    for k = 1:2
        for i = 1:nb
            for j = 1:nb
                Mp{k}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*phi2d(:,:,j).*gpsi{k},'All');
            end
        end
    end

    for i = 1:nb
        for j = 1:nb
            Sp{1}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*dxphi2d(:,:,j).*psi,'All');
            Sp{2}(i,j,e) = sum(Jac.*w2d.*phi2d(:,:,i).*dyphi2d(:,:,j).*psi,'All');
        end
    end
    
    for k = 1:2
        for j = 1:nb
            T1{k}(j,e) = hpsi{k}(j);%sum(Jac.*w2d.*phi2d(:,:,j).*hpsi{k},'All');
            T2{k}(j,e) = 0;
        end
    end
    
end