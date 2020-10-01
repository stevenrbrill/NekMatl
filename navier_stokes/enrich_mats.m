function[Mp,Sp,T1,T2,z,w] =  enrich_mats(Xe,Ye,E,N,fpsi,fgpsi,fhpsi)
%
%                                                 ^
%     Compute the single element 1D SEM Stiffness Mass, and Convection
%     matrices, as well as the points and weights for a polynomial
%     of degree N
%

[Ah,Bh,Ch,Dh,z,w] =  semhat(N);
[z,w] = zwgll(N);
basis = get_nodal_basis_coeffs(z);
phi = zeros(N+1,length(z));
phi_p = zeros(N+1,length(z));
nb = (N+1)*(N+1);
phi2d = zeros(length(z),length(z),nb);
dxphi2d = zeros(length(z),length(z),nb);
dyphi2d = zeros(length(z),length(z),nb);
w2d = w*w';

psi = fpsi(z,z');
gpsi{1} = fgpsi{1}(z',z);
gpsi{2} = fgpsi{2}(z',z);
hpsi{1} = fhpsi{1}(z',z);
hpsi{2} = fhpsi{2}(z',z);

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
    X=Xe(:,:,e); Y=Ye(:,:,e); 
    Xr  = Dh*X; Xs = X*Dh';
    Yr  = Dh*Y; Ys = Y*Dh';
    Jac = Xr.*Ys - Xs.*Yr;
    rx  =  Ys./Jac; ry  = -Xs./Jac;
    sx  = -Yr./Jac; sy  =  Xr./Jac;
    
    for k = 1:2
        for i = 1:nb
            for j = 1:nb
                Mp{k}(i,j,e) = sum(Jac*w2d.*phi2d(:,:,i).*phi2d(:,:,j).*gpsi{k},'All');
            end
        end
    end
    
    for i = 1:nb
        for j = 1:nb
            Sp{1}(i,j,e) = sum(Jac*w2d.*phi2d(:,:,i).*dxphi2d(:,:,j).*psi,'All');
            Sp{2}(i,j,e) = sum(Jac*w2d.*phi2d(:,:,i).*dyphi2d(:,:,j).*psi,'All');
        end
    end
    
    for k = 1:2
        for j = 1:nb
            T1{k}(j,e) = sum(Jac*w2d.*phi2d(:,:,j).*hpsi{k},'All');
            T2{k}(j,e) = 0;
        end
    end
    
end