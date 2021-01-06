function[L_x,L_y,Jac,gpsi,Jac_flat,gpsi_flat] =  precomp_en(Xe,Ye,E,N,fpsi,fgpsi,fhpsi)

N_quad = N;
[z,w] = zwgll(N_quad);

L_x = zeros(E,1);
L_y = zeros(E,1);
gpsi{1} = zeros(length(z),length(z),E);
gpsi{2} = zeros(length(z),length(z),E);
gpsi{3} = zeros(length(z),length(z),E);
gpsi{4} = zeros(length(z),length(z),E);

for e = 1:E
    X_min = min(min(Xe(:,:,e)));
    X_max = max(max(Xe(:,:,e)));
    Y_min = min(min(Ye(:,:,e)));
    Y_max = max(max(Ye(:,:,e)));
    L_x(e) = X_max-X_min;
    L_y(e) = Y_max-Y_min;
    w2d = w*w';
    
    z_x = L_x(e)/2*(z--1)+X_min;
    z_y = L_y(e)/2*(z--1)+Y_min;
    
    gpsi{1}(:,:,e) = fgpsi{1}(z_x,z_y');
    gpsi{2}(:,:,e) = fgpsi{2}(z_x,z_y');
    gpsi{3}(:,:,e) = fgpsi{3}(z_x,z_y');
    gpsi{4}(:,:,e) = fgpsi{4}(z_x,z_y');

    Jac(:,:,e) = L_x(e)/2*L_y(e)/2*ones(size(w2d));
end

Jac_flat = reshape(Jac,[(N+1)*(N+1),E]);
gpsi_flat{1} = reshape(gpsi{1},[(N+1)*(N+1),E]);
gpsi_flat{2} = reshape(gpsi{2},[(N+1)*(N+1),E]);
gpsi_flat{3} = reshape(gpsi{3},[(N+1)*(N+1),E]);
gpsi_flat{4} = reshape(gpsi{4},[(N+1)*(N+1),E]);
