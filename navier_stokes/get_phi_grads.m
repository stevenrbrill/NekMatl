function [dphi_dxi, dphi_deta, dphi_dx, dphi_dy, dpdx_dpdx, dpdy_dpdy, dpdx_dpdy] ...
    = get_phi_grads(N1,Dh,E,J_x,J_y)

dphi_dxi = zeros(N1,N1,N1*N1);
dphi_deta = zeros(N1,N1,N1*N1);
ind=1;
for i=1:N1
    for j=1:N1
        x=zeros(N1);
        x(j,i)=1;
        
        dphi_dxi(:,:,ind) = Dh*x;
        dphi_deta(:,:,ind) = x*Dh';
        ind = ind +1;
    end
end

dphi_dx = zeros(N1,N1,N1*N1,E);
dphi_dy = zeros(N1,N1,N1*N1,E);
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

dpdx_dpdx = zeros(N1,N1,N1*N1,N1*N1,E);
dpdy_dpdy = zeros(N1,N1,N1*N1,N1*N1,E);
dpdx_dpdy = zeros(N1,N1,N1*N1,N1*N1,E);
for ie=1:E
for i=1:N1*N1
    for j=1:N1*N1
        dpdx_dpdx(:,:,i,j,ie) = dphi_dx(:,:,i,ie).*dphi_dx(:,:,j,ie);
        dpdy_dpdy(:,:,i,j,ie) = dphi_dy(:,:,i,ie).*dphi_dy(:,:,j,ie);
        dpdx_dpdy(:,:,i,j,ie) = dphi_dx(:,:,i,ie).*dphi_dy(:,:,j,ie);
    end
end   
end
