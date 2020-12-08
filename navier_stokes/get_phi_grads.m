function [dphi_dxi, dphi_deta, dpdx_dpdx, dpdy_dpdy] = get_phi_grads(N1,Dh)

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

dpdx_dpdx = zeros(N1,N1,N1*N1,N1*N1);
dpdy_dpdy = zeros(N1,N1,N1*N1,N1*N1);
for i=1:N1*N1
    for j=1:N1*N1
        dpdx_dpdx(:,:,i,j) = dphi_dxi(:,:,i).*dphi_dxi(:,:,j);
        dpdy_dpdy(:,:,i,j) = dphi_deta(:,:,i).*dphi_deta(:,:,j);
    end
end   
