function[T1] =  form_T1_psi(E,N,w1d,Jac_flat,dyphi2d_flat,gpsi_flat,Re,N_en_y,phi_2d_flat)

nb = (N+1)*(N+1);
T1 = zeros(nb,E);

Re_flat = reshape(1./Re,[],E);

k=1;
for e = 1:E
    if (e <= N_en_y) || E-e < N_en_y
        Re_quad = phi_2d_flat*Re_flat(:,e); % Interpolate Re at quadrature points. Has some aliasing
        T1(:,e) = -((w1d'.*Re_quad.*Jac_flat(:,e).*gpsi_flat{2*(k-1)+2}(:,e))'*dyphi2d_flat(:,:,e))';
    end
end
