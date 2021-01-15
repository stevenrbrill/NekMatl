function[T1] =  form_T1_psi(E,N,w1d,Jac_flat,dyphi2d_flat,gpsi_flat,Re,N_en_y)

nb = (N+1)*(N+1);
T1 = zeros(nb,E);

Re_flat = reshape(Re,[],E);

k=1;
for e = 1:E
    if (e <= N_en_y) || E-e < N_en_y
        T1(:,e) = -((w1d'.*Re_flat(:,e).*Jac_flat(:,e).*gpsi_flat{2*(k-1)+2}(:,e))'*dyphi2d_flat(:,:,e))';
    end
end