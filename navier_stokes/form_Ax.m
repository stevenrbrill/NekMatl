function A_x = form_Ax(N1,E,w2d,J_x,J,dpdx_dpdx,Re)

A_x = zeros(N1*N1*E);
for ie=1:E
    for i=1:N1*N1
        for j=1:N1*N1
            A_x((ie-1)*N1*N1+i,(ie-1)*N1*N1+j) = ... 
                sum(J_x(ie)*J_x(ie)*J(:,:,ie).*1./Re(:,:,ie).*w2d.*dpdx_dpdx(:,:,i,j),'All');           
        end
    end    
end