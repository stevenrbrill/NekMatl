function A_xy = form_Axy(N1,E,w2d,J_x,J_y,J,dpdx_dpdy,Re)

A_xy = zeros(N1*N1*E);
for ie=1:E
    for i=1:N1*N1
        for j=1:N1*N1
            A_xy((ie-1)*N1*N1+i,(ie-1)*N1*N1+j) = ... 
                sum(J_x(ie)*J_y(ie)*J(:,:,ie).*1./Re(:,:,ie).*w2d.*dpdx_dpdy(:,:,i,j),'All');           
        end
    end    
end