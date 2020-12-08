function A_y = form_Ay(N1,E,w2d,J_y,J,dpdy_dpdy,Re)

A_y = zeros(N1*N1*E);
for ie=1:E
    for i=1:N1*N1
        for j=1:N1*N1
            A_y((ie-1)*N1*N1+i,(ie-1)*N1*N1+j) = ... 
                sum(J_y(ie)*J_y(ie)*J(:,:,ie).*1./Re(:,:,ie).*w2d.*dpdy_dpdy(:,:,i,j),'All');           
        end
    end    
end