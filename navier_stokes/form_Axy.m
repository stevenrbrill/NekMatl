function A_xy = form_Axy(N1,E,w2d,J_x,J_y,J,dpdx_dpdy,Re)

A_xy = zeros(N1*N1*E);
JRe = J.*1./Re;
for ie=1:E
    wJRe = JRe(:,:,ie).*w2d;
    for i=1:N1*N1
        for j=1:N1*N1
            A_xy((ie-1)*N1*N1+i,(ie-1)*N1*N1+j) = ... 
                sum(wJRe.*dpdx_dpdy(:,:,i,j,ie),'All');           
        end
    end    
end