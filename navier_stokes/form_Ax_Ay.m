function [A_x,A_y] = form_Ax_Ay(N1,E,w2d,J_x,J_y,J,dpdx_dpdx,dpdy_dpdy,Re)

A_x = zeros(N1*N1*E);
A_y = zeros(N1*N1*E);
JRe = J.*1./Re;
for ie=1:E
    wJRe = JRe(:,:,ie).*w2d;
    for i=1:N1*N1
        for j=1:N1*N1
            A_x((ie-1)*N1*N1+i,(ie-1)*N1*N1+j) = ... 
                sum(wJRe.*dpdx_dpdx(:,:,i,j,ie),'All');     
            A_y((ie-1)*N1*N1+i,(ie-1)*N1*N1+j) = ... 
                sum(wJRe.*dpdy_dpdy(:,:,i,j,ie),'All');     
        end
    end    
end