function A_x = form_Ax(N1,E,w2d,J_x,J,dpdx_dpdx,Re)

A_x = zeros(N1*N1*E);
JRe = J.*1./Re;
for ie=1:E
    wJRe = JRe(:,:,ie).*w2d;
    for i=1:N1*N1
        for j=1:N1*N1
            A_x((ie-1)*N1*N1+i,(ie-1)*N1*N1+j) = ... 
                sum(wJRe.*dpdx_dpdx(:,:,i,j,ie),'All');           
        end
    end    
end