function A_y = form_Ay(N1,E,w2d,J_y,J,dpdy_dpdy,Re)

A_y = zeros(N1*N1*E);
JRe = J.*1./Re;
for ie=1:E
    wJRe = JRe(:,:,ie).*w2d;
    for i=1:N1*N1
        for j=1:N1*N1
            A_y((ie-1)*N1*N1+i,(ie-1)*N1*N1+j) = ... 
                sum(wJRe.*dpdy_dpdy(:,:,i,j,ie),'All');           
        end
    end    
end