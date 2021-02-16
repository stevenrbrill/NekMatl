function [u_new, v_new] = project_to_enrich(X,Y,E,Ex,Ey,N,N_en_y,u,v,psi)

    u_new = u;
    v_new = v;
    for e = 1:E
        for i = 1:N+1
            for j = 1:N+1
                if ((e <= Ex*N_en_y) || (e > (Ey-N_en_y)*Ex))
                    u_new(i,j,e) = u(i,j,e) - psi{1}(X(i,j,e),Y(i,j,e));
                    v_new(i,j,e) = v(i,j,e) - psi{2}(X(i,j,e),Y(i,j,e));
                end
            end
        end
    end
end