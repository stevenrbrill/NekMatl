function u = apply_en_cont_soln(Ey,N_en_y,nodes,u,psi)

if N_en_y < Ey/2
    u(nodes) = u(nodes) + psi;    
end