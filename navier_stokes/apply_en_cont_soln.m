function u = apply_en_cont_soln(Ex,Ey,N_en_y,u,psi)

if N_en_y < Ey/2
    % Top
    Estart = N_en_y*Ex;
    u(:,1,Estart+1:Estart+Ex) = u(:,1,Estart+1:Estart+Ex) + psi;
    
    % Bottom
    Estart = (Ey-N_en_y-1)*Ex;
    u(:,end,Estart+1:Estart+Ex) = u(:,end,Estart+1:Estart+Ex) + psi;
end