function [psi_xy,psi_xy_act,gpsi_xy_act,Sp_all,Sp_all_Q,Jac_e_flat,gpsi_e_flat,Mp_uv,Sp_uv,Mp_all_c,Sp_all_c] = assemble_enrichment(X,Y,Ex,Ey,E,N,N1,nn,nL,J,Q,R,N_en_y,psi,gpsi,hpsi,en_b_nodes,psi_p)

%% Assemble enrichment matrices
psi_xy{1} = psi{1}(X,Y);
psi_xy{2} = psi{2}(X,Y);
gpsi_xy{1} = gpsi{1}(X,Y);
gpsi_xy{2} = gpsi{2}(X,Y);
gpsi_xy{3} = gpsi{3}(X,Y);
gpsi_xy{4} = gpsi{4}(X,Y);
psi_xy_act{1} = zeros(size(X));
psi_xy_act{2} = zeros(size(X));
gpsi_xy_act{1} = zeros(size(X));
gpsi_xy_act{2} = zeros(size(X));
gpsi_xy_act{3} = zeros(size(X));
gpsi_xy_act{4} = zeros(size(X));
    for iy = 1:Ey
        for ix = 1:Ex
            if ((iy <= N_en_y) || (iy > Ey-N_en_y))
                i = (iy-1)*Ex+ix;
                psi_xy_act{1}(:,:,i) = psi_xy{1}(:,:,i);
                psi_xy_act{2}(:,:,i) = psi_xy{2}(:,:,i);
                gpsi_xy_act{1}(:,:,i) = gpsi_xy{1}(:,:,i);
                gpsi_xy_act{2}(:,:,i) = gpsi_xy{2}(:,:,i);
                gpsi_xy_act{3}(:,:,i) = gpsi_xy{3}(:,:,i);
                gpsi_xy_act{4}(:,:,i) = gpsi_xy{4}(:,:,i);
            end
        end
    end

    disp("Computing enrichment")
    [Mp,Sp,T1,T2,T1_alt,T1_alt2,Mp_alt,Sp_alt,z_en,w_en] = enrich_mats(X,Y,E,N,psi,gpsi,hpsi,J);
    nb = N1*N1;
    
    [L_x_e,L_y_e,Jac_e,gpsi_e,Jac_e_flat,gpsi_e_flat] =  precomp_en(X,Y,E,N,psi,gpsi,hpsi);  
        
    for is=1:2
        Mp_all{is} = zeros(nb*E,nb*E);
        Sp_all{is} = zeros(nb*E,nb*E);
        T1_all{is} = zeros(N+1,N+1,E);
        T2_all{is} = zeros(N+1,N+1,E);
        T1_alt_all{is} = zeros(N+1,N+1,E);
        T1_alt2_all{is} = zeros(N+1,N+1,E);
        Mp_alt_all{is} = zeros(nb*E,nb*E);
        Sp_alt_all{is} = zeros(nb*E,nb*E);
        
        T1_rs{is} = reshape(T1{is},[N+1,N+1,E]);
        T2_rs{is} = reshape(T2{is},[N+1,N+1,E]);
        T1_alt_rs{is} = reshape(T1_alt{is},[N+1,N+1,E]);
        T1_alt2_rs{is} = reshape(T1_alt2{is},[N+1,N+1,E]);
        
        % TODO: Assemble differently for different elements
        for iy = 1:Ey
            for ix = 1:Ex
                i = (iy-1)*Ex+ix;
                if ((iy <= N_en_y) || (iy > Ey-N_en_y))
                    Mp_all{is}((i-1)*nb+1:i*nb,(i-1)*nb+1:i*nb) = Mp{is}(:,:,i);
                    Sp_all{is}((i-1)*nb+1:i*nb,(i-1)*nb+1:i*nb) = Sp{is}(:,:,i);
                    T1_all{is}(:,:,i) = T1_rs{is}(:,:,i); 
                    T2_all{is}(:,:,i) = T2_rs{is}(:,:,i);
                    T1_alt_all{is}(:,:,i) = T1_alt_rs{is}(:,:,i); 
                    T1_alt2_all{is}(:,:,i) = T1_alt2_rs{is}(:,:,i); 
                end                
            end
        end
        %%
        if is == 1
            Mp_check_1 = full(Mp_all{is});
            Sp_check_1 = full(Sp_all{is});
            T1_check_1 = full(T1_all{is});
            T2_check_1 = full(T2_all{is});
            Mp_check_2 = zeros(size(Mp_all{is}));
            Sp_check_2 = zeros(size(Mp_all{is}));
            T1_check_2 = zeros(size(Mp_all{is}));
            T2_check_2 = zeros(size(Mp_all{is}));
        elseif is == 2
            Mp_check_2 = full(Mp_all{is});
            Sp_check_2 = full(Sp_all{is});
            T1_check_2 = full(T1_all{is});
            T2_check_2 = full(T2_all{is});
        end
        %%
        Mp_all{is} = sparse(Mp_all{is});        
        Mp_full{is} = Mp_all{is};
        Mp_all_c{is} = R*Q'*apply_en_cont(Mp_all{is},en_b_nodes,psi_p);
        Mp_all{is} = R*Q'*Mp_all{is}*Q*R';
        Sp_all{is} = sparse(Sp_all{is});
        Sp_full{is} = Sp_all{is};
        Sp_all_c{is} = R*Q'*apply_en_cont(Sp_all{is},en_b_nodes,psi_p);
        Sp_all_Q{is} = Q'*Sp_all{is}*Q;
        Sp_all{is} = R*Q'*Sp_all{is}*Q*R';
        
        Mp_alt_all{is} = sparse(Mp_alt_all{is});        
        Mp_alt_all_c{is} = R*Q'*apply_en_cont(Mp_alt_all{is},en_b_nodes,psi_p);
        Mp_alt_all{is} = R*Q'*Mp_alt_all{is}*Q*R';   
    end
    
    Mp_uv = zeros(2*nn);
    Sp_uv = zeros(2*nn);
    Mp_uv(1:nn,1:nn) = Mp_all{1};
    Mp_uv(1:nn,nn+1:2*nn) = Mp_all{2};
    Sp_uv(1:nn,1:nn) = Sp_all{1};
    Sp_uv(nn+1:2*nn,nn+1:2*nn) = Sp_all{1};
    
    %%
    Mp_uv_check = zeros(2*nL);
    Mp_uv_check(1:nL,1:nL) = Mp_check_1;
    Mp_uv_check(1:nL,nL+1:2*nL) = Mp_check_2;
    Sp_uv_check(1:nL,1:nL) = Sp_check_1;
    Sp_uv_check(nL+1:2*nL,nL+1:2*nL) = Sp_check_1;
    Mp_q_1 = full(Mp_all{1});
    Sp_q_1 = full(Sp_all{1});
    
    Mp_q_2 = full(Mp_all{2});
    Sp_q_2 = full(Sp_all{2});
    %%

    Mp_uv = sparse(Mp_uv);
    Sp_uv = sparse(Sp_uv);   


end