function [mu_t,gam_k,gam_omg,G_k,G_omg,Y_k,Y_omg,S_k,S_omg] ...
        = get_rans_coeffs(rho,mu,k,omg,SS,OS,dkdx,dkdy,domgdx,domgdy)

    sigma_k = 2;
    sigma_omg = 2;

    R_beta = 8;
    R_k = 6;
    R_omg = 2.95;
    
    beta_inf_star = 0.09;
    beta_0 = 0.072;
    
    alpha_inf = 0.52;
    alpha_inf_star = 1;
    alpha_0 = 1/9;
    alpha_0_star = beta_0/3;
    
    Re_t = rho.*k./(mu.*omg);
    
    alpha_star = alpha_inf_star.*(alpha_0_star+Re_t./R_k)./(1+Re_t./R_k);
    % alpha_star can b 1
    alpha = alpha_inf./alpha_star.*(alpha_0+Re_t./R_omg)./(1+Re_t./R_omg);
    % alpha can be 13/25
        
    
    xk = 1./(omg.^3).*(dkdx.*domgdx+dkdy.*domgdy); 
    f_beta_star = 1.*(xk < 0) + (1+680*xk.^2)./(1+400*xk.^2).*(xk>0);
    
    x_omg = 0.*k; 
    for i=1:2
        for j=1:2
            for ik=1:2
                x_omg = x_omg + OS(:,:,:,i,j).*OS(:,:,:,j,ik).*SS(:,:,:,ik,i);
            end
        end
    end
    x_omg = abs(x_omg./(beta_inf_star.*omg).^3);
    S2 = 2*(SS(:,:,:,1,1).^2+SS(:,:,:,1,2).^2+SS(:,:,:,2,1).^2+SS(:,:,:,2,2).^2); 
    
    f_beta = (1+70*x_omg)./(1+80*x_omg);
    beta = beta_0.*f_beta;

    
    beta_star = beta_inf_star.*f_beta_star.*((4/15+(Re_t./R_beta).^4)./(1+(Re_t./R_beta).^4));
    % Last term is only in some formulations
    
    mu_t = alpha_star.*rho.*k./omg;
    gam_k = mu + mu_t./sigma_k;
    gam_omg = mu + mu_t./sigma_omg;
    G_k = mu_t.*S2;
    G_omg = alpha.*omg./k.*G_k;
    Y_k = rho.*beta_star.*k.*omg;
    Y_omg = rho.*beta.*omg.^2;
    S_k = 0.*k;
    S_omg = 0.*k;

end
