function [mu_t,gam_k,gam_omg,G_k,G_omg,Y_k,Y_omg,S_k,S_omg] ...
        = get_rans_coeffs(rho,mu,k,omg)

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
    alpha = alpha_inf./alpha_star.*(alpha_0+Re_t./R_omg)./(1+Re_t./R_omg);
        
    
    xk = 0.*k; %1./(omg.^3).*grad_k.*grad_omg; % TODO: Check this
    f_beta_star = 1.*(xk < 0) + (1+680*xk.^2)./(1+400*xk.^2).*(xk>0);
    
    x_omg = 0.*k; % TODO: Add this. it's complicated
    S2 = 0.*k; % TODO: Add this.
    
    f_beta = (1+70*x_omg)./(1+80*x_omg);
    beta = beta_0.*f_beta;

    
    beta_star = beta_inf_star.*f_beta_star.*((4/15+(Re_t./R_beta).^4)./(1+(Re_t./R_beta).^4));
    
    mu_t = alpha_star.*rho.*k./omg;
    gam_k = mu + mu_t./sigma_k;
    gam_omg = mu + mu_t./sigma_omg;
    G_k = mu_t.*S2;
    G_omg = alpha.*omg./k.*G_k;
    Y_k = rho.*beta_star.*k.*omg;
    Y_omg = rho.*beta.*f_beta.*omg.^2;
    S_k = 0.*k;
    S_omg = 0.*k;

end
