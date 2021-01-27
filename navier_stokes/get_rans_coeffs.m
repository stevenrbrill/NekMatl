function [mu_t,gam_k,gam_omg,G_k,G_omg,Y_k,Y_omg,S_k,S_omg,R1,R2,R3,omg_w] ...
        = get_rans_coeffs(rho,mu,k,omg_prime,SS,OS,dkdx,dkdy,domg_primedx,domg_primedy,y,u,v)

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
    
    a=6;
    y_w = 1-abs(y);
    dy_wdx = 0;
    dy_wdy = -1+2.*(y>0);
    hess_y_w = 0;
    omg_w = a.*mu./rho./(beta_0.*y_w.*y_w);
    domg_wdx = -2*dy_wdx./y_w.*omg_w;
    domg_wdy = -2*dy_wdy./y_w.*omg_w;
    omg = omg_prime + omg_w;
    domgdx = domg_primedx + domg_wdx;
    domgdy = domg_primedy + domg_wdy;
    
    Re_t = rho.*k./(mu.*omg+eps);
    
    alpha_star = alpha_inf_star.*(alpha_0_star+Re_t./R_k)./(1+Re_t./R_k);
    % alpha_star can b 1
    alpha = alpha_inf./alpha_star.*(alpha_0+Re_t./R_omg)./(1+Re_t./R_omg);
    % alpha can be 13/25
            
    xk = 1./(omg.^3 + eps).*(dkdx.*domgdx+dkdy.*domgdy); 
    f_beta_star = 1.*(xk < 0) + (1+680*xk.^2)./(1+400*xk.^2).*(xk>=0);
    
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
    
    mu_t = alpha_star.*rho.*k./(omg+eps);
    gam_k = mu + mu_t./sigma_k;
    gam_omg = mu + mu_t./sigma_omg;
    G_k = mu_t.*S2;
    G_omg = alpha.*omg./(k+eps).*G_k;
    Y_k = rho.*beta_star.*k.*omg;
    Y_omg = rho.*beta.*omg.^2;
    S_k = 0.*k;
    S_omg = 0.*k;
    
    
    % Regularlized terms
    hess_omg_w = omg_w.*(6*(dy_wdx.^2+dy_wdy.^2)./(y_w.^2)-2*(hess_y_w)./(y_w));
    R1 = (mu+mu_t./sigma_omg).*hess_omg_w;
    
    r_tilde = Re_t./R_k; % typo in paper?
    f_r_tilde = (1-alpha_0_star)./((alpha_0_star+r_tilde).*(1+r_tilde));
    R2 = -2*rho.*alpha_star./(sigma_omg).*(1+r_tilde.*f_r_tilde).*(omg_w./omg)...
        .*((dy_wdx.*dkdx + dy_wdy.*dkdy)./y_w ...
        + 2*(omg_w./omg).*(dy_wdx.^2+dy_wdy.^2)./(y_w.^2).*k ...
        -(omg_w./omg).*k.*(dy_wdx.*domg_primedx + dy_wdy.*domg_primedy)./(y_w.*omg_w));
    
    R3 = 2*omg_w.*rho.*(u.*dy_wdx + v.*dy_wdy)./y_w;
    

end
