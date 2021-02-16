function [Yp,Up] = get_Up(Y,U,u_tau,mu,lotw,Ey,N,psi,en_on)

    U_poly = U-psi{1}(0,Y)*en_on;
    N_sample = 1000;
    for i = 1:Ey
        basis = get_nodal_basis_coeffs(Y((N+1)*(i-1)+1:(N+1)*(i)));
        y_sample((N_sample)*(i-1)+1:(N_sample)*(i)) = linspace(Y((N+1)*(i-1)+1),Y((N+1)*(i)),N_sample);
        for j = 1:N+1 
            phi_sample(:,j) = polyval(basis(j,:),y_sample((N_sample)*(i-1)+1:(N_sample)*(i)));
        end
        u_sample((N_sample)*(i-1)+1:(N_sample)*(i)) = phi_sample*U_poly((N+1)*(i-1)+1:(N+1)*(i));
    end
    u_sample = u_sample + psi{1}(0,y_sample)*en_on;

    Yp = (Y+1)*u_tau/mu;
    Up = U/u_tau;
    
    
    
    Yp_sample = (y_sample+1)*u_tau/mu;
    Up_sample = u_sample/u_tau;
    
    figure
    semilogx(Yp_sample,lotw(Yp_sample),'k-')
    hold on
    semilogx(Yp_sample,Up_sample,'r--')
    semilogx(Yp,Up,'bo')
    legend('Law of the Wall','Subsampled Solution','Nodal Values')
    

end