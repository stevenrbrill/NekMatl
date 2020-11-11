function[nodes] =  get_en_bound_nodes2(Ex,Ey,N,N_en_y)

% Probably doesn't work for Ex>1?

E = Ex*Ey;

m = Ex*N;  % Number of x locations (counting periodicity)
n = Ey*N+1;  % Number of y locations

m_orig = Ex*(N+1);
n_orig = Ey*(N+1);
nb = (N+1)*(N+1);

if N_en_y >= Ey/2
    nodes = zeros(0);
else
    count = 1;
    % Bottom
    y_orig = (N_en_y*(N+1))-1;
    for i = 1:Ex*N
        nodes(count,1) = (y_orig-1)*Ex*N+i;
        nodes(count,2) = (y_orig)*Ex*N+i;
        count = count + 1;
    end
    
    
    % Top
    y_orig = ((Ey-N_en_y)*(N+1))-1;
    for i = 1:Ex*N
        nodes(count,2) = (y_orig-1)*Ex*N+i;
        nodes(count,1) = (y_orig)*Ex*N+i;
        count = count + 1;
    end
end