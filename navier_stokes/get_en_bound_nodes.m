function[nodes] =  get_en_bound_nodes(Ex,Ey,N,N_en_y)

E = Ex*Ey;

m = Ex*N;  % Number of x locations (counting periodicity)
n = Ey*N+1;  % Number of y locations

m_orig = Ex*(N+1);
n_orig = Ey*(N+1);

count = 1;
% Bottom
y_orig = (N_en_y*(N+1));
for j = 1:Ex
    for i = 1:N+1
        nodes(count) = y_orig*m_orig+(j-1)*(N+1)+i;
        count = count + 1;
    end
end

% Top
y_orig = ((Ey-N_en_y)*(N+1))-1;
for j = 1:Ex
    for i = 1:N+1
        nodes(count) = y_orig*m_orig+(j-1)*(N+1)+i;
        count = count + 1;
    end
end
