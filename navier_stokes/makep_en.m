function[P] =  makep_en(Ex,Ey,N,N_en_y)

E = Ex*Ey;

m = Ex*N;  % Number of x locations (counting periodicity)
n = Ey*N+1;  % Number of y locations

m_orig = Ex*(N+1);
n_orig = Ey*(N+1);

P = zeros(n_orig*m_orig,n*m);

% Bottom
y_orig = (N_en_y*(N+1));
y_aglo = (N_en_y*(N)-1);
for j = 1:Ex-1
    for i = 1:N+1
        P(y_orig*m_orig+(j-1)*(N+1)+i,y_aglo*m+(j-1)*N+i) = 1;
    end
end
j=Ex;
for i = 1:N
    P(y_orig*m_orig+(j-1)*(N+1)+i,y_aglo*m+(j-1)*N+i) = 1;
end
P((y_orig+1)*m_orig,y_aglo*m+1) = 1;

% Top
y_orig = ((Ey-N_en_y)*(N+1))+1;
y_aglo = ((Ey-N_en_y)*(N));
for j = 1:Ex-1
    for i = 1:N+1
        P(y_orig*m_orig+(j-1)*(N+1)+i,y_aglo*m+(j-1)*N+i) = -1;
    end
end
j=Ex;
for i = 1:N
    P(y_orig*m_orig+(j-1)*(N+1)+i,y_aglo*m+(j-1)*N+i) = -1;
end
P((y_orig+1)*m_orig,y_aglo*m+1) = -1;
