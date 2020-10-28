function[Q] =  makeqp(Ex,Ey,N)

%
%  Make the Q array for an [ Ex x Ey ] array of elements of order N
%
%  Assumes periodicity in x direction
%
%  Numbers top and bottom rows last to make Dirichlet BCs easy
%
%
%        +--- numbered 2nd to last
%        |
%        v
%    +---------+---------+---------+---------+---------P
%    |         |         |         |         |         P <-- peiodic
%    |         |         |         |         |         P
%    |  e=1    |  e=2    |  e=3    |  ...    |  e=E    P
%    |         |         |         |         |         P
%    +---------+---------+---------+---------+---------P
%        ^
%        |
%        +--- numbered last
%
%
E = Ex*Ey;

m = Ex*N;  % Number of x locations (counting periodicity)
n = Ey*N+1;  % Number of y locations

q = zeros(m+1,n);


ig0 = 0;
for j=2:n
    ig = ig0+[1:m]; 
    q(1:end-1,j) = ig;
    q(end,j)=q(1,j); 
    ig0=ig0+m;
end
j=1; 
ig = ig0+[1:m]; 
q(1:end-1,j) = ig; 
q(end,j)=q(1,j); 

% q', pause

N1=N+1; Q=zeros(N1*N1,E);
j0=1; j1=j0+N;
for ey=1:Ey
    i0=1; i1=i0+N;
    for ex=1:Ex
        Q(:,(ey-1)*Ex+ex)=reshape(q(i0:i1,j0:j1),N1*N1,1);
        i0=i0+N; i1=i0+N;
        %   qq=reshape(Q(:,e),N1,N1)
    end
    j0=j0+N; j1=j0+N;
end
%pause
q=reshape(Q,N1*N1*E,1);
nL = E*N1*N1; nb=max(max(q)); % nb = nbar
Q=spalloc(nL,nb,nL);       % sparse Q matrix
for il=1:nL
    ig=q(il); 
    Q(il,ig)=1; 
end

% Make P fo

