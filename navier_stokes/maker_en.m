function[R] =  maker_en(Q,Ex,N,N_en_y)

%
%  Make the restriction matrix
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
%    |  e=1    |  e=2    |  e=3    |  ...    |  e=Ex   P
%    |         |         |         |         |         P
%    +---------+---------+---------+---------+---------P
%    |         |         |         |         |         P <-- peiodic
%    |         |         |         |         |         P
%    |  e=Ex+1 |  e=Ex+2 |  e=Ex+3 |  ...    |  e=2*Ex P
%    |         |         |         |         |         P
%    +---------+---------+---------+---------+---------P
%        ^
%        |
%        +--- numbered last
%


[nL,nb]=size(Q); 
n_dirichlet = 2*Ex*N;
R=speye(nb); 
R=R(1:(nb-n_dirichlet),:);

