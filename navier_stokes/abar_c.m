function[w] = abar_c(u,Dh,G,Q)

E=size(G,4); N1=size(G,1); N=N1-1;

uL=reshape(u,N1,N1,E);
wL=aloc(uL,Dh,G);
wL=reshape(wL,N1*N1*E,1);
w=wL;
