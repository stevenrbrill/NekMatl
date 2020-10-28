function[w,w_p] = abar_en(u,Dh,G,Q,P)

E=size(G,4); N1=size(G,1); N=N1-1;

uL=reshape(Q*u,N1,N1,E);
wL=aloc(uL,Dh,G);
wL=reshape(wL,N1*N1*E,1);
w=Q'*wL;


uL_p=reshape(P*u,N1,N1,E);
wL_p=aloc(uL_p,Dh,G);
wL_p=reshape(wL_p,N1*N1*E,1);
w_p=P'*wL_p;