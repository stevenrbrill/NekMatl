function[X,Y] =  make_geom_channel(Ex,Ey,N)

E = Ex*Ey;

[z,w]=zwgll(N);

N1=N+1;

[X0,Y0]=ndgrid(z,z); X=zeros(N1,N1,E); Y=X;

y0=0; y1=y0+1;
for ey=1:Ey
    x0=0; x1=x0+1; 
    for ex=1:Ex;  % Domain of unit [0,1] squares
        X(:,:,ex+Ex*(ey-1)) = x0+1*(x1-x0)*(X0+1)/2;
        x0=x0+1; x1=x0+1;
        Y(:,:,ex+Ex*(ey-1)) = y0+1*(y1-y0)*(Y0+1)/2;
    end;
    y0=y0+1; y1=y0+1;
end;


xM=glmax(X); 
xm=glmin(X); 
s=0.375/(xM-xm); %2*pi/(xM-xm); 
X=s*X; %  X on [0,2pi]

yM=glmax(Y); 
ym=glmin(Y); 
s=2/(yM-ym); 
Y=s*Y-1; %  Y on [-1,1]



% %  Map to circle, starting at theta=pi/2
% xmax=glmax(X); xmin=glmin(X);
% ymax=glmax(Y); ymin=glmin(Y);
% theta = (pi/2)-(2*pi)*(X-xmin)/(xmax-xmin);
% rad   = 0.5 + (Y-ymin)/(ymax-ymin);
% X=rad.*cos(theta); Y=rad.*sin(theta);

%% 
%% Below just for diagnostics
%% 
%% F=X.*X-Y.*Y; fmx=max(max(max(F))); fmn=min(min(min(F)));
%% 
%% hold off
%% for e=1:E; 
%%    mesh(X(:,:,e),Y(:,:,e),F(:,:,e)); hold on
%%    contour(X(:,:,e),Y(:,:,e),F(:,:,e),[fmn:.1:fmx]); hold on
%% end;

