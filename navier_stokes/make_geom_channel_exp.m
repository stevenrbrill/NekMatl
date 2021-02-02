function[X,Y] =  make_geom_channel_exp(Ex,Ey,N)

E = Ex*Ey;

[z,w]=zwgll(N);

N1=N+1;

[X0,Y0]=ndgrid(z,z); X=zeros(N1,N1,E); Y=X;

N_side = Ey/2;
max_side = 0;
if mod(Ey,2) == 1
    N_side = (Ey-1)/2;
    max_side = N_side*(2/Ey);
end

Pr = 0.55;
a = 1./((1-Pr.^N_side)./(1-Pr));
f_geo = @(n) a*(Pr).^n;
dist_geo = f_geo(linspace(0,N_side-1,N_side));
nodes_y = zeros(Ey+1,1);
nodes_y(1) = -1;
pt_geo = zeros(size(dist_geo));
for i = 1:N_side
    for j = 1:i
        pt_geo(i) = pt_geo(i) + dist_geo(j);
    end
end
% for ey = 1:Ey
%     nodes_y(ey+1) = nodes_y(ey) + 2/Ey;
% end
if mod(Ey,2) == 0
    nodes_y(1) = -1;
    for ey = 1:N_side
        nodes_y(ey) = -pt_geo(N_side-ey+1);
    end
%     nodes_y(N_side+1) = 0;
    for ey = 1:N_side-1
        nodes_y(N_side+ey+1) = pt_geo(ey);
    end
    nodes_y(end) = 1;
else
    error("Need even Ey");
end


y0=0; y1=y0+1;
for ey=1:Ey
    x0=0; x1=x0+1; 
    for ex=1:Ex;  % Domain of unit [0,1] squares
        X(:,:,ex+Ex*(ey-1)) = x0+1*(x1-x0)*(X0+1)/2;
        x0=x0+1; x1=x0+1;
        Y(:,:,ex+Ex*(ey-1)) = nodes_y(ey)+1*(nodes_y(ey+1)-nodes_y(ey))*(Y0+1)/2;
    end;
    y0=y0+1; y1=y0+1;
end;


xM=glmax(X); 
xm=glmin(X); 
s=0.375/(xM-xm); %2*pi/(xM-xm); 
X=s*X; %  X on [0,2pi]

% yM=glmax(Y); 
% ym=glmin(Y); 
% s=2/(yM-ym); 
% Y=s*Y-1; %  Y on [-1,1]

