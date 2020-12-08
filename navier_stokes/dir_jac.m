function [L_x,L_y,J_x,J_y] = dir_jac(E,X,Y)

L_x = zeros(E,1);
L_y = zeros(E,1);
for ie=1:E
    X_min = min(min(X(:,:,ie)));
    X_max = max(max(X(:,:,ie)));
    Y_min = min(min(Y(:,:,ie)));
    Y_max = max(max(Y(:,:,ie)));
    L_x(ie) = X_max-X_min;
    L_y(ie) = Y_max-Y_min; 
end
J_x = 2./L_x;
J_y = 2./L_y;