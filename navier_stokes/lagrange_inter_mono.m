function m = lagrange_inter_mono(x,y)
    % function m = lagrange_inter_mono(x,y)
    % Computes the monomial coefficients for the Lagrange interpolating
    % polynomial that interpolates the values y at the points x
    % m(1)x^p + m(2)x^(p-1) + ...
    % polyval(m,pts) to evaluate the polynomial
    % x - node locations
    % y - function values at nodes

    n = length(x);
    
    mat = zeros(n);
    for i = 1:n
        for j = 1:n
            mat(i,n-j+1) = x(i)^(j-1);
        end
    end
    
    m = mat\y;    
end