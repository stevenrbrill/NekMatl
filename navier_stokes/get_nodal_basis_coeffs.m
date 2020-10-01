function basis_coeffs = get_nodal_basis_coeffs(soln_pts)
    % basis_coeffs = get_nodal_basis_coeffs(soln_pts)
    % Generates a nodal basis at nodes soln_pts(P+1,N)
    % soln_pts - Nodes to create the basis around

    [n_bases,N] = size(soln_pts);

    basis_coeffs = zeros(n_bases,n_bases,N);
    for n = 1:N
        e = eye(n_bases);
        for i = 1:n_bases
            basis_coeffs(i,:,n) = lagrange_inter_mono(soln_pts(:,n),e(:,i))';
        end
    end
end