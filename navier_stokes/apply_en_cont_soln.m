function u = apply_en_cont_soln(u,nodes,psi)
% M is the full size matrix

for i = 1:length(nodes)
    u(nodes(i)) = u(nodes(i)) + psi;
end