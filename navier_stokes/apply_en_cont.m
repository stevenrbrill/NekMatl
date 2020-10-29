function out = apply_en_cont(M,nodes)
% M is the full size matrix

out = zeros(size(M(:,1)));

for i = 1:length(nodes)
    out = out + M(:,nodes(i));
end