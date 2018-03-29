function [indices] = max_indx_row(A)

M = max(A,[],2);
indices = zeros(size(A,1),1);
for i = 1:size(A,1)
    indices(i) = find(A(i,:) == M(i));
end

end