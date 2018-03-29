function [DKL_sum] = DKL2(p,q)

p = p(:);
q = q(:);

DKL_sum = sum(p.*log2(p./q));