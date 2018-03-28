function [ pC_I, pX_C, PC, XC ] = Geo_all_iteration( pIX, beta, p0X_C,p0C,X,maxiter, tol )
%Input: 
%       pIX: p(i,x)
%       beta: coefficient
%       p0X_C: p0(x|c)
%       p0C: p0(c)
%       X: input matrix
%       maxiter: maximum iteration
%       tolerance
%Output:
%       pC_I: p(c|i)
%       pX_C: p(x|c)
%       pC: p(c)
%       XC: final cluster point
[pC_I, pX_C, PC,XC,L] = Geo_IB_per_iteration(pIX, beta, p0X_C,p0C,X);
for iter = 1:maxiter
    pX_C_pre = pX_C;
    pC_pre = PC;
    L_pre = L;
    [pC_I, pX_C, PC, XC,L] = Geo_IB_per_iteration(pIX, beta, pX_C_pre,pC_pre,X);
    difference = L_pre - L;
    fprintf('iter: %3i  difference = %6i\n',iter,difference);
    if difference < tol, break; end
end

    
    



end

