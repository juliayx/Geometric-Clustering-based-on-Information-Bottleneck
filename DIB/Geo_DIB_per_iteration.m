function [pC_I, pX_C, pC, XC,L] = Geo_IB_per_iteration(pIX, beta, p0X_C,p0C,X)
%Input:
%       pIX: Joint distribution p(i,x)
%       beta: coefficient of model
%       p0X_C: p0(x|c), initial value
%       p0C: p0(c), initial value
%       X: Location of testing points
%Output:
%       pC_I: p(c|i)
%       pX_C: p(x|c)
%       pC: p(c)
%       XC: location of cluster points

iDim = size(pIX,1);
xDim = size(pIX,2);
cDim = size(p0C,1);
pI = sum(pIX, 2);
pX_I = (pIX ./ repmat(pI, [1, xDim]))';


%Bayes:
%p0I_C = ((p0C_I .* repmat(pI', [cDim 1])) ./ (repmat(p0C, [1 iDim])))';

%pX_I = (pIX./repmat(pI, [1 xDim]))';
%p0X_C = pX_I * p0I_C;

% distance_C_I = zeros(cDim, iDim);
% for i = 1:cDim
%     for j = 1:iDim
%         distance_C_I(i,j) = mydistance(pX_I,p0X_C);
%     end
% end
%d_I_C : d(i,c) = DKL(p(x|i) || p(x|c))
d_I_C = repmat(sum(pX_I.*log2(pX_I),1)', [1 cDim]) - pX_I' * log2(p0X_C);
%L_I_C: l(i,c) = log(p0(c)) - beta * d(i,c)
L_I_C = repmat(p0C', [iDim 1]) - beta * d_I_C;
%c_I: c(i) = argmax_c l(i,c)
c_I = max_indx_row(L_I_C);
pC_I = zeros(cDim,iDim);
for j = 1:iDim
    pC_I(c_I(j),j) = 1;
end


% unnorm_p_C_I = repmat(p0C, [1,iDim]) .* exp(1*beta * (log2(p0X_C))' * pX_I);
% Z_I_beta = repmat(sum(unnorm_p_C_I), [cDim 1]);
% pC_I = unnorm_p_C_I ./ Z_I_beta;

pC = pC_I*pI;

pI_C = ((pC_I .* repmat(pI', [cDim 1])) ./ (repmat(pC, [1 iDim])))';
pX_C = pX_I * pI_C;

XC = X * pX_C;

[~, IXC, HC] = cal_information(pI,pC,pC_I,pIX,pX_C);
L = HC - beta * IXC;


end


