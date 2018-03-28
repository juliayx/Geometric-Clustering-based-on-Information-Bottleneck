function [ICI, IXC] = cal_information(pI,pC,pC_I,pIX,pX_C)
%Input:
%       pI: p(i)
%       pC: p(C)
%       pC_I: p(c|i)
%       pIX: p(i,x)
%       pX_C: p(x|c)
%Output:
%       ICI: I(c,i)
%       IXC: I(x,c)

iDim = size(pIX,1);
xDim = size(pIX,2);
cDim = size(pC,1);

pCI = pC_I .* repmat(pI', [cDim 1]);
pX = (sum(pIX, 1))';
pXC = pX_C .* repmat(pC', [xDim 1]);
ICI = DKL2(pCI,pC * pI');
IXC = DKL2(pXC,pX * pC');

end

