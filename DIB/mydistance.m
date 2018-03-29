function [distance] = mydistance(X , Y)

%Input:
%   row vector X
%   row vector Y
%Output:
%   distance of (X,Y)
distance = sum((X-Y).^2);
end
