function Px = Project_1s(x)

% syntax: Px = Project_1s(x);
%
% This function takes any vector input x and projects it to a space
% orthogonal to the all 1's vector.  It is to be used in the Power Method
% for the Negative Graph Laplacian to solve for the second largest
% eigenvector because the largest eigenvector is known to be the all 1's
% vector.

v = ones(size(x));                   % the 1's vector
P = eye(numel(x)) - (v*v')/(v'*v);   % projection matrix orthogonal to v       
Px = P*x;                  % x vector projected orthogonal to 1's vector

% If you want to check that Px is orthogonal to v, uncomment this line.
%dot(Px,v) % should be 0.