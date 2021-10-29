% A script to test projection matrices and how to project a vector onto a
% vector (the all 1's vector) and orthogonal to it.  It works.



dim = 10               % dimension / size of vectors & matrices

v = ones(dim,1)         % the 1's vector
P_on = (v*v')/(v'*v)   % projection matrix onto v

x = rand(dim,1)  % some vector we want to project         

P_off = eye(dim) - P_on % projection matrix perpendicular to 1's (= x -P*x)

x_on = P_on*x    % x vector projected onto 1's vector
x_off = P_off*x  % x vector projected orthogonal to 1's vector

dot(x_on,v)  % should be 1.
dot(x_off,v) % should be 0.