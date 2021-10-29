X=rand(101,101);

mask = ones(size(X));
mask(60,2:7) = 0;

% processing step
Y = imfilter(X,ones(3,3)/9,'replicate');
Y(find(~mask)) = X(find(~mask));