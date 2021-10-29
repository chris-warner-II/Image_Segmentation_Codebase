function [dprime] = compute_dprime_metric(V,s)

% syntax: [dprime] = compute_dprime_metric(V,s);
%
% This function ...
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% INPUT: (V,s)
%    V - a vector of length N that is the eigenvector of one method'smatrix.
%    s - a vector of length N that has k unique values. It represents the 
%        segmentation of partition of an N node network into k communities.
%
% OUTPUT: [dprime]
%    dprime - dprime or sensitivity metric.
%





%%

% Compute sensitivity index (d') with mu and sig for 2 classes from eigenvector.
k = numel(unique(s));

if(k~=2)
    disp('The dprime metric function is not equipped to handle anything other than bipartitions.')
    dprime=0;
    return
end

for i = 1:k

    ind = find(s==i);
    %
    mu0(i) = mean(V(ind));
    sig0(i) = std(V(ind));

end

dprime = abs( ( mu0(1) - mu0(2) )./sqrt( 0.5.*( sig0(1).^2 +  sig0(2).^2 ) ) );