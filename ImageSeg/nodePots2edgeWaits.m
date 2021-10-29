function [Wt] = nodePots2edgeWaits(Pt)

% syntax: Wt = nodePots2edgeWaits(Pt);
%
%   Input: Pt = a vector of node potentials calculated iteratively in findMLE
%   Outpt: Wt = the maximum entropy distribution of edge weights on graph  
%  
%       Wt_{ij} = 1 / ( 1 + exp( Pt_i + Pt_j ) )
%
%   These Edge Weights (Wt) will be used as the null model for Modularity.

% calculate equation above
adtn = repmat(Pt,1,numel(Pt));
adtn = adtn + adtn';
Wt = 1./(exp(adtn)+1);

% zero diagonals (no self-weights)
Wt = Wt.*(1-eye(size(Wt)));
