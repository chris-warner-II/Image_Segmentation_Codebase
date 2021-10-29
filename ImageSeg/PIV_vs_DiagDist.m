function [ Bm ] = PIV_vs_DiagDist(W,W_conn)

% syntax: [ B ] = PIV_vs_Dist(W,W_conn)
%
% INPUT: [W & W_conn]
% This function takes in an Adjacency Matrix and a "Hardware" Connectivity
% Matrix (which is equivalent to an Adjacency Matrix calculated on a 
% uniform/null image).  
%
% OUTPUT: [B]
% It outputs a vector B which is the Expected Weight Value (averaged over
% the entire image) of pixels separated by a distance and direction given
% by the distance from the Identity Diagonal of the matrix. The first entry
% in the B vector should be zero because the Adjacency Matrix is
% constrained to have no self-connectivity and therefore the diagonal is
% full of zero weights.  Each B entry is calculated by the ratio of the sum 
% of entries in a diagonal of the Adjacency Matrix and the sum of entries
% in the corresponding diagonal of the Connectivity Matrix.  Each B entry
% take a continuous value between 0 and 1.  The B vector is   

%% Allocate Sparse Matrices for data
imsz = size(W,1); % number of pixels in image.
pctNN = max(sum(logical(W_conn)))/(imsz-1);   % ratio of number of connections by number of possible connections

B=spalloc(imsz,1,ceil(pctNN*(imsz-1)));
Bn=spalloc(imsz,1,ceil(pctNN*(imsz-1)));
Bm = spalloc(imsz,imsz,ceil(pctNN*numel(W)));

%% Build up diagonal sums - connectivity of nodes separated by given distance and direction.
for i=1:imsz
    B(i) = sum(diag(W,i-1));       % sum of each diagonal in Adjacency matrix (for this image)
    Bn(i) = sum(diag(W_conn,i-1)); % maximum possible diagonal sum from adjacency in uniform image
end
B = B./Bn;     % The expected weight between pixels with given distance / direction of separation
B(isnan(B))=0; % Because you probably divided by zero in previous step.

for i=1:numel(B)
    Bm = Bm + diag(B(i)*ones(imsz+1-i,1),i-1) + diag(B(i)*ones(imsz+1-i,1),-(i-1));
end