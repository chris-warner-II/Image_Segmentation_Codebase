function [ Bm ] = PIV_vs_DiagDistB(W,Wdist,Mask)

% syntax: [ B ] = PIV_vs_Dist(W,Wdist,Mask)
%
% INPUT: [W & Wdist]
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

disp(['Inside PIV_vs_DiagDist function : Computing Normalized Diagonal Sums : '])
tic

imsz = size(W,1); % number of pixels in image.
pctNN = max(sum(logical(Wdist)))/(imsz-1);   % ratio of number of connections by number of possible connections

B  = spalloc(imsz,1,ceil(pctNN*(imsz-1)));    % Vector of diagonal sums of W matrix (normalized by Wdist to make it avg.)
%B2  = spalloc(imsz,1,ceil(pctNN*(imsz-1))); 
Bm = spalloc(imsz,imsz,ceil(pctNN*numel(W)));

%% Build up diagonal sums - connectivity of nodes separated by given distance and direction.
disp(['Iteration # out of ',num2str(imsz)])
for i=1:imsz
    
    fprintf('%s',[num2str(i),' '])
    
    if ( sum(diag(Wdist,i-1)) ~= 0 ) % whenever I wont be dividing by zero (ie. whenever 2 nodes are possibly connected)
        B(i) = sum(diag(W,i-1)) ./ sum(diag(Wdist,i-1));       % sum of each diagonal in Adjacency matrix (for this image) ./
                                                                % maximum possible diagonal sum from adjacency in uniform image
        %B2(i) = mean(diag(W,i-1));                                                        
    end
end
fprintf('\n')


% Bm is the B vector packed into a matrix where each element in matrix tells average connectivity of pixels separated by 
% distance that separates pixels corresponding to matrix element.
nonzeroBs = find(B);
%
disp(['Converting Vector (B) to Matrix (Bm): Number nonzero elements = ',num2str(numel(nonzeroBs))])
disp(['Should be 1/2 number of nearest neighbors in Rmax. If Rmax = 4, then 48/2 = 24 or 48, not sure.'])
%
if( numel(nonzeroBs) > ceil(pctNN*(imsz-1)) )
    disp('Something wrong with number of non-zero elements in Diagonal Sums computation.')
    keyboard
end
%
for i=1:numel(nonzeroBs)
    disp(['Iteration # ',num2str(i),' / ',num2str(numel(nonzeroBs)),' : ',num2str(nonzeroBs(i))])
    Bm = Bm + diag(B(nonzeroBs(i))*ones(imsz+1-nonzeroBs(i),1),   nonzeroBs(i)-1) ...
            + diag(B(nonzeroBs(i))*ones(imsz+1-nonzeroBs(i),1), -(nonzeroBs(i)-1));
end



% NOTE: I should use Mask or at least compare it to Bm, they should be the same!




toc


% Check memory usage.
disp('Memory Usage in PIV_vs_DiagDistB:')
% whos
check_memory_usage