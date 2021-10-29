function [Q] = compute_ModularityB(W,Wdist,Mask,imagin,topo,maskFlg,rmax)

% syntax: [Q] = compute_ModularityB(W,Wdist,Mask,imagin,topo,maskFlg,rmax);
%
% This function ... 
%
% INPUT: (W,Wdist,Mask,ximg,yimg,topo,maskFlg,rmax)
%    W       - N^2 x N^2 Adjacency Matrix using both distance-dependent connectivity constraints and pixel differences.
%    Wdist   - N^2 x N^2 "Hardware" Connectivity constraints imposed by programmer - Distance dependence of network.
%    Mask    - structure which contains a mask of all pixel pairs in adjacency separated by a each distance less than rmax.
%              [ Mask(i).distance = the distance (scalar value).                                                           ].
%              [ Mask(i).mask     = the binary mask (is N^2 x N^2 like W matrix).                                          ].
%    imagin - input image that was used to compute W, etc.
%    topo    - a flag whether to do topographic modularity (topo=1) or nontopographic (topo-0)
%    maskFlg - a flag to choose which mask to compute. (0=No Mask. 1=Mask of pairwise pixel distances in image. 2=Not Working/Dont Use.)
%    rmax    - maximum distance across which pixels/nodes in network can be connected. Implemented by heaviside/step function.
%
% OUTPUT: [Q]
%    Q       - Modularity Matrix. Compute eigenvectors or plug into Kuramoto coupled oscillator simulation to get segmentation.





disp('In compute_ModularityB : Before computing Null Model : ')


ximg = size(imagin,1); 
yimg = size(imagin,2);

if isinf(rmax)
    nNeighbors = ximg*yimg;
else
    nNeighbors = sum(2*(1:rmax)*4); % this is an overestimate right now (takes everything in a box)
    % pctNN = max(sum(logical(Wdist)))/(ximg*yimg-1);   % ratio of number of connections by number of possible connections
	% ceil(pctNN*numel(Wdist))
    % TODO: For more efficient code, we can pare it down from 80 to 48 for Rmax = 4.
end

% preallocate modularity matrix (will be mostly zeros because of nearest neighbor connectivity)
Q = spalloc( ximg*yimg , ximg*yimg, nNeighbors*ximg*yimg );


[NM] = compute_null_modelB(W,Wdist,Mask,imagin,topo,maskFlg); % Compute Null Model
    
disp('Back out in compute_ModularityB : After Null Model has been computed')

% whos
check_memory_usage

disp('In compute_ModularityB : Here I do a full matrix subtraction.')
tic
Q = W - NM; % Compute Modularity from Adjacency and Null Model
toc




%% Some Analysis (Looking at Weights, Null Model & Modularity Matrices vs Image)
if(0)
    figure, 
    subplot(131), imagesc(W), colorbar, title('A'), axis square
    subplot(132), imagesc(NM), colorbar, title('N'), axis square
    subplot(133), imagesc(Q), colorbar, title('Q'), axis square

    Qf = full(Q);
    Wf = full(W);
    Nf = full(NM);

    Wfn = Wf.*(Qf<0);
    Nfn = Nf.*(Qf<0);

    figure,
    subplot(221), imagesc(Qf.*(Qf>0)), colorbar, title('Q>0')
    subplot(222), imagesc(Qf.*(Qf<0)), colorbar, title('Q<0')
    subplot(223), imagesc(Nfn, [min(min(Nfn(:)),min(Wfn(:))),max(max(Nfn(:)),max(Wfn(:)))]), colorbar, title('N(Q<0)')
    subplot(224), imagesc(Wfn, [min(min(Nfn(:)),min(Wfn(:))),max(max(Nfn(:)),max(Wfn(:)))]), colorbar, title('W(Q<0)')

    keyboard
end


