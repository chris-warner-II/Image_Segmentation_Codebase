function [q, qc, ind, clust_size] = compute_Modularity_scalar(W,s)


disp('NOTE: This isnt working.  Using modularity function in graph dir in CodeDownloads')


% syntax: [q] = compute_Modularity_scalar(W,s);
%
% This function quantifies how modular a segmentation (s) of a network (W)
% is or how well a segmentation captures the community structure within a
% network. 
%
% Scalar modularity metric (q) "measures the fraction of within-community 
% edges minus the expected value of the same quantity in a network with the 
% same community divisions but random connections between the vertices." 
%                             -Newman Phys Rev E 2004 
%
% It is the total weight of within-community edges minus the total weight
% of across-community edges, with communities defined by the segmentation (s).
%
% It can be used to rate and compare different segmentations of the
% same network.
%      q = 1 for a network with all weights within-community.
%      q = 0 for completely random network where within-community weights 
%            are equal to across-community weights.
%      q < 0 if weights across-community are larger than weights
%            within-communities. This is possible.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% INPUT: (W,s)
%    W - N^2 x N^2 Adjacency Matrix using both distance-dependent connectivity 
%        constraints and pixel differences.
%    s - a vector of length N that has k unique values. It represents the 
%        segmentation of partition of an N node network into k communities.
%
% OUTPUT: [q]
%    q - Modularity scalar value. 
%

%% Communities defined from segmentation vector (s)
N = numel(s);
communities = unique(s);
k = numel(communities);



% %% If there is only a single community, set modularity to zero (because this is the uninteresting segmentation).
% if(k==1)
%     disp('only one cluster. uninteresting.')
%     q = 0;
%     clust_size = [N, 0];
%     ind = [1:N];
%     return
% end




%% Build the "e matrix" - a k x k matrix where entries denote fraction of edges in the network that connect community i to community j.
% e = zeros(k); 
% for i = 1:k
%     comm_i = find(s==communities(i));
%     for j = i:k
%         comm_j = find(s==communities(j));
%         e(i,j) = sum(sum(W(comm_i,comm_j))) ./ sum(W(:));
%         if i~=j
%             e(j,i) = e(i,j);
%         end
%         
%         keyboard
%         
%     end
% end
% 
% 
% %% Modularity - Fraction of edges in whole network within communities minus fraction of edges crossing communities.
% q = trace(e) - sum(sum(triu(e-diag(diag(e)))));




m = sum(W(:));
m_wt = numel(W) - numel(s);
%
e = zeros(1,k); 
a = zeros(1,k); 
for i = 1:k
    comm_i = find(s==communities(i));
    not_comm_i = find(s~=communities(i));
    e(i) = sum(sum( W(comm_i,comm_i) ));
    e_wt(i) = numel(comm_i).^2 - numel(comm_i);
    a(i) = 2.*sum(sum( W(comm_i,:) ));
    a_wt(i) = 2.*numel(comm_i).*numel(s);
       % keyboard
end

% NOTE:  THESE ARE MY OWN VERSIONS OF SOMETHING LIKE MODULARITY...
q = sum( e - a ) ./ m

qb = sum( e - a.^2 ) ./ (m)

qc = sum( (e./e_wt) - (a./a_wt) )


%% Record the number of nodes in each cluster (clust_size) and the identity of the nodes in each cluster (ind)
clust_size = zeros(1,k);
ind = zeros(1,N);
clust = unique(s);
st = 1;
for i = 1:numel(clust)
    members = find(s==clust(i));
    nd = st-1+numel(members);
    ind(st:nd) = find(s==clust(i));
    clust_size(i) = numel(members);
    st = nd+1;
end
%
clust_size = sort(clust_size,'descend');



