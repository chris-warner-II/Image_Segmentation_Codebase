%% Structure of parameters used to construct a HandCooked Network. 
% Parameters to build Coupling Matrix K1     





%% This section is from the 
if ~exist('netParams','var')

    %
    for i = 1:size(GT,1)
        netParams.gndTruth{i} = GT(i,:);           % mask of which cluster oscillators belong to - can be hierarchical
    end
    netParams.Rmax = inf;                          % max distance between oscillators over which they can have a connection.
                                                   % Note: set rmax to inf (or N) if there are no topological connectivity constraints  
    netParams.PconnFar = 0;                        % Probability of a connection beyond rmax (could make distance dependent too).
    netParams.Ndims = size(netParams.gndTruth{1}); % Number of oscillators in x, y (,and maybe z) dimension.
    netParams.N = N;                               % Number of oscillators is product of number in each dimension.                                
    netParams.C = num_clust;                       % Number of clusters.   
    netParamsCsize = sizeC_GT;                     % Size of each cluster in each level of hierarchy
    %
    
    
    % FOR HANDCRAFTED HIERARCHICAL CLUSTER NETWORKS, NEED TO CHANGE THESE !!
    
    
    netParams.Weak = embed_clusters_params.mu(2);        % Mean Coupling weights (or p if Bin) across Clusters. 
    netParams.Strng = embed_clusters_params.mu(1);       % Mean Coupling Weights (or p if Bin) within Clusters.
    %
    try
        netParams.sigWk = embed_clusters_params.sig(2);  % Gaussian Std of Coupling Weights across Clusters.  
        netParams.sigSg = embed_clusters_params.sig(1);  % GaussianStd of Coupling Weights within Clusters. 
    catch
        netParams.sigWk = 0;                             % Binomial network has no std. 
        netParams.sigSg = 0;                             % Binomial network has no std.
    end
    %
    
else
    
    %disp('Not creating: HCN - netParams already exists.')
    %netParams
    
end
