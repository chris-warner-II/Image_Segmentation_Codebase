function [MC] = metaClusterAnalysis(theta, netParams, kurParams, circ)


% syntax: [MC] = metaClusterAnalysis(theta, netParams, kurParams, circ);
%
% This function does a meta analysis on how clusterable oscillators are in
% phase space (or how clusterable pixels are in original image if strawman).
% Using the ground truth, we calculate Divisive Margine - ratio of average
% pairwise distance between elements within cluster vs. across cluster
% boundaries.
%
% Input:  theta     - NxT matrix that contains position of all N
%                     oscillators at all T timesteps.  If strawman, T=1.
%         netParams - parameters that construct network from image.
%         kurParams - parameters that determines kuramoto simulation.
%         circ      - a flag whether theta is a circular variable (phase)
%                     or not (pixels with strawman)
%
% Output: MC        - a structure that contains for each ground truth
%                     segmentation, DistAvgPW - this is a 2xT vector. It
%                     has avg pairwise distance between oscillators in the
%                     same cluster & across cluster boundaries for each
%                     timepoint.





% Extract data from kurParams
vars = fieldnames(kurParams);
for i=1:numel(vars)
    eval([vars{i},'=kurParams.',vars{i},';'])
end

% Extract data from netParams
vars = fieldnames(netParams);
for i=1:numel(vars)
    eval([vars{i},'=netParams.',vars{i},';'])
end


colors = 'rgbkcmyrgbkcmyrgbkcmyrgbkcmyrgbkcmyrgbkcmyrgbkcmyrgbkcmyrgbkcmyrgbkcmyrgbkcmy';

if(size(theta,1)==T)
    theta = theta'; % flip theta because all these functions want it the other way..
else
    T = 1;
    spp = 1;
    
    % The T=1 statement is for when you are running on eigenvector.  Not sure
    % how to scale eigenvector when I calculate these thing to make it
    % comparable to phase in kuramoto.  Positive and negative should have same
    % meaning but wont know about scale.  Can compare Evec to self and Kur to
    % self, but not across to one another.
    %
    % Actually, T can be 1 when dealing with image as theta (for strawman
    % model) too.

end


if ~circ
    theta = theta.*pi; % spread out "oscillators" so that 0 -> 0 and 1 -> pi.
                       % a first cut, but it seems to be working.
end



%% NOTE : Need to Loop Thru Ground Truths for Image Patches.
numSegmentations = numel(gT);

for S = 1:numSegmentations

    %% Find size of each cluster (number of nodes/pixels/oscillators).  Can use to do weighted average of cluster extent below.
    Clusters = unique(gT{S}(:));
    
    % C = numel(Clusters);
    % C will be different for each Segmentation.  Can we figure that out earlier on?
    
    ClusterSize = zeros(1,C(S));
    for i = 1:C(S); % loop through clusters
        ClusterSize(i) = numel(find(gT{S}==Clusters(i)));
    end




    %% Calculate Subtractive Relative Margin - Compare average distance between oscillator pairs inside clusters vs across clusters.
    %disp('Meta Cluster Analysis: Computing average pairwise distance within and across clusters')
    %tic
    DistAvgPW = zeros(ceil(T/spp)+1,2);
    AvgPairwiseDist = zeros(C(S),ceil(T/spp)+1,2);

    for i = 1:C(S) % loop through clusters

        ind_in = find(gT{S}==Clusters(i));
        ind_out = find(gT{S}~=Clusters(i));
        
        
        % Pairwise Distances at Phase Initialization
        dists_in = circ_dist(repmat(theta(ind_in,1),1,numel(ind_in)),repmat(theta(ind_in,1)',numel(ind_in),1));
        x = abs(dists_in(dists_in~=0));
        if isempty(x)
            x=0;
        end
        %
        AvgPairwiseDist(i,1,1) = circ_mean(x(:)); % avg distance between oscillators within same cluster
        %
        dists_out = circ_dist(repmat(theta(ind_in,1),1,numel(ind_out)),repmat(theta(ind_out,1)',numel(ind_in),1));
        y = abs(dists_out(dists_out~=0));
        if isempty(y)
            y=0;
        end
        %
        AvgPairwiseDist(i,1,2) = circ_mean(y(:)); % avg distance between oscillators in different clusters
        
        % Pairwise Distance as Time Advances
        for k = 1:T/spp
            %
            dists_in = circ_dist(repmat(theta(ind_in,spp*k),1,numel(ind_in)),repmat(theta(ind_in,spp*k)',numel(ind_in),1));
            x = abs(dists_in(dists_in~=0));
            if isempty(x)
                x=0;
            end
            %
            AvgPairwiseDist(i,k+1,1) = circ_mean(x(:)); % avg distance between oscillators within same cluster
            %
            dists_out = circ_dist(repmat(theta(ind_in,spp*k),1,numel(ind_out)),repmat(theta(ind_out,spp*k)',numel(ind_in),1));
            y = abs(dists_out(dists_out~=0));
            if isempty(y)
                y=0;
            end
            %
            AvgPairwiseDist(i,k+1,2) = circ_mean(y(:)); % avg distance between oscillators in different clusters
            %
        end

        DistAvgPW(:,1) = DistAvgPW(:,1) + AvgPairwiseDist(i,:,1)'.*ClusterSize(i); % weighting by cluster size for within cluster distance
        DistAvgPW(:,2) = DistAvgPW(:,2) + AvgPairwiseDist(i,:,2)'.*ClusterSize(i); % weighting by cluster size for across cluster distance

    end

    DistAvgPW = DistAvgPW./N; % finish weighting Pairwise distance between clusters by size of cluster.
    %toc

    % NOTE:  I Will Have To Normalize This Properly And Determine if Divisive or Subtractive Comparison is Better.
    % NOTE:  I Will Have To Write How To Do This For Eigenvectors (Non-Circular Calculations)

    if(0)

        figure, 
        subplot(211), hold on, 
        plot(squeeze(AvgPairwiseDist(1,:,1)),'b'),  
        plot(squeeze(AvgPairwiseDist(2,:,1)),'r')
        plot(DistAvgPW(:,1),'g')


        subplot(212), hold on, 
        plot(squeeze(AvgPairwiseDist(1,:,2)),'b'),  
        plot(squeeze(AvgPairwiseDist(2,:,2)),'r')
        plot(DistAvgPW(:,2),'g--')

        figure, plot(DistAvgPW)

    end


    %% Package it all into a data structure to be sent to calling function
    MC{S}.DistAvgPW = DistAvgPW;

end

