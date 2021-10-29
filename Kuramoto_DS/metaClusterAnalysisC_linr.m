function [MC] = metaClusterAnalysisC_linr(theta, netParams, kurParams)


% syntax: [MC] = metaClusterAnalysisB(theta, netParams, kurParams, circ);
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


disp(['Inside metaClusterAnalysis (Linear Variable) function : '])

tic


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


T = 0;   % no concept of time when running on eigenvectors or on image pixels (strawman)
spp = 1; % 'samples per period' not useful for eigenvectors or on image pixels (strawman)
theta = theta'; % theta structure should be (# units) x (# dims per unit)



thresholdNodes2DoBatch = 2.5e7; % if a matrix has beyond this many elements, do batch processing on it.
batchSize = 1000;               % number of oscillators to process at one time in single batch.


%% NOTE : Need to Loop Thru Ground Truths for Image Patches.
numSegmentations = numel(gT);

for S = 1:numSegmentations
    
    disp(['Segmentation # ',num2str(S),' / ',num2str(numSegmentations)])

    %% Find size of each cluster (number of nodes/pixels/oscillators).  Can use to do weighted average of cluster extent below.
    Clusters = unique(gT{S}(:));
    
    % C = numel(Clusters);
    % C will be different for each Segmentation.  Can we figure that out earlier on?
    
    ClusterSize = zeros(1,C(S));
    for i = 1:C(S); % loop through clusters
        ClusterSize(i) = numel(find(gT{S}==Clusters(i)));
    end




    %% Calculate Divisive Relative Margin - Compare average distance between oscillator pairs inside clusters vs across clusters.
    %disp('Meta Cluster Analysis: Computing average pairwise distance within and across clusters')
    %tic
    DistAvgPW = zeros(ceil(T/spp)+1,2);
    AvgPairwiseDist = zeros(C(S),ceil(T/spp)+1,2);

    for i = 1:C(S) % loop through clusters
        
        disp(['Cluster # ',num2str(i),' / ',num2str( C(S) )])

        ind_in = find(gT{S}==Clusters(i));  % Nodes inside this cluster
        ind_out = find(gT{S}~=Clusters(i)); % Nodes outside this cluster.
        

        
        % Nodes within the same cluster - Average Pairwise Distance of (Numerator of Divisive Margin)
        if( numel(ind_in)*numel(ind_in) < thresholdNodes2DoBatch )

            % Pairwise Euclidian Distances (of linear variable) of either image pixels or Eigenvector intensities.
            dists_in = dist( theta(:,ind_in)', theta(:,ind_in) );
            idx = eye(size(dists_in)); 
            x = abs(dists_in(~idx)); % get rid of diagonal (distance between unit and itself)
            if isempty(x)
                x=0;
            end
            %
            AvgPairwiseDist(i,1,1) = mean(x(:)); % avg distance between oscillators within same cluster

        else
           
            disp('Number of node pairs within cluster too many.  Do batch mode to compute Avg Pairwise Distance.')
            
            % Run in Batch Mode - Doing a Subset of oscillators...
            batchBeg = [1:batchSize:numel(ind_in)];
            batchEnd = [batchSize:batchSize:numel(ind_in)];
            if ( numel(batchEnd) < numel(batchBeg) )
                batchEnd = [batchEnd,numel(ind_in)];
            end
            

            % Pairwise Distances within Cluster at Phase Initialization
            disp(['Batch # out of ',num2str(numel(batchBeg))])
            for j = 1:numel(batchBeg)

                fprintf('%s',[num2str(j),' '])
        
                st = batchBeg(j);
                nd = batchEnd(j);
                num = nd-st+1;

                dists_in = dist( theta(:,ind_in(st:nd))', theta(:,ind_in(st:nd)) );
                idx = eye(size(dists_in)); 
                x = abs(dists_in(~idx)); % get rid of diagonal (distance between unit and itself)
                if isempty(x)
                    x=0;
                end
                
                batchMN(j) = mean(x(:)); % build up these stats over batches and maybe I can combine them after.
                batchWT(j) = numel(x);
                
            end
            fprintf('\n')
            
            AvgPairwiseDist(i,1,1) = sum(batchMN.*batchWT)./sum(batchWT); % straight mean (weighted)
            
            clear batchMN batchWT

            % VIP: AT SOME POINT I HAVE TO CHECK THAT THE WEIGHTED CIRC_MEAN OF
            % CIRC_MEANS IS THE SAME AS THE CIRC_MEAN ON WHOLE BIG MATRIX.
                        
            
            
        end % check if number of node pairs in cluster enough to warrant batch processing
        
        
        
        
        
        
        
        
        
        
        % Nodes across cluster boundary (1 node in this cluster and other node in different cluster) - 
        % Average Pairwise Distance of (Denominator of Divisive Margin)
        if( numel(ind_in)*numel(ind_out) < thresholdNodes2DoBatch )

            % Pairwise Distances across cluster boundaries at Phase Initialization in One Shot
            dists_out = dist( theta(:,ind_in)', theta(:,ind_out) );
            y = abs(dists_out); % do not need to get rid of diagonal (distance to self) on this one.
            if isempty(y)
                y=0;
            end
            %
            AvgPairwiseDist(i,1,2) = mean(y(:)); % avg distance between oscillators in different clusters
        
        else
            
            disp('Number of node pairs across cluster boundary too many.  Do batch mode to compute Avg Pairwise Distance.')
            
            % Run in Batch Mode - Doing a Subset of oscillators...
            batchBeg = [1:batchSize:numel(ind_out)];
            batchEnd = [batchSize:batchSize:numel(ind_out)];
            if ( numel(batchEnd) < numel(batchBeg) )
                batchEnd = [batchEnd,numel(ind_out)];
            end
            

            % Pairwise Distances across cluster boundaries at Phase Initialization in Batches
            disp(['Batch # out of ',num2str(numel(batchBeg))])
            for j = 1:numel(batchBeg)

                fprintf('%s',[num2str(j),' '])
        
                st = batchBeg(j);
                nd = batchEnd(j);
                num = nd-st+1;

                dists_out = dist( theta(:,ind_in)', theta(:,ind_out(st:nd)) );
                x = abs(dists_out);
                if isempty(x)
                    x=0;
                end
                
                batchMN(j) = mean(x(:)); % build up these stats over batches and maybe I can combine them after.
                batchWT(j) = numel(x);
                
            end
            fprintf('\n')
            
            AvgPairwiseDist(i,1,2) = sum(batchMN.*batchWT)./sum(batchWT); % straight mean (weighted)
            
            clear batchMN batchWT
                        
            
        end % check if number of node pairs across cluster boundary enough to warrant batch processing

        
        
    end % Loop Over Clusters in this Segmentation

    DistAvgPW(:,1) = sum(AvgPairwiseDist(:,:,1)'.*ClusterSize)./N; % weighting by cluster size for within cluster distance
    DistAvgPW(:,2) = sum(AvgPairwiseDist(:,:,2)'.*ClusterSize)./N; % weighting by cluster size for across cluster distance

        
    % DistAvgPW = DistAvgPW./N; % finish weighting Pairwise distance between clusters by size of cluster.
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

end % Loop over Human Segmentations


toc

disp(['Inside MetaClusterAnalysisB, Memory Usage Check : '])
%whos
check_memory_usage
