function [MC] = metaClusterAnalysisB(theta, netParams, kurParams, circ)


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


disp(['Inside metaClusterAnalysisB function : '])

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

if(size(theta,1)==T)
    theta = theta'; % flip theta because all these functions want it the other way..
else
    T = size(theta,2);
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
    % a first cut at calculating linear variable mean & dists using circular functions , but it seems to be working.
end


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

            % Pairwise Distances at Phase Initialization
            dists_in = circ_dist(repmat(theta(ind_in,1),1,numel(ind_in)),repmat(theta(ind_in,1)',numel(ind_in),1));
            x = abs(dists_in(dists_in~=0));
            if isempty(x)
                x=0;
            end
            %
            AvgPairwiseDist(i,1,1) = circ_mean(x(:)); % avg distance between oscillators within same cluster
            
            
            % Pairwise Distance within Cluster as Time Advances
            for k = 1:(T/spp)

                disp(['Within Cluster: Time Step # ',num2str(k),' / ',num2str(T/spp)])

                dists_in = circ_dist(repmat(theta(ind_in,spp*k),1,numel(ind_in)),repmat(theta(ind_in,spp*k)',numel(ind_in),1));
                x = abs(dists_in(dists_in~=0));
                if isempty(x)
                    x=0;
                end
                %
                AvgPairwiseDist(i,k+1,1) = circ_mean(x(:)); % avg distance between oscillators within same cluster

            end
            
        
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

                dists_in = circ_dist(repmat(theta(ind_in,1),1,num),repmat(theta(ind_in(st:nd),1)',numel(ind_in),1));
                x = abs(dists_in(dists_in~=0));
                if isempty(x)
                    x=0;
                end
                
                batchMN(j) = circ_mean(x(:)); % build up these stats over batches and maybe I can combine them after.
                batchWT(j) = numel(x);
                
            end
            fprintf('\n')
            
            AvgPairwiseDist(i,1,1) = circ_mean(batchMN',batchWT');        % circular mean (weighted)
            
            clear batchMN batchWT
            
            % or
            % sum(batchMN.*batchWT)./sum(batchWT) % straight mean (weighted)
            
            % VIP: AT SOME POINT I HAVE TO CHECK THAT THE WEIGHTED CIRC_MEAN OF
            % CIRC_MEANS IS THE SAME AS THE CIRC_MEAN ON WHOLE BIG MATRIX.
                        
            
            
            
            % Pairwise Distance within Cluster as Time Advances
            for k = 1:(T/spp)

                disp(['Within Cluster: Time Step # ',num2str(k),' / ',num2str(T/spp)])
                disp(['Batch # out of ',num2str(numel(batchBeg))])
                for j = 1:numel(batchBeg)

                    fprintf('%s',[num2str(j),' '])

                    st = batchBeg(j);
                    nd = batchEnd(j);
                    num = nd-st+1;

                    % Pairwise Distances at Phase Initialization
                    dists_in = circ_dist(repmat(theta(ind_in,spp*k),1,num),repmat(theta(ind_in(st:nd),spp*k)',numel(ind_in),1));
                    x = abs(dists_in(dists_in~=0));
                    if isempty(x)
                        x=0;
                    end

                    batchMN(j) = circ_mean(x(:)); % build up these stats over batches and maybe I can combine them after.
                    batchWT(j) = numel(x);

                end
                fprintf('\n')
                
                AvgPairwiseDist(i,k+1,1) = circ_mean(batchMN',batchWT');        % circular mean (weighted)
                
                clear batchMN batchWT
                
            end
            
            
        end % check if number of node pairs in cluster enough to warrant batch processing
        
        
        
        
        
        
        
        
        
        
        % Nodes across cluster boundary (1 node in this cluster and other node in different cluster) - 
        % Average Pairwise Distance of (Denominator of Divisive Margin)
        if( numel(ind_in)*numel(ind_out) < thresholdNodes2DoBatch )

            % Pairwise Distances across cluster boundaries at Phase Initialization in One Shot
            dists_out = circ_dist(repmat(theta(ind_in,1),1,numel(ind_out)),repmat(theta(ind_out,1)',numel(ind_in),1));
            y = abs(dists_out(dists_out~=0));
            if isempty(y)
                y=0;
            end
            %
            AvgPairwiseDist(i,1,2) = circ_mean(y(:)); % avg distance between oscillators in different clusters
            
            
            
            % Pairwise Distance across Cluster Boundary as Time Advances in One Shot
            for k = 1:(T/spp)

                disp(['Across Clusters : Time Step # ',num2str(k),' / ',num2str(T/spp)])

                dists_out = circ_dist(repmat(theta(ind_in,spp*k),1,numel(ind_out)),repmat(theta(ind_out,spp*k)',numel(ind_in),1));
                y = abs(dists_out(dists_out~=0));
                if isempty(y)
                    y=0;
                end
                %
                AvgPairwiseDist(i,k+1,2) = circ_mean(y(:)); % avg distance between oscillators in different clusters
                %
            end
        
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

                dists_out = circ_dist(repmat(theta(ind_in,1),1,num),repmat(theta(ind_out(st:nd),1)',numel(ind_in),1));
                x = abs(dists_out(dists_out~=0));
                if isempty(x)
                    x=0;
                end
                
                batchMN(j) = circ_mean(x(:)); % build up these stats over batches and maybe I can combine them after.
                batchWT(j) = numel(x);
                
            end
            fprintf('\n')
            
            AvgPairwiseDist(i,1,2) = circ_mean(batchMN',batchWT');        % circular mean (weighted)
            
            clear batchMN batchWT
                        
            
            % Pairwise Distance across Cluster boundaries as Time Advances in Batches
            for k = 1:(T/spp)

                disp(['Across Cluster Boundary: Time Step # ',num2str(k),' / ',num2str(T/spp)])
                
                disp(['Batch # out of ',num2str(numel(batchBeg))])
                for j = 1:numel(batchBeg)

                    fprintf('%s',[num2str(j),' '])

                    st = batchBeg(j);
                    nd = batchEnd(j);
                    num = nd-st+1;

                    dists_out = circ_dist(repmat(theta(ind_in,spp*k),1,num),repmat(theta(ind_out(st:nd),spp*k)',numel(ind_in),1));
                    x = abs(dists_out(dists_out~=0));
                    if isempty(x)
                        x=0;
                    end

                    batchMN(j) = circ_mean(x(:)); % build up these stats over batches and maybe I can combine them after.
                    batchWT(j) = numel(x);

                end
                fprintf('\n')
                
                AvgPairwiseDist(i,k+1,2) = circ_mean(batchMN',batchWT');        % circular mean (weighted)
                
                clear batchMN batchWT
                
            end
            
            
        end % check if number of node pairs across cluster boundary enough to warrant batch processing
        
        
        

        DistAvgPW(:,1) = DistAvgPW(:,1) + AvgPairwiseDist(i,:,1)'.*ClusterSize(i); % weighting by cluster size for within cluster distance
        DistAvgPW(:,2) = DistAvgPW(:,2) + AvgPairwiseDist(i,:,2)'.*ClusterSize(i); % weighting by cluster size for across cluster distance

        
        
    end % Loop Over Clusters in this Segmentation

    
    
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

end % Loop over Human Segmentations


toc

disp(['Inside MetaClusterAnalysisB, Memory Usage Check : '])
whos
check_memory_usage
