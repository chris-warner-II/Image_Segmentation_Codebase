function[D] = compute_SensitivityIndex(Cmns, Cstds, cluster_id, cluster_size, RateDistStd, circ)

% ,Dmin]

% This function takes in a vector of Cluster Means & Cluster Stds and
% compute the d' (Sensitivity Index) metric.  
%
% This is working reasonably.  But I think there is still some things to
% make sense of before I can be convinced it is working well to compare one
% phase distribution evolution vs another and esp circular vs linear vars.


directed_ntwrk_flg = 0; % FLAG that = 1 if directed network, = 0 if undirected

if ~exist('RateDistStd','var')
    RateDistStd = 0;
end

% if ~exist('weightClustersFlag','var')
%     weightClustersFlag = 1; 
% end


if(circ)
    RateDistStd = RateDistStd.*(2.*pi); % Does this make the 2 comparable? I dont think it necessarily would, esp with >1 eigenvector
end


for i = 1:numel(Cmns) % Loop through each segmentation (with varying numbers of segments)
    
    i
    
    % Preallocate memory
    numPairs = size(Cmns{i},1);   % number of segments in this segmentation
    pts = size(Cmns{i},2);        % number of timepoints in simulation (or eigenvectors for spectral seg)
    Csize = cluster_size{i};
    D{i} = zeros( numPairs, pts, numel(RateDistStd) ); % d' tensor of size(pairs x pts X numRD) 

    
    
    % COME BACK TO THIS. UNCOMMENT AT SOME POINT!!
%     % Weight each d' term with the proportion of pairs in those 2 clusters
%     if weightClustersFlag
%         N = sum(Csize);
%         S = (repmat(Csize',1,numSegs).*repmat(Csize,numSegs,1)) - diag(Csize); % this S = number of edges in DIRECTED network without self-weights
%         
%         if ~directed_ntwrk_flg
%             S = S./2; % this S is number of edges in UNDIRECTED network
%         end
% 
%         % Wt is percentage of node pairs in each cluster pair
%         Wt = S./( N.*(N-1) );    % this Wt is percent of weights in DIRECTED network  
% 
%         % WHY?? What does this mean? -> Check is that: sum(sum((1-eye(size(Wt))).*Wt)) should be 1. (it is)
% 
%     else 
%         Wt = 1./(numSegs.^2 - numSegs); % Again, want sum(sum((1-eye(size(Wt))).*Wt)) should be 1. (it is) - WHY??
%         % WHAT IS THIS TOO?
%     end
%     
%     
%     
%     
%     if ~directed_ntwrk_flg
%         Wt = Wt./2; % this Wt is percent of weights in UNDIRECTED network    
%     end

    Wt = 1; % NOT SURE WHAT WT is doing right now.  A short circuit line here.
    
    
    
    
    
    for t = 1:pts % Loop through time points
        
        Cmn = squeeze(Cmns{i}(:,t,:));
        Cstd = squeeze(Cstds{i}(:,t,:));
        
        % compute distance between all pairs of cluster means.
        if(circ)
            Num = abs( circ_dist(repmat(Cmn,1,numSegs),repmat(Cmn',numSegs,1)) ); % if circular variable.
            NumUB = 1;
            NumLB = 1;
        else
            Num = abs(Cmn(:,1) - Cmn(:,2)); % if linear variable.
            NumUB = 1;
            NumLB = 1;
        end



        for j = 1:numel(RateDistStd)
            % Compute average standard deviation of cluster pairs and include rate distortion term 
            % (which keeps d' metric from going to infinity. Rate distortion adds a bit of std onto all points.
            Den = sqrt( 0.5*sum(Cstd.^2,2) ) + RateDistStd(j);
            DenUB = 1;
            DenLB = 1;
            
            d = ( Num ./ Den ).*Wt; % assumes 1 rateDistortion & 1 segmentation.

            D{i}(:,t,j) = d;

        end % Loop over different rate distortions.

        t

    end % Loop over t time points
    
end % Loop over i # of segmentations




% Plot d' vs eigenvector for each cluster pair as well as mean across cluster pairs.
if(0)
    figure, hold on,
    plot(D{1}')
    plot(mean(D{1}),'k--','LineWidth',2)
    title('Sensitivity (d'') for each cluster and \mu')
    xlabel('Eigenvector #')
    ylabel('d''')
end


% Visualization / Error Checking plot of oscillator phases alongsize mean & stds.
if(0)
    
    skip_t = time_pts/100;
    figure, 
    subplot(311), hold on
    errorbar(repmat([1:skip_t:time_pts]',1,numSegs),Cmns{1}(:,1:skip_t:time_pts)',Cstds{1}(:,1:skip_t:time_pts)')
    plot(repmat([1:skip_t:time_pts]',1,numSegs),Cmns{1}(:,1:skip_t:time_pts)','LineWidth',2)
    xlabel('Simulation Time','FontSize',18,'FontWeight','Bold')
    ylabel('Phase','FontSize',18,'FontWeight','Bold')
    title('Cluster \mu''s & \sigma''s','FontSize',20,'FontWeight','Bold')
    set(gca,'ytick',[-pi:pi/2:pi],'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},'FontSize',16,'FontWeight','Bold','DefaultTextInterpreter','latex')
    xlim([0 time_pts])
    ylim([-pi pi])
    legend({'1','2','3','4'},'Location','EastOutside')
    %
    subplot(312)
    plot(D)
    legend({num2str(RateDistStd(1),2),num2str(RateDistStd(2),2),num2str(RateDistStd(3),2),num2str(RateDistStd(4),2),...
        num2str(RateDistStd(5),2),num2str(RateDistStd(6),2),num2str(RateDistStd(7),2)},'Location','EastOutside')
    ylabel('d'' mean pairwise','FontSize',18,'FontWeight','Bold')
    set(gca,'FontSize',16,'FontWeight','Bold')
    
    %
    subplot(313)
    plot(Dmin)
    legend({num2str(RateDistStd(1),2),num2str(RateDistStd(2),2),num2str(RateDistStd(3),2),num2str(RateDistStd(4),2),...
        num2str(RateDistStd(5),2),num2str(RateDistStd(6),2),num2str(RateDistStd(7),2)},'Location','EastOutside')
    ylabel('d'' min pairwise','FontSize',18,'FontWeight','Bold')
    set(gca,'FontSize',16,'FontWeight','Bold')
    
    keyboard
    
end

