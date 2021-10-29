function [AUC_ROC_1D, AUC_ROC_GS] = calc_ClusterPairROC_B(gT,inpt,circ)

% What to output actually?
%
% For circular:
% sum_hits, sum_fas, AUC1 - with thresholds tied.
% hits_accum, fas_accum, AUC2 - with thresholds separate.
%
% For linear: 
% 
%

% Syntax: [AUC_ROC] = calc_ClusterPairROC(gT,inpt,circ)
%
% This function takes in a cell of ground truths (gT) and an input (inpt) 
% which can be an input image, the phases of oscillators or values of an
% eigenvector. (circ) is a flag that tells whether this is a circular
% variable (phase) or not (eigen or img).
%
% The function will output a cell of matrices which tell hit rate - 
% probability of correct classification of units in a pair of clusters
% using an optimal threshold (Pcorrect).  Each different cell corresponds
% to a cell in the gT.  And each i,j pair in the matrix is the probability
% of correct classification using that cluster pair.  (threshs) is a
% similarly organized cell array of matrices of the optimal threshold values
% used.  In the linear case, both the outputs will be symmetric matrices.  
% In the circular case, the threshs matrix will contain THETA1 in the i,j
% position and THETA2 in the j,i position.  Sometimes these values can be
% equal when you get the best results by classifying all units into one
% cluster.



plt_flg = 0;

bins = 100; % # of bins to chop interval [0,1] or [0,2pi] into.
% Note: This value says something about the temporal/phase resolution of
% the downstream readout neuron. Whole phase from 0-2pi has a temporal
% duration of ~17ms if we think it is the gamma (60Hz) oscillation.
% Splitting this into 100 bins means that the downstream neuron can
% differentiate between spikes time of arrival on the order of 0.2ms.


if(circ)
    xbins = linspace(0,2*pi,bins); % this is for circular
else
    xbins = linspace(0,1,bins); % this is for non-circular
end



% Now, computing mean & std for truncated distr
for i = 1:numel(gT) % loop thru ground truths

    C = unique(gT{i});
    for j = 1:numel(C)
        
        for k = (j+1):numel(C)

            c1 = inpt(gT{i}==C(j));        % all units in cluster 1.
            c2 = inpt(gT{i}==C(k));        % all units in cluster 2.

            N(1) = numel(c1);              % number of samples in cluster 1.
            N(2) = numel(c2);              % number of samples in cluster 2.
           
            n1 = hist(c1,xbins) ./ N(1); % histogram (PDF) of cluster 1.
            n2 = hist(c2,xbins) ./ N(2); % histogram (PDF) of cluster 2.
            % normalize away the source properties to only look at channel.
            
            ncs1 = cumsum(n1);             % cumulative histogram (CDF) of c1.    
            ncs2 = cumsum(n2);             % cumulative histogram (CDF) of c2.
            

            if(circ)
               

                % STEP 1:  Find peak of distribution of larger cluster
                n = [n1;n2];
                
                [ori,th0a] = find(n==max(n(:)));
                %
                if (numel(th0a)>1)
                    th0a = th0a(round(numel(th0a)/2)); % take threshold in middle if there are many values with minimum intersct value
                    ori = ori(round(numel(ori)/2));
                end
                
                
                % STEP 2: Compute cumulative right- & left-going sums from THETA0. (These will be used when tieing thresholds together)
                n1_rightward = [n1(th0a:end),n1(1:th0a-1)]; 
                n2_rightward = [n2(th0a:end),n2(1:th0a-1)];
                %th_rightward = [th0a:bins,1:th0a-1];
                %
                n1_leftward = [n1(th0a:-1:1),n1(end:-1:th0a+1)]; 
                n2_leftward = [n2(th0a:-1:1),n2(end:-1:th0a+1)];
                %th_leftward = [th0a:-1:1, bins:-1:th0a+1];
                %
                nRcs1 = cumsum(n1_rightward(2:end));            % cumulative histogram (CDF) of c1 going rightward from THETA0.    
                nRcs2 = cumsum(n2_rightward(2:end));            % cumulative histogram (CDF) of c2 going rightward from THETA0.
                %
                nLcs1 = cumsum(n1_leftward(2:end));             % cumulative histogram (CDF) of c1 going leftward from THETA0.    
                nLcs2 = cumsum(n2_leftward(2:end));             % cumulative histogram (CDF) of c2 going leftward from THETA0.
                
               
                
                
                
                % STEP 3: Add up cumulative sums (assuming each step we move both thresholds out from THETA0 by same amount) to
                % get Hits & FA's for each threshold setting.
                
                % {nRcs1, nLcs1} = hits & {nRcs2, nLcs2} = false alarms.
                if(ori==1) 
                    sum_hits = n1(th0a) + nRcs1(1:ceil(bins/2)) + nLcs1(1:ceil(bins/2));
                    sum_fas = n2(th0a) + nRcs2(1:ceil(bins/2)) + nLcs2(1:ceil(bins/2)); % method1 that ties both THETAs together.
                    %
                    hits = n1;
                    fas = n2; % using these for greedy optimization for ROC curve moving thresholds separately.
                    
                % {nRcs2, nLcs2} = hits & {nRcs1, nLcs1} = false alarms.                    
                elseif(ori==2)
                    sum_hits = n2(th0a) + nRcs2(1:ceil(bins/2)) + nLcs2(1:ceil(bins/2));
                    sum_fas = n1(th0a) + nRcs1(1:ceil(bins/2)) + nLcs1(1:ceil(bins/2)); % method1 that ties both THETAs together.
                    %
                    hits = n2;
                    fas = n1; % using these for greedy optimization for ROC curve moving thresholds separately.
                   
                else
                    disp('Uh oh!  Well this is unexpected and embarrassing.')
                end
                
                
                
                %
                sum_hits = [0,sum_hits,1];
                sum_fas = [0,sum_fas,1];   % explicitly add endpoints to avoid Kur_AUC_1D > Kur_AUC_GS.
                
                
                AUC1 = trapz(sum_fas,sum_hits); % area underneath ROC curve tying THETA1 & THETA2 together
                AUC1 = max(AUC1,1-AUC1);
                
                
                % STEP 4: Move thresholds THETA1 & THETA2 greedily-optimally so that each move increases hits as
                % much as possible while increasing fa's as little as possible.
                greedy_2d_search = 0;
                if(greedy_2d_search)

                    fitness = hits./(hits+fas);    
                    fitness(isnan(fitness)) = inf; % if there are no hits or fa's by moving threshold, move preferentially there.

                    fit_wind = 0; % window of extent to look out in calculation of fitness when deciding which threshold to move.
                                  % if I set this to zero, it is like doing greedy optimization looking only one unit out to move threshold.

                    hits_accum = zeros(1,bins);
                    fas_accum = zeros(1,bins);

                    th1 = th0a; % leftward moving threshold.
                    th2 = th0a; % rightward moving threshold.

                    hits_accum(1) = hits(th0a);
                    fas_accum(1) = fas(th0a);

                    for T = 2:bins

                        % FOR NONZERO FITNESS WINDOW...
    %                     stL = th1;
    %                     ndL = 0
    %                     
    %                     stR = th2;
    %                     enR = 0
    %                     
    %                     fitness_leftlooking = 1
    %                     fitness_rightlooking = 1
    %                     %
    %                     hits_leftlooking = 1
    %                     hits_rightlooking = 1

                        % decide with threshold to move by fitness, then hits, then flip a coin.
                        if (fitness(th2+fit_wind) > fitness(th1-fit_wind)) % NOTE: fitness(th) is greedy and only looks one step away.

                            th2 = mod(th2+1,bins); % move THETA2 rightward
                            if(th2==0); th2=1; end
                            hits_by_move = hits(th2);
                            fas_by_move = fas(th2);

                        elseif(fitness(th2+fit_wind) < fitness(th1-fit_wind))

                            th1 = mod(th1-1,bins); % move THETA1 leftward
                            if(th1==0); th1=bins; end
                            hits_by_move = hits(th1);
                            fas_by_move = fas(th1);

                        else % fitnesses must be equal. Look at #hits in each

                            if (hits(th2+fit_wind) > hits(th1-fit_wind))

                                th2 = mod(th2+1,bins); % move THETA2 rightward
                                if(th2==0); th2=1; end
                                hits_by_move = hits(th2);
                                fas_by_move = fas(th2);

                            elseif(hits(th2+fit_wind) < hits(th1-fit_wind))

                                th1 = mod(th1-1,bins); % move THETA1 leftward
                                if(th1==0); th1=bins; end
                                hits_by_move = hits(th1);
                                fas_by_move = fas(th1);

                            else % #hits must be equal. flip a coin.

                                if(rand > 0.5)
                                    th2 = mod(th2+1,bins); % move THETA2 rightward
                                    if(th2==0); th2=1; end
                                    hits_by_move = hits(th2);
                                    fas_by_move = fas(th2);
                                else
                                    th1 = mod(th1-1,bins); % move THETA1 leftward
                                    if(th1==0); th1=bins; end
                                    hits_by_move = hits(th1);
                                    fas_by_move = fas(th1);
                                end

                            end

                        end

                        % calculate hits & FAs with new threshold settings.
                        hits_accum(T) = hits_accum(T-1) + hits_by_move;
                        fas_accum(T) = fas_accum(T-1) + fas_by_move;

                        [T,th1,th2]


                    end


                    AUC2 = trapz(fas_accum,hits_accum); % area underneath ROC curve tying THETA1 & THETA2 together






%                     % This plot shows ROC curves for 2 different Circular Variable methods.
%                     if(plt_flg)
%                         % DOING THIS NOW IN LARGER FIGURE BELOW.
%                         figure, hold on,
%                         plot(sum_fas,sum_hits,'b')
%                         plot(fas_accum,hits_accum,'r')
%                         plot([0 1], [0 1],'k--')
%                         title('ROC Curves : Circular Variables')
%                         xlabel('false alarms')
%                         ylabel('hits')
%                         axis([0 1 0 1])
%                         text(0.75, 0.25, ['\color{blue}{AUC = ',num2str(AUC1,2),'}'])
%                         text(0.75, 0.22, ['\color{red}{AUC = ',num2str(AUC2,2),'}'])
%                         legend({'tying \theta1 & \theta2 together','moving \theta1 & \theta2 separately'},'Location','NorthWest')
%                     end


                
                end
                
                % DO GRID SEARCH SETTING THRESHOLD1 AND VARYING THRESHOLD2 TO GET A BUNCH OF ROC CURVES. THE OPTIMAL 
                % THING TO DO WOULD BE TO TAKE THE CONVEX HULL OF ALL POINTS TO GET OPTIMAL ROC CURVE.
                convhull_gridsearch = 1;
                if(convhull_gridsearch)
                    
                    TPg_R = zeros(bins,bins);
                    TNg_R = zeros(bins,bins);
                    %
                    TPg_L = zeros(bins,bins); % separate out left and right going analyses from each THETA1
                    TNg_L = zeros(bins,bins);

                    for B = 1:bins

                        %B

                        hits_temp_Rit = [hits(B:end),hits(1:B-1)]; 
                        fas_temp_Rit = [fas(B:end),fas(1:B-1)];
                        %
                        hits_temp_Lft = [hits(B-1:-1:1),hits(end:-1:B)]; 
                        fas_temp_Lft = [fas(B-1:-1:1),fas(end:-1:B)];
                        
                        % %
                        
                        hits_cs_R = cumsum(hits_temp_Rit);    % rightward going cumulative histogram (CDF) of c1.    
                        fas_cs_R = cumsum(fas_temp_Rit);      % rightward going cumulative histogram (CDF) of c2.
                        %
                        hits_cs_L = cumsum(hits_temp_Lft);    % leftward going cumulative histogram (CDF) of c1.    
                        fas_cs_L = cumsum(fas_temp_Lft);      % leftward going cumulative histogram (CDF) of c2.
                        
                        % %
                        
                        TPg_R(B,:) = hits_cs_R;             % number correctly classified to be inside (rightgoing)
                        FPg_R(B,:) = fas_cs_R;              % number incorrectly classified to be inside (rightgoing)
                        %
                        TPg_L(B,:) = hits_cs_L;             % number correctly classified to be inside (leftgoing)
                        FPg_L(B,:) = fas_cs_L;              % number incorrectly classified to be inside (leftgoing)
                        
                        % %
                        
                        x1_R(B) = trapz(FPg_R(B,:),TPg_R(B,:)); % ROC AUC with one threshold set at B (rightgoing)
                        x1_L(B) = trapz(FPg_L(B,:),TPg_L(B,:)); % ROC AUC with one threshold set at B (leftgoing)

                    end
                    
                    
                    
                    % Recompute the one threshold value that maximizes Area under ROC curve using this grid search method.
                    th1L = find(x1_L==max(x1_L)); th1L = th1L(1);
                    th1R = find(x1_R==max(x1_R)); th1R = th1R(1);
                    %
                    [xbins(th1L), x1_L(th1L)]
                    [xbins(th1R), x1_L(th1R)]
                    
                    
                    % Compute the convex hull or enevlope of all the points from the grid search to get optimal ROC curve.
                    try
                        K1 = flipud(convhull(FPg_L(:),TPg_L(:)));
                        K2 = flipud(convhull(FPg_R(:),TPg_R(:))); % flipud will traverse convex hull of grid search from origin in counterclockwise fashion. Up and then to right.

                        AUC3 = max(abs(cumtrapz([0;FPg_L(K1)],[0;TPg_L(K1)])));
                        AUC4 = max(abs(cumtrapz([0;FPg_R(K2)],[0;TPg_R(K2)]))); % note: again, explicitly add (fp,tp)=(0,0) point to avoid 1D>GS problem
                        AUCgs = max([AUC3,AUC4]);
                    catch
                        AUCgs = nan; % something failed with convex hull calculation.
                    end
                    
                    
                    
                    
                    % Save Area Under ROC Curve into data structure to be passed out to calling function.
                    if(numel(C)==2)

                        AUC_ROC_1D{i} = AUC1;
                        AUC_ROC_GS{i} = AUCgs;

                    else

                        AUC_ROC_1D{i}(j,k) = AUC1;
                        AUC_ROC_1D{i}(k,j) = AUC1;
                        %
                        AUC_ROC_GS{i}(j,k) = AUCgs;
                        AUC_ROC_GS{i}(k,j) = AUCgs;

                    end
                    
                    if(1)
                    % Why does Grid Search with Convex Hull Upperbound sometimes give smaller AUC than 1D search?
                    if(AUC1 - AUCgs > 1e-6)
                        disp('Weird that Upper Bound is smaller than 1D AUC for Kur (beyond machine precision)')

                        N

                        [AUC1, AUCgs]

                       % figure, hold on, plot(n1,'b'),plot(n2,'g')

                        plt_flg=1;
                    end
                    end
                    
                    

                    % Plot Cluster Distributions and Area Under ROC Curve for Circular Variables
                    if(plt_flg)
                        figure
                        %
                        subplot(261); hold on, axis([0 1 0 1])
                        for B = 1:bins
                            plot(FPg_R(B,:),TPg_R(B,:))
                        end
                        %
                        plot(FPg_R(th1R,:),TPg_R(th1R,:),'r--','LineWidth',2)
                        plot( FPg_R(K2), TPg_R(K2) ,'m')
                        plot([0 1],[0 1],'w--')
                        xlabel('FP')
                        ylabel('TP')
                        title('Moving Right From Tallest Peak')
                        axis square
                        %
                        subplot(2,6,2:5), hold on, plot(xbins,n1,'b'), plot(xbins,n2,'g')
                        title('Cluster Probability Distributions')
                        xlabel('Phase Angle (\Theta)')
                        ylabel('p(C_i=\Theta)')
                        xlim([0 2*pi])
                        legend({'C1','C2'},'Location','Best')
                        %
                        subplot(266); hold on, axis([0 1 0 1])
                        for B = 1:bins
                            plot(FPg_L(B,:),TPg_L(B,:))
                        end
                        %
                        plot(FPg_L(th1L,:),TPg_L(th1L,:),'r--','LineWidth',2)
                        plot( FPg_L(K1), TPg_L(K1) ,'m')
                        plot([0 1],[0 1],'w--')
                        xlabel('FP')
                        ylabel('TP')
                        title('Moving Left From Tallest Peak')
                        axis square
                        %
                        subplot(2,5,6:7), hold on, 
                        plot(xbins,abs(x1_R),'b')
                        plot(xbins,repmat(AUC1,size(xbins)),'c--')
                        %plot(xbins,repmat(AUC2,size(xbins)),'r--')
                        plot(xbins,repmat(AUCgs,size(xbins)),'k--')
                        title('GridSearch 1: AUC Setting \Theta_1 and moving \Theta_2 Rightward')
                        xlabel('\Theta_1 Fixed (note: \Theta_2 is Varied rightward from \Theta_1)')
                        ylabel('Area Under ROC')
                        legend({'TP & FP (grid)','AUC tying \Theta1 & \Theta2','convex hull of grid search'},'Location','Best') % 'AUC independent \Theta1 & \Theta2',
                        %
                        scatter(xbins(th1R),abs(x1_R(th1R)),'bx','LineWidth',2)
                        axis([0 2*pi 0 1])
                        %
                        subplot(2,5,8), hold on
                        plot(sum_fas,sum_hits,'c')
                        %plot(fas_accum,hits_accum,'r')
                        plot( FPg_L(K1), TPg_L(K1) ,'k')
                        plot( FPg_R(K2), TPg_R(K2) ,'k')
                        plot([0 1], [0 1],'k--')
                        title('Heuristic ROC Curves')
                        xlabel('FP')
                        ylabel('TP')
                        axis([0 1 0 1])
                        text(0.75, 0.25, ['\color{cyan}{AUC = ',num2str(AUC1,2),'}'])
                        %text(0.75, 0.20, ['\color{red}{AUC = ',num2str(AUC2,2),'}'])
                        text(0.75, 0.15, ['\color{black}{AUC = ',num2str(AUCgs,2),'}'])
                        legend({'tying \theta1 & \theta2 together','convex hull of grid search'},'Location','NorthWest') % 'moving \theta1 & \theta2 separately',
                        axis square
                        %
                        subplot(2,5,9:10), hold on, 
                        plot(xbins,abs(x1_L),'b')
                        plot(xbins,repmat(AUC1,size(xbins)),'c--')
                        %plot(xbins,repmat(AUC2,size(xbins)),'r--')
                        plot(xbins,repmat(AUCgs,size(xbins)),'k--')
                        title('GridSearch 2: AUC Setting \Theta_1 and moving \Theta_2 Leftward')
                        xlabel('\Theta_1 Fixed (note: \Theta_2 is Varied leftward from \Theta_1)')
                        ylabel('Area Under ROC')
                        legend({'TP & FP (grid)','AUC tying \Theta1 & \Theta2','convex hull of grid search'},'Location','Best') % 'AUC independent \Theta1 & \Theta2',
                        %
                        scatter(xbins(th1L),abs(x1_L(th1L)),'bx','LineWidth',2)
                        axis([0 2*pi 0 1])
                        
                        
                        
                        
                        figure, hold on, axis([0 1 0 1])
                        scatter(FPg_L(:),TPg_L(:),'r.'),
                        scatter(FPg_R(:),TPg_R(:),'b.'),
                        plot([0 1],[0 1],'k')
                        plot( FPg_R(K2), TPg_R(K2) ,'m')
                        plot( FPg_L(K1), TPg_L(K1) ,'c')
                        plot(FPg_L(th1L,:),TPg_L(th1L,:),'r--','LineWidth',2)
                        
                        
                        keyboard
                        
                    end

                    
                end
                
                
                
                
                
                
                
            else % if linear variables (not circular)
                
                

                
                % Find larger cluster and call it "hits". Call other "false alarms".
                n = [n1;n2];
                
                [ori,th0a] = find(n==max(n(:)));
                %
                if (numel(th0a)>1)
                    th0a = th0a(round(numel(th0a)/2)); % take threshold in middle if there are many values with minimum intersct value
                    ori = ori(round(numel(ori)/2));
                end

                
                % {nRcs1, nLcs1} = hits & {nRcs2, nLcs2} = false alarms.
                if(ori==1) 
                    hits = n1;
                    fas = n2; % using these for greedy optimization for ROC curve moving thresholds separately.
                    
                % {nRcs2, nLcs2} = hits & {nRcs1, nLcs1} = false alarms.                    
                elseif(ori==2)
                    hits = n2;
                    fas = n1; % using these for greedy optimization for ROC curve moving thresholds separately.
                    
                else
                    disp('Uh oh!  Well this is unexpected and embarrassing.')
                end
                
                

                % NOW USE THIS TO GET HITS VS. FA's FOR EACH THRESHOLD SETTING.
                for B = 1:bins

                    %B

                    TPg(B) = sum(hits(1:B-1));             % number correctly classified to be inside (leftgoing)
                    FPg(B) = sum(fas(1:B-1));              % number incorrectly classified to be inside (leftgoing)

                end
                
                AUC1 = trapz(FPg,TPg); % ROC AUC with one threshold set at B (leftgoing)
                AUC1 = max(AUC1,1-AUC1);
                AUC_ROC_GS = 0;        % where we do a gridsearch across pairs of thresholds for circular variables
                
                
                
                % Save Area Under ROC Curve into data structure to be passed out to calling function.
                if(numel(C)==2)
            
                    AUC_ROC_1D{i} = AUC1;
                    
                else

                    AUC_ROC_1D{i}(j,k) = AUC1;
                    AUC_ROC_1D{i}(k,j) = AUC1;

                end
                
                
                
                
                
                
                % Plot ROC Curve and cluster distributions for linear case.
                if(plt_flg)
                    figure, 
                    subplot(211), hold on
                    plot(hits,'b'), hold on, plot(fas,'g')
                    title('Pair of Cluster Distributions in LINEAR')
                    %
                    subplot(212), hold on
                    plot(FPg,TPg,'r','LineWidth',2)
                    plot([0 1],[0 1],'k--')
                    text(0.75, 0.25, ['\color{black}{AUC = ',num2str(AUC1,2),'}'])
                    xlabel('FP')
                    ylabel('TP')
                    title('ROC Curve')
                    axis([0 1 0 1])
                end

                
            end % if circular or linear variable
            
            
            
            
            
        end % Loop over Cluster 1

    end % Loop over Cluster 2
    
end % Loop over Ground Truths


if(plt_flg)
    keyboard
    close all
end
