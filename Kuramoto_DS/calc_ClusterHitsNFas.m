function [Pcorrect,threshs] = calc_ClusterHitsNFas(gT,inpt,circ)

% Syntax: [Pcorrect,threshs] = calc_ClusterHitsNFas(gT,inpt,circ)
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

bins = 1000; % # of bins to chop interval [0,1] or [0,2pi] into.
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
           
            n1 = hist(c1,xbins) ./ sum(N); % histogram (PDF) of cluster 1.
            n2 = hist(c2,xbins) ./ sum(N); % histogram (PDF) of cluster 2.
            
            ncs1 = cumsum(n1);             % cumulative histogram (CDF) of c1.    
            ncs2 = cumsum(n2);             % cumulative histogram (CDF) of c2.
            

            if(circ)
                
                
                % Step 0: Assume that peak of larger distribution is going
                % to be inside region defined by thresholds {THETA1,THETA2}
                % that contains that cluster. I think it is a safe assumption.
                
                
                % Step 1:  Set temporary threshold THETA0 at peak of larger distribution
                n = [n1;n2];
                [ori,th0] = find(n==max(n(:)));
                %
                if (numel(th0)>1)
                    th0 = th0(round(numel(th0)/2)); % take threshold in middle if there are many values with minimum intersct value
                    ori = ori(round(numel(ori)/2));
                end
                
                
                % Step 2: Compute cumulative right-going sums from THETA0.
                n1_temp = [n1(th0:end),n1(1:th0-1)]; 
                n2_temp = [n2(th0:end),n2(1:th0-1)];
                %
                ncs1 = cumsum(n1_temp);             % cumulative histogram (CDF) of c1.    
                ncs2 = cumsum(n2_temp);             % cumulative histogram (CDF) of c2.
                
                
                % Step 3: Compute False Alarms as you move threshold THETA1 right 
                % (find optimal spot, where total false alarms is minimized)
                if(ori==1)
                    sum_fa = N(1)./sum(N)-ncs1+ncs2; % this assumes cluster1 is on left
                elseif(ori==2)
                    sum_fa = N(2)./sum(N)-ncs2+ncs1; % this assumes cluster2 is on left
                    % NOTE: sum_fa2 = ( 1 - sum_fa1 )
                else
                    disp('Uh oh')
                end
                %
                th1 = find(sum_fa==min(sum_fa(:)));
                %
                if (numel(th1)>1)
                    th1 = th1(round(numel(th1)/2)); % take threshold in middle if there are many values with minimum intersct value
                end
                %
                th1 = mod(th1+th0,bins); % transform THETA1 to live back in original (not temporary) distribution space.
                
                
                % Step 4: Redo this from THETA0 going leftward to find THETA2.
                n1_temp = fliplr(n1_temp); 
                n2_temp = fliplr(n2_temp);
                %
                ncs1 = cumsum(n1_temp);             % cumulative histogram (CDF) of c1.    
                ncs2 = cumsum(n2_temp);             % cumulative histogram (CDF) of c2.
                %
                if(ori==1)
                    sum_fa1 = N(1)./sum(N)-ncs1+ncs2; % this assumes cluster1 is on left
                elseif(ori==2)
                    sum_fa1 = N(2)./sum(N)-ncs2+ncs1; % this assumes cluster2 is on left
                else
                    disp('Uh oh')
                end
                % NOTE: sum_fa2 = ( 1 - sum_fa1 )
                %       sum_fa2 = 1 - sum_fa1;
                %
                th2 = find(sum_fa1==min(sum_fa1(:)));
                %
                if (numel(th2)>1)
                    th2 = th2(round(numel(th2)/2)); % take threshold in middle if there are many values with minimum intersct value
                end
                %
                th2 = mod(th0-th2,bins);
                
                
                % Step 5: Compute Hits & False Alarms for both Clusters using THETA1 & THETA2
                if(th1==th2)
                    
                    hit = max(N) ./ sum(N);
                    fa = min(N) ./ sum(N);
                    %
                    inClust1 = hit;
                    outClust1 = 0;
                    inClust2 = 0;
                    outClust2 = fa;
                    
                    
                else
                
                    th_a = min([th1,th2]);
                    th_b = max([th1,th2]);                
                    %
                    LHS1 = sum(n1(th_a:th_b));
                    RHS1 = sum([n1(1:(th_a-1)),n1((th_b+1):end)]);
                    LHS2 = sum(n2(th_a:th_b));
                    RHS2 = sum([n2(1:(th_a-1)),n2((th_b+1):end)]);
                    %
                    class1 = (LHS1+RHS2); % either hits or fa's.
                    class2 = (LHS2+RHS1); % either hits or fa's.
                    %
                    if( class1 > class2 )
                        inClust1 = LHS1;
                        outClust1 = RHS1;
                        inClust2 = RHS2;
                        outClust2 = LHS2;
                    else
                        inClust1 = RHS1;
                        outClust1 = LHS1;
                        inClust2 = LHS2;
                        outClust2 = RHS2;
                    end
                    %
                    % Hits & False Alarms of Pair of Clusters using the Optimal Threshold.
                    hit = sum(inClust1) + sum(inClust2);
                    fa = sum(outClust1) + sum(outClust2);
                
                end
                
                
                % STEP 6:  HOW TO COMPUTE ROC CURVE (TP,FP,TN,FN) FOR THE
                % CIRCULAR CASE WITH 2 THRESHOLDS WITHOUT DOING GRID SEARCH
                
                
                
                % Step 6-1:  Set temporary threshold THETA0 at peak of larger distribution
                n = [n1;n2];
                [ori,th0] = find(n==max(n(:)));
                %
                if (numel(th0)>1)
                    th0 = th0(round(numel(th0)/2)); % take threshold in middle if there are many values with minimum intersct value
                    ori = ori(round(numel(ori)/2));
                end
                
                
                % Step 6-2: Compute cumulative right-going sums from THETA0.
                n1_temp = [n1(th0:end),n1(1:th0-1)]; 
                n2_temp = [n2(th0:end),n2(1:th0-1)];
                %
                ncs1 = cumsum(n1_temp);             % cumulative histogram (CDF) of c1.    
                ncs2 = cumsum(n2_temp);             % cumulative histogram (CDF) of c2.
                
                % Step 6-3: Compute False Alarms as you move threshold THETA1 right  (find optimal spot, where total false alarms is minimized)
                if(ori==1)
                    sum_fa = N(1)./sum(N)-ncs1+ncs2; % this assumes cluster1 is on left
                elseif(ori==2)
                    sum_fa = N(2)./sum(N)-ncs2+ncs1; % this assumes cluster2 is on left
                    % NOTE: sum_fa2 = ( 1 - sum_fa1 )
                else
                    disp('Uh oh')
                end
                %
                th1 = find(sum_fa==min(sum_fa(:)));
                %
                if (numel(th1)>1)
                    th1 = th1(round(numel(th1)/2)); % take threshold in middle if there are many values with minimum intersct value
                end
                %
                th1 = mod(th1+th0,bins); % transform THETA1 to live back in original (not temporary) distribution space.

                % Step 6-4: Starting from THRESH1, do line search of
                % THRESH2 values and compute associated (TP,FP,TN,FN)
                % values. Like Linear case.
                
                n1_temp2 = [n1(th1:end),n1(1:th1-1)]; 
                n2_temp2 = [n2(th1:end),n2(1:th1-1)];
                %
                ncs1 = cumsum(n1_temp2);             % cumulative histogram (CDF) of c1.    
                ncs2 = cumsum(n2_temp2);             % cumulative histogram (CDF) of c2.
                
                TP = ncs1;              % number correctly classified to be inside cluster 1 
                FP = ncs2;              % number incorrectly classified to be outside cluster 2
                TN = N(2)./sum(N)-ncs2; % number correctly classified to be inside cluster 2
                FN = N(1)./sum(N)-ncs1; % number incorrectly classified to be outside cluster 1
                
                [abs(trapz(FP,TP)),abs(trapz(FN,TN))] % ROC_AUC

                if(1)
                    figure, 
                    subplot(121), hold on,
                    plot(TP,'r--')
                    plot(FP,'r-')
                    plot(TN,'b--')
                    plot(FN,'b-')
                    legend('TP','FP','TN','FN')
                    title('Fixing \Theta1 to right of peak of larger cluster')

                    % Plot ROC Curve
                    subplot(122), hold on
                    plot(FP,TP,'b') % assumes cluster 1 is on the left of threshold
                    plot(FN,TN,'r') % assumes cluster 1 is on the right of threshold
                    xlabel('FP')
                    ylabel('TP')
                    title('ROC Curve')
                    axis([0 1 0 1])
                end
                
                % Step 6-5: Redo this from THETA0 going leftward to find THETA2.
                n1_temp = fliplr(n1_temp); 
                n2_temp = fliplr(n2_temp);
                %
                ncs1 = cumsum(n1_temp);             % cumulative histogram (CDF) of c1.    
                ncs2 = cumsum(n2_temp);             % cumulative histogram (CDF) of c2.
                %
                if(ori==1)
                    sum_fa1 = N(1)./sum(N)-ncs1+ncs2; % this assumes cluster1 is on left
                elseif(ori==2)
                    sum_fa1 = N(2)./sum(N)-ncs2+ncs1; % this assumes cluster2 is on left
                else
                    disp('Uh oh')
                end
                % NOTE: sum_fa2 = ( 1 - sum_fa1 )
                %       sum_fa2 = 1 - sum_fa1;
                %
                th2 = find(sum_fa1==min(sum_fa1(:)));
                %
                if (numel(th2)>1)
                    th2 = th2(round(numel(th2)/2)); % take threshold in middle if there are many values with minimum intersct value
                end
                %
                th2 = mod(th0-th2,bins);
                
                n1_temp2 = [n1(th2:end),n1(1:th2-1)]; 
                n2_temp2 = [n2(th2:end),n2(1:th2-1)];
                %
                ncs1 = cumsum(n1_temp2);             % cumulative histogram (CDF) of c1.    
                ncs2 = cumsum(n2_temp2);             % cumulative histogram (CDF) of c2.
                
                TP = ncs1;              % number correctly classified to be inside cluster 1 
                FP = ncs2;              % number incorrectly classified to be outside cluster 2
                TN = N(2)./sum(N)-ncs2; % number correctly classified to be inside cluster 2
                FN = N(1)./sum(N)-ncs1; % number incorrectly classified to be outside cluster 1
                
                [abs(trapz(FP,TP)),abs(trapz(FN,TN))] % ROC_AUC

                if(1)
                    figure, 
                    subplot(121), hold on,
                    plot(TP,'r--')
                    plot(FP,'r-')
                    plot(TN,'b--')
                    plot(FN,'b-')
                    legend('TP','FP','TN','FN')
                    title('Fixing \Theta2 to left of peak of larger cluster')

                    % Plot ROC Curve
                    subplot(122), hold on
                    plot(FP,TP,'b') % assumes cluster 1 is on the left of threshold
                    plot(FN,TN,'r') % assumes cluster 1 is on the right of threshold
                    xlabel('FP')
                    ylabel('TP')
                    title('ROC Curve')
                    axis([0 1 0 1])
                end
                
                
                % STEP 7:  COMPUTE ROC CURVE (TP,FP,TN,FN) FOR THE
                % CIRCULAR CASE WITH 2 THRESHOLDS BY DOING GRID SEARCH
                
                
                
                
                
                keyboard
                
                
            else % if linear variables (not circular)
                
                sum_fa1 = N(1)./sum(N)-ncs1+ncs2; % this assumes cluster1 is on left
                sum_fa2 = N(2)./sum(N)-ncs2+ncs1; % this assumes cluster2 is on left
                % NOTE: sum_fa2 = ( 1 - sum_fa1 )
                
                
                
                % Calculate ROC curve area and or Precision-Recall Curve Area
                TP = ncs1;              % number correctly classified to be inside cluster 1 
                FP = ncs2;              % number incorrectly classified to be outside cluster 2
                TN = N(2)./sum(N)-ncs2; % number correctly classified to be inside cluster 2
                FN = N(1)./sum(N)-ncs1; % number incorrectly classified to be outside cluster 1

                ROC_AUC = max([abs(trapz(FP,TP)),abs(trapz(FN,TN))]);

                if(plt_flg)
                    figure, 
                    subplot(121), hold on,
                    plot(TP,'r--')
                    plot(FP,'r-')
                    plot(TN,'b--')
                    plot(FN,'b-')
                    legend('TP','FP','TN','FN')

                    % Plot ROC Curve
                    subplot(122), hold on
                    plot(FP,TP,'b') % assumes cluster 1 is on the left of threshold
                    plot(FN,TN,'r') % assumes cluster 1 is on the right of threshold
                    xlabel('FP')
                    ylabel('TP')
                    title('ROC Curve')
                    axis([0 1 0 1])
                end

                % Plot Precision-Recall Curve maybe too

                
                
                

                sum_fa = [sum_fa1;sum_fa2];
                [ori,th1] = find(sum_fa==min(sum_fa(:)));

                if (numel(th1)>1)
                    th1 = th1(round(numel(th1)/2)); % take threshold in middle if there are many values with minimum intersct value
                    ori = ori(round(numel(ori)/2));
                end

                hit = 1-sum_fa(ori,th1); % probability of correct classification

            end
            

            
            
            
            
            
            % To plot the distribution pair and threshold and CDFs etc.
            if(plt_flg)
                if(circ) % CIRCULAR VARIABLE
                    
                    
                    figure, hold on
                    plot(xbins,n1,'b')
                    plot(xbins,n2,'g')
                    plot( [xbins(th1),xbins(th1)],[0,max([n1,n2])],'r--')
                    text(xbins(th1),max([n1,n2]),'\color{red}{\theta_1}')
                    title('Circular PDF')
                    xlim([0,2*pi])
                    plot( [xbins(th2),xbins(th2)],[0,max([n1,n2])],'r--')
                    text(xbins(th2),max([n1,n2]),'\color{red}{\theta_2}')
                    text(0,max([n1,n2]),['\color{blue}{N=',num2str(N(1)),'}'])
                    text(0.9*2*pi,max([n1,n2]),['\color{green}{N=',num2str(N(2)),'}'])
                    %
                    text( xbins(round(th1/2)), 0.7.*max([n1,n2]), ['\color{blue}{P_{hit} = ',num2str( inClust1, 2 ),'}'] )
                    text( xbins(round(th1/2)), 0.5.*max([n1,n2]), ['\color{green}{P_{fa} = ',num2str(  outClust2, 2),'}'] )
                    text( xbins(round((th1+bins)/2)), 0.7.*max([n1,n2]), ['\color{green}{P_{hit} = ',num2str( inClust2, 2),'}'] )
                    text( xbins(round((th1+bins)/2)), 0.5.*max([n1,n2]), ['\color{blue}{P_{fa} = ',num2str( outClust1, 2 ),'}'] )
                    %
                    plot( [xbins(th2),xbins(th2)],[0,max([n1,n2])],'r--')
                    
                    
                    
                else      % LINEAR VARIABLE
                
                
                    figure, 
                    subplot(311), hold on
                    plot(xbins,n1,'b')
                    plot(xbins,n2,'g')
                    plot( [xbins(th1),xbins(th1)],[0,max([n1,n2])],'r--')
                    text(xbins(th1),max([n1,n2]),'\color{red}{\theta_1}')
                    title('PDF')
                    xlim([0,1])
                    text(0,max([n1,n2]),['\color{blue}{N=',num2str(N(1)),'}'])
                    text(0.9,max([n1,n2]),['\color{green}{N=',num2str(N(2)),'}'])
                    %
                    subplot(312), hold on
                    plot(xbins,N(1)./sum(N)-ncs1,'b')
                    plot(xbins,ncs2,'g')
                    plot( [xbins(th1),xbins(th1)],[0,1],'r--')
                    text(xbins(th1),1,'\color{red}{\theta_1}')
                    text( xbins(round(th1/2)), 0.7, ['\color{blue}{P_{hit} = ',num2str( ncs1(th1), 2 ),'}'] )
                    text( xbins(round(th1/2)), 0.5, ['\color{green}{P_{fa} = ',num2str(  ncs2(th1), 2),'}'] )
                    text( xbins(round((th1+bins)/2)), 0.7, ['\color{green}{P_{hit} = ',num2str( N(2)./sum(N)-ncs2(th1), 2),'}'] )
                    text( xbins(round((th1+bins)/2)), 0.5, ['\color{blue}{P_{fa} = ',num2str( N(1)./sum(N)-ncs1(th1), 2 ),'}'] )
                    title('CDF')
                    xlim([0,1])
                    ylim([0,1])
                    %
                    subplot(313), hold on
                    plot(xbins,sum_fa1,'b'),
                    xlabel('\theta')
                    ylabel('P_{fa}')
                    plot( [xbins(th1),xbins(th1)],[0,1],'r--')
                    plot(xbins,sum_fa2,'g'),
                    

                end
            end
            
            
            
            Pcorrect{i}(j,k) = ROC_AUC;
            Pcorrect{i}(k,j) = ROC_AUC;
            
            
            
            threshs{i}(j,k) = xbins(th1);
            if(circ)
                threshs{i}(k,j) = xbins(th2); % for circular variables, can store 2 thresholds in a non-symmetric matrix.
            else
                threshs{i}(k,j) = xbins(th1);
            end
            
            
            
        end % Loop over Cluster 1

    end % Loop over Cluster 2
    
end % Loop over Ground Truths




if(plt_flg)
    keyboard
    close all
end
