function [mn, sig,  d_prime_1D_multEv, d_prime_1D_singEv, cluster_id, cluster_size, notTrunc] = ...
          calc_ClusterMnNVars(gT, inpt, RateDistStd, circ, do_plots, method_str)
% mnCI,sigCI, were also input.
%
%
% Syntax: [mn,sig,notTrunc] = calc_ClusterMnNVars(gT,inpt,circ)
%
% This function takes the input (image, phase resulting from kuramoto,
% eigenvector, whatever) and ground truth provided by human segmentations
% and it will return vectors of mean cluster location (mn), standard
% deviation of cluster (sig) and whether the std results from fitting
% parameters on a truncated distribution on a bounded range vs. the std
% measured empirically because the MLE distribution fitting failed
% (notTrunc).
%
% inpt_T = N x t in size.

notTrunc = [];

a_95 = 0.05;        % Note: Going with 95% Confidence Interval.

pts = size(inpt,2); % for kuramoto, number of timepoints in input.  
                    % For spectral, number of eigenvectors. 

% Now, computing mean & std for truncated distr
for k = 1:numel(gT) % loop thru ground truths

    N = numel(gT{k});
    C = unique(gT{k});
    num_pairs = numel(C).*(numel(C)-1)./2; % N(N-1)/2

    % Gather info on cluster pairs: id, # of nodes, weightings
    cluster_id{k} = zeros(num_pairs,2);       % id's of each cluster pair.
    cluster_size{k} = zeros(num_pairs,2);     % sizes each cluster pair.
    cluster_weights{k} = zeros(num_pairs,2);  % weights of each cluster pair.
    %
    c_pair_cnt = 0;                       % initialize counter of cluster pairs
    for i = 1:numel(C)
        for j = i+1:numel(C)
            c_pair_cnt = c_pair_cnt+1;    % increment counter to index into mn, sig, etc...
            %
            clusterA = find(gT{k}==C(i)); % indices to all units in clusterA.
            clusterB = find(gT{k}==C(j)); % indices to all units in clusterB.
            %
            cluster_id{k}(c_pair_cnt,1) = i;
            cluster_id{k}(c_pair_cnt,2) = j;
            %
            cluster_size{k}(c_pair_cnt,1) = numel(clusterA); % number of nodes in clusterA
            cluster_size{k}(c_pair_cnt,2) = numel(clusterB); % number of nodes in clusterB
        end
    end
    %
    cluster_weights{k} = prod(cluster_size{k},2) ./ sum( prod(cluster_size{k},2) );
    
    if(circ)
        RateDistStd = RateDistStd.*pi;
        disp('How to convert RateDistStd to circ? Not sure this is best or principled even.')
    end

    % preallocate memory
    mn{k}   = zeros(num_pairs,pts,2);
    cstd{k} = zeros(num_pairs,pts,2);
    sig{k}  = zeros(num_pairs,pts,2);
    %mnCI{k} = zeros(numel(C),pts,4); 
    %sigCI{k} = zeros(numel(C),pts,4);
    %
    d_prime_1D_multEv = zeros(num_pairs,pts);
    d_prime_kD_multEv = zeros(num_pairs,pts); 
    d_prime_1D_singEv = zeros(num_pairs,pts);
    %
    ang_1D_projK         = zeros(num_pairs,pts);
    ang_cluster_vecs_kD  = zeros(num_pairs,pts);
    ang_cluster_vecs_1D  = zeros(num_pairs,pts);
    %
    normA1   = zeros(num_pairs,pts); 
    normB1   = zeros(num_pairs,pts); 
    normV1   = zeros(num_pairs,pts);
    %
    normAk   = zeros(num_pairs,pts); 
    normBk   = zeros(num_pairs,pts); 
    normVk   = zeros(num_pairs,pts);
    %
    normAk_1D = zeros(num_pairs,pts); 
    normBk_1D = zeros(num_pairs,pts); 
    normVk_1D = zeros(num_pairs,pts); 
    %
    AB_dist_kD      = zeros(num_pairs,pts); 
    AB_dist_1D      = zeros(num_pairs,pts); 
    ABk_dist_proj1D = zeros(num_pairs,pts);
    %
    AB_std_kD      = zeros(num_pairs,pts); 
    AB_std_1D      = zeros(num_pairs,pts); 
    ABk_std_proj1D = zeros(num_pairs,pts);     % Maybe all these can be 6 x 100 2D matrices
    %
    % normV  = zeros(1,pts); % not sure this is right yet...
    normVb = zeros(N,pts); % THINKING ABOUT NORMS FOR KURAMOTO PHASES...
    
    

    for t = 1:pts
        
        c_pair_cnt = 0;                        % initialize counter of cluster pairs

        inptT = inpt(:,1:t);
        t

        for i = 1:numel(C)
            for j = (i+1):numel(C)
                
                c_pair_cnt = c_pair_cnt+1;    % increment counter to index into mn, sig, etc...

                clusterA = find(gT{k}==C(i)); % indices to all units in clusterA.
                clusterB = find(gT{k}==C(j)); % indices to all units in clusterB.
                
                nA = numel(clusterA);         % number of nodes in clusterA
                nB = numel(clusterB);         % number of nodes in clusterB
                
                %if( size(inptT,2)>1 ) 
                if(circ) % circular variables (kuramoto phase)
                    
                    % Project from kD down to 1D
                    A = inptT(clusterA,:);
                    B = inptT(clusterB,:);

                    Amn = circ_mean(A);  % cluster centers in each of k eigenvector dimensions individually
                    Bmn = circ_mean(B);
                    Astd = circ_std(A);  % standard deviations of clusters in each of k eigenvector dimensions individually
                    Bstd = circ_std(B); 

                    AB_dist_kD(c_pair_cnt,t) = sum( abs( circ_dist(Bmn,Amn) ) );      % distance between cluster centers in k-dims
                    AB_std_kD(c_pair_cnt,t)  = circ_mean( sqrt( 0.5.*(Astd.^2 + Bstd.^2 ) )' + RateDistStd ); % avg std of cluster pair averaged across eigenvectors 1:k

                    AB_dist_1D(c_pair_cnt,t) = sum( abs( circ_dist( Bmn(t), Amn(t) ) ) );      % distance between cluster centers in k-dims
                    AB_std_1D(c_pair_cnt,t)  = circ_mean( sqrt( 0.5.*(Astd(t).^2 + Bstd(t).^2 ) )' + RateDistStd ); % avg std of cluster pair averaged across eigenvectors 1:k


                    % Check for 1D case that sees if projection changes d'. Seems to change 
                    % mn & std values, but hopefully scales them by same amount
                    
                    % disp('Note: Want to do LDA instead of just projecting onto vector connecting cluster means.')
                    % keyboard
                    
                    direcK = circ_dist( Bmn, Amn) ./ ( Astd.^2 + Bstd.^2 + RateDistStd );           % normalized 1D vector between cluster means in multiple K-dims. 
                    direcK = direcK ./ norm(direcK);
                    %
                    if(all(isnan(direcK)))
                        direcK = ones(1,t);
                        direcK = direcK ./ norm(direcK);
                    end
                    %
                    Ak_1D = dot(A,repmat(direcK,nA,1),2);    % project all points in cluster A in k-dims down onto 1D direcK vector.
                    Bk_1D = dot(B,repmat(direcK,nB,1),2);    % project all points in cluster B in k-dims down onto 1D direcK vector.
                    Vk_1D = dot(inptT,repmat(direcK,N,1),2); % project all points in V in k-dims down onto 1D direcK vector.
                    %
                    normAk(c_pair_cnt,t) = norm(A,'fro');      % frobenius matrix norm tells distance in high dim space.
                    normBk(c_pair_cnt,t) = norm(B,'fro');
                    normVk(c_pair_cnt,t) = norm(inptT,'fro');
                    %
                    normA1(c_pair_cnt,t) = norm(A(:,t),'fro'); % not sure it makes sense to compute norm of clusters (not = 1)
                    normB1(c_pair_cnt,t) = norm(B(:,t),'fro'); 
                    normV1(c_pair_cnt,t) = norm(inptT(:,t),'fro');
                    %
                    for h = 1:N
                        normVb(h,t) = norm(inptT(h,:),'fro'); % norm/length of vector of all single nodes node.  With all eigenvectors (h=N), it converges to 1.
                    end
                    %
                    normAk_1D(c_pair_cnt,t) = norm(Ak_1D,'fro'); % not sure it makes sense to compute norm of clusters (not = 1)
                    normBk_1D(c_pair_cnt,t) = norm(Bk_1D,'fro'); 
                    normVk_1D(c_pair_cnt,t) = norm(Vk_1D,'fro'); 
                     
                    
                    

%                         Akmn_1D = circ_mean(Ak_1D); % cluster centers in 1D projection from multiple k-dims.
%                         Bkmn_1D = circ_mean(Bk_1D);
%                         Akstd_1D = circ_std(Ak_1D); % standard deviations of distributions in 1D projection from multiple k-dims.
%                         Bkstd_1D = circ_std(Bk_1D);



                    %NOTE : I should first, i think, still try to compute these with a truncated gaussian, no? (below)
                    mn{k}(c_pair_cnt,t,1) = circ_mean(Ak_1D);             % cluster centers in 1D projection from multiple k-dims.
                    mn{k}(c_pair_cnt,t,2) = circ_mean(Bk_1D);             % empirical mean for clusters
                    sig{k}(c_pair_cnt,t,1) = circ_std(Ak_1D);             % standard deviations of distributions in 1D projection from multiple k-dims.
                    sig{k}(c_pair_cnt,t,2) = circ_std(Bk_1D);             % empirical std (if truncated dist fitting fails)

                    %
                    ABk_dist_proj1D(c_pair_cnt,t) = abs( circ_dist( mn{k}(c_pair_cnt,t,1), mn{k}(c_pair_cnt,t,2) ) );
                    ABk_std_proj1D(c_pair_cnt,t)  =  sqrt( 0.5.*(sig{k}(c_pair_cnt,t,1).^2 + sig{k}(c_pair_cnt,t,2).^2 ) ) + RateDistStd;

                    one = [1,zeros(1,t-1)];
                    aMn_tmpk = circ_mean(A(:,1:t));
                    bMn_tmpk = circ_mean(B(:,1:t));
                    aMn_tmp1 = circ_mean(A(:,t));
                    bMn_tmp1 = circ_mean(B(:,t));
                    ang_1D_projK(c_pair_cnt,t) = acosd( dot(one,direcK) ./ ( norm(one) .* norm(direcK) ) );                      % angle between direc vector (btwn cluster means) and the vector [1,0,0,0,0,...0]
                    ang_cluster_vecs_kD(c_pair_cnt,t) = acosd( dot(aMn_tmpk,bMn_tmpk) ./ ( norm(aMn_tmpk) .* norm(bMn_tmpk) ) ); % angle between cluster means in k-dim space.
                    ang_cluster_vecs_1D(c_pair_cnt,t) = acosd( dot(aMn_tmp1,bMn_tmp1) ./ ( norm(aMn_tmp1) .* norm(bMn_tmp1) ) ); % angle between cluster means in 1-dim space.

                    
                    
                    
                    % Plots for error-checking
                    if(0)
                    if( mod(t,100)==0 && c_pair_cnt == size( cluster_id{k},1) )
                        keyboard
                        
                        figure,
                        subplot(211), plot(mean(AB_dist_1D(:,1:t)))
                        subplot(212), plot(std(AB_dist_1D(:,1:t)))
                        
                        figure,
                        subplot(211), plot(mean(ABk_dist_proj1D(:,1:t)))
                        subplot(212), plot(std(ABk_dist_proj1D(:,1:t)))
                        
                        figure, imagesc([A;B]), colormap(hsv), colorbar
                        
                    end
                    end
                    

                else % linear variables (eigenvectors)

                    % Project from kD down to 1D
                    A = inptT(clusterA,:);
                    B = inptT(clusterB,:);

                    Amn = mean(A);  % cluster centers in each of k eigenvector dimensions individually
                    Bmn = mean(B);
                    Astd = std(A);  % standard deviations of clusters in each of k eigenvector dimensions individually
                    Bstd = std(B); 

                    AB_dist_kD(c_pair_cnt,t) = sum( abs(Bmn - Amn) );      % distance between cluster centers in k-dims
                    AB_std_kD(c_pair_cnt,t)  = mean( sqrt( 0.5.*(Astd.^2 + Bstd.^2 ) ) + RateDistStd ); % avg std of cluster pair averaged across eigenvectors 1:k

                    AB_dist_1D(c_pair_cnt,t) = sum( abs(Bmn(t) - Amn(t)) );      % distance between cluster centers in k-dims
                    AB_std_1D(c_pair_cnt,t)  = mean( sqrt( 0.5.*(Astd(t).^2 + Bstd(t).^2 ) ) + RateDistStd ); % avg std of cluster pair averaged across eigenvectors 1:k


                    % Check for 1D case that sees if projection changes d'. Seems to change 
                    % mn & std values, but hopefully scales them by same amount
                    
%                     disp('Note: Want to do LDA instead of just projecting onto vector connecting cluster means.')
%                     keyboard
                    
                    direcK = (Bmn - Amn) ./ (Bstd.^2 + Astd.^2 + RateDistStd);              % normalized 1D vector between cluster means in multiple K-dims. 
                    direcK = direcK ./ norm(direcK);
                    %
                    if(all(isnan(direcK)))
                        direcK = ones(1,t);
                        direcK = direcK ./ norm(direcK);
                    end
                    %
                    
                    
                    % NOTE: TRY LDA / QDA. AM I PROJECTING DOWN TO 1 DIM FROM K DIM PROPERLY? STRANGE THAT D' 
                    %  FOR K DIMS < D' FOR 1 DIM I THINK.HAVE TO RECOMPUTE DIRECK_LDA AND OR DIRECK_QDA, ETC.
                    
                    
                    %
                    Ak_1D = dot(A,repmat(direcK,nA,1),2);                  % project all points in cluster A in k-dims down onto 1D direcK vector.
                    Bk_1D = dot(B,repmat(direcK,nB,1),2);                  % project all points in cluster B in k-dims down onto 1D direcK vector.
                    Vk_1D = dot(inptT,repmat(direcK,N,1),2);               % project all points in V in k-dims down onto 1D direcK vector.
                    %
                    normA1(c_pair_cnt,t) = norm(A(:,t),'fro');             % not sure it makes sense to compute norm of clusters (not = 1)
                    normB1(c_pair_cnt,t) = norm(B(:,t),'fro'); 
                    normV1(c_pair_cnt,t) = norm(inptT(:,t),'fro');
                    %
                    normAk(c_pair_cnt,t) = norm(A,'fro');                  % frobenius matrix norm tells distance in high dim space.
                    normBk(c_pair_cnt,t) = norm(B,'fro');
                    normVk(c_pair_cnt,t) = norm(inptT,'fro');
                    %
                    for h = 1:N
                        normVb(h,t) = norm(inptT(h,:),'fro');              % norm/length of vector of all single nodes node.  With all eigenvectors (h=N), it converges to 1.
                    end
                    %
                    normAk_1D(c_pair_cnt,t) = norm(Ak_1D,'fro');           % not sure it makes sense to compute norm of clusters (not = 1)
                    normBk_1D(c_pair_cnt,t) = norm(Bk_1D,'fro'); 
                    normVk_1D(c_pair_cnt,t) = norm(Vk_1D,'fro'); 
                    %
                    normA1_1D(c_pair_cnt,t) = norm(A(:,t),'fro');          % not sure it makes sense to compute norm of clusters (not = 1)
                    normB1_1D(c_pair_cnt,t) = norm(B(:,t),'fro'); 
                    normV1_1D(c_pair_cnt,t) = norm(inptT(:,t),'fro');

%                         Akmn_1D = mean(Ak_1D); % cluster centers in 1D projection from multiple k-dims.
%                         Bkmn_1D = mean(Bk_1D);
%                         Akstd_1D = std(Ak_1D); % standard deviations of distributions in 1D projection from multiple k-dims.
%                         Bkstd_1D = std(Bk_1D);



                    if(any (isnan([mean(Ak_1D), mean(Bk_1D), std(Ak_1D), std(Bk_1D)]) ) & t>1 )
                        keyboard
                    end

                    %NOTE : I should first, i think, still try to compute these with a truncated gaussian, no? (below)
                    mn{k}(c_pair_cnt,t,1) = mean(Ak_1D);             % cluster centers in 1D projection from multiple k-dims.
                    mn{k}(c_pair_cnt,t,2) = mean(Bk_1D);             % empirical mean for clusters
                    sig{k}(c_pair_cnt,t,1) = std(Ak_1D);             % standard deviations of distributions in 1D projection from multiple k-dims.
                    sig{k}(c_pair_cnt,t,2) = std(Bk_1D);             % empirical std (if truncated dist fitting fails)

                    %
                    ABk_dist_proj1D(c_pair_cnt,t) = abs( mn{k}(c_pair_cnt,t,1) - mn{k}(c_pair_cnt,t,2) );
                    ABk_std_proj1D(c_pair_cnt,t)  = sqrt( 0.5.*(sig{k}(c_pair_cnt,t,1).^2 + sig{k}(c_pair_cnt,t,2).^2 ) ) + RateDistStd;

                    one = [1,zeros(1,t-1)];
                    aMn_tmpk = mean(A(:,1:t));
                    bMn_tmpk = mean(B(:,1:t));
                    aMn_tmp1 = mean(A(:,t));
                    bMn_tmp1 = mean(B(:,t));
                    ang_1D_projK(c_pair_cnt,t) = acosd( dot(one,direcK) ./ ( norm(one) .* norm(direcK) ) );                      % angle between direc vector (btwn cluster means) and the vector [1,0,0,0,0,...0]
                    ang_cluster_vecs_kD(c_pair_cnt,t) = acosd( dot(aMn_tmpk,bMn_tmpk) ./ ( norm(aMn_tmpk) .* norm(bMn_tmpk) ) ); % angle between cluster means in k-dim space.
                    ang_cluster_vecs_1D(c_pair_cnt,t) = acosd( dot(aMn_tmp1,bMn_tmp1) ./ ( norm(aMn_tmp1) .* norm(bMn_tmp1) ) ); % angle between cluster means in 1-dim space.

                    
                    
                end % if(circ)
                
                
                %end % End analysis for clustering in multiple dimensions.
                
                

%                 % MAYBE REINSTATE THESE DIFFERENT COMPUTATIONS FOR KAPPA AND JUNK...
%                 if(circ) % circular variables (kuramoto phase)
% 
% %                     [mn{k}(c_pair_cnt,t,1),ubA,lbA] =  circ_mean(inptT(clusterA));
% %                     [mn{k}(c_pair_cnt,t,2),ubB,lbB] =  circ_mean(inptT(clusterB));
% % 
% %                     [sig1{k}(c_pair_cnt,t,1),cstd{k}(c_pair_cnt,t,1)] = circ_std(inptT(clusterA));
% %                     [sig1{k}(c_pair_cnt,t,2),cstd{k}(c_pair_cnt,t,2)] = circ_std(inptT(clusterB));
% %                     % sig1 is "angular deviation":           s = sqrt(2(1-R))
% %                     % cstd is "circular standard deviation": s = sqrt(-2lnR)
% % 
% %                     [kapA, ubA, lbA] = circ_kappa(inptT(clusterA)); % ML estimate of concentration parameter
% %                     [kapB, ubB, lbB] = circ_kappa(inptT(clusterB)); % ML estimate of concentration parameter
% % 
% %                     sig{k}(c_pair_cnt,t,1) = sqrt(1./kapA); % std is "analogous" to inverse of concentration parameter
% %                     sig{k}(c_pair_cnt,t,2) = sqrt(1./kapB); % 
% % 
% %                     notTrunc = 0;
% % 
% %                     % plot distribution for this cluster to compare inverse kappa vs angular deviation vs circular std.
% %                     if(0)
% %                         figure, 
% % 
% %                         ctrs = linspace(0,2*pi,100);
% % 
% %                         subplot(311), hist(inptT(cluster),ctrs),xlim([0 2*pi])
% %                         title(['Von Mises Distribution w/ \mu = ',num2str(mn{k}(j,t),2),' [ ',num2str(mnCI{k}(j,t,1),2),' , ',num2str(mnCI{k}(j,t,2),2),'] & 1/\kappa = \sigma = ',num2str(sig{k}(j,t),2),' [ ',num2str(sigCI{k}(j,t,1),2),' , ',num2str(sigCI{k}(j,t,2),2),' ] '])
% %                         %
% %                         X = normrnd(mn{k}(j,t),sig{k}(j,t),size(cluster));
% %                         nx = hist(X,ctrs);
% %                         Y = normrnd(mn{k}(j,t),sigCI{k}(j,t,2),size(cluster));
% %                         ny = hist(Y,ctrs);
% %                         subplot(312), bar(ctrs,ny,'r'), hold on,
% %                         bar(ctrs,nx,'b'), xlim([0 2*pi])
% %                         title(['Normal Distribution w/ \mu = ',num2str(mn{k}(j,t),2),' & 1/\kappa = \sigma = ',num2str(sig{k}(j,t),2),' & ',num2str(sigCI{k}(j,t,2),2)])
% %                         %
% %                         Z = normrnd(mn{k}(j,t),sig1{k}(j,t),size(cluster));
% %                         subplot(313), hist(Z,ctrs), xlim([0 2*pi])
% %                         title(['Normal Distribution w/ \mu = ',num2str(mn{k}(j,t),2),' & sqrt{2(1-R)} = \sigma = ',num2str(sig1{k}(j,t),2)])
% %                         %
% %                         A = normrnd(mn{k}(j,t),sig1{k}(j,t),size(cluster));
% %                         subplot(313), hist(A,ctrs), xlim([0 2*pi])
% %                         title(['Normal Distribution w/ \mu = ',num2str(mn{k}(j,t),2),' & sqrt{-2ln(R)} = \sigma = ',num2str(sig1{k}(j,t),2)])
% % 
% %                         keyboard
% %                     end
% 
% 
% 
% 
% 
%                 else % linear variable - eigenvectors
% 
% %                     t_95A = tinv([a_95/2  1-a_95/2],nA-1);           % T-Score for clusterA
% %                     t_95B = tinv([a_95/2  1-a_95/2],nB-1);           % T-Score for clusterB
% %                     
% %
% %                     try
% %                         [paramEstsA,paramCIsA] = truncNormStats_MLEest(inptT(clusterA),0,1);
% %                         [paramEstsB,paramCIsB] = truncNormStats_MLEest(inptT(clusterB),0,1);
% %                         sig{k}(c_pair_cnt,t,1) = paramEstsA(2);                                     % fitted std to truncated normal dist for clusterA
% %                         sig{k}(c_pair_cnt,t,2) = paramEstsB(2);                                     % fitted std to truncated normal dist for clusterB
% %     
% %         %                 sig_UB_distA = abs(paramEstsA(2) - paramCIsA(1,2));
% %         %                 sig_LB_distA = abs(paramEstsA(2) - paramCIsA(2,2));
% %         %                 sig_UB_distB = abs(paramEstsB(2) - paramCIsB(1,2));
% %         %                 sig_LB_distB = abs(paramEstsB(2) - paramCIsB(2,2));
% %         %                 if( sig_UB_dist - sig_LB_dist > 1e-6 )
% %         %                     disp('But I thought the Confidence Intervals were symmetric for STD')
% %         %                     keyboard
% %         %                 end
% %         %                 sigCI{k}(c_pair_cnt,t,1) = mean([sig_UB_distA, sig_LB_distA]);               % 95% Confidence interval for Standard Deviation of fit Normal Distribution
% %         %                 sigCI{k}(c_pair_cnt,t,2) = mean([sig_UB_distB, sig_LB_distB]);               % 95% Confidence interval for Standard Deviation of fit Normal Distribution
% %     
% %                         sigCI{k}(c_pair_cnt,t,1) = paramCIsA(1,2);                                  % lower bound (95% Confidence interval) for Standard Deviation of fit Normal Distribution
% %                         sigCI{k}(c_pair_cnt,t,2) = paramCIsA(2,2);                                  % upper bound (95% Confidence interval) for Standard Deviation of fit Normal Distribution
% %                         sigCI{k}(c_pair_cnt,t,3) = paramCIsB(1,2);                                  % lower bound (95% Confidence interval) for Standard Deviation of fit Normal Distribution
% %                         sigCI{k}(c_pair_cnt,t,4) = paramCIsB(2,2);                                  % upper bound (95% Confidence interval) for Standard Deviation of fit Normal Distribution
% %     
% %     
% %                     catch
% % 
% %                         mn{k}(c_pair_cnt,t,1) = mean(Ak_1D);             % cluster centers in 1D projection from multiple k-dims.
% %                         mn{k}(c_pair_cnt,t,2) = mean(Bk_1D);             % empirical mean for clusters
% %                         sig{k}(c_pair_cnt,t,1) = std(Ak_1D);             % standard deviations of distributions in 1D projection from multiple k-dims.
% %                         sig{k}(c_pair_cnt,t,2) = std(Bk_1D);             % empirical std (if truncated dist fitting fails)
% %                         
% %                         notTrunc(end+1,:) = [k,j];
% 
%     %                 end
% 
%                 end % switch for circ or linear case

            end % Loop over j - C (number of unique clusters in gt)
        end % Loop over i - C (number of unique clusters in gt)

    end % Loop over t (time points with unique phase distributions) (t=1 for spectral)

end % Loop over gt (number of ground truths - multiple for segmentations) 

% Compute weighted d' measures by distance between cluster means and avg cluster stds. 
% Note: weighting cluster pairs by relative number of nodes in that pair vs other pairs.
d_prime_kD_multEv = repmat( cluster_weights{1}, 1, pts)  .* AB_dist_kD      ./  AB_std_kD;
d_prime_1D_multEv = repmat( cluster_weights{1}, 1, pts)  .* ABk_dist_proj1D ./  ABk_std_proj1D;
d_prime_1D_singEv = repmat( cluster_weights{1}, 1, pts)  .* AB_dist_1D      ./  AB_std_1D;



%% Various plots for error and/or sanity checking in multi-dimensional space.
%if( size(inptT,2)>1 )

% if(circ)
%     do_plots = 1;
% else
%     do_plots=1;
% end




% (1). Plot difference in cluster means and cluster std's (averaged across all clusters & pairs)
       % for 2 cases: (a). each eigenvector used individually. (b). multiple eigenvectors in k-dim
       % space projected onto 1D line between cluster means.
if(do_plots & 0)
   figure,
   subplot(211), hold on,
   plot( mean(AB_dist_kD), 'r','Linewidth',2);
   plot( mean(ABk_dist_proj1D), 'b','Linewidth',2);
   plot( mean(AB_dist_1D), 'g','Linewidth',2);
   ylabel('\Delta \mu of Cluster pairs','Fontsize',18,'Fontweight','Bold')
   set(gca,'Fontsize',16,'Fontweight','Bold')
   title(['Cluster means (\Delta \mu) and STD''s (\sigma) accumulating Dimensions - ',method_str],'Fontsize',20,'Fontweight','Bold')
   grid on
   axis tight
    %
   subplot(212), hold on,
   plot( mean(AB_std_kD), 'r','Linewidth',2);
   plot( mean(ABk_std_proj1D), 'b','Linewidth',2);
   plot( mean(AB_std_1D), 'g','Linewidth',2); 
   ylabel('avg \sigma of Clusters','Fontsize',18,'Fontweight','Bold')
   xlabel('Dimension number (k)','Fontsize',18,'Fontweight','Bold')
   set(gca,'Fontsize',16,'Fontweight','Bold')
   legend({'k dims','k dims projected down to 1D','single dim'},'Location','East')
   grid on
   axis tight
end


% (2). Plot angle between angle between line connecting cluster centers
       % in k dimensions and vector [1,0,...,0] for each cluster pair
       % Should settle down to constant after first couple eigenvectors.
if(do_plots & 0)
   %
   figure, 
   subplot(411),hold on
   plot(1:pts, ang_1D_projK,'LineWidth',2)
   title({['\Delta \Theta settles when accumulating dimensions  - ',method_str],...
          'Is very different from [1,0,...,0] for accumulating.'},'Fontsize',20,'Fontweight','Bold')
   ylabel({'\Theta (in deg)'},'Fontsize',18,'Fontweight','Bold')
   text(0.5*pts,90,'Colors denote different cluster pairs','Fontsize',16,'Fontweight','Bold')
   set(gca,'Fontsize',16,'Fontweight','Bold') %,'YTick',[0,90,180,270,360])
   %
   subplot(412),hold on
   plot(1:pts, ang_cluster_vecs_kD ,'LineWidth',2)
   title({'\Theta between cluster mean vectors in k-dim as k increases.'},'Fontsize',20,'Fontweight','Bold')
   ylabel({'\Theta'},'Fontsize',18,'Fontweight','Bold')
   set(gca,'Fontsize',16,'Fontweight','Bold') %,'YTick',[0,90,180,270,360])
   %
   subplot(413),hold on
   plot(1:pts, ang_cluster_vecs_1D ,'LineWidth',2)
   title({'\Theta between cluster mean vectors in 1D.'},'Fontsize',20,'Fontweight','Bold')
   ylabel({'\Theta'},'Fontsize',18,'Fontweight','Bold')
   set(gca,'Fontsize',16,'Fontweight','Bold') %,'YTick',[0,90,180,270,360])
   %
   subplot(414), imagesc(inpt), 
   xlabel('Dimension number (k)','Fontsize',18,'Fontweight','Bold')
   ylabel('Node #','Fontsize',18,'Fontweight','Bold')
   title('Eigenvector or Osc Phase.','Fontsize',20,'Fontweight','Bold')
   set(gca,'Fontsize',16,'Fontweight','Bold')
   if(circ)
        colormap(hsv)
        colorbar('Location','West')
   else
        colormap(jet)
        colorbar('Location','East')
   end
   %
end


% (3). Plot d' measured from individual eigenvectors and from
       % accumilating eigenvectors from 1:k in k-dims and projecting the
       % down to 1D which is vector connecting cluster means.
if(do_plots)
   figure, 
   subplot(3,1,[1,2]), hold on,
   plot( mean(d_prime_1D_multEv) ,'b','Linewidth',2);
   %plot( mean(d_prime_1D_multEv)./sqrt(1:pts),'b--','Linewidth',2) % this is same as d_prime_1D_multEv_divK.
   plot( mean(d_prime_kD_multEv) ,'r','Linewidth',2);
   %plot( mean(d_prime_kD_multEv)./sqrt(1:pts),'r--','Linewidth',2) % this is same as d_prime_1D_multEv_divK.
   plot( mean(d_prime_1D_singEv) ,'g','Linewidth',2);
   title(['Avg d'' across all cluster pairs  - ',method_str],'Fontsize',20,'Fontweight','Bold')
   %xlabel('Dimension number (k)','Fontsize',18,'Fontweight','Bold')
   ylabel('< d'' >','Fontsize',18,'Fontweight','Bold')
   legend({'multiple k dims (projected to 1D)',...%'multiple k dims (projected to 1D) ./ \sqrt{k}',...
       'mean of k single dims',...%'mean of k single dims ./ \sqrt{k}',...
       'single 1-dim'},'Location','North')
   set(gca,'Fontsize',16,'Fontweight','Bold')
   grid on
   axis tight
   %
   subplot(313), imagesc(inpt), 
   xlabel('Dimension number (k)','Fontsize',18,'Fontweight','Bold')
   ylabel('Node #','Fontsize',18,'Fontweight','Bold')
   title('Eigenvector or Osc Phase.','Fontsize',20,'Fontweight','Bold')
   set(gca,'Fontsize',16,'Fontweight','Bold')
   if(circ)
        colormap(hsv)
        colorbar('Location','West','Fontsize',16,'Fontweight','Bold')
   else
        colormap(jet)
        colorbar('Location','East','Fontsize',16,'Fontweight','Bold')
   end
end



% (4). Plot Norm of A vs A_1D and B vs. B_1D for each i,j combo for all k's
if(do_plots & 0)
   figure, hold on,
   for i = 1:num_pairs
       plot(normVk(i,:),'r-','LineWidth',2)
       plot(normVk(i,:) ./ sqrt([1:pts]),'g-','LineWidth',2)
       plot(normVk_1D(i,:),'b-','LineWidth',2)
       plot(normV1(i,:),'k-','LineWidth',2)
       %
       plot(normAk(i,:),'r--','Linewidth',2)
       plot(normAk(i,:) ./ sqrt([1:pts]) ,'g--','Linewidth',2) %  normalizes for added distance that adding k'th dimension gets you.
       plot(normAk_1D(i,:),'b--','Linewidth',2)
       plot(normA1(i,:),'k--','Linewidth',2)
       %
       plot(normBk(i,:),'r--','Linewidth',2)
       plot(normBk(i,:) ./ sqrt([1:pts]) ,'g--','Linewidth',2) %  normalizes for added distance that adding k'th dimension gets you.
       plot(normBk_1D(i,:),'b--','Linewidth',2)
       plot(normB1(i,:),'k--','Linewidth',2)
   end
   %
   xlabel('Dimension number (k)','Fontsize',18,'Fontweight','Bold')
   ylabel('Frobenius Norm','Fontsize',18,'Fontweight','Bold')
   legend({'Multiple k dims',...
       'Multiple k dims /sqrt(k)',...
       'Multiple k dims 1D proj',...
       'Single 1 dim'},'Location','NorthWest')
   text(0,0.7*sqrt(pts),{'solid lines - all ntwrk nodes','dashed - nodes in clusters'},'Fontsize',16,...
       'Fontweight','Bold','HorizontalAlignment','Left','VerticalAlignment','Top')
   title(['Spatial Extent (Norm) of clusters & whole network - ',method_str],'Fontsize',20,'Fontweight','Bold')
   set(gca,'Fontsize',16,'Fontweight','Bold')
   grid on
   axis tight
end


% Meh, not necessary. Was looking at norm of each node for increasing dimensions to see if we could 
% segment just based on angle and not magnitude.  Turns out, kinda, we cant.
% % (5). Plot norm of vector pertaining to every node in k-dimensional subspace of k eigenvectors vs k.  
% if(do_plots)
%    mVb = mean(normVb);
%    sVb = std(normVb);
%    plt_inds = [1:N];
%    figure, 
%    subplot(211), imagesc(normVb), h=colorbar('Location','East');
%    %xlabel('Eigenvector number (k)','Fontsize',18,'Fontweight','Bold')
%    ylabel('Each Individual Node','Fontsize',18,'Fontweight','Bold')
%    title(h,'Norm','Fontsize',18,'Fontweight','Bold')
%    set(gca,'Fontsize',16,'Fontweight','Bold')
%    title('Note: Can consider ANGLE CLUSTERING (ignoring magnitude) when Node Norm   \sigma -> 0','Fontsize',20,'Fontweight','Bold')
%    %
%    subplot(212), hold on
%    plot(D./max(D),'g','Linewidth',2)
%    plot(plt_inds,mVb(plt_inds),'r--','Linewidth',2)
%    plot(plt_inds,mVb(plt_inds)+sVb(plt_inds),'r.','Linewidth',2)
%    plot(plt_inds,mVb(plt_inds)-sVb(plt_inds),'r.','Linewidth',2)
%    plot(ones(1,N),'k--','Linewidth',2)
%    text(N/10,0.92,'on the k-dim unit hypersphere','Fontsize',16,'Fontweight','Bold')
%    legend({'Normalized Eigenvalues','Node Norms (\mu & \sigma )'},'Location','East')
%    xlabel('Eigenvector number (k)','Fontsize',18,'Fontweight','Bold')
%    ylabel(['Norm of Vector (\mu & \sigma )'],'Fontsize',18,'Fontweight','Bold')
%    set(gca,'Fontsize',16,'Fontweight','Bold')
% end



if(do_plots & 1)
    keyboard
end




% (7). Plot all Eigenvectors in subplots of single figure (for trouble shooting)
if(0)
   figure,
   for i = 1:N
       for j = 1:num_clust
           ind = find(s==j);
           subplot(10,10,i), hold on
           scatter(ind,V(ind,i), 'MarkerEdgeColor', colors(j,:) )
       end
   end
end



% (8). Visualization / Error Checking plot of oscillator phases alongsize mean & stds.
if(0)
    skip_t = pts/100;
    figure, 
    subplot(211), hold on
    errorbar(repmat([1:skip_t:pts]',1,numel(C)),mn{1}(:,1:skip_t:pts)',sig{1}(:,1:skip_t:pts)')
    plot(repmat([1:skip_t:pts]',1,numel(C)),mn{1}(:,1:skip_t:pts)','LineWidth',2)
    xlabel('Simulation Time','FontSize',18,'FontWeight','Bold')
    ylabel('Phase','FontSize',18,'FontWeight','Bold')
    title('Cluster \mu''s & \sigma''s','FontSize',20,'FontWeight','Bold')
    set(gca,'ytick',[-pi:pi/2:pi],'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},'FontSize',16,'FontWeight','Bold','DefaultTextInterpreter','latex')
    xlim([0 pts])
    legend({'1','2','3','4'},'Location','EastOutside')
    %
    subplot(212)
    imagesc(inpt), colormap(hsv), colorbar
    set(gca,'FontSize',16,'FontWeight','Bold')
    
    keyboard
end



 
% x = randi(50, 1, 100);                      % Create Data
% SEM = std(x)/sqrt(length(x));               % Standard Error
% ts = tinv([0.025  0.975],length(x)-1);      % T-Score
% CI = mean(x) + ts*SEM;                      % Confidence Intervals