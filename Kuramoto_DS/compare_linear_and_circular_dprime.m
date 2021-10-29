% This script serves as a sanity check of our d' measure.  We will
% construct a pair of distributions (denoting 2 different clusters) in the
% linear space contained between [0,1].  We can convert this pair of
% distributions into the equivalent circular pair by scaling values by 2pi.
% (For simplicity, at least for now, we only consider distributions fully
% contained within the [0,1]-interval and with no wrap-around the circular
% interval).  We can then compute means & STDs of the 2 clusters in both
% linear and circular cases.  From those, we can compute d' in both cases
% and expect to find similar d' values.  If not, perhaps, we have
% discovered a bias between linear & circular conceptions that we need to
% correct for.
%
% This also looks at d' in both circular & linear regimes as a function of
% hits and false alarms at the ideal segmentation threshold. So, we find
% the threshold where the 2 clusters' distributions intersect.  We then
% construct the cumulative distributions (CDF) of each cluster and count up
% hits (points in each cluster on the correct side of the threshold) and
% false alarms (points on the wrong side of the threshold).  


del_mu_loop = 0.3; %[0.1:0.1:1];
sig_loop = 0.1; % [0.01, 0.1:0.1:2];

num_samples = 1; %10;

plt_flg = 1; % flag to make plots of PDF's & CDF's

for i = 1:numel(del_mu_loop)
    for j = 1:numel(sig_loop)
        for k = 1:num_samples
        
        disp([num2str(i),' / ',num2str(numel(del_mu_loop)),' : ',num2str(j),' / ',num2str(numel(sig_loop)),' : ',num2str(k),' / ',num2str(num_samples)])

        %% (0). User Input to Construct Distributions:

        bins = 1000; % # of bins to chop interval [0,1] or [0,2pi] into.

        xbins = linspace(0,1,bins);

        N = [1000, 20000];                      % number of units in cluster (for now same in both)
        del_mu = del_mu_loop(i);
        mu = [ 0.5-del_mu./2, 0.5+del_mu./2 ];   % means of 2 clusters (between 0 & 1) Note: mean1 <= mean2.  
        sig = [ sig_loop(j), sig_loop(j) ];      % std's of 2 cluster

        %% (1). Construct 2 distributions in linear space between [0,1], find optimal 
        % threshold & compute hit rate & false alarm rate at that threshold.
        c1 = normrnd(mu(1),sig(1),1,N(1));
        c2 = normrnd(mu(2),sig(2),1,N(2));

        n1 = hist(c1,xbins) ./ sum(N);
        n2 = hist(c2,xbins) ./ sum(N);
        
        ncs1 = cumsum(n1);
        ncs2 = cumsum(n2);

        intersct = abs(N(1)./sum(N)-ncs1-ncs2);
        xovr = find(intersct==min(intersct));

        if (numel(xovr)>1)
            xovr = xovr(round(numel(xovr)/2)); % take threshold in middle if there are many values with minimum intersct value
        end

        hit = sum(n1(1:xovr)) + sum(n2((xovr+1):end)); % probability of correct classification
        fa = sum(n2(1:xovr)) + sum(n1((xovr+1):end)); % probability of incorrect classification

        d_prime1(i,j,k) = fa; % hit - fa;
        
        if(plt_flg)
            figure, 
            subplot(211), hold on
            plot(xbins,n1,'b')
            plot(xbins,n2,'g')
            plot( [xbins(xovr),xbins(xovr)],[0,max([n1,n2])],'r--')
            title('UnTruncated Linear PDF')
            xlim([0,1])
            %
            subplot(212), hold on
            plot(xbins,N(1)./sum(N)-ncs1,'b')
            plot(xbins,ncs2,'g')
            plot( [xbins(xovr),xbins(xovr)],[0,1],'r--')
            text( xbins(round(xovr/2)), 0.5, ['\color{blue}{P_{hit} = ',num2str( sum(n1(1:xovr)), 2 ),'}'] )
            text( xbins(round(xovr/2)), 0.3, ['\color{green}{P_{fa} = ',num2str(  sum(n2(1:xovr)), 2),'}'] )
            text( xbins(round((xovr+bins)/2)), 0.5, ['\color{green}{P_{hit} = ',num2str( sum(n2((xovr+1):end)), 2),'}'] )
            text( xbins(round((xovr+bins)/2)), 0.3, ['\color{blue}{P_{fa} = ',num2str( sum(n1((xovr+1):end)), 2 ),'}'] )
            title('UnTruncated Linear CDF')
            xlim([0,1])
            ylim([0,1])
            
%             % find crossover point of PDFs (not CDFs) - because not sure they are the same thing.
%             figure, plot(abs(n1 - n2))
        end
        


        %% (2). Compute d' in linear for untruncated normal distribution parameters
        d_prime2(i,j,k) = abs(mu(1) - mu(2)) ./ (0.5*sqrt( sig(1).^2 + sig(2).^2 ));
        % computed from the parameters used to construct the two distributions




        %% (3). Compute d' in linear for truncated normal distribution parameters (Here using Hits & False Alarms)

        xTruncL = 0;
        xTruncR = 1;

        ct1 = c1( (c1<=xTruncR) & (c1>=xTruncL) );
        Nt(1) = length(ct1);

        ct2 = c2( (c2<=xTruncR) & (c2>=xTruncL) );
        Nt(2) = length(ct2);

        truncation_loss(i,j,k) = ( mean(N) - mean(Nt) ) ./ mean(N);

        nt1 = hist(ct1,xbins)./Nt(1);
        nt2 = hist(ct2,xbins)./Nt(2);

        ntcs1 = cumsum(nt1);
        ntcs2 = cumsum(nt2);

        intersctt = abs(1-ntcs1-ntcs2);
        xovrt = find(intersctt==min(intersctt));

        if (numel(xovrt)>1)
            xovrt = xovrt(round(numel(xovrt)/2)); % take threshold in middle if there are many values with minimum intersct value
        end

        hit = ( sum(n1(1:xovrt)) + sum(n2((xovrt+1):end)) ) ./ 2;  % probability of correct classification
        fa = ( sum(n2(1:xovrt)) + sum(n1((xovrt+1):end)) ) ./ 2; % probability of incorrect classification

        d_prime1t(i,j,k) = fa; %hit - fa;
        
        
        if(plt_flg)
            figure, 
            subplot(211), hold on
            plot(xbins,nt1,'b')
            plot(xbins,nt2,'g')
            plot( [xbins(xovrt),xbins(xovrt)],[0,max([nt1,nt2])],'r--')
            title('Truncated Linear PDF')
            xlim([0,1])
            %
            subplot(212), hold on
            plot(xbins,1-ntcs1,'b')
            plot(xbins,ntcs2,'g')
            plot( [xbins(xovrt),xbins(xovrt)],[0,1],'r--')
            text( xbins(round(xovrt/2)), 0.5, ['\color{blue}{P_{hit} = ',num2str( sum(nt1(1:xovrt)), 2 ),'}'] )
            text( xbins(round(xovrt/2)), 0.3, ['\color{green}{P_{fa} = ',num2str(  sum(nt2(1:xovrt)), 2),'}'] )
            text( xbins(round((xovrt+bins)/2)), 0.5, ['\color{green}{P_{hit} = ',num2str( sum(nt2((xovrt+1):end)), 2),'}'] )
            text( xbins(round((xovrt+bins)/2)), 0.3, ['\color{blue}{P_{fa} = ',num2str( sum(nt1((xovrt+1):end)), 2 ),'}'] )
            title('Truncated Linear CDF')
            xlim([0,1])
            ylim([0,1])
        end



        %% (4). Compute d' from mean & std computed using ML fit estimates with truncated distribution.

        [params_ct1] = truncNormStats_MLEest(ct1,xTruncL,xTruncR); % ML estimates of mn & std for truncated cluster #1
        [params_ct2] = truncNormStats_MLEest(ct2,xTruncL,xTruncR); % ML estimates of mn & std for truncated cluster #1

        mn_ct1 = mean(ct1); % empirical mean of truncated cluster #1
        mn_ct2 = mean(ct2); % empirical mean of truncated cluster #2


        d_prime2t(i,j,k) = abs(mn_ct1 - mn_ct2) ./ (0.5*sqrt( params_ct1(2).^2 + params_ct2(2).^2 ));


        %% (5). Convert distributions over to circular space between [0,2pi], find
        % optimal threshold & compute hit rate & false alarm rate at that threshold.

        cc1 = wrapTo2Pi(2*pi*c1); % these are the untruncated linear distributions
        cc2 = wrapTo2Pi(2*pi*c2); % and so will have wrap-around.

        xbins = linspace(0,2*pi,bins);

        nC1 = hist(cc1,xbins)./N(1);
        nC2 = hist(cc2,xbins)./N(2);

        nCcs1 = cumsum(nC1);
        nCcs2 = cumsum(nC2);

        

        intersct = abs(1-nCcs1-nCcs2);
        xovrC1 = find(intersct==min(intersct));

        if (numel(xovrC1)>1)
            xovrC1 = xovrC1(round(numel(xovrC1)/2)); % take threshold in middle if there are many values with minimum intersct value
        end


        % Need to set another threshold between the two distributions where the circular space wraps around.
        % How: Reorder distributions to begin from the first xovr point.  Then do
        % same analysis of minimum intersection of CDFs
        nr1 = [nC1(xovrC1+1:end),nC1(1:xovrC1)];
        nr2 = [nC2(xovrC1+1:end),nC2(1:xovrC1)];
        xrbins = [xbins(xovrC1+1:end),xbins(1:xovrC1)];

        nrcs1 = cumsum(nr1);
        nrcs2 = cumsum(nr2);

        intersctr = abs(1-nrcs1-nrcs2);
        xovrr = find(intersctr==min(intersctr));

        if (numel(xovrr)>1)
            xovrr = xovrr(round(numel(xovrr)/2)); % take threshold in middle if there are many values with minimum intersct value
        end

        xovrC2 = find(xbins == xrbins(xovrr));

%         if(plt_flg)
%             figure, 
%             subplot(211), hold on
%             plot(nr1,'b')
%             plot(nr2,'g')
%             title('Circular PDF #2')
%             xlim([1,bins])
%             % set(gca,'XTick',[1,bins])
%             %
%             subplot(212), hold on
%             plot(1-nrcs2,'g')
%             plot(nrcs1,'b')
%             title('Circular CDF #2')
%             xlim([0,bins])
%         end
        
        

        if (xovrC1 > xovrC2) % if 2nd threshold is between 0 and 1st one
            inClust1 = nC1((xovrC2+1):xovrC1);
            outClust1 = [nC1(1:xovrC2),nC1((xovrC1+1):end)];
            inClust2 = [nC2(1:xovrC2),nC2((xovrC1+1):end)];
            outClust2 = nC2((xovrC2+1):xovrC1);
        else                 % if 2nd threshold is between 1st one and 2pi
            inClust1 = [nC1(1:xovrC1),nC1((xovrC2+1):end)];
            outClust1 = nC1((xovrC1+1):xovrC2);
            inClust2 = nC2((xovrC1+1):xovrC2);
            outClust2 = [nC2(1:xovrC1),nC2((xovrC2+1):end)];
        end


        % Need to still use 2nd threshold to determine hits & false alarms.  Not done yet.
        hit = ( sum(inClust1) + sum(inClust2) ) ./ 2;
        fa = ( sum(outClust1) + sum(outClust2) ) ./ 2;

        d_prime3(i,j,k) = fa; % hit - fa;
        
        
        
        if(plt_flg)
            figure, 
            subplot(211), hold on
            plot(xbins,nC1,'b')
            plot(xbins,nC2,'g')
            plot( [xbins(xovrC1),xbins(xovrC1)],[0,max([nC1,nC2])],'r--')
            plot( [xbins(xovrC2),xbins(xovrC2)],[0,max([nC1,nC2])],'r--')
            title('Circular PDF')
            xlim([0,2*pi])
            %
            subplot(212), hold on
            plot(xbins,1-nCcs1,'b')
            plot(xbins,nCcs2,'g')
            plot( [xbins(xovrC1),xbins(xovrC1)],[0,1],'r--')
            plot( [xbins(xovrC2),xbins(xovrC2)],[0,1],'r--')
            text( xbins(round(xovrC1/2)), 0.5, ['\color{blue}{P_{hit} = ',num2str( sum( inClust1 ), 2 ),'}'] )
            text( xbins(round(xovrC1/2)), 0.3, ['\color{green}{P_{fa} = ',num2str(  sum( outClust1 ), 2),'}'] )
            text( xbins(round((xovrC1+bins)/2)), 0.5, ['\color{green}{P_{hit} = ',num2str( sum( inClust2 ), 2),'}'] )
            text( xbins(round((xovrC1+bins)/2)), 0.3, ['\color{blue}{P_{fa} = ',num2str( sum( outClust2 ), 2 ),'}'] )
            
            title('Circular CDF')
            xlim([0,2*pi])
            ylim([0,1])
        end
        
        
        


        %% (6). Find mean & stds of circular distributions & compute d'

        r1 = circ_r(cc1');
        r2 = circ_r(cc2');

        kappa1 = circ_kappa(cc1);
        kappa2 = circ_kappa(cc2);

        mn1 = circ_mean(cc1');
        mn2 = circ_mean(cc2');

        [s01, s11] = circ_std(cc1');
        [s02, s12] = circ_std(cc2');

        s21 = 1./sqrt(kappa1);
        s22 = 1./sqrt(kappa2);

        d_prime4_0(i,j,k) = abs(circ_dist(mn1,mn2)) ./ ( 0.5*sqrt(s01.^2 + s02.^2) );
        d_prime4_1(i,j,k) = abs(circ_dist(mn1,mn2)) ./ ( 0.5*sqrt(s11.^2 + s12.^2) );
        d_prime4_2(i,j,k) = abs(circ_dist(mn1,mn2)) ./ ( 0.5*sqrt(s21.^2 + s22.^2) );
        
        % keyboard

        end % loop over number of samples
        
    end % loop over sigma
    
end % loop over delta mu




% Compute Mean & STD values for different d_prime measurements averaged
% over the number of samples taken to get errorbars for plots.
mn_dp1 = mean(d_prime1,3);
std_dp1 = std(d_prime1,[],3);
%
mn_dp1t = mean(d_prime1t,3);
std_dp1t = std(d_prime1t,[],3);
%
mn_dp2 = mean(d_prime2,3);
std_dp2 = std(d_prime2,[],3);
%
mn_dp2t = mean(d_prime2t,3);
std_dp2t = std(d_prime2t,[],3);
%
mn_dp1 = mean(d_prime1,3);
std_dp1 = std(d_prime1,[],3);
%
mn_dp3 = mean(d_prime3,3);
std_dp3 = std(d_prime3,[],3);
%
mn_dp4_0 = mean(d_prime4_0,3);
std_dp4_0 = std(d_prime4_0,[],3);
%
mn_dp4_1 = mean(d_prime4_1,3);
std_dp4_1 = std(d_prime4_1,[],3);
%
mn_dp4_2 = mean(d_prime4_2,3);
std_dp4_2 = std(d_prime4_2,[],3);
%


%% Plot del_mu vs. sig vs. d' computed the 5 different ways (also the 2 ways which are: hit - false alarm)

if(0)
    for i = 1:numel(sig_loop)

        figure, hold on,
        errorbar(1:numel(del_mu_loop),mn_dp2(:,i),sqrt(std_dp2(:,i)),'k-','LineWidth',2)
        errorbar(1:numel(del_mu_loop),mn_dp2t(:,i),sqrt(std_dp2t(:,i)),'b-','LineWidth',2)
        errorbar(1:numel(del_mu_loop),mn_dp4_0(:,i),sqrt(std_dp4_0(:,i)),'r-','LineWidth',2)
        errorbar(1:numel(del_mu_loop),mn_dp4_1(:,i),sqrt(std_dp4_1(:,i)),'g-','LineWidth',2)
        errorbar(1:numel(del_mu_loop),mn_dp4_2(:,i),sqrt(std_dp4_2(:,i)),'c-','LineWidth',2)

        % indicated percentage of points lost in truncation as a measure of how
        % different the circular is from the linear since ciruclar has
        % wrap-around and linear has strict truncation or lopping off.
        for j = 1:numel(del_mu_loop)
            text(j,d_prime2t(j,i),[num2str(truncation_loss(j,i),2)],'FontSize',16,'FontWeight','Bold')
        end

        legend({'Linear UnTruncated','Linear Truncated','Circ w/ std = \sqrt{2(1-r)}','Circ w/ std = \sqrt{-Ln(r)}','Circ w/ std = 1/\sqrt{\kappa}'},'Location','NorthWest')

        title(['Cluster spread or std (\sigma = ',num2str(sig_loop(i)),')'],'FontSize',20,'FontWeight','Bold')
        xlabel(['Distance between Cluster Centers (\Delta \mu)'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' - metric'],'FontSize',18,'FontWeight','Bold')

        set(gca,'XTick',1:numel(del_mu_loop),'XTickLabel',del_mu_loop,'FontSize',16,'FontWeight','Bold')

    end
end

% Now plot, sig on x axis for del_mu set.
if(0)
    for i = 1:5 %numel(del_mu_loop)

        figure, hold on,
        errorbar(1:numel(sig_loop),mn_dp2(i,:),std_dp2(i,:),'k-','LineWidth',2,'Markers',18)
        errorbar(1:numel(sig_loop),mn_dp2t(i,:),std_dp2t(i,:),'b-','LineWidth',2,'Markers',18)
        errorbar(1:numel(sig_loop),mn_dp4_0(i,:),std_dp4_0(i,:),'rx','LineWidth',2,'Markers',18)
        errorbar(1:numel(sig_loop),mn_dp4_1(i,:),std_dp4_1(i,:),'gx','LineWidth',2,'Markers',18)
        errorbar(1:numel(sig_loop),mn_dp4_2(i,:),std_dp4_2(i,:),'cx','LineWidth',2,'Markers',18)

        % indicated percentage of points lost in truncation as a measure of how
        % different the circular is from the linear since ciruclar has
        % wrap-around and linear has strict truncation or lopping off.
        for j = 1:numel(sig_loop)
            text(j,d_prime2t(i,j),[num2str(truncation_loss(i,j),2)],'FontSize',16,'FontWeight','Bold')
        end

        legend({'Linear UnTruncated','Linear Truncated','Circ w/ std = \sqrt{2(1-r)}','Circ w/ std = \sqrt{-Ln(r)}','Circ w/ std = 1/\sqrt{\kappa}'},'Location','NorthEast')

        title(['Distance Between Cluster Centers (\Delta \mu = ',num2str(del_mu_loop(i)),')'],'FontSize',20,'FontWeight','Bold')
        xlabel(['Cluster Spread or std (\sigma)'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' - metric'],'FontSize',18,'FontWeight','Bold')

        set(gca,'XTick',1:numel(sig_loop),'XTickLabel',sig_loop,'FontSize',16,'FontWeight','Bold')

    end
end


% do same plots for d' determined by Z(hit) - Z(F.A.) where Z(x) is the CDF

if(1)
    for i = 1:numel(sig_loop)

        figure, hold on,
        
        errorbar(1:numel(del_mu_loop),mn_dp1(:,i),sqrt(std_dp1(:,i)),'ko-','LineWidth',2)
        errorbar(1:numel(del_mu_loop),mn_dp1t(:,i),sqrt(std_dp1t(:,i)),'bx-','LineWidth',2)
        errorbar(1:numel(del_mu_loop),mn_dp3(:,i),sqrt(std_dp3(:,i)),'rs-','LineWidth',2)
        
%         plot(d_prime1(:,i),'ko-','LineWidth',2)
%         plot(d_prime1t(:,i),'bx-','LineWidth',2)
%         plot(d_prime3(:,i),'rs-','LineWidth',2)

        % indicated percentage of points lost in truncation as a measure of how
        % different the circular is from the linear since ciruclar has
        % wrap-around and linear has strict truncation or lopping off.
        for j = 1:numel(del_mu_loop)
            text(j,mn_dp1(j,i),[num2str(truncation_loss(j,i),2)],'FontSize',16,'FontWeight','Bold')
        end

        legend({'Linear UnTruncated','Linear Truncated','Circular'},'Location','NorthWest')

        title(['Cluster spread or std (\sigma = ',num2str(sig_loop(i)),')'],'FontSize',20,'FontWeight','Bold')
        xlabel(['Distance between Cluster Centers (\Delta \mu)'],'FontSize',18,'FontWeight','Bold')
        ylabel(['Prob Misclassification'],'FontSize',18,'FontWeight','Bold')

        set(gca,'XTick',1:numel(del_mu_loop),'XTickLabel',del_mu_loop,'FontSize',16,'FontWeight','Bold')

        ylim([0,0.5])

    end
end
