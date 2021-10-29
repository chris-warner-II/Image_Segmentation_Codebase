% Script to play around with finding the parameters (mu, sig) for a normal
% distribution on a finite interval using MLE.

[dirPre,sizeGoodIm] = onCluster;


dirSave = [dirPre,'output/truncatedNormalFittingPlots/'];


if ~exist('dirSave','dir')
    mkdir(dirSave)
end



xTruncR=1;
xTruncL=0;

Ne_tgt = [500 400 300 200 100 50]; % Number of points in truncated distribution will be same always (set here)

mu_t  = [-0.1:0.1:1.1]'; % true distribution's parameter values to range over
sig_t = [0.01 0.10 0.20 0.30 0.40 0.50];


numIter=5;

for L = 1:numel(Ne_tgt)

    for k = 1:numIter
        
        [L,k]

        for i = 1:numel(mu_t)
            for j = 1:numel(sig_t)

                %[i,j]

                nogood = 0; % flag that becomes 1 if true distribution has nothing in truncated region.
                Ne = 0; % starting measure of how man points are in empirical/truncated distribution
                xe =[]; % actual values in truncated distribution
                n=0;    % counter thru while loop below


                while (Ne < Ne_tgt(L))

                    n=n+1;
                    xt(n) = normrnd(mu_t(i),sig_t(j));
                    xe = xt( (xt<=xTruncR) & (xt>=xTruncL) );
                    Ne = length(xe);


                    if (n>1000 & Ne==0)
                        nogood = 1;
                        break

                    end


                end

                if nogood
                    Nt(i,j) = nan;
                    mu_e(i,j) = nan;
                    sig_e(i,j) = nan;
                    mu_f(i,j) = nan;
                    sig_f(i,j) = nan;

                    continue
                end

                Nt(i,j) = length(xt);
                clear xt




                mu_e(i,j) = mean(xe);
                sig_e(i,j) = std(xe);



                if(0)
                    figure,
                    subplot(211),hist(xt); xlim([0 1])
                    subplot(212),hist(xe); xlim([0 1])
                end

                % Now find the fit parameters for the truncated distribution
                %try
                    [paramEsts1,paramCIs,acov,stderr] = truncNormStats_MLEest(xe,xTruncL,xTruncR);
                    
%                     [paramEsts2,paramCIs,acov,stderr] = truncNormStats_MLEest(xe,xTruncL,xTruncR);
%                     
%                     [paramEsts3,paramCIs,acov,stderr] = truncNormStats_MLEest(xe,xTruncL,xTruncR);
                %catch
                %   keyboard
                %end
                
                
%                 [paramEsts1(1),paramEsts2(1),paramEsts3(1), diff([paramEsts1(1),paramEsts2(1),paramEsts3(1)])]
%                 [paramEsts1(2),paramEsts2(2),paramEsts3(2), diff([paramEsts1(2),paramEsts2(2),paramEsts3(2)])]
                



                mu_f(i,j) = paramEsts1(1);
                sig_f(i,j) = paramEsts1(2);

                
                muLB_f(i,j) = paramCIs(1,1);
                muUB_f(i,j) = paramCIs(2,1);
                sigLB_f(i,j) = paramCIs(1,2);
                sigUB_f(i,j) = paramCIs(2,2);
                


                %disp('True / Empirical / Fit mu & sig')
                %[ [mu_t(i);sig_t(j)], [mu_e(i,j);sig_e(i,j)], [mu_f(i,j);sig_f(i,j)] ]
                %
                %disp('# measurements in true & empirical distribution')
                %[Nt(i,j), Ne]


                

            end
        end
        
        
        
        
        
        % Check if upper bound and lower bound are symetrically displaced from estimate.  
        % This stuff below assumes they are.
        mn_UB_dist = muUB_f - mu_f;
        mn_LB_dist = mu_f - muLB_f;
        %
        sig_UB_dist = sigUB_f - sig_f;
        sig_LB_dist = sig_f - sigLB_f;
        %
        if any( mn_UB_dist(:) - mn_LB_dist(:) > 1e-6 )
            disp('But I thought the Confidence Intervals were symmetric for Mean')
            keyboard
        end
        %
        if any( sig_UB_dist(:) - sig_LB_dist(:) > 1e-6 )
            disp('Confidence Intervals should not be symmetric for STD. But it is.')
            keyboard
        end


        sigCI_f = [sigUB_f - sig_f]; % this just looks at distance from middle to upper bound
        muCI_f =  [muUB_f - mu_f];   % but the distance from middle to lower is the same if symmetric. 
                
        
        
        
        
        
        
        
        



        % Find where true distribution was completely outside of truncation region.
        [x,y] = find(isnan(Nt));






        % Plot 1: Compare Empirical & Fitted STD's (how close they are to the true distribution.
        diff_sig_e = abs( repmat(sig_t,numel(mu_t),1) - sig_e );
        diff_sig_f = abs( repmat(sig_t,numel(mu_t),1) - sig_f );
        rat_N = Nt./Ne;
        [xx,yy] = find(rat_N==1);

        cbmax = max([diff_sig_e(:)]); % ;diff_sig_f(:)
        %cbmin = min([diff_sig_e(:)]); % ;diff_sig_f(:)
           
        [xxx, yyy] = find(diff_sig_f > cbmax); % find where diff_sig_f > cbmax and annotate those pointson imagesc plot.
        

        H=figure; 
        subplot(221), imagesc(diff_sig_e'), title('\Delta \sigma empirical','FontSize',18,'FontWeight','Bold'), caxis([0 cbmax]),
        colorbar('SouthOutside')
        hold on, scatter(x,y,'wx')
        set(gca,'YTick',[1:numel(sig_t)],'YTickLabel',sig_t,'XTick',[1:numel(mu_t)],'XTickLabel',mu_t,'FontSize',16,'FontWeight','Bold')
        ylabel('\sigma true'), xlabel('\mu true')
        %
        subplot(222), imagesc(diff_sig_f'), title('\Delta \sigma fitted','FontSize',18,'FontWeight','Bold'), caxis([0 cbmax]), 
        colorbar('SouthOutside')
        hold on, scatter(x,y,'wx')
        for i = 1:numel(xxx)
            text(xxx(i),yyy(i),['\color{white}{',num2str(diff_sig_f(xxx(i),yyy(i)),2),'}'],'HorizontalAlignment','Center')
        end
        set(gca,'YTick',[1:numel(sig_t)],'YTickLabel',sig_t,'XTick',[1:numel(mu_t)],'XTickLabel',mu_t,'FontSize',16,'FontWeight','Bold')
        ylabel('\sigma true'), xlabel('\mu true')
        %
        subplot(223), imagesc(rat_N'), title('ratio #pts','FontSize',18,'FontWeight','Bold'), 
        colorbar('SouthOutside')
        hold on, scatter(x,y,'wx'), scatter(xx,yy,'wo')
        set(gca,'YTick',[1:numel(sig_t)],'YTickLabel',sig_t,'XTick',[1:numel(mu_t)],'XTickLabel',mu_t,'FontSize',16,'FontWeight','Bold')
        ylabel('\sigma true'), xlabel('\mu true')
        %
        subplot(224), imagesc(sigCI_f'), title('\sigma fitted C.I.','FontSize',18,'FontWeight','Bold')
        colorbar('SouthOutside')
        hold on, scatter(x,y,'wx')
        set(gca,'YTick',[1:numel(sig_t)],'YTickLabel',sig_t,'XTick',[1:numel(mu_t)],'XTickLabel',mu_t,'FontSize',16,'FontWeight','Bold')
        ylabel('\sigma true'), xlabel('\mu true')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
                    ['# Pts Empirical = ',num2str(Ne_tgt(L))], ...
                    'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')


        saveGoodImg(H,[dirSave,'Sig_Ne',num2str(Ne_tgt(L)),' ',num2str(k)],sizeGoodIm)
        close(H)



        % Plot 2: Compare Empirical & Fitted MEAN's (how close they are to the true distribution.
        diff_mu_e = abs( repmat(mu_t,1,numel(sig_t)) - mu_e );
        diff_mu_f = abs( repmat(mu_t,1,numel(sig_t)) - mu_f );
        rat_N = Nt./Ne;
        [xx,yy] = find(rat_N==1);

        cbmax = max([diff_mu_e(:)]); % ;diff_mu_f(:)
        %cbmin = min([diff_mu_e(:)]); % ;diff_mu_f(:)
        
        [xxx, yyy] = find(diff_mu_f > cbmax); % find where diff_sig_f > cbmax and annotate those pointson imagesc plot.

        H=figure; 
        subplot(221), imagesc(diff_mu_e'),  title('\Delta \mu empirical','FontSize',18,'FontWeight','Bold'), caxis([0 cbmax]),
        colorbar('SouthOutside')
        hold on, scatter(x,y,'wx')
        set(gca,'YTick',[1:numel(sig_t)],'YTickLabel',sig_t,'XTick',[1:numel(mu_t)],'XTickLabel',mu_t,'FontSize',16,'FontWeight','Bold')
        ylabel('\sigma true'), xlabel('\mu true')
        %
        subplot(222), imagesc(diff_mu_f'), title('\Delta \mu fitted','FontSize',18,'FontWeight','Bold'), caxis([0 cbmax]),
        colorbar('SouthOutside')
        hold on, scatter(x,y,'wx')
        for i = 1:numel(xxx)
            text(xxx(i),yyy(i),['\color{white}{',num2str(diff_mu_f(xxx(i),yyy(i)),2),'}'],'HorizontalAlignment','Center')
        end
        set(gca,'YTick',[1:numel(sig_t)],'YTickLabel',sig_t,'XTick',[1:numel(mu_t)],'XTickLabel',mu_t,'FontSize',16,'FontWeight','Bold')
        ylabel('\sigma true'), xlabel('\mu true')
        %
        subplot(223), imagesc(rat_N'), title('ratio #pts','FontSize',18,'FontWeight','Bold'), 
        colorbar('SouthOutside')
        hold on, scatter(x,y,'wx'), scatter(xx,yy,'wo')
        set(gca,'YTick',[1:numel(sig_t)],'YTickLabel',sig_t,'XTick',[1:numel(mu_t)],'XTickLabel',mu_t,'FontSize',16,'FontWeight','Bold')
        ylabel('\sigma true'), xlabel('\mu true')
        %
        subplot(224), imagesc(muCI_f'), title('\mu fitted C.I.','FontSize',18,'FontWeight','Bold')
        colorbar('SouthOutside')
        hold on, scatter(x,y,'wx')
        set(gca,'YTick',[1:numel(sig_t)],'YTickLabel',sig_t,'XTick',[1:numel(mu_t)],'XTickLabel',mu_t,'FontSize',16,'FontWeight','Bold')
        ylabel('\sigma true'), xlabel('\mu true')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
                    ['# Pts Empirical = ',num2str(Ne_tgt(L))], ...
                    'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')


        saveGoodImg(H,[dirSave,'Mu_Ne',num2str(Ne_tgt(L)),' ',num2str(k)],sizeGoodIm)
        close(H);

    end

end


