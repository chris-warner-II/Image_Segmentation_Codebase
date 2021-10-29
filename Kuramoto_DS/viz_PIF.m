function viz_PIF(PIFlg, PIFparams, saveDir)

% % For now, this is just a script to visualize the Phase Interaction
% % Functions for different types and parameter settings.

if ~exist([saveDir,'/PIF_plots/',PIFlg,'_params',num2str(PIFparams)],'file')
    
    x = linspace(-pi, pi, 100);            
    gam = pick_PIF(x, PIFlg, PIFparams);            % choose PIF here.
    normalization = max(abs(gam));                  % compute normalization 
    gam = gam ./normalization;                      % normalize PIF


    if ~exist([saveDir,'PIF_plots'],'dir')
        mkdir([saveDir,'PIF_plots'])
    end

    h=figure; hold on, grid on
    plot([-pi pi],[0 0],'r--','LineWidth',2)
    plot([0 0], [min(gam) max(gam)],'r--','LineWidth',2)
    axis([-pi pi min(gam) max(gam)])
    %
    plot(x,gam,'LineWidth',4)
    %
    xlabel(['Phase Difference (x = \theta_j - \theta_i)'],'FontSize',18,'FontWeight','Bold')
    ylabel({'Strength of PIF ( \Gamma_{ij}(x) )', '< - - - [ Delay Oscillator i Phase ] - - - [ Advance Oscillator i Phase ] - - - >'  },'FontSize',18,'FontWeight','Bold')
    title(['Coupled Oscillator PIF - ',PIFlg,' w/ params = (',num2str(PIFparams),')'],'FontSize',20,'FontWeight','Bold')
    set(gca,'FontSize',16,'FontWeight','Bold','Xtick',[-pi -pi/2 0 pi/2 pi],'XtickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
    text(pi/2,-0.8,'Osc i phase lags','FontSize',18,'FontWeight','Bold','HorizontalAlignment','Center')
    text(-pi/2,0.8,'Osc i phase leads','FontSize',18,'FontWeight','Bold','HorizontalAlignment','Center')

    saveGoodImg(h,[saveDir,'/PIF_plots/',PIFlg,'_params',num2str(PIFparams)],[0 0 1 1])
    close(h)

else
    disp('Phase Interaction Function image file exists already.  Not saving.')
end