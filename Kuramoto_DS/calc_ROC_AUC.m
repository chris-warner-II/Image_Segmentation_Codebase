function [ROCparams, best, fast, stab, AUCm2, chng2] = calc_ROC_AUC()


% This function takes in the output from main_Kuramoto function which
% includes results from segmentation at each time step of the Kuramoto
% simulation, for an array of cosine distance threshold values and a number
% of runs.  This uses results (true pos, false pos, true neg, false neg) to
% compute Sensitivity and Selectivity metrics and uses them to compute a
% Receiver Operating Characteristic (ROC) Curve at each time point of the
% simulation.  This yields a curve of AUC vs. t for each parameter ensemble 
% telling how accurate a segmentation from the current phase clustering of 
% oscillators would be.  From this curve, we extract two bits of
% information: (1). The maximum AUC and (2). The time it takes for the
% network to reach 95% of that maximum value.  These two values tell us
% about the quality and speed of segmentation.  We save these results along
% with the parameters for each run that produced them. On top of this we
% save the time course of AUC measure and chng measure (Hamming distance
% between segmentation at previous time using threshold and current
% segmentation).



ROCmovFlg = 0;        % Flag to plot a movie of the ROC curve as it changes with time.
AUCvT_singParams = 0; % Flag to plot images of AUC vs T for single parameters settings.


grey = [0.6 0.6 0.6]; % 3x1 color vector - grey line for plotting.


% Directory to save output plots into...
imgsDir = ['/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/',...
    'Projects/output/Kuramoto/HandCookedNetwork/imgs/SegStatsParamSearch/'];
if ~exist(imgsDir,'dir') 
    mkdir(imgsDir)
end


% Directory to save Single Parameter AUC vs Time Plots
if ( AUCvT_singParams && ~exist([imgsDir,'AUCvT_singParams/'],'dir') )
    mkdir([imgsDir,'AUCvT_singParams/'])
end

% Directory to look for data output from Kuramoto main to analyze here
dataKurDir = ['/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/',...
    'Projects/output/Kuramoto/HandCookedNetwork/data/SegResKur/'];
%
dataEigDir = ['/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/',...
    'Projects/output/Kuramoto/HandCookedNetwork/data/SegResEig/'];


k=0; % a counter to keep track of number of files.

subsetSpec = 'N1x48'; % to look at only certain files or directories

dirsK = dir([dataKurDir,'*',subsetSpec,'*']);








% 1st Loop through directories containing different clusterings of parameters we searched over.
for hh = 1:numel(dirsK)
    
    % Load data file saved from main_Kuramoto_HandCookedNetwork
    if ~isempty(strfind(dirsK(hh).name,'PIF_'))
        files = dir([dataKurDir,'/',dirsK(hh).name]);
    else
        files = [];
        continue
    end
    
    

    % 2nd Loop through output files saved from main_Kuramoto code and plot statistics
    for i = 1:numel(files)


        % Load data file saved from main_Kuramoto_HandCookedNetwork
        if strfind(files(i).name,'.mat')
            load([dataKurDir,'/',dirsK(hh).name,'/',files(i).name]);
            k=k+1
        else
            continue
        end


        %% 
        tp = segres.tp; % true positive: 
        fp = segres.fp; % false positive:
        tn = segres.tn; % true negative:
        fn = segres.fn; % false negative:
        chng = segres.chng;
        chng(:,1,:) = 0;


%         % Get rid of last time point (something weird.  no big deal)
%         tp = tp(:,(1:end-1),:);
%         fp = fp(:,(1:end-1),:);
%         tn = tn(:,(1:end-1),:);
%         fn = fn(:,(1:end-1),:);
%         chng = chng(:,(1:end-1),:);



%         if ~isfield(runParams,'PconnFar')
%             runParams.PconnFar = 0; % if not other wise stated probability of a connection out beyond Rmax is zero.
%         end

        %
        thresh_cosDist_seg = linspace(-1,1,runParams.thresh_seg);
        for j = 1:runParams.thresh_seg
            threshCell{j} = num2str(thresh_cosDist_seg(j),2);
        end
        
        


        %
        runParamsTag = {['PIF ',runParams.PIFlg,' N',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),' - C',num2str(runParams.C)],...
          [' - rmax',num2str(runParams.rmax),' - NF ',num2str(runParams.muW),' ',num2str(runParams.sigW)],[' - Win ',num2str(runParams.Strng),' ',...
                    num2str(runParams.sigSg),' - Wout ',num2str(runParams.Weak),' ',num2str(runParams.sigWk)]};


        %% ROC Curve using Sensitivity & Specificity !!

        sens = tp ./ (tp+fn); % 
        spec = tn ./ (fp+tn); % 


        % calculating AUC of ROC individually on each run to get error bars
        for r = 1:runParams.runs

            xr = sens(:,:,r); 
            yr = spec(:,:,r);

            if(ROCmovFlg)
                ROCwriterObj = VideoWriter([imgsDir,'ROC_time_lapse_',num2str(r),'.avi']);
                ROCwriterObj.FrameRate = 2;
                open(ROCwriterObj);
                h=figure;
            end
            %
            for j = 1:size(xr,2)
                xpr = [1, xr(:,j)', 0];
                ypr = 1 - [0, yr(:,j)', 1];
                AUCr(j,r) = abs( trapz(ypr,xpr) );

                % to plot ROC curve that changes with time.
                if(ROCmovFlg) % can make a movie out of it.
                    
                    figure(h), subplot(3,1,1:2), hold on
                    plot(ypr,xpr,'b','LineWidth',2) 
                    plot([0 1],[0 1],'r--','LineWidth',2)
                    text(0.02, 1, 'Best','Color','red','VerticalAlignment','Top','FontSize',16,'FontWeight','Bold')
                    text(0.5, 0.5, 'Uninformative','Color','red','VerticalAlignment','Top','FontSize',16,'FontWeight','Bold')
                    text(0.6, 0.2, runParamsTag,'VerticalAlignment','Top','FontSize',16,'FontWeight','Bold')
                    title(['Run ',num2str(r),' Time ',num2str(j),' ROC Curve'],'FontSize',18,'FontWeight','Bold')
                    ylabel('Sensitivity','FontSize',18,'FontWeight','Bold'), 
                    xlabel('1 - Specificity','FontSize',18,'FontWeight','Bold')
                    set(gca,'FontSize',18,'FontWeight','Bold')
                    set(h,'Units','Normalized','Position',[0.1 0.1 0.5 0.9])
                    %
                    subplot(3,1,3), hold on, cla
                    plot([0.5 size(xr,2)], [0.5 0.5], 'r--')
                    text(size(xr,2)./2, 0.5, 'Uninformative','Color','red','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold')
                    text(size(xr,2)./2, 1, 'Best','Color','red','VerticalAlignment','Top','FontSize',16,'FontWeight','Bold')
                    plot(AUCr(1:j,r),'b','LineWidth',2), 
                    xlabel('Time - Period of 60Hz osc.','FontSize',18,'FontWeight','Bold')
                    ylabel('Area Under ROC Curve','FontSize',18,'FontWeight','Bold')
                    axis([0 size(xr,2) 0.4 1])
                    set(gca,'FontSize',16,'FontWeight','Bold')
                    %
                    mov = getframe(h);
                    writeVideo(ROCwriterObj,mov);
                    %
                    subplot(3,1,1:2), plot(ypr,xpr,'Color',grey,'LineWidth',2)
                    
                end
            end
            
            
            if(ROCmovFlg)
                close(ROCwriterObj);
            end
            
%             keyboard
            
            
        end
        
        if(1)
%             h=figure; hold on, shadedErrorBar(1:size(AUCr,1), mean(AUCr'), std(AUCr'))
        end
        
        
%         % HOW IT WAS.  THIS IS FINE IF I DONT WANT ERRORBARS        
%         x = mean(sens,3); 
%         y = mean(spec,3);
% 
%         %figure,
%         for j = 1:size(x,2)
%             xp = [1, x(:,j)', 0];
%             yp = 1 - [0, y(:,j)', 1];
%             AUC(j) = abs( trapz(yp,xp) );
%             
%             % to plot ROC curve that changes with time.
%             if(0) % can make a movie out of it.
%                 plot(yp,xp), title(num2str(j))
%                 pause
%             end
%             
%         end


        % Parameters attached to performance (may need to grab more).
        %
        ROCparams.Strng(k) = runParams.Strng;
        ROCparams.Weak(k) = runParams.Weak;
        ROCparams.C(k) = runParams.C;
        ROCparams.sigW(k) = runParams.sigW;
        ROCparams.rmax(k) = runParams.rmax;
        ROCparams.pfar(k) = runParams.PconnFar;
        
        AUC = mean(AUCr');

        AUCm2{k} = mean(AUCr'); % Curve of AUC vs T averaged over simulations for each parameter setting
        AUCs2{k} = std(AUCr'); % Curve of AUC vs T averaged over simulations for each parameter setting
        chng2{k} = mean(mean(chng,3)); % avg first across runs, then across thresholds
        %
        best(k) = max(AUC);
        bestestT = find(AUC == max(AUC) );
        bestT(k) = bestestT(1);
        %
        fastest = find(AUC > 0.98*max(AUC));
        try
            fast(k) = fastest(1);
        catch
            fast(k) = find(AUC == max(AUC));
        end
        %
        stab(k) = chng2{k}(fast(k)-1);
        %
        %
         dif = diff(AUCm2{k});
         cut = find(dif < 0);
         try
            tFlat(k) = cut(1);
         catch
            tFlat(k) = numel(dif);
         end
         qFlat(k) = AUCm2{k}(tFlat(k));



        
        
        
        

        % plot for single parameter setting, time course of AUC
        if(AUCvT_singParams)
            h=figure; hold on, 
            shadedErrorBar(1:size(AUCr,1), mean(AUCr'), std(AUCr'))
            plot(AUC,'k','LineWidth',2),
            scatter(fast(k),AUC(fast(k)),100,'gx','LineWidth',2)
            text(fast(k),0.98*AUC(fast(k)),'t_{98%}','Color','g','FontSize',16,'FontWeight','Bold','VerticalAlignment','Top')
            scatter(bestT(k),best(k),100,'ro','LineWidth',2)
            text(bestT(k),1.02*best(k),'Best','Color','r','FontSize',16,'FontWeight','Bold','VerticalAlignment','Bottom')
            title([runParamsTag],'FontSize',20,'FontWeight','Bold') % ['File# ',num2str(i)],
            xlabel('Time - period of 60Hz Osc (60=1sec)','FontSize',18,'FontWeight','Bold')
            ylabel('Area under ROC Curve (100 runs)','FontSize',18,'FontWeight','Bold')
            axis([1 numel(AUC) 0.5 1])
            set(gca,'FontSize',16,'FontWeight','Bold')
            grid on
            %
            % Save image
            rpt = ['PIF_',runParams.PIFlg,'_N',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),'_C',...
                num2str(runParams.C),'_rmax',num2str(runParams.rmax),'_NF_',num2str(runParams.muW),'_',...
                num2str(runParams.sigW),'_Win ',num2str(runParams.Strng),'_',num2str(runParams.sigSg),...
                '_Wout ',num2str(runParams.Weak),'_',num2str(runParams.sigWk)];
            saveGoodImg(h,[imgsDir,'AUCvT_singParams/',rpt],[0 0 1 1])
            close(h) 
        end
        

%         keyboard
        clear AUC % get rid of this too.
        
        disp([num2str(hh),'/',num2str(numel(dirsK))])
        disp([num2str(i),'/',num2str(numel(files))])
        
        

    end % Loop over files in a given directory
    

end % Loop over directories in data dir




%% save an output mat file.
save([dataKurDir,'AUC_ROC_stats'],'ROCparams','runParams','best','bestT','tFlat','qFlat','fast','stab','AUCm2','AUCs2','chng2','subsetSpec'); 
