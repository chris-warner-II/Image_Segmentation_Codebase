%
% In loop_main_Kuramoto_HandCookedNetwork, I have a flag called
% saveEigMCdata.  If you set it to 1, then a data file called
% EigMC_ParamSearch_data.mat will appear in the
% Kuramoto/HandCookedNetwork/data directory.  It contains values of run
% parameters {R,C,Win,Wout}.  This script will load in that data file and
% produce a series of plots.  Each plot will be for fixed C.  A grid of
% subplots for different {Win, Wout} combinations.  Each subplot will
% contain a number of plots of the dominant eigenvector in different
% colors.  Each color is for a different R value.

SeparationInFlg = 'location'; % 'location' or 'phase'

colorsJ = colormap('jet');
colorsB = flipud(colormap('bone'));

cmapRWB = rd_plotColorbar('redwhiteblue',256);
cmapWR = rd_plotColorbar('whitered',256);
cmapWR = [[0 0 0]; cmapWR];      % Set zero = black because

dirPre = onCluster;

% Directory to look for data output from Kuramoto main to analyze here
dataDir = [dirPre,'output/Kuramoto/HandCookedNetwork/data/'];
matfile = 'EigMC_ParamSearch_data_wNegs';
load([dataDir,matfile,'.mat'])



% Directory to save output plots into...
imgsDir = [dirPre,'output/Kuramoto/HandCookedNetwork/imgs/MetaCluster/Evec_MetaCluster/'];
if ~exist(imgsDir,'dir') 
    mkdir(imgsDir)
end



% (0). extract parameters from ROCparams structure
C = EigMCdata.C;
N = EigMCdata.N;
Win = EigMCdata.Win;
Wout = EigMCdata.Wout;
Rmax = EigMCdata.Rmax;
Pfar = EigMCdata.Pfar;
%
gndTruth = EigMCdata.gndTruth;
meanCVar = EigMCdata.meanCVar;
stdCVar = EigMCdata.stdCVar;
meanCExt = EigMCdata.meanCExt;
stdCExt = EigMCdata.stdCExt;
%
switch SeparationInFlg

    case 'location'
        meanCSep = EigMCdata.location.meanCSep;
        stdCSep = EigMCdata.location.stdCSep;
        minCSep = EigMCdata.location.minCSep;
        minCSepID = EigMCdata.location.minCSepID;
        meanCDist = EigMCdata.location.meanCDist;
        stdCDist = EigMCdata.location.stdCDist;

    case 'phase'
        meanCSep = EigMCdata.phase.meanCSep;
        stdCSep = EigMCdata.phase.stdCSep;
        minCSep = EigMCdata.phase.minCSep;
        minCSepID = EigMCdata.phase.minCSepID;
        meanCDist = EigMCdata.phase.meanCDist;
        stdCDist = EigMCdata.phase.stdCDist;
        
end

   

uC = unique(C);
uST = unique(Win);
uWK = unique(Wout);
uRM = unique(Rmax);
uPF = unique(Pfar);


% to Label axis of plot
for i = 1:numel(uWK)
    uWKc{i} = num2str(uWK(i));
end
%
for i = 1:numel(uST)
    uSTc{i} = num2str(uST(i));
end


% colorScat(1,:) = [1 0 0];       % red 
% colorScat(2,:) = [0 1 0];       % green
% colorScat(3,:) = [0 0 1];       % blue
% colorScat(4,:) = [0 0 0];       % black
% colorScat(5,:) = [0 1 1];       % cyan
% colorScat(6,:) = [1 0 1];       % magenta
% colorScat(7,:) = [1 1 0];       % yellow
% colorScat(8,:) = [0.6 0.6 0.6]; % grey


colorScat = colormap(jet);
colorScat = downsample(colorScat, round(size(colorScat,1)./numel(uWK)) );


colorScat2 = {'red','green','blue','black','cyan','magenta','yellow','black'};

mrkrScat = 'oxs^+';


ltGray = [0.4 0.4 0.4]; % light gray color


%% Plot eigenvectors for different parameter combinations (Win,Wout,R)
if(1) 
    
    for BB=1:numel(uC)
    
        h=figure; hold on
    
        for i = 1:numel(uST)
            for j = 1:numel(uRM)   
                
                subplot( numel(uST), numel(uRM), j + numel(uRM)*(i-1) ), hold on 
                
%                 figure, hold on

                %minsofar = 0;
                
                for k = 1:numel(uWK)
                
                    ind = find( Win == uST(i) & Wout == uWK(k) & Rmax == uRM(j) & C == uC(BB));
                    
                    if k==1
                        template = EigMCdata.Vdom(:,ind);
                    end
                    
                    x = EigMCdata.Vdom(:,ind);
                    
                    
                    if (sum(template.*x) > sum(template.*(-x)) )
                        plot(x + 0.0*(k-1), 'Color', colorScat(k,:),'LineWidth',2)
                    else
                        plot(-x + 0.0*(k-1), 'Color', colorScat(k,:),'LineWidth',2)
                    end
                    
                    %minsofar = min(minsofar, min(x));
                    
                    % [EigMCdata.Rmax(ind), EigMCdata.AUC(ind)]
                    
%                     keyboard
                
                end
                
                
                 
                if(j==1)
                    ylabel(['Win=',num2str(uST(i))],'FontSize',16,'FontWeight','Bold')
                    %xlabel('max AUC')
                    %ylabel('% time @ max')
                end
                
                if(i==1)
                    title(['R=',num2str(uRM(j))],'FontSize',16,'FontWeight','Bold')
                end
                
                set(gca,'XTick',[],'YTick',[])
                axis tight
                
                
            end
        end
        
        
        % label axes on one plot
%         subplot(numel(uST), numel(uWK), numel(uST).*numel(uWK))
%         ylabel('Wout','FontSize',12,'FontWeight','Bold')
%         xlabel('Win','FontSize',12,'FontWeight','Bold')
%         set(gca,'YTick',[1:numel(uWK)] ,'XTick',[1:numel(uST)], 'YTickLabel',uWKc, 'XTickLabel',uSTc)
%         rotateXLabels(gca(),90)
        %
        
        
        TitleVec = ['Wout ', ...
                '\color{',colorScat2{1},'}(=',num2str(uWK(1)),')','\color{',colorScat2{2},'}(=',num2str(uWK(2)),')',...
                '\color{',colorScat2{3},'}(=',num2str(uWK(3)),')','\color{',colorScat2{4},'}(=',num2str(uWK(4)),')',...
                '\color{',colorScat2{5},'}(=',num2str(uWK(5)),')','\color{',colorScat2{6},'}(=',num2str(uWK(6)),')',...
                '\color{',colorScat2{7},'}(=',num2str(uWK(7)),')','\color{',colorScat2{8},'}(=',num2str(uWK(8)),' gray)'];
            
            
        TitleVec = ['Wout ', ...
                '\color{blue}(=',num2str(uWK(1)),')','\color{red}(=',num2str(uWK(end)),')'];    
        
        % title containing other information that is fixed.
        annotation('textbox', [0 0.9 1 0.1],'String', ...
                    ['Eigenvectors for C = ',num2str(uC(BB)),' : with  ',TitleVec],...
                    'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')

        xlabel('Oscillators N')
        ylabel('a.u')
                        
        %
        % Save image
        saveGoodImg(h,[imgsDir,'Evector_vs_Rmax_Win_Wout_N',num2str(N(1)),'_C_',num2str(BB)],[0 0 1 1])
        close(h)
        
%         keyboard
        
        
    end
    
end






%% Look at results of MetaCluster Analysis (Separation measure) using Eigenvectors.
%  Analogous analysis for Kuramoto done in metaClusterVisualize & parseSeparationStruct.

data2plot = meanCSep;
%           this can be mean, median, max, etc.
datamax = max(data2plot);
datamin = min(data2plot);
datalarger = max( datamax, abs(datamin) );




% HOW TO CENTER RED WHITE BLUE COLORMAP WHITE AT ZERO?
mid = size(cmapRWB,1)/2;
% redtop = mid + round((datamax/datalarger)*size(cmapRWB,1)/2);
% blubot = mid + round((datamin/datalarger)*size(cmapRWB,1)/2);





h=figure; hold on

for i = 1:numel(uC)
    for j = 1:numel(uRM)   

        subplot( numel(uC), numel(uRM), j + numel(uRM)*(i-1) ), hold on 

        ind = find( C == uC(i) & Rmax == uRM(j) );

        for I = 1:numel(ind)
            % position of scatter points on grid.
            x = find( uWK == Wout(ind(I)) );
            y = find( uST == Win(ind(I)) );

            % use colorscale to indicate Separation
            % cindx = max( round( ( data2plot(ind(I)) ./ datalarger ) .* size(cmapRWB,1) ), 1);
            cindx = max( mid + round( ( data2plot(ind(I)) ./ datalarger ) .* size(cmapRWB,1)/2 ), 1) ;
            scatter(y,x, 60, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmapRWB(cindx,:), 'LineWidth',1)
            
%             keyboard
            
        end

        axis([0.5 numel(uST)+0.5 0.5 numel(uWK)+0.5])
        pbaspect([numel(uST),numel(uWK),1])
        set(gca,'XTick',[],'YTick',[])

        if(j==1)
            ylabel(['C=',num2str(uC(i))],'FontSize',16,'FontWeight','Bold')
        end

        if(i==1)
            title(['Rmax=',num2str(uRM(j))],'FontSize',16,'FontWeight','Bold')
        end

    end
end


% label axes on one plot
subplot(numel(uC), numel(uRM), numel(uC).*numel(uRM))
ylabel('Wout','FontSize',12,'FontWeight','Bold')
xlabel('Win','FontSize',12,'FontWeight','Bold')
set(gca,'YTick',[1:numel(uWK)] ,'XTick',[1:numel(uST)], 'YTickLabel',uWKc, 'XTickLabel',uSTc)
rotateXLabels(gca(),90)
%

% title containing other information that is fixed.
annotation('textbox', [0 0.9 1 0.1],'String', ...
          ['Cluster Separation Measure Eigenvector  (extremes = [\color{blue}{',num2str(datamin,3),'},\color{red}{',num2str(datamax,3),'}\color{black}{])}'],...
                    'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')

%
% Save image
saveGoodImg(h,[imgsDir,'ClusterSepMean_Eig_',SeparationInFlg],[0 0 0.5 0.8])
close(h)
        






%% Plot Eigenvectors when ratio of Win/Wout = (1):(-1)  then (1):(-2)  then (2):(-1) then (1):(-10) then (10):(-1)

% (1):(-1)  = 100,    50,      20,        10
% (1):(-2)  = 5:-10,  10:-20,  50:-100
% (2):(-1)  = 20:-10, 100:-50
% (1):(-10) = 10:-100, 5:-50
% (10):(-1)  = 100:-10, 50:-5, 10:-1
