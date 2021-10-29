% This script will go through the Clustering_Results_allPatches mat files and will
% plot rate distortion curves for combinations of {netMethod, segMethod} and 
% parameter values {rM,sP}.  Different types of plots:
%
% (1). Plot 4 different figures for 4 sP values. With 4x11 subplots for the
% different netMethods and different rM values.  Plot the 5 different rate
% distortion curves corresponding to the different segMethods {im, kur,
% ev1, ev2, ev3}
%
% (2). Make 4 subplots (one for each netMethod) with 5 different rate
% distortion curves in each (for different segMethods) but only take the
% parameter value that performed best {rM*,sP*} for each.


WeightedClusters = 1; % flag to use the D measure that summed up d'_ij values and weighted them by the size of clusters


[dirPre,sizeGoodIm] = onCluster;


dirSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/imgs/RateDistPlots/'];

if ~exist('dirSave','dir')
    mkdir(dirSave)
end



netMethods = {'GLnrm','AAnrm','Mod_SKHAdj','Mod_N&G'};
netMethC = {'GL','AA','Mod SKH','Mod N&G'};
rM = [1:10,inf];
sP = [0.1,0.2,0.3,0.4];
sPc = {'0p1','0p2','0p3','0p4'};

ymax = 2.2;



% This is for GLnrm (unweighted or equally weighted clusters regardless of their size)
runningBestTrapz1 = zeros(1,5); % size of segMethods {im, kur, ev1, ev2, ev3}
runningBestRm1 = zeros(1,5); 
runningBestSp1 = zeros(1,5); 
runningBestRDcurve1 = zeros(5,6); % numel(segMethods) x numel(RateDistSig)-1
%
% This is for AAnrm (unweighted clusters)
runningBestTrapz2 = zeros(1,5); % size of segMethods {im, kur, ev1, ev2, ev3}
runningBestRm2 = zeros(1,5); 
runningBestSp2 = zeros(1,5); 
runningBestRDcurve2 = zeros(5,6); % numel(segMethods) x numel(RateDistSig)-1
%
% This is for Mod_SKHAdj (unweighted clusters)
runningBestTrapz3 = zeros(1,5); % size of segMethods {im, kur, ev1, ev2, ev3}
runningBestRm3 = zeros(1,5); 
runningBestSp3 = zeros(1,5); 
runningBestRDcurve3 = zeros(5,6); % numel(segMethods) x numel(RateDistSig)-1
%
% This is for Mod_N&G (unweighted clusters)
runningBestTrapz4 = zeros(1,5); % size of segMethods {im, kur, ev1, ev2, ev3}
runningBestRm4 = zeros(1,5); 
runningBestSp4 = zeros(1,5); 
runningBestRDcurve4 = zeros(5,6); % numel(segMethods) x numel(RateDistSig)-1

for K = 1:numel(sP)

    H=figure; 

    for j = 1:numel(rM)

        for i = 1:numel(netMethods)

            try

                load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/RateDistPlots/',netMethods{i},'/RDdata_',netMethods{i},'_sP',sPc{K},'_rM',num2str(rM(j)),'_sDInf_sW0_Ks300_Ts1.mat']);

                if(WeightedClusters)
                    x = meanFins(:,2:end);
                else
                    x = meanFinsWt(:,2:end);
                end


                if(K==2 & i==3 & j==1)
                    xpcrd % call a variable that doesnt exist to fall into the catch. (this file is weird)
                end


                subplot(numel(netMethods),numel(rM),j+(i-1)*numel(rM)),
                plot( RateDistSig(2:end),  x' ,'LineWidth',3), hold on
                area(RateDistSig(2:end),x(2,:)','FaceColor',[0.6 0.6 0.6]), alpha(0.5)
                xlabel('\sigma_ ','FontSize',14,'FontWeight','Bold')
                ylabel('d''','FontSize',14,'FontWeight','Bold')
                axis([RateDistSig(2) RateDistSig(end) 0 ymax])
                axis square
                set(gca,'FontSize',12,'FontWeight','Bold','XTick',[RateDistSig(2) RateDistSig(end)])
                %
                if(i==1)
                    title([['R=',num2str(rM(j))]],'FontSize',16,'FontWeight','Bold')
                end
                %clr
                if(j==1)
                    ylabel({netMethC{i},'d'''},'FontSize',16,'FontWeight','Bold')
                end
                %
                if( i==numel(rM) & j==numel(netMethods) )
                    legend(segMethods)
                end

                % compute max d' value and area under curve.

                currentTrapz = trapz(meanFins(:,2:end)');
                
                switch i
                    
                    case 1

                        for t = 1:numel(segMethods)
                            if(currentTrapz(t) > runningBestTrapz1(t))
                                runningBestRm1(t) = rM(j); 
                                runningBestSp1(t) = sP(K);
                                runningBestRDcurve1(t,:) = x(t,:);
                            end
                        end

                        runningBestTrapz1 = max( [runningBestTrapz1; currentTrapz] );
                        
                    case 2
                        
                        for t = 1:numel(segMethods)
                            if(currentTrapz(t) > runningBestTrapz2(t))
                                runningBestRm2(t) = rM(j); 
                                runningBestSp2(t) = sP(K);
                                runningBestRDcurve2(t,:) = x(t,:);
                            end
                        end

                        runningBestTrapz2 = max( [runningBestTrapz2; currentTrapz] );
                        
                        
                    case 3
                        
                        for t = 1:numel(segMethods)
                            if(currentTrapz(t) > runningBestTrapz3(t))
                                runningBestRm3(t) = rM(j); 
                                runningBestSp3(t) = sP(K);
                                runningBestRDcurve3(t,:) = x(t,:);
                            end
                        end

                        runningBestTrapz3 = max( [runningBestTrapz3; currentTrapz] );
                        
                    case 4
                        
                        for t = 1:numel(segMethods)
                            if(currentTrapz(t) > runningBestTrapz4(t))
                                runningBestRm4(t) = rM(j); 
                                runningBestSp4(t) = sP(K);
                                runningBestRDcurve4(t,:) = x(t,:);
                            end
                        end

                        runningBestTrapz4 = max( [runningBestTrapz4; currentTrapz] );
                        
                end

                % % Another thing I could do is look at max value or RD curve instead of area underneath it.
                % max(meanFins(:,2:end)')

            catch
                
                % Leave blank. Leave Subplot blank if data does not exist.
                % Just put label, title or legend.
                subplot(numel(netMethods),numel(rM),j+(i-1)*numel(rM)),
                axis off
                if(i==1)
                    title([['R=',num2str(rM(j))]],'FontSize',16,'FontWeight','Bold')
                end
                %
                if(j==1)
                    ylabel({netMethC{i},'d'''},'FontSize',16,'FontWeight','Bold')
                end
                %
                if( i==numel(rM) & j==numel(netMethods) )
                    legend(segMethods)
                end

            end

        end

    end

    if(WeightedClusters)
        titlestr = ['\sigmaP = ',num2str(sP(K)),' (d'' weighted by cluster size)'];
        filenamestr = [dirSave,'RDplots_sP',sPc{K},'_wtd'];
    else
        titlestr = ['\sigmaP = ',num2str(sP(K)),' (d'' unweighted)'];
        filenamestr = [dirSave,'RDplots_sP',sPc{K}];
    end
    
    annotation('textbox', [0 0.9 1 0.1],'String', titlestr, ...
         'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')


    % save each different sP figure.
    saveGoodImg(H,filenamestr,sizeGoodIm)
    close(H)


end








% Now, Plot RD-Curves for each segMethod for best Performing Parameter {rM,sP} settings
for i=1:numel(segMethods)
    bestLegend1{i} = [segMethods{i},' (rM',num2str(runningBestRm1(i)),' , \sigmaP',num2str(runningBestSp1(i)),')'];
    bestLegend2{i} = [segMethods{i},' (rM',num2str(runningBestRm2(i)),' , \sigmaP',num2str(runningBestSp2(i)),')'];
    bestLegend3{i} = [segMethods{i},' (rM',num2str(runningBestRm3(i)),' , \sigmaP',num2str(runningBestSp3(i)),')'];
    bestLegend4{i} = [segMethods{i},' (rM',num2str(runningBestRm4(i)),' , \sigmaP',num2str(runningBestSp4(i)),')'];
end



H=figure;

subplot(2,2,1)
plot( RateDistSig(2:end),  runningBestRDcurve1' ,'LineWidth',3), hold on
area(RateDistSig(2:end),runningBestRDcurve1(2,:)','FaceColor',[0.6 0.6 0.6]), alpha(0.5)
xlabel('\sigma_ ','FontSize',14,'FontWeight','Bold')
ylabel('d''','FontSize',14,'FontWeight','Bold')
axis([RateDistSig(2) RateDistSig(end) 0 ymax])
set(gca,'FontSize',12,'FontWeight','Bold')
title([netMethC{1}],'FontSize',16,'FontWeight','Bold')
legend(bestLegend1)

subplot(2,2,2)
plot( RateDistSig(2:end),  runningBestRDcurve2' ,'LineWidth',3), hold on
area(RateDistSig(2:end),runningBestRDcurve2(2,:)','FaceColor',[0.6 0.6 0.6]), alpha(0.5)
xlabel('\sigma_ ','FontSize',14,'FontWeight','Bold')
ylabel('d''','FontSize',14,'FontWeight','Bold')
axis([RateDistSig(2) RateDistSig(end) 0 ymax])
set(gca,'FontSize',12,'FontWeight','Bold')
title([netMethC{2}],'FontSize',16,'FontWeight','Bold')
legend(bestLegend2)

subplot(2,2,3)
plot( RateDistSig(2:end),  runningBestRDcurve3' ,'LineWidth',3), hold on
area(RateDistSig(2:end),runningBestRDcurve3(2,:)','FaceColor',[0.6 0.6 0.6]), alpha(0.5)
xlabel('\sigma_ ','FontSize',14,'FontWeight','Bold')
ylabel('d''','FontSize',14,'FontWeight','Bold')
axis([RateDistSig(2) RateDistSig(end) 0 ymax])
set(gca,'FontSize',12,'FontWeight','Bold')
title([netMethC{3}],'FontSize',16,'FontWeight','Bold')
legend(bestLegend3)

subplot(2,2,4)
plot( RateDistSig(2:end),  runningBestRDcurve4' ,'LineWidth',3), hold on
area(RateDistSig(2:end),runningBestRDcurve4(2,:)','FaceColor',[0.6 0.6 0.6]), alpha(0.5)
xlabel('\sigma_ ','FontSize',14,'FontWeight','Bold')
ylabel('d''','FontSize',14,'FontWeight','Bold')
axis([RateDistSig(2) RateDistSig(end) 0 ymax])
set(gca,'FontSize',12,'FontWeight','Bold')
title([netMethC{4}],'FontSize',16,'FontWeight','Bold')
legend(bestLegend4)

if(WeightedClusters)
    titlestr = ['Segmentation Performance with Optimized Parameters (d'' weighted by cluster size)'];
    filenamestr = [dirSave,'RDplots_optimal_wtd'];
else
    titlestr = ['Segmentation Performance with Optimized Parameters (d'' unweighted)'];
    filenamestr = [dirSave,'RDplots_optimal'];
end

annotation('textbox', [0 0.9 1 0.1],'String', titlestr, ...
         'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')


% save each different sP figure.
saveGoodImg(H,filenamestr,sizeGoodIm)
close(H)
