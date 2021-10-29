function [] = plot_bench_stats(tag,comp_meths,comp_params)

foci = dir(['./output/ImgSeg/simpleExamples/BSDS/data/*',tag,'*']); % <---- can use a * here.

% comp_meths = 1;  % hold parameter values constant and compare different methods
% comp_params = 0; % hold method constant and compare different parameter settings

colour = 'brgkcmy';

h1=figure; hold on
h2=figure; hold on
h3= hgload('BSDS_bench_jailbreak/benchmarks/isoF.fig'); hold on % load iso-contours for precision-recall

for j = 1:numel(foci)
    
    % copy and paste directory name with output files that you want to feed through this benchmark code.
    focus = foci(j).name

    txtDir = ['./output/ImgSeg/simpleExamples/BSDS/data/',focus,'/eval'];


    if exist(fullfile(txtDir,'eval_cover.txt'),'file'),
    
        cover = dlmread(fullfile(txtDir,'eval_cover.txt'));
        PrecisionBest(j) = cover(2);
        RecallBest(j) = cover(3);
        %
        RI_VOI = dlmread(fullfile(txtDir,'eval_RI_VOI.txt'));
        RIsum(j) = RI_VOI(3);
        VOIsum(j) = RI_VOI(6);
        %
        Info = dlmread(fullfile(txtDir,'eval_Info.txt'));
        InfoMeanSum(j) = Info(2);
        InfoStdSum(j) = Info(3);
        
        % calculate legend entries for plots
        x = strfind(focus,'sP');
        if(comp_meths)
            legos{j} = focus(1:x-2);
            legos{j}(legos{j}=='_')=' ';
        end
        if(comp_params)
            legos{j} = focus(x:end);
        end
        
        % plot info histogram
        figure(h1), plot([0:0.1:1],Info(3:end),colour(j),'LineWidth',2)
        
        % plot precision-recall
%         figure(h3),hold on, scatter(PrecisionBest(j), RecallBest(j), 3, colour(j))

    end
    
    
    
end

% calculate legend entries for plots
x = strfind(focus,'sP');
if(comp_meths)
    foctit = focus(x:end);
    plttit = 'Methods Comparison';
end
if(comp_params)
    foctit = focus(1:x-2);
    plttit = 'Parameters Comparison';
end
foctit(foctit=='_')=' ';

figure(h1), legend(legos')
title(foctit,'Fontsize',20,'Fontweight','Bold')
xlabel('Information value bin','Fontsize',18,'Fontweight','Bold')
ylabel('Percentage of Images falling in Bin','Fontsize',18,'Fontweight','Bold')
grid on

% Plot bar plot of different Information measures
Y = [InfoMeanSum', InfoStdSum', RIsum', VOIsum', PrecisionBest', RecallBest'];
X = 1:numel(foci);
figure(h2), bar(X,Y), grid on
legend('InfoMean', 'InfoStd', 'RandIndex', 'VariationOfInformation', 'Best Precision', 'Best Recall')
title(foctit,'Fontsize',20,'Fontweight','Bold')
xlabel(legos)
ylabel('Value of Quality Metric','Fontsize',18,'Fontweight','Bold')

% Plot Precision-Recall curves (or points)
PR_legos = cell(1,numel(legos)+2);
PR_legos{1} = 'NCut Full';
PR_legos{2} = 'IsoContours';
figure(h3), grid on
for i = 1:numel(legos)
    PR_legos{i+2} = legos{i};
    scatter(PrecisionBest(i), RecallBest(i), 100, colour(i)) % 
end
legend(PR_legos)
title(foctit,'Fontsize',20,'Fontweight','Bold')
xlabel('Precision','Fontsize',18,'Fontweight','Bold')
ylabel('Recall','Fontsize',18,'Fontweight','Bold')



% Save Plots
saveas(h1,['./output/ImgSeg/simpleExamples/BSDS/',foctit,' InfoHist'],'tif'); %name is a string 
                            
saveas(h2,['./output/ImgSeg/simpleExamples/BSDS/',foctit,' Bars'],'tif'); %name is a string 

saveas(h3,['./output/ImgSeg/simpleExamples/BSDS/',foctit,' PrecionRecall'],'tif'); %name is a string 