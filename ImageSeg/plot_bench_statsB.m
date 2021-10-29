function [] = plot_bench_statsB(tag,comp_meths,comp_params)

foci = dir(['./output/ImgSeg/simpleExamples/gradBox/data/*',tag,'*']); % <---- can use a * here.

% comp_meths = 1;  % hold parameter values constant and compare different methods
% comp_params = 0; % hold method constant and compare different parameter settings


%% Organize the directories by method and parameter values
focus = foci(4).name; % starts at 4 because dir finds '.' '..' & '.Dsstore'
x = strfind(focus,'sP');
km = 1; kp = 1;
method{km} = focus(1:x-2);
params{kp} = focus(x:end);

for j = 5:numel(foci)
    
    focus = foci(j).name;
    x = strfind(focus,'sP');
    temp_method = focus(1:x-2);
    temp_params = focus(x:end);
    
    % if this file is for a new method, record it.
    tagm=1;
    for i = 1:numel(method)
        if ~isempty(strmatch(method{i},temp_method))
            tagm = 0;
        end
    end
    %
    if(tagm)
        km=km+1;
        method{km} = temp_method;
    end
    
    
    % if this file is for a new method, record it.
    tagp=1;
    for i = 1:numel(params)
        if ~isempty(strmatch(params{i},temp_params))
            tagp = 0;
        end
    end
    %
    if(tagp)
        kp=kp+1;
        params{kp} = temp_params;
    end
    
end


%% Extract quality metric data from txt files and store it in matrices
for j = 1:numel(method)
    for k = 1:numel(params)
    
        focus = [method{j},' ',params{k}];

        txtDir = ['./output/ImgSeg/simpleExamples/BSDS/data/',focus,'/eval'];

        if exist(fullfile(txtDir,'eval_cover.txt'),'file'),
            cover = dlmread(fullfile(txtDir,'eval_cover.txt'));
            PrecisionBest(j,k) = cover(2);
            RecallBest(j,k) = cover(3);
            %
            RI_VOI = dlmread(fullfile(txtDir,'eval_RI_VOI.txt'));
            RIsum(j,k) = RI_VOI(3);
            VOIsum(j,k) = RI_VOI(6);
            %
            Info = dlmread(fullfile(txtDir,'eval_Info.txt'));
            InfoMeanSum(j,k) = Info(2);
            InfoStdSum(j,k) = Info(3);
            InfoHist(j,k,:) = Info(3:end); % Use This Later:  reshape(InfoHist(1,1,:),1,11)
        end
    
    end
    
end


%% Make various plots (means across parameters for each method) (best results & parameter settings for each method)

% Bar plot (Mean & Std) of quality metric (averaged across parameter settings) for each method
if(1) 

    
    
    
    
end



keyboard













%% Make Some Plots of different Quality Metrics for Different Methods and Parameter Values
colour = 'brgkcmy';

h1=figure; hold on % plot histogram of our Information Measure
h2=figure; hold on % bar plot of all quality metrics 
h3= hgload('BSDS_bench_jailbreak/benchmarks/isoF.fig'); hold on % load iso-contours for precision-recall


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