% This script is me trying to find ideal distribution of units in linear
% space given the number and size of clusters.  The idea is when Din is
% zero and Dout is maximal.
%
%
function [Dout_BestC,Dout_BestL,Dout_DiffC,Dout_DiffL] = compute_Dout_ideal(fileSt,fileNd,numIter,plotPic,outFileCheck)

    % Load an image patch file (which includes ground truth)
    [dirPre, sizeGoodIm] = onCluster;
    patchesDir = [dirPre,'images/BSDS_patch/51x51_ds1/'];
    
    % Directory to save patches and result from this optimization into
    saveDir = [dirPre,'images/BSDS_patch/51x51_ds1_wDopt/'];
    if ~exist(saveDir,'dir')
        mkdir(saveDir)
    end
    

    files = dir([patchesDir,'*.mat']);
    
    for F = fileSt:fileNd
    
    
        % Construct simple example with equal sized clusters for sanity check.
        if(0)

            ximg=50;
            yimg=50;
            gT{1} = zeros(ximg,yimg);

            factor(ximg*yimg)

            numC = 2; 
            s = repmat( (ximg*yimg./numC), 1, numC);
            for i = 1:numC
                gT{1}((i-1)*s(i) + 1:i*s(i))=i;
            end
            im=gT{1};

            % What would optimal Dout for circular case be?
            x = linspace(0,2*pi,numC+1);
            x = x(1:end-1);
            numpairs = sum(s).^2 - sum(s.^2); % note double counting here and in counting up pairwise distances.

            d_optimal_C = wtd_dist_C(x,s);
            d_optimal_C = -d_optimal_C./numpairs;

        end
        
        % Check if the image patch file already exists in the _wDopt directory
        if(outFileCheck & exist([saveDir,files(F).name],'file'))
            disp(['File with Dout Optimal Computed Already Exists.  Moving on...'])
            [saveDir,files(F).name]
            continue
        end
        
        
        
        
        load([patchesDir,files(F).name])

        Dout_BestC = zeros(1,numel(gT));
        Dout_DiffC = zeros(1,numel(gT));
        Dout_BestL = zeros(1,numel(gT));
        Dout_DiffL = zeros(1,numel(gT));


        if(plotPic)
            HE=figure; 
            subplot(6+numIter,1,1:5), imagesc(im),axis off square, colormap('bone'), freezeColors, title('Linear Eigenvector')
            for i = 1:numel(gT)
                subplot(6+numIter, numel(gT), 5*numel(gT)+i), imagesc(gT{i}),colormap('jet'), axis off square, freezeColors
            end
            
            %
            HL=figure; 
            subplot(6+numIter,1,1:5), imagesc(im),axis off square, colormap('bone'), freezeColors, title('Linear Strawman')
            for i = 1:numel(gT)
                subplot(6+numIter, numel(gT), 5*numel(gT)+i), imagesc(gT{i}),colormap('jet'), axis off square, freezeColors
            end
            
            %
            HC=figure; 
            subplot(6+numIter,1,1:5), imagesc(im),axis off square, colormap('bone'), freezeColors, title('Circular Kuramoto')
            for i = 1:numel(gT)
                subplot(6+numIter, numel(gT), 5*numel(gT)+i), imagesc(gT{i}),colormap('jet'), axis off square, freezeColors
            end
            
        end


        for i = 1:numel(gT) % Loop thru each human segmentation.

               disp(['GT#',num2str(i),' / ',num2str(numel(gT))]) 

               x = unique(gT{i}); % number of different clusters in this segmentation.

               for j = 1:numel(x)
                   s(j) = numel(find(gT{i}==x(j))); % size of each cluster
               end
               
               numpairs = sum(s).^2 - sum(s.^2); % number of pairs of units across clusters (normalize D by this)
                                                 % note double counting here and in counting up pairwise distances.

               if(plotPic)
                    figure(HE), subplot(6+numIter, numel(gT), 5*numel(gT)+i),title(['#Clusters=',num2str(numel(s))])
                    figure(HL), subplot(6+numIter, numel(gT), 5*numel(gT)+i),title(['#Clusters=',num2str(numel(s))])
                    figure(HC), subplot(6+numIter, numel(gT), 5*numel(gT)+i),title(['#Clusters=',num2str(numel(s))])
               end


           if(1)

               for k = 1:numIter 

                   disp(['Optimization Iter#',num2str(k)])

                    % Use Matlab's fmincon to solve this constrained optimization problem.
                    x0 = rand(size(s)); % initialize all clusters at zero.
                    % x0 = s; % initialize clusters location based on size.
                    % x0 = 1/s; % initialize at inverse of cluster size and stagger pos and neg?
                    % x0 = [0 -1 1];

                    options = optimset('fmincon');
                    options.MaxFunEvals = 99999;
                    options.MaxIter = 4000;

                    % Note:  This is for the Linear Variable Eigenvector Optimization Stuff
                    [xE,dvalE] = fmincon(@(x)wtd_dist_L(x,s),x0,[],[],[],[],[],[],@(x)evecNorm(x,s),options);
                    if( sum(s.*(xE).^2) - 1 > 1e-3 )
                        disp('Error: Eigenvector not normalized')
                        keyboard
                    end

                    % Note: This is for the Linear Variable Strawman Optimization Stuff
                    [xL,dvalL] = fmincon(@(x)wtd_dist_C(x,s),x0,[],[],[],[],zeros(size(x)),ones(size(x)),[],options);

                    % Note: This is for the Circular Variable Kuramoto Optimization Stuff
                    x0c = 2*pi*x0;
                    [xC,dvalC] = fmincon(@(x)wtd_dist_C(x,s),x0c,[],[],[],[],zeros(size(x)),2*pi*ones(size(x)),[],options);
                    xC = wrapTo2Pi(xC-xC(1)); % make 1st cluster at zero (just for plotting)


        %             % check against spreading largest clusters most and putting smaller one in between.
        %             s
        %             x
        %             -dval % dout is essentially -dval.

                    if(plotPic)
                        figure(HE)
                        subplot(6+numIter, numel(gT), (6+k-1)*numel(gT) + i)
                        scatter(xE,ones(1,numel(xE)),s/5)
                        xlabel(['D=',num2str(-dvalE./numpairs,2)])
                        axis tight
                        set(gca,'YTick',[])
                        %
                        figure(HL)
                        subplot(6+numIter, numel(gT), (6+k-1)*numel(gT) + i)
                        scatter(xL,ones(1,numel(xL)),s/5)
                        xlabel(['D=',num2str(-dvalL./numpairs,2)])
                        xlim([0 1])
                        set(gca,'YTick',[],'XTick',[0 1])
                        
                        %
                        figure(HC)
                        subplot(6+numIter, numel(gT), (6+k-1)*numel(gT) + i)
                        scatter(xC,ones(1,numel(xC)),s/5)
                        xlabel(['D=',num2str(-dvalC./numpairs,4)])
                        xlim([0 2*pi])
                        set(gca,'YTick',[],'XTick',[0 2*pi],'XTickLabel',{'0','2Pi'})
                        
                    end

                    DE(k) = -dvalE;
                    DL(k) = -dvalL;
                    DC(k) = -dvalC;

               end

               DE = DE./numpairs;
               DL = DL./numpairs;
               DC = DC./numpairs;
               %
               Dbest_E = max(DE);
               Dbest_L = max(DL);
               Dbest_C = max(DC);
               %
               Ddiff_E = mean(DE-Dbest_E);
               Ddiff_L = mean(DL-Dbest_L);
               Ddiff_C = mean(DC-Dbest_C);
               

               if(1)
                    Hh=figure; 
                    subplot(311), hist(DE-Dbest_E),title(['Eigenvector Max D = ',num2str(Dbest_L),' : Mean Diff ',num2str(Ddiff_L),' : #C=',num2str(numel(s))]);
                    subplot(312), hist(DL-Dbest_L),title(['Strawman Max D = ',num2str(Dbest_L),' : Mean Diff ',num2str(Ddiff_L)]);
                    subplot(313), hist(DC-Dbest_C),title(['Kuramoto Max D = ',num2str(Dbest_C),' : Mean Diff ',num2str(Ddiff_C)]);
               end

           end

           
           %
           Dout_BestE(i) = Dbest_E;
           Dout_DiffE(i) = Ddiff_E;
           %
           Dout_BestL(i) = Dbest_L;
           Dout_DiffL(i) = Ddiff_L;
           %
           Dout_BestC(i) = Dbest_C;
           Dout_DiffC(i) = Ddiff_C;
           
           
           clusterSize{i} = s;

           if exist('d_optimal_C','var')
                d_optimal_C
           end

           % keyboard
           
           clear s

        end


        % Save it somewhere (in patch file?) (in other Dout_ideal_Opt file?)
        save([saveDir,files(F).name], 'Dout_BestC', 'Dout_DiffC', 'Dout_BestE', 'Dout_DiffE', 'Dout_BestL', 'Dout_DiffL', 'clusterSize', ...
            'numIter', 'options', 'yimg', 'ds_fctr', 'gT', 'gTfull', 'im', 'imFull', 'pach', 'xFimg', 'ximg', 'yFimg')

    
    end
    
    

end


% define the weighted distance function for Linear (Eivenvector) Variable that
% I am trying to maximize (I am actually trying to minimize the negative of it)
function [d] = wtd_dist_L(x,s)

    d = 0;
    for i = 1:numel(s)
        for j = 1:numel(s)
            d = d - abs(x(i) - x(j)).*s(i).*s(j); % note double counting here and in numpairs
        end
    end
    
    
    
end


% define Eigenvector normalization which is a nonlinear constraint.
function [c,ceq] = evecNorm(x,s)

    c = [];                    % nonlinear inequality constraints
    ceq = sum( s.*(x.^2) )  -1; % nonlinear equality constraints (=0)
    
end




% define the weighted distance function for Circular (Kuramoto) Variable that
% I am trying to maximize (I am actually trying to minimize the negative of it)
function [d] = wtd_dist_C(x,s)

    d = 0;
    for i = 1:numel(s)
        for j = 1:numel(s)
            d = d - abs(circ_dist(x(i),x(j))).*s(i).*s(j); % note double counting here and in numpairs
        end
    end
    
end