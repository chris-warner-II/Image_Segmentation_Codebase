% This script / function will loop thru the Ground Truth files for all
% image patches and determine how similar each ground truth is to each
% other ground truth by computing Precision, Recall & F-measure of one
% ground truther's boundary using each other ground truther's boundary as
% ground truth.


% maxDist = 1e-3; % 0.0075; % this is the default value used in my code in bench_bsds500

[dirPre,sizeGoodIm] = onCluster;
addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])

gTdir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];

gTfiles = dir([gTdir,'*.mat']);

outImgDir = [gTdir,'AgreementPics/'];
if ~exist(outImgDir,'dir')
    mkdir(outImgDir)
end



for F = 1:numel(gTfiles)

    
    disp(['Ground Truth file # ',num2str(F),' / ',num2str(numel(gTfiles))])
    gTfiles(F).name

    load([gTdir,gTfiles(F).name])
    
    % Break out and dont rerun if the variables I am about to create already exist.
    if exist('N_edge','var')
        clear N_edge
        continue
    end
    
    
    % preallocate memory
    N_edge = zeros(1,numel(groundTruth));    % number of edge pixels in each gT
    %
    F_cross = zeros(numel(groundTruth));     % pairwise measures of (dis-)similarity between gTs
    F_cross1a = zeros(numel(groundTruth));
    F_cross2a = zeros(numel(groundTruth));
    F_cross1b = zeros(numel(groundTruth));
    F_cross2b = zeros(numel(groundTruth));
    %
    F_1vrest = zeros(2,numel(groundTruth));  % for each gT, a measure of (dis-)similarity between itself and all other gTs
    F_1vrest1 = zeros(2,numel(groundTruth));
    F_1vrest2 = zeros(2,numel(groundTruth));
    %
    F_tot = [];
    F_tot1 = [];
    F_tot2 = [];
    
    
    

    for i = 1:numel(groundTruth) % loop thru ground truths
        N_edge(i) = numel(find(groundTruth{i}.Boundaries)); % # of pixels that are edges
    end
    
    
    

    for i = 1:numel(groundTruth) % loop thru ground truths

        for j = 1:numel(groundTruth) % loop thru other ground truths
            
            % CW: 1st way: computes direct overlap between pixels of gT boundaries.
            match0 = groundTruth{i}.Boundaries&groundTruth{j}.Boundaries;
            F_cross(i,j) = numel(find(match0)) ./ sqrt(N_edge(i)*N_edge(j));

            % A 2nd way to compute pixel "overlap" (allowing proximity with small maxDist).
            maxDist = 1e-3;
            [match1a,match1b] = correspondPixels(double(groundTruth{i}.Boundaries), double(groundTruth{j}.Boundaries), maxDist);
            F_cross1a(i,j) = numel(find(match1a)) ./ sqrt(N_edge(i)*N_edge(j));
            F_cross1b(i,j) = numel(find(match1b)) ./ sqrt(N_edge(i)*N_edge(j));
            
            
            % A 3nd way to compute pixel "overlap" (allowing proximity with larger maxDist).
            maxDist = 0.0075;
            [match2a,match2b] = correspondPixels(double(groundTruth{i}.Boundaries), double(groundTruth{j}.Boundaries), maxDist);
            F_cross2a(i,j) = numel(find(match2a)) ./ sqrt(N_edge(i)*N_edge(j));
            F_cross2b(i,j) = numel(find(match2b)) ./ sqrt(N_edge(i)*N_edge(j));

        end
        
        % Compute Mean values for each GT relative to others and one mass mean
        % to tell how similar all ground truths are to one another.
        F_1vrest(1,i) = mean(F_cross(i,[1:i-1,i+1:end]));
        F_1vrest(2,i) = std(F_cross(i,[1:i-1,i+1:end]));
        %
        F_1vrest1(1,i) = mean(F_cross1a(i,[1:i-1,i+1:end]));
        F_1vrest1(2,i) = std(F_cross1a(i,[1:i-1,i+1:end]));
        %
        F_1vrest2(1,i) = mean(F_cross2a(i,[1:i-1,i+1:end]));
        F_1vrest2(2,i) = std(F_cross2a(i,[1:i-1,i+1:end]));
        
        F_tot = [ F_tot, F_cross(i,[1:i-1,i+1:end]) ];
        F_tot1 = [ F_tot1, F_cross1a(i,[1:i-1,i+1:end]) ];
        F_tot2 = [ F_tot2, F_cross2a(i,[1:i-1,i+1:end]) ];
         
    end
    
    
    % Save mean & std of (dis-)similarity of all pairs of ground truths
    F_tot_stats = [mean(F_tot), std(F_tot)];
    F_tot_stats1 = [mean(F_tot1), std(F_tot1)];
    F_tot_stats2 = [mean(F_tot2), std(F_tot2)];

    
    
    
    
    
    
    
    
    
    
    
    
    if(0)
        % check that ...
        F_cross1a - F_cross1b
        F_cross2a - F_cross2b
        F_cross1a - F_cross
        F_cross2a - F_cross

        % check if each F_cross matrix is symmetric (same using one gT as truth vs using other)
        F_cross - F_cross'
        F_cross1a - F_cross1a'
        F_cross2a - F_cross2a' 
        
        if ~any(any(F_cross1a - F_cross1b))
            clear F_cross1b 
        end

        if ~any(any(F_cross2a - F_cross2b))
            clear F_cross2b 
        end
    
    end
    

    
    
    
    
    
    %% Plot these different F_cross matrices along with gT images
    H0=figure;
    subplot(231), imagesc(F_cross),   caxis([0 1]), axis square, title('Direct Pixel Overlap'), freezeColors
    subplot(232), imagesc(F_cross1a), caxis([0 1]), axis square, title('Correspond Pixels (maxDist=0.0001)'), freezeColors
    subplot(233), imagesc(F_cross2a), caxis([0 1]), axis square, title('Correspond Pixels (maxDist=0.0075)'), freezeColors
    %
    for i = 1:numel(groundTruth) % loop thru ground truths
        subplot(2,numel(groundTruth),numel(groundTruth)+i),
        imagesc(groundTruth{i}.Boundaries), axis square
        set(gca,'XTick',[],'YTick',[])
        title(['gT',num2str(i)])
        colormap('bone'), freezeColors
    end
    
    saveGoodImg(H0,[outImgDir,'gT_agreement0_',gTfiles(F).name(1:end-4)],sizeGoodIm)
    close(H0)
    
    
    
    %% Plot Pairwise Pixel Correspondence/Overlap along with numbers quantifying this.
    H1 = figure;
    ha = tight_subplot(numel(groundTruth),numel(groundTruth),[.01 .01],[.01 .05],[.01 .01]);
    
    H2 = figure;
    hb = tight_subplot(numel(groundTruth),numel(groundTruth),[.01 .01],[.01 .05],[.01 .01]);
    
    
    for i = 1:numel(groundTruth) % loop thru ground truths

        for j = 1:numel(groundTruth) % loop thru other ground truths
            
            
            % Plot overlap of each pair of gT's (match0)
            figure(H1)
            axes(ha( (i-1)*numel(groundTruth)+j ) )
            match0 = groundTruth{i}.Boundaries&groundTruth{j}.Boundaries;
            imagesc(match0), axis square, colormap('bone'), set(gca,'XTick',[],'YTick',[])
            text(size(groundTruth{1}.Boundaries),size(groundTruth{1}.Boundaries),['\color{red}{',num2str(F_cross(i,j),2),'}'],'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16,'FontWeight','Bold')
            %
            if(i==j)
                text(1,1,{['\color{red}{gT#',num2str(i),'}'], ['\color{red}{N=',num2str(N_edge(i)),'}']},'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            end
            %
            if(j==1)
                ylabel(['[',num2str(F_1vrest(1,i),2), ' +/- ', num2str(F_1vrest(2,i),2),']'],'FontSize',16,'FontWeight','Bold')
            end
            
            % Plot overlap of each pair of gT's (match0)
            figure(H2)
            axes(hb( (i-1)*numel(groundTruth)+j ) )
            maxDist = 0.0075;
            [match2a,match2b] = correspondPixels(double(groundTruth{i}.Boundaries), double(groundTruth{j}.Boundaries), maxDist);
            imagesc(match2a>0), axis square, colormap('bone'), set(gca,'XTick',[],'YTick',[])
            text(size(groundTruth{1}.Boundaries),size(groundTruth{1}.Boundaries),['\color{red}{',num2str(F_cross2a(i,j),2),'}'],'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16,'FontWeight','Bold')
            %
            if(i==j)
                text(1,1,{['\color{red}{gT#',num2str(i),'}'], ['\color{red}{N=',num2str(N_edge(i)),'}']},'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            end
            %
            if(j==1)
                ylabel(['[',num2str(F_1vrest2(1,i),2), ' +/- ', num2str(F_1vrest2(2,i),2),']'],'FontSize',16,'FontWeight','Bold')
            end
            
            
        end
        
    end

    figure(H1), annotation('textbox', [0 0.9 1 0.1],'String',['Pairwise GT Direct Boundary Overlap : ',gTfiles(F).name(1:end-4),' : [',num2str(mean(F_tot),2), ' +/- ', num2str(std(F_tot),2),']'], 'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')
    figure(H2), annotation('textbox', [0 0.9 1 0.1],'String',['Pairwise GT Boundary Proximity (d = 0.0075) : ',gTfiles(F).name(1:end-4),' : [',num2str(mean(F_tot2),2), ' +/- ', num2str(std(F_tot2),2),']'], 'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')
    %
    saveGoodImg(H1,[outImgDir,'gT_agreement1_',gTfiles(F).name(1:end-4)],sizeGoodIm)
    close(H1)
    %
    saveGoodImg(H2,[outImgDir,'gT_agreement2_',gTfiles(F).name(1:end-4)],sizeGoodIm)
    close(H2)
    
    

    
    %% Save gT mat file back in place of the old one.
    save([gTdir,gTfiles(F).name],'groundTruth','pach','N_edge','F_cross','F_cross1a','F_cross2a',...
        'F_1vrest','F_1vrest1','F_1vrest2','F_tot_stats','F_tot_stats1','F_tot_stats2')
    
    clear N_edge % doing this to get around break out check above.
    
    
end