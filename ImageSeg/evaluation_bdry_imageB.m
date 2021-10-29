function [thresh,cntR,sumR,cntP,sumP] = evaluation_bdry_imageB(inFile, gtFile, prFile, nthresh, maxDist)
% [thresh,cntR,sumR,cntP,sumP] = evaluation_bdry_imageB(inFile,gtFile, prFile, nthresh, maxDist)
%
% Calculate precision/recall curve for an image patch.  Written by
% C.Warner. We are now using correspondPixelsB - which allows us to compute
% match between thresholded pb image and groundtruth for different distance
% tolerances from 0 to maxDist. Precision & Recall values will be
% calculated from these overlaps at different distance tolerances and will
% be saved into *_ev1.txt files that can be parsed and plotted later.
%
% INPUT
%	inFile  : Can be one of the following:
%             - a soft or hard boundary map in image format.
%             - a collection of segmentations in a cell 'segs' stored in a mat file
%             - an ultrametric contour map in 'doubleSize' format, 'ucm2'
%               stored in a mat file with values in [0 1].
%
%	gtFile	: File containing a cell of ground truth boundaries
%   prFile  : Temporary output for this image.
%	nthresh	: Number of points in PR curve.
%   MaxDist : For computing Precision / Recall.
%   thinpb  : option to apply morphological thinning on segmentation
%             boundaries. (NO LONGER AVAILABLE!)
%
% OUTPUT
%	thresh		Vector of threshold values.
%	cntR,sumR	Ratio gives recall.
%	cntP,sumP	Ratio gives precision.
%
% What I am actually now saving in the *_d#_ev2.txt files.  Chris Warner
% I am saving one output txt file for each distance between 0 & maxDist.
% In it, are:
%
%  (1). 'thresh'          - threshold used on pb_png image for that row/column of results. 
%  (2). 'd' E [0,maxDist] - distance allowed for correspondPixels function (0 = direct pixel overlap)
%  (3). mean Recall (across gT's at a given d value)
%  (4). std  Recall (across gT's at a given d value)
%  (5). max  Recall (for single best gT at a given d value)    - Rather it is R that leads to max F-measure
%  (6). mean Precision (across gT's at a given d value)
%  (7). std  Precision (across gT's at a given d value)
%  (8). max  Precision (for single best gT at a given d value) - Rather it is P that leads to max F-measure
%  (9). mean F-measure (across gT's at a given d value)
% (10). std  F-measure (across gT's at a given d value)
% (11). max  F-measure (for single best gT at a given d value)
% (12). Number of groundTruths (gT's) for this image patch.
% (13). Fmax_whichGT      - which gT matched up best to make max F-measure (#11)
%

plot_th_vs_gt_diagnose = 0; % flag to plot gt vs thresholded pb file and overlap.


[p,n,e]=fileparts(inFile);
if strcmp(e,'.mat'),
    load(inFile);
end

if exist('ucm2', 'var'),
    pb = double(ucm2(3:2:end, 3:2:end));
    clear ucm2;
elseif ~exist('segs', 'var')
    pb = double(imread(inFile))/255;
end


load(gtFile);
if isempty(groundTruth),
    error(' bad gtFile !');
end

if ~exist('segs', 'var')
    thresh = linspace(1/(nthresh+1),1-1/(nthresh+1),nthresh)';
else
    if nthresh ~= numel(segs)
        warning('Setting nthresh to number of segmentations');
        nthresh = numel(segs);
    end
    thresh = 1:nthresh; thresh=thresh';
end

% zero all counts
cntR = zeros(numel(thresh),4);
sumR = zeros(size(thresh));
cntP = zeros(size(thresh));
sumP = zeros(size(thresh));

cntR0 = zeros(size(thresh)); % added by CW
cntP0 = zeros(size(thresh)); % added by CW

sumR0 = zeros(size(thresh)); % added by CW

N = zeros(numel(thresh),numel(groundTruth)); % added by CW
% N1 = zeros(numel(thresh),numel(groundTruth)); % added by CW
% N2 = zeros(numel(thresh),numel(groundTruth)); % added by CW
N3 = zeros(numel(thresh),numel(groundTruth),4); % added by CW (4 = numel(cpB_dists)
Dp = zeros(numel(thresh),numel(groundTruth)); % added by CW
Dr = zeros(numel(thresh),numel(groundTruth)); % added by CW

if(plot_th_vs_gt_diagnose)
    H1=figure;
    H2=figure;
    H3=figure;
    H4=figure;
end




for t = 1:nthresh, % Loop thru different threshold values in the probabalistic boundaries (pb) image.

    
    if(0)
       disp(['PB Threshold = ',num2str(t),' / ',num2str(nthresh)]) 
    end
    

    if ~exist('segs', 'var')
        bmap = (pb>=thresh(t));
    else
        bmap = logical(seg2bdry(segs{t},'imageSize'));
    end


    if(plot_th_vs_gt_diagnose)
        % plot thresholded pb_png images
        for d=1:4
            eval(['figure(H',num2str(d),');']) 
            subplot(numel(groundTruth)+1,nthresh+1,t+1)
            imagesc(bmap), axis square, colormap(bone)
            set(gca,'XTick',[],'YTick',[])
            title(['th=',num2str(thresh(t),2)])
            xlabel(['N=',num2str(numel(find(bmap)))])
        end
    end



    % accumulate machine matches, since the machine pixels are
    % allowed to match with any segmentation
%     accP = zeros(size(bmap));
    accP0 = zeros(size(bmap)); % added by CW.
    accGT = zeros(size(bmap)); % added by CW.

    % compare to each seg in turn
    for i = 1:numel(groundTruth),

        % CW: This just computes direct overlap or match or correspondence between gT boundaries and bmap.
        match0 = bmap&groundTruth{i}.Boundaries;
        % accumulate machine matches
        accP0 = accP0 | match0;

        accGT = accGT | groundTruth{i}.Boundaries; % the union of all gT boundaries (CW)

        % Writing my own correspondPixels. (Theirs is a nonsensical black box)        
        [match3,cpB_dists] = correspondPixelsB(bmap, double(groundTruth{i}.Boundaries), maxDist);
        
        
        
        
        % Note: match3(:,:,1) must = match0!
        if any(any(match0 - match3(:,:,1)))
            disp('Something is not right!')
            keyboard
        end
        
        
        
        
        % compute recall
        sumR(t) = sumR(t) + sum(groundTruth{i}.Boundaries(:));

        cntR0(t) = cntR0(t) + sum(match0(:)>0); % added by CW

        N(t,i) = sum(match0(:)>0);              % count pixels in overlap of thresholded pb image at t and groundtruth (numerator of P&R) : True Positives.

        Dp(t,i) = sum(bmap(:));                       % count pixels in thresholded pb image at t (does not vary with gT(i)) (denominator of P) : Items Selected.
        Dr(t,i) = sum(groundTruth{i}.Boundaries(:));  % count pixels in ground truth i (does not vary with pb threshold t)   (denominator of R) : Relevant Items.

        for d = 1:numel(cpB_dists) % Loop over distances allowed in correspondPixels computation.

            if(plot_th_vs_gt_diagnose)

                % plot ground truths
                eval(['figure(H',num2str(d),');']) 
                subplot(numel(groundTruth)+1,nthresh+1,i*(nthresh+1)+1)
                imagesc(groundTruth{i}.Boundaries), axis square, colormap(bone)
                set(gca,'XTick',[],'YTick',[])
                ylabel(['gt#',num2str(i,2)])
                if(i==1)
                    xlabel(['\color{green}P | \color{blue}R | \color{red}F'])
                end

                % plot overlap or correspondence (match0, match1, match2) between gt and thresholded pb_png
                eval(['figure(H',num2str(d),');']) 
                subplot(numel(groundTruth)+1,nthresh+1,i*(nthresh+1)+1+t)
                imagesc(logical(match3(:,:,d))), axis square, colormap(bone)
                set(gca,'XTick',[],'YTick',[])

            end
            
            cntR(t,d) = cntR(t,d) + sum(sum(match3(:,:,d)>0));

            N3(t,i,d) = sum(sum(match3(:,:,d)>0));

        end

        clear match0 match3
        
        

    end % finish looping over groundTruth boundaries



    % save little images for figure 2 of NIPS paper. Save image patch,
    % blurred image patch, phase map, pb, thresholded pb's, gT's, match
    % of one bb-gT pair for different cPd values.
    if(0) % CW

        keyboard

        % (1). plot probabalistic boundary map
        imwrite(pb,'./pb_map.png')
        %
        % (2). plot ground truths
        for g = 1:numel(groundTruth)
            imwrite(1-groundTruth{g}.Boundaries,['./gT',num2str(g),'.png'])
        end
        %
        % (3). plot binary thresholded pb maps
        for t = 1:2:nthresh
            bb = (pb>thresh(t));
            imwrite(1-bb,['./b_pb',num2str(t),'.png'])
        end
        %
        % (4). plot match between gT and bb for different cpD values
        [match3,cpB_dists] = correspondPixelsB(bmap, double(groundTruth{3}.Boundaries), maxDist);
        for d = 1:numel(cpB_dists)
            imwrite(1-match3(:,:,d),['./match_d',num2str(d),'.png'])
        end

    end




    % compute precision     % Note: These arent being computed in a loop so they dont need this (self +) term
    sumP(t) = sum(bmap(:)); % sumP(t) +    
%    cntP(t) = sum(accP(:)); % cntP(t) + 
    cntP0(t) = sum(accP0(:)); % CW. % cntP0(t) + 

    sumR0(t) = sum(accGT(:)); % CW (note: this is the denominator for Recall 
    % that would be paired with Precision as the benchmark code is originally implemented)
    % (This is the number of edge pixels in the union of all groundtruths)

    if(0) % CW
        [sumR(t) cntR(t) cntR0(t) sumP(t) cntP(t) cntP0(t)]
        [cntP./sumP, cntP0./sumP] % Precision.
        [cntR./sumR, cntR0./sumR] % Recall.

        keyboard
    end

end % looping over t (different values of threshold)




% Compute Precision, Recall & F-measure using direct pixel overlap between thresholded pb & groundtruth i.
% (This is used to compute maxGT & meanGT performance metric calculations)
R_overlap = N./Dr;
P_overlap = N./Dp;
F_overlap = 2*P_overlap.*R_overlap./(P_overlap+R_overlap+((P_overlap+R_overlap)==0));
%
[Fmax, Fmax_whichGT] = max(F_overlap,[],2);
Rmax = zeros(size(Fmax));
Pmax = zeros(size(Fmax));
for i = 1:nthresh
    Rmax(i) = R_overlap(i,Fmax_whichGT(i));
    Pmax(i) = P_overlap(i,Fmax_whichGT(i));
end



% Compute Precision, Recall & F-measure using correspondPixelsB match3 between thresholded pb & groundtruth i at different distances (cpB_dists).
% (This is used to compute maxGT & meanGT performance metric calculations)
R_cP3 = N3./repmat(Dr,[1,1,numel(cpB_dists)]);
P_cP3 = N3./repmat(Dp,[1,1,numel(cpB_dists)]);
F_cP3 = 2*P_cP3.*R_cP3./(P_cP3+R_cP3+((P_cP3+R_cP3)==0));
%
[Fmax_cP3, Fmax_whichGT_cP3] = max(F_cP3,[],2);
Fmax_cP3 = squeeze(Fmax_cP3);
Fmax_whichGT_cP3 = squeeze(Fmax_whichGT_cP3);

Rmax_cP3 = zeros(size(Fmax_cP3));
Pmax_cP3 = zeros(size(Fmax_cP3));
for d = 1:numel(cpB_dists)
    for i = 1:nthresh
        Rmax_cP3(i,d) = R_cP3(i,Fmax_whichGT_cP3(i),d);
        Pmax_cP3(i,d) = P_cP3(i,Fmax_whichGT_cP3(i),d);
    end
end




% Record Precision, Recall, F-measure on the various plots
if(plot_th_vs_gt_diagnose)
    
    for d = 1:numel(cpB_dists)
        d
        for i = 1:numel(groundTruth)
            for t = 1:nthresh
                if(F_cP3(t,i,d)>0)
                    eval(['figure(H',num2str(d),');']) 
                    subplot(numel(groundTruth)+1,nthresh+1,i*(nthresh+1)+1+t)
                    xlabel(['\color{green}',num2str(P_cP3(t,i,d),2),' | \color{blue}',num2str(R_cP3(t,i,d),2),' | \color{red}',num2str(F_cP3(t,i,d),2)])
                end
            end
        end
        %
        annotation('textbox', [0 0.9 1 0.1],'String',['Pixel Correspondence between bmap & gT allowing for distance of ',num2str(cpB_dists(d),2),' pixels.'], ...
                   'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')
    end
          
end





% To displace Precision, Recall, F-measure (mean, std, max) data.
if(0)
    % Precision and Recall (computed new way using N, Dp & Dr - 12/31/15)
    disp('Precision [old, new_overlap, new_cP1, new_cP2]')
    [cntP0./sumP, mean(P_overlap,2)] %, mean(P_cP1,2), mean(P_cP2,2)]

    disp('Precision New Way Direct Pixel Overlap -0- (mean, std, max)')
    [P_overlap, zeros(size(thresh)), mean(P_overlap,2), std(P_overlap,[],2), Pmax]
    
    disp('Precision New Way correspondPixelsB (match3) -0- (mean, std, max)')
    for d = 1:numel(cpB_dists)
        d
        [P_cP3(:,:,d), zeros(size(thresh)), mean(P_cP3(:,:,d),2), std(P_cP3(:,:,d),[],2), Pmax_cP3(:,d), Fmax_whichGT_cP3(:,d)]
    end
    
    disp('%')
    disp('%%%%%')
    disp('%')

    disp('Recall [old, new_overlap, new_cP1, new_cP2]')
    [cntR0./sumR, mean(R_overlap,2)] % , mean(R_cP1,2), mean(R_cP2,2)]

    disp('Recall New Way Direct Pixel Overlap -0- (mean, std, max)')
    [R_overlap, zeros(size(thresh)), mean(R_overlap,2), std(R_overlap,[],2), Rmax]
    
    disp('Recall New Way correspondPixelsB (match3) -0- (mean, std, max)')
    for d = 1:numel(cpB_dists)
        d
        [R_cP3(:,:,d), zeros(size(thresh)), mean(R_cP3(:,:,d),2), std(R_cP3(:,:,d),[],2), Rmax_cP3(:,d), Fmax_whichGT_cP3(:,d)]
    end
    
    disp('%')
    disp('%%%%%')
    disp('%')

    disp('F-measure New Way Direct Pixel Overlap -0- (mean, std, max)')
    [F_overlap, zeros(size(thresh)), mean(F_overlap,2), std(F_overlap,[],2), Fmax, Fmax_whichGT]
    
    disp('F-measure New Way correspondPixelsB (match3) -0- (mean, std, max)')
    for d = 1:numel(cpB_dists)
        d
        [F_cP3(:,:,d), zeros(size(thresh)), mean(F_cP3(:,:,d),2), std(F_cP3(:,:,d),[],2), Fmax_cP3(:,d), Fmax_whichGT_cP3(:,d)]
    end
    
end


% Save output the _ev2.txt files for each image patch.
for d = 1:numel(cpB_dists)
    
    prdFile = [prFile(1:end-7),'d',num2str(d),'_',prFile(end-6:end)];
    disp(['Saving output file : ',prdFile])

    fid = fopen(prdFile,'w');
    if fid==-1,
        error('Could not open file %s for writing.', prdFile);
    end


    fprintf(fid,'%10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g\n',...
        [thresh, repmat(cpB_dists(d),numel(thresh),1), ...
        mean(R_cP3(:,:,d),2), std(R_cP3(:,:,d),[],2), Rmax_cP3(:,d), ...
        mean(P_cP3(:,:,d),2), std(P_cP3(:,:,d),[],2), Pmax_cP3(:,d), ...
        mean(F_cP3(:,:,d),2), std(F_cP3(:,:,d),[],2), Fmax_cP3(:,d), ...
        repmat(numel(groundTruth),numel(thresh),1), Fmax_whichGT_cP3(:,d)]');
        % [thresh, dt, mnR, stdR, Rmax, meanP, stdP, Pmax, meanF, stdF, Fmax, #GT, BestGT]
    fclose(fid);

end


% keyboard