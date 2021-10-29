function [patch,gTpatch,xpbeg,xpfin,ypbeg,ypfin,pnum] = extractPatch(image,xp,yp,groundTruth,fname,dirOut,pnum)

% [patch,xpbeg,xpfin,ypbeg,ypfin] = extractPatch(image,xp,yp,groundTruth,dirOut)
%
% This function will take in an image and patch dimensions and will output
% a patch randomly chosen from that image of the proper dimensions.

%% Setup
xi = size(image,2); % size of input image
yi = size(image,1);
%
xpmax = xi -(1+xp/2); % randomly choose center point for patch
xpmin = 1+xp/2;
%
ypmax = yi -(1+yp/2);
ypmin = 1+yp/2;

% disp(['Finding a Good Patch'])

keepGoing=1; % while loop breakout flag
i=0; % while loop counter variable


% Look through groundTruth Segmentations for Pathalogical cases with one segment.
goodGts = [];
for i = 1:numel(groundTruth)
    x = hist(double(groundTruth{i}.Segmentation(:)), double(unique(groundTruth{i}.Segmentation)));
    y = sort(x,'descend');
    y = y./max(y);
    y = y(2:end);
    if(  y(1) > 0.01 ) % if there are fewer than 5 segments and 2nd largest one is 100x smaller than largest
        goodGts = [goodGts,i]; % numel(y)>4 && -- gettting rid of that 5 segment thing...
    end
end

    
%% Looping through random patches to find good ones (ones where they stradle a segment in truth files)

while(keepGoing)
    
    % randomly pick a center point for patch
    xpmid = xpmin + round((xpmax-xpmin)*rand);
    ypmid = ypmin + round((ypmax-ypmin)*rand);

    xpbeg = xpmid - xp/2;
    xpfin = xpmid + xp/2;
    ypbeg = ypmid - yp/2; 
    ypfin = ypmid + yp/2; 

    

    % look at ground truth and make sure there is some overlap with a region edge
    goodPatch=0;
    for j = 1:numel(goodGts)
        gTpatch{j} = groundTruth{goodGts(j)}.Segmentation(ypbeg:ypfin,xpbeg:xpfin);
        x = hist(double(gTpatch{j}(:)));
        if (numel(find(x))>1)
            y=sort(x,'descend');
            %
            if(y(2)/y(1) > 0.5)
                goodPatch=goodPatch+1;
            end
            %
            if (goodPatch == numel(goodGts))
                keepGoing=0; % if patch overlaps significantly with region boundary in all good Gts
                patch = image(ypbeg:ypfin,xpbeg:xpfin); % take the patch
            end
        end
    end
            
    i = i+1;
    if (i>10000)
        disp('Cant find a good patch.')
        patch=0;
        % keyboard
        keepGoing=0; % Giving up.
    end
            
end

patch = patch - min(patch(:));
patch = patch ./ max(patch(:)); % normalize patch to spread pixel values between 0 & 1.

% keyboard

%% Plot full image and image patch if they dont already exist in the directory
if( (0) && ~exist(['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/patches/',fname,'_patch',num2str(pnum),'_x',num2str(xpmid),'_y',num2str(ypmid),'.tif'],'file') )
   
    xpl = ceil(sqrt(numel(groundTruth)+2));
    ypl = round(sqrt(numel(groundTruth)+2));
    h=figure; hold on, colormap('bone')
    
    subplot(ypl,xpl,1), imagesc(image)
    title([num2str(xp),'x',num2str(yp),' Patch @ (',num2str(xpmid),',',num2str(ypmid),')'],'FontSize',20,'FontWeight','Bold')
    subplot(ypl,xpl,2),imagesc(patch),box on
    set(gca,'xtick',[],'ytick',[]);
    
    for i=1:(numel(groundTruth))
        gTpatch = groundTruth{i}.Segmentation(ypbeg:ypfin,xpbeg:xpfin);
        subplot(ypl,xpl,i+2), imagesc(gTpatch),box on,
        set(gca,'xtick',[],'ytick',[]);
        title(['GT',num2str(i)],'FontSize',18,'FontWeight','Bold')
    end
    
    if ~exist(['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/patches/'],'dir')
        mkdir(['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/patches/']);
    end
    
    saveas(h,['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/patches/',fname,'_patch',num2str(pnum),'_x',num2str(xpmid),'_y',num2str(ypmid)],'tif');
    close
end