function [] = GradientBoxLoopGen(boxDims)

% Create an ENSEMBLE (?) OF IMAGES that have a box inside a larger box.
% Both boxes have pixel gradients, but they go in opposite directions.  The
% mean pixel value of both boxes are the same so that segmenting by mean
% pixel value can not work.  


% OTHER IDEAS FOR TESTS::
% INSTEAD OF GRADIENTS GOING FROM 0 TO 1, THEY CAN BE SMALLER...
% INSIDE BOX MOVES AROUND INSIDE OUTSIDE BOX TO INVESTIGATE EDGE EFFECTS...
% MORE COMPLICATED EDGES AND CONTOURS THAN JUST BOX... (Curvature effects segmentation?)
 
% Check if you are using cluster to adjust home directory
dirPre = onCluster;

ds_fctr = 1;

for i = 1:numel(boxDims)
    
    ximg = boxDims(i); % Horizontal
    yimg = boxDims(i); % Vertical

    % Make directories to house stacks of input images generated
    if ~exist([dirPre,'images/GradientBox/',num2str(ximg),'x',num2str(yimg),'_ds1'],'dir')
        mkdir([dirPre,'images/GradientBox/',num2str(ximg),'x',num2str(yimg),'_ds1'])
    end

    x2 = round(0.4*ximg); % start and end of x&y of inner rectangle
    x3 = round(0.7*ximg);
    y2 = round(0.4*yimg);
    y3 = round(0.7*yimg);

    disp(['Input Image is ',num2str(ximg),'x',num2str(yimg),' box with greyscale gradient inside and opposing gradient outside.']);
    iids = struct('name','gradBox    '); % note: spaces after are important
    outbeg = 0;
    outfin = 1;
    inbeg = 0;
    infin = 1;

    im = GradientBox(ximg,yimg,x2,x3,y2,y3,outbeg,outfin,inbeg,infin,0);

    fOut = ['GradientBox_',num2str(outbeg),'',num2str(outfin),num2str(inbeg),'',num2str(infin)]; % Note: ds1 means downsampled by 1 (ie. not downsampled)
    fOut = fOut(fOut~='.');

    gT{1} = binary_box(ximg,yimg,x2,x3,y2,y3,0,0); % binary ground truth of where box is.
    gT{1}(gT{1}==-1) = 2;                          % so regions in ground truth are labeled 1 & 2.

    % save each image stack to a mat file to be used as input later.
    save([dirPre,'images/GradientBox/',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr),'/',fOut],'im','gT','outbeg','outfin','inbeg','infin','ximg','yimg','x2','x3','y2','y3','ds_fctr');

    if(0)
        figure, imagesc(im)
        keyboard
    end
    
    clear im

end
