%% This script saves the flags that go into constructing a graph / network from an image into a data structure.

if ~exist('netflags','var')
    netflags.AvgAssoc = AvgAssoc;
    netflags.Modularity = Modularity;
    netflags.maxEnt = maxEnt;
    netflags.maskFlg = maskFlg;
    netflags.GraphLap = GraphLap;
    netflags.NormCut = NormCut;
    netflags.shiftEvals = shiftEvals;
    netflags.normlze = normlze;
    netflags.neg = neg;         % Note: All but neg and proj are input into this SegmentMethod function
    netflags.proj = proj;
    netflags.topo = topo;
    netflags.maxMEiter = maxMEiter;
    %netflags.evectors = evectors;
    %
    %netflags.save_Evecs = save_Evecs;
    % netflags.vizSlope = vizSlope; 
    netflags.save_Evecs_Mat = save_Evecs_Mat;
    
    if(blur_image | DoG_image)
        netflags.imageIn = [InputImages,'_',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr),'_blur_sig',sB];  
    else
        netflags.imageIn = [InputImages,'_',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr)];
    end

    netflags.fname =  fname(1:end-4);
    
else
    
    %disp('Not creating: netflags already exists.')
    %netflags
    
end