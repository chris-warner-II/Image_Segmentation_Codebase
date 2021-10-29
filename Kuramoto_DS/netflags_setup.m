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
    netflags.evectors = evectors;
    %
    netflags.save_Evecs = save_Evecs;
    netflags.vizSlope = vizSlope; 
    netflags.save_Mat = save_Mat;
    
else
    
    disp('Not creating: netflags already exists.')
    netflags
    
end