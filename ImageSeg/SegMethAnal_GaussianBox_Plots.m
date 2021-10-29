% Script to loop through plotting code and make plots.  This will probably
% take a while to run ....
%


SegAt2 = {'Optimal','MeanPix'};

for pp = 1:numel(SegAt2)
    
    SegAtIn = SegAt2{pp};
    
    disp(['SegDataViz1 - Seg At ',SegAtIn])
    SegDataViz1
    %
    disp(['SegDataViz2 - Seg At ',SegAtIn])
    SegDataViz2
    %
    disp(['SegDataViz3 - Seg At ',SegAtIn])
    SegDataViz3
    %
    disp(['SegDataViz3b - Seg At ',SegAtIn])
    SegDataViz3b
    %
    disp(['SegDataViz4 - Seg At ',SegAtIn])
    SegDataViz4
    %
    disp(['SegDataViz4b - Seg At ',SegAtIn])
    SegDataViz4b
    %
    disp(['SubplotEvecNSeg - Seg At ',SegAtIn])
    SubplotEvecNSeg
    
end

