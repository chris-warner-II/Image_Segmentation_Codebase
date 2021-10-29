function out = TopDownPyramid(x1,y1,x2,x3,y2,y3,outbeg,outfin,inbeg,infin,disp)

% This function makes an rectangular image of size (x1,y1).  The image
% contains a rectangle spanning from (x2-x3) and (y2-y3) within the
% rectangular image.  
%
% The external rectangle has a horizontal greyscale
% gradient going from dark on the left to light on the right.  The internal
% rectangle has a horizontal greyscale gradient going from light on the 
% left to dark on the right.  This second test will be a slightly harder
% image to segment.


% Resize image and rectangle to make number of pixels divisible by 4.
% (useful in building pyramid)
x3 = x3 - mod(x3-x2,4);
y3 = y3 - mod(y3-y2,4);
x1 = x1 - mod(x1,4);
y1 = y1 - mod(y1,4);


%% make internal rectangle be a gradient coming down from center
PeakIn = [(x3-x2)/2 , (y3-y2)/2]; % can decide to put peak point elsewhere too.

TentInX = [linspace(inbeg,infin,PeakIn(1)),linspace(infin,inbeg,PeakIn(1))];
TentInX = repmat(TentInX,2*PeakIn(2),1);
TentInY = [linspace(inbeg,infin,PeakIn(2)),linspace(infin,inbeg,PeakIn(2))]';
TentInY = repmat(TentInY,1,2*PeakIn(1));

RectIn = TentInX.*TentInY; % build inside rectangle

%% make external rectangle be a gradient coming down from center
PeakOut = [x1/2 , x1/2]; % can decide to put peak point elsewhere too.

TentOutX = [linspace(outbeg,outfin,PeakOut(1)),linspace(outfin,outbeg,PeakOut(1))];
TentOutX = repmat(TentOutX,2*PeakOut(2),1);
TentOutY = [linspace(outbeg,outfin,PeakOut(2)),linspace(outfin,outbeg,PeakOut(2))]';
TentOutY = repmat(TentOutY,1,2*PeakOut(1));

RectOut = TentOutX.*TentOutY; % build outside rectangle
RectOut(x2+1:x3,y2+1:y3) = RectIn; % put inside rectangle inside outside one

%% Display image to be segmented
if(1)
    figure, imagesc(RectOut),colorbar
end