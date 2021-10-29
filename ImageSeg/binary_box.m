function out = binary_box(x1,y1,x2,x3,y2,y3,color,display)

% syntax: out = binary_box(x1,y1,x2,x3,y2,y3,color,display)
% This function makes a binary image of size (x1,y1).  It places within that
% image a rectangle spanning from (x2-x3) and (y2-y3). Color is a flag that
% if 1, puts a white rectangle in a black space and if 0, puts a black
% rectangle in a white space.

% as a default make test1 image to be while with black rectangle
out = 1*ones(y1,x1); % was just ones (not 3*)
out(y2:y3,x2:x3) = -1*ones(1+y3-y2,1+x3-x2); % was zeros()

% if color flag is one, flip pixels
if color
    out = abs(1-out);
end

if display
   % plot something if you wanna 
end