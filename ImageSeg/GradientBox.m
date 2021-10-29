function out = GradientBox(x1,y1,x2,x3,y2,y3,outbeg,outfin,inbeg,infin,disp)

% This function makes an rectangular image of size (x1,y1).  The image
% contains a rectangle spanning from (x2-x3) and (y2-y3) within the
% rectangular image.  The external rectangle has a horizontal greyscale
% gradient going from dark on the left to light on the right.  The internal
% rectangle has a horizontal greyscale gradient going from light on the 
% left to dark on the right.  This second test will be a slightly harder
% image to segment.

rec_out = repmat( linspace(outbeg,outfin,x1) ,y1 ,1 );
rec_in = fliplr( repmat( linspace(inbeg,infin,x3-x2+1) ,y3-y2+1 ,1 ) );

out = rec_out;
out(y2:y3,x2:x3) = rec_in;

out = out./max(max(out)); % Normalize output to be from 0 to 1.

if disp
   % plot something if you wanna 
end