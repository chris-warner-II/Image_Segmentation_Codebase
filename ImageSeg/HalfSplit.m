function out = HalfSplit(x, y, mu1, sig1, mu2, sig2, plt)

% syntax: out = HalfSplit(x, y, mu1, sig1, mu2, sig2, plt);
%
% This function creates an image imDim x imDim that is split into two
% halves (one light and one dark).
%
% Inputs:
%   x - size of image along horizontal
%   y - size of image along vertical
%   mu1 - mean pixel value in 1st part of image
%   sig1 - spread of pixel values in 1st part of image
%   mu2 - mean pixel value in 2nd part of image
%   sig2 - spread of pixel values in 2nd part of image
%   disp - a flag to display image (or not).
%
% Outputs:
%   out - the image (a matrix of numbers representing the image)
%

%% Make image
rec1 = mu1 + sig1.*randn(x, floor(y/2)); % out(:,1:floor(y/2))
rec2 = mu2 + sig2.*randn(x, floor(y/2)); % out(:,floor(y/2)+1:end)


%% Iterate with while loop bc some may be < 0 or > 1 depending on mu & sigma
redo1 = 1;
while any(any(redo1))
    redo1 = (rec1 < 0 | rec1 > 1);
    rec1(redo1) = 0;
    rec1 = rec1 + ( mu1 + sig1.*randn(x, floor(y/2)) ) .* redo1;
end
%
redo2 = 1;
while any(any(redo2))
    redo2 = (rec2 < 0 | rec2 > 1);
    rec2(redo2) = 0;
    rec2 = rec2 + ( mu2 + sig2.*randn(x, floor(y/2)) ) .* redo2;
end


out = [rec1, rec2];


%% Plot Image
if(plt)
    
    figure, imagesc(out), colormap('bone'), colorbar
    
    
end