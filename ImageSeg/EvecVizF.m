function [xViz] = EvecVizF(x,a)

% syntax: xViz = EvecVizF(x,a)
% 
% Function to display Eigenvector that is better than just taking log of it.
%                   It is a transformation to the eigenvector so that its
% structure is more easily seen by eye.  It is a monotonically increasing
% function that does not saturate and has an arbitrarily steep slope at
% (and narrowly around) one point with the slope controlled by the a
% parameter.  Small a means steeper slope.

xViz = log( x./sqrt(a) + sqrt(1 + (x.^2)./a ) );