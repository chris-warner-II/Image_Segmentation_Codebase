function [In] = computeDiscriminability(In)


% NOTE: I AM NOT SURE I AM EVER USING THIS FUNCTION.  IT SEEMS THAT FIX
% FILES FUNCTIONS AND SEGMENTMETHOD AND MAIN_KURAMOTO ALL CALL THE
% COMPUTE_BOUDNARY



In.boundaryDiscriminability = ( In.meanGradientOnBoundary - In.meanGradientOffBoundary ) ...
    ./ sqrt( 0.5*( In.stdGradientOnBoundary.^2 + In.stdGradientOffBoundary.^2 ) );

% BoundaryDiscriminability was not right because previously because I was squaring the standard deviations.
%
% Actually, DAMMIT, squaring standard deviations WAS the right thing to do.
% Not squaring them causes this linear bias where d' changes if you
% multiply the input (before computing spatial gradients) with a scalar
% value). 

