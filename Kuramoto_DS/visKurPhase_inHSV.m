function phaseFinal = visKurPhase_inHSV(img,phase)

% Convert final phase solution back into linear variable to compare to pixels (just for visualization).
[one, two] = find(img == min(img(:))); % find brightest pixel.
pixMax(1) = one(1);
pixMax(2) = two(1);

phaseFinal = phase;
phaseOffset = pi - phaseFinal(pixMax(1),pixMax(2));           % global phase to add to all phases to make phase at brightest pixel = pi
phaseFinal = mod( phaseFinal + phaseOffset , 2*pi);           % rotate all phases by phase offset
% phaseFinal = abs(phaseFinal - pi) ./ pi;                    % fold circular variables down to [0 1] number line - ie. pixel space


    