function [a] = calculate_Pmetric()

% This will take in a single DistsPW_Results_allPatches mat file and the
% DoutIdeal_allPatches mat file.  It will loop through entries in the
% DistsPW file and match them up with entries in the DoutIdeal file by the
% ImgPtchID and ImgGtID.  Then it will calculate the value of the P-metric
% using numbers from both.

% ... maybe later