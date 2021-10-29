function [dirPre, sizeGoodIm] = onCluster

% syntax:  [dirPre] = onCluster
%
% This will return current directory to prepend to save and load paths if
% working on a machine other than the cluster.  Will prepend save and load
% paths with my scratch directory if hostname call returns signature from
% cortex cluster.

[x,comp] = system('hostname');

if ~isempty(strfind(comp,'cortex'))
    %disp('On Cluster')
    dirPre = '/clusterfs/cortex/scratch/cwarner/'; % my scratch directory on Cortex
    sizeGoodIm = [0 0 0.5 1];
% elseif ~isempty(strfind(comp,'edison'))
%     %disp('On Cluster')
%     dirPre = '/scratch1/scratchdirs/cwarner/'; % my scratch directory on NERSC
%     sizeGoodIm = [0 0 1 1]; % was [0 0 1 0.5] but plots all elongated horizontally.
    
else
    %disp('On Local Machine')
    dirPre = './';
    sizeGoodIm = [0 0 1 1];
end