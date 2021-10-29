function [] = saveGoodImg(h,fname,fig_pos)

% I found this code online & implemented it in this function because Matlab's
% saveas() functions sucks when you want save jpg's.  Use this instead.
%
% This takes a very long time to save images tho.  But what other options
% do I have?

if ~exist('fig_pos','var')
    fig_pos = [0 0 1 1];
end

drawnow;
set(h,'units','normalized','position', fig_pos)
set(h, 'Color', 'white');   % white bckgr
export_fig( h, ...          % figure handle
    fname, ...              % name of output file without extension
    '-painters', ...        % renderer
    '-jpg', ...             % file format
    '-r100' );               % resolution in dpi (150 or 300 for better res)
