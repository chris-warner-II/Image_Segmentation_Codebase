function stdErr = sem(data) %,which_errbars)

% function to compute the standard error of the mean (instead of std) for
% smaller errorbars.


% switch which_errbars
%     case 'sem'
%         %disp('sem')
stdErr = std(data)./sqrt(numel(data));
%     case 'std'
%         %disp('std')
%         stdErr = std(data);
% end
