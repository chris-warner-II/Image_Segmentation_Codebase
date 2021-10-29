% This script is like a probe that you stick in your script or function to
% determine how much memory you are using at that point.
%
% But I have to think of how to use this.  WHere in function does it make
% sense to use this?


memUsed = 0;
xxx = whos;

for bbb = 1:numel(xxx)
    memUsed = memUsed + xxx(bbb).bytes;

end

% disp(['Matlab variables currently using ',num2str(memUsed),' bytes.'])
disp(['All Matlab variables currently using ',num2str(memUsed/1000000),' MB.'])