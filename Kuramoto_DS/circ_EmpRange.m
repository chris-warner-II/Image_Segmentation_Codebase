function [EmpRange_best] = circ_EmpRange(vars)

% This function will take in an array of circular (phase) variables and
% compute their Empirical Range using circ_dist.  It is not the most
% efficient thing to do, but what the hell?

EmpRange_best = 0;

for i = 1:numel(vars)

    EmpRange_curr = max(circ_dist(vars(i),vars));
    EmpRange_best = max(EmpRange_best,EmpRange_curr);
    
    if(EmpRange_best==pi) % if range = max possible range.
        break
    end
    
end

%disp(['In circ_EmpRange, did ',num2str(i),' / ',num2str(numel(vars)),' iterations.'])