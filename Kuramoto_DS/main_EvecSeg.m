function [Vdom] = main_EvecSeg(K)

% This script will setup inputs to Kuramoto.m function, call it, make
% diagnostic plots for output, and save results.



if ~exist('K','var')
    K = build_K2(netParams); % Build coupling matrix for hand crafted network using Win, Wout, C, Rmax, etc in kurParams.
end



%% Dominant Eigenvector of Coupling Matrix K - To compare with results of Kuramoto DS simulation.

[V,D] = eig(K);
dom = find(diag(D)==max(diag(D)));
dom = dom(1);
Vdom = V(:,dom);

% Vdom_pair = abs(repmat(Vdom,1,numel(Vdom)) - repmat(Vdom',numel(Vdom),1));




