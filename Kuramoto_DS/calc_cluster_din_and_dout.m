function [ dm, ds ] = calc_cluster_din_and_dout(num_tg, num_clust, GT, Y)

% Syntax: [ dm, ds ] = calc_cluster_din_and_dout(num_tg, num_clust, GT, Y)
%
% Compute mean and stardard deviation of pairwise distance for all cluster
% pairs at different points in time. For k clusters and t time points, dm
% and ds are k x k x t. The diagonal values indicate mn & std within each
% cluster. The strict upper triangular values indicate mn & std across each
% combination of cluster pairs. The value in the extreme lower part of
% matrix of each time point indicates the mean and std for all oscillator
% pairs in the whole network regardless of their cluster designation.
%
% Input: (num_tg, num_clust, GT, Y)
%
%   num_tg = number of time points for which we have pairwise phase
%            differences (note a distance at single time point may be 
%            computed over a phase trajectory).
%   num_clust = number of clusters in network
%   GT = ground truth indicating which cluster each node belongs to
%   Y = a 3D tensor containing phase difference for each oscillator pair at
%       each time point (num_clust x num_clust x num_tg)
%
% Output: [dm, ds]
%
%   dm = matrix of mean pairwise phase difference between osc in a given
%        cluster pair. (num_clust x num_clust x num_tg)
%
%   ds = matrix of std of pairwise phase difference between osc in a given
%        cluster pair. (num_clust x num_clust x num_tg)




% , dinn, doutt, dtott, drat_IO, ddiff_IO, drat_IT, ddiff_IT
%
% Compute <d_in> & <d_out> & <d_rat> OR normalized <d_rat> for each of the using phase & phase
% velocity for each of the 5 methods of grouping time (different y's).
%
% TURN THIS INTO A FUNCTION TO CALL
% inputs:  num_tg, num_clust, GT, Y = y1 or yv1, 
% outputs: dm, ds, dinn, doutt, dtott, drat_IO, ddiff_IO, drat_IT, ddiff_IT


dm = zeros(num_clust, num_clust, num_tg);
ds = zeros(num_clust, num_clust, num_tg);
% dinn = zeros(num_tg,4);
% doutt = zeros(num_tg,4);
% dtott = zeros(num_tg,4);

for j = 1:num_tg
    
    for c1 = 1:num_clust
        %
        ind1 = find(GT==c1);
        thing = triu(squeeze(Y(ind1,ind1,j)),1); % d_in: take upper triangle to not bias things with n_c zero values.
        dm(c1,c1,j) = mean(thing(thing>0));
        ds(c1,c1,j) = std( thing(thing>0));
        %
        for c2 = c1+1:num_clust
            ind2 = find(GT==c2);
            thing = squeeze(Y(ind1,ind2,j));     % d_out: pairwise distance between oscillator pairs in 2 different clusters.
            dm(c1,c2,j) = mean(thing(:));
            ds(c1,c2,j) = std( thing(:));
        end
    end
    %
    dtot = triu( squeeze(Y(:,:,j)) ,1);         % d_tot: pairwise distance between all oscillator pairs in the network.
    dtot = dtot(dtot > 10e-15);
    if isempty(dtot)
        dtot = nan;
    end
    dm(num_clust,1,j) = mean(dtot);
    ds(num_clust,1,j) = std(dtot);
    
    
    
    
    
%     %
%     % Get Stats on
%     thing = squeeze(dm(:,:,j));
%     din = diag(thing);
%     dinn(j,:) = [ max(din), mean(din), std(din), min(din) ];
%     %
%     dout = triu(thing,1);
%     dout = dout(dout>0);
%     if isempty(dout)
%         dout = nan;
%     end
%     doutt(j,:) = [ max(dout), mean(dout), std(dout), min(dout) ];
%     %
%     dtot = triu( squeeze(Y(:,:,j)) ,1);
%     dtot = dtot(dtot > 10e-15);
%     if isempty(dtot)
%         dtot = nan;
%     end
%     dtott(j,:) = [ max(dtot), mean(dtot), std(dtot), min(dtot) ];

end


% drat_IO = doutt ./ dinn;
% drat_IT = dtott ./ dinn;
% ddiff_IO = doutt - dinn;
% ddiff_IT = dtott - dinn;