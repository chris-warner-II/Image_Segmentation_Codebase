% A script to directly compare a method (such as Avg Assoc. or Graph Lap.)
% and its normalized version.

%% Load MAT files
TODO = 'GL'; %'GL' or 'AA'
AA = load(TODO);
AAn = load([TODO,'n']);
load(TODO)

%% Plot Auxiliary Matrices.
figure
subplot(121), imagesc(AA.Q), colorbar, title(TODO)
subplot(122), imagesc(AAn.Q), colorbar, title([TODO,'n'])

%% Plot Eigenvectors.
figure
for i = 1:numel(evectors)
    subplot(numel(evectors),2,2*i-1), imagesc(AA.EvecsML{evectors(i)}), colorbar, xlabel(['Evec ',num2str(evectors(i)),' (\lambda=',num2str(AA.EvalsML(i)),')'])
    subplot(numel(evectors),2,2*i), imagesc(AAn.EvecsML{i}), colorbar, xlabel(['Evec ',num2str(evectors(i)),' (\lambda=',num2str(AAn.EvalsML(i)),')'])
end

subplot(numel(evectors),2,1), title(['Eigenvectors:  ', TODO])
subplot(numel(evectors),2,2), title([TODO,'n'])

%% Plot Segmentations (when thresholding Eigenvector at Mean Value.
figure
for i = 1:numel(evectors)
    subplot(numel(evectors),2,2*i-1), imagesc(AA.seg{evectors(i)}), colorbar, xlabel(['Evec ',num2str(evectors(i)),' (\lambda=',num2str(AA.EvalsML(i)),')'])
    subplot(numel(evectors),2,2*i), imagesc(AAn.seg{i}), colorbar, xlabel(['Evec ',num2str(evectors(i)),' (\lambda=',num2str(AAn.EvalsML(i)),')'])
end

subplot(numel(evectors),2,1), title(['Segmentation at mean:  ', TODO])
subplot(numel(evectors),2,2), title([TODO,'n'])