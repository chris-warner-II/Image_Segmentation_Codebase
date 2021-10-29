iev = 1; % which eigenvector to look at.

kurParams_setup
vizNonLins = logspace(0,-16,17);
figure
for bb = 1:numel(vizNonLins)

    [MC] = metaClusterAnalysis(reshape(EvecVizF(EvecsML{iev},vizNonLins(bb)),numel(EvecsML{1}),1), netParams, kurParams, 0);

    subplot(5,4,bb),imagesc( EvecVizF(EvecsML{iev},vizNonLins(bb)) ), colorbar, axis off, title(num2str(vizNonLins(bb),2))
    vizVsSep(bb) = meanCSep1;

end

subplot(5,4,[bb+1:5*4]), semilogx(vizNonLins,vizVsSep), 
title([netflags.imageIn,' - ',netflags.method,' - sD ',netflags.sD,' - rM ',netflags.rM,' - sP ',netflags.sP])
xlabel('Nonlinearity Strength')
ylabel('Separation')

