function plot_Ebars_KurHCN(indx, means, stds, I_Wout, ytag) 

% Define Colors for Points on plot.
colors = colormap('jet');
num_colors = size(colors,1);
num_Wout = numel(I_Wout);
cstep = floor(num_colors/num_Wout);

% Loop through all data points and plot mean and std error bars colorcoded for Wout parameter
b=0;
for L = indx
    b = b+1;
    errorbar(L, means(L), stds(L), 'Color', colors(cstep*mod(b,num_Wout)+1,:),'LineStyle','None')
    scatter(L, means(L), 'Filled','MarkerFaceColor', colors(cstep*mod(b,num_Wout)+1,:), 'Marker','o')
    text(L, means(L)+stds(L),num2str(means(L),2),...
        'VerticalAlignment','Bottom','HorizontalAlignment','Center','Color',colors(cstep*mod(b,num_Wout)+1,:))
end
ylabel(ytag,'FontSize',14,'FontWeight','Bold')
axis([0 numel(indx)+1 0 1.2* max( means+stds )])
for p = 1:round( numel(indx) / numel(I_Wout) )
    plot([p*numel(I_Wout)+0.5 p*numel(I_Wout)+0.5],[0 1.1*max( means+stds )], 'k--')
end
set(gca,'XTick',[],'XTickLabel',[],'FontSize',16,'FontWeight','Bold')