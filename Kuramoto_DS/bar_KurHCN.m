function bar_KurHCN(indx, data, I_Wout, ytag, params)

% Define Colors for Points on plot.
colors = colormap('jet');
num_colors = size(colors,1);
num_Wout = numel(I_Wout);
cstep = floor(num_colors/num_Wout);

% plot
b = 0;
for L = indx
    b = b+1;
    bar(L,data(L),'FaceColor',colors(cstep*mod(b,num_Wout)+1,:))
    text(L, data(L), num2str(data(L),2),...
        'VerticalAlignment','Bottom','HorizontalAlignment','Center','Color',colors(cstep*mod(b,num_Wout)+1,:))
end
axis([0 numel(indx)+1 0 1.2* max( data )])
for p = 1:round( numel(indx) / numel(I_Wout) )
    plot([p*numel(I_Wout)+0.5 p*numel(I_Wout)+0.5],[0 1.2*max( data )], 'k--')
end
ylabel(ytag,'FontSize',14,'FontWeight','Bold')
set(gca,'XTick',1:numel(data),'XTickLabel',params,'FontSize',16,'FontWeight','Bold')
rotateXLabels(gca(),90)