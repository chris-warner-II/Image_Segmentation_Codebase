function plot_DivMarg(MC,kurParams,kurflags,netflags)


DivMarg = MC{1}.DistAvgPW(:,1) ./ MC{1}.DistAvgPW(:,2);

plot(linspace(0,kurParams.Tsec,kurParams.T./kurParams.spp+1), DivMarg'), 
xlim([0 kurParams.Tsec]), ylim([0 1.1])
hold on, plot([0 kurParams.Tsec], [0.385 0.385], 'k--','LineWidth',2)
%title([netflags.imageIn,' ',fileSubset,' ',netflags.method,' ',kurflags.KurTitleTag,' Kscale ',num2str(kurParams.Kscale)],'FontSize',14,'FontWeight','Bold')
xlabel('Time')
ylabel('Divisive Margin')
