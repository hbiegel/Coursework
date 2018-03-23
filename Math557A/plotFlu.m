% plot influenza data to exhibit seasonality

load('wILIvec.mat') %load weighted ILI data downloaded from the CDC FluView
wILI = wILI(747:955);
% loads wILI
T = [0:length(wILI)-1]/52 + 2012;

flufig = figure();
hold on
plot([2012,2012],[0,7],'--','Color',[169,169,169]/255)
plot([2013,2013],[0,7],'--','Color',[169,169,169]/255)
plot([2014,2014],[0,7],'--','Color',[169,169,169]/255)
plot([2015,2015],[0,7],'--','Color',[169,169,169]/255)
plot([2016,2016],[0,7],'--','Color',[169,169,169]/255)
plot(T,wILI,'LineWidth',4,'Color',[51,153,255]/255)
set(gca,'FontSize',16)
title('Influenza-Like Illness Incidence, 2012-2015')
xlim([2012,2016])
xticks([2012,2013,2014,2015,2016])
ylim([0,7])
xlabel('Time (year)')
ylabel('Incidence (% of total hospital cases)')
hold off

saveas(flufig,'seasonalflu.png')