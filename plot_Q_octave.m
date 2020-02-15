function plot_Q_octave()

fs = 15;
load mean_Q_int_hs_ca_wt_37.mat
figure
plot(lon,mean_Q_int,'linewidth',1.5,'color','k');
xlabel('Longitude','fontsize',fs)
ylabel('W/m^2','fontsize',fs)
title('a) Q','fontsize',fs)
xlim([0,360])

m = mean(mean_Q_int);
hold all
plot([0,360],[m,m],'linewidth',1.0,'color','k');
set(gca,'fontsize',fs)

end