function plotQ0figure()

lon = 0:1:360;
Q0_ref = 3*(1-((lon-180)/50).^2); 
Q0_ref(Q0_ref<0) = nan;
Q0_ct = 100/160*3*(1-((lon-180)/80).^2); 
Q0_ca = 3*(1-((lon-180)/80).^2); 
Q0_ct(Q0_ct<0) = nan;
Q0_ca(Q0_ca<0) = nan;

figure;
plot(lon,Q0_ref,'color','k','linewidth',2.0)
hold all;
plot(lon,Q0_ct,'k--')
plot(lon,Q0_ca,'color','k')
fs = 15;
set(gca,'fontsize',fs-1)
ylabel('$\widehat{Q}$ [K/day]','fontsize',fs,'interpreter','latex')
xlabel('Longitude','fontsize',fs,'interpreter','latex')
%legend('Original','Domain-proportional','Constant total')
xlim([80,280])


end