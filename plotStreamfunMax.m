function plotStreamfunMax()

h = (1:0.2:20)*1E3; l = (1800:100:0.5E4)*1E3;
r = 1/(15*86400);
eps = 1/(5*86400);
N2 = 1E-4;
g = 9.81;
T0 = 250;
Q0 = 6*g/T0/86400; % m/s^3

[L,H] = meshgrid(l,h);
Q = zeros(size(L)) + Q0;
strfun_max = computeStreamfunMax(Q,L,H,sqrt(N2),r,eps);

figure;
surf(L/1000,H/1E3,strfun_max,'edgecolor','none');
fs = 15;
set(gca,'fontsize',fs)
xlabel('L [km]','fontsize',fs)
ylabel('H^z [km]','fontsize',fs)
title('Constant heating amplitude', 'fontsize',fs)
c1=colorbar;
clim=caxis;
view(2)
axis tight
hold all
contour3(L/1000,H/1000,strfun_max,'k')

Q = Q.*(5E3./H).*(2000E3./L);
strfun_max = (computeStreamfunMax(Q,L,H,sqrt(N2),r,eps));

figure;
surf(L/1000,H/1E3,strfun_max,'edgecolor','none');
fs = 15;
set(gca,'fontsize',fs)
c2=colorbar;
caxis(clim)
xlabel('L [km]','fontsize',fs)
ylabel('H^z [km]','fontsize',fs)
title('Constant total heating', 'fontsize',fs)
% mintick = 1E4*ceil(min(10.^strfun_max(:))/1E4);
% maxtick = 1E4*floor(max(10.^strfun_max(:))/1E4);
% ticks = linspace(log10(mintick),log10(maxtick),10);
% set(cb,'Ticks',ticks)
% set(cb,'TickLabels',1E4*round(10.^ticks/1E4)/1E6)
view(2)
axis tight
hold all
contour3(L/1000,H/1000,strfun_max,'k')

end