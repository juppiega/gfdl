function plotStreamfunMax(r_inv_days, eps_inv_hours)

h = (1:0.1:7.5)*1E3; l = (2000:100:1E4)*1E3;
r = 1/(r_inv_days*86400);
eps = 1/(eps_inv_hours*3600);
N2 = 1.5E-4;
g = 9.81;
T0 = 250;
Q0 = 3*g/T0/86400; % m/s^3
rho0 = 0.67;
dy = 2*pi*6371E3*40/360;
lambda = r*eps/N2

[L,H] = meshgrid(l,h);
Q = zeros(size(L)) + Q0;
%r = 1./(eps * (pi*L./(sqrt(N2)*15E3)).^2);
strfun_max = rho0*dy*computeStreamfunMax(Q,L,H,sqrt(N2),r,eps)/1E9;

figure;
colormap(parula(23))
contourf(L/1000,H/1E3,strfun_max,'edgecolor','none');
fs = 15;
set(gca,'fontsize',fs)
xlabel('L [km]','fontsize',fs)
ylabel('H^z [km]','fontsize',fs)
title('\textbf{a) Domain-prop., $\Lambda = 3.6 \cdot 10^{-8}$}', 'fontsize',fs, 'interpreter','latex')
c1=colorbar;
clim=[5,120];
caxis(clim)
view(2)
axis tight
%hold all
%contour3(L/1000,H/1000,strfun_max,'k')

Q = Q.*(2.5E3./H).*(50/360*2*pi*6371E3./L);
strfun_max = rho0*dy*(computeStreamfunMax(Q,L,H,sqrt(N2),r,eps))/1E9;
% H_s = sqrt(N2)*H./sqrt(r*eps);
% sqrt_part = 1./sqrt(2*L.^2+3*L.*H_s+H_s.^2);

figure;
colormap(parula(23))
contourf(L/1000,H/1E3,strfun_max,'edgecolor','none');
fs = 15;
set(gca,'fontsize',fs)
c2=colorbar;
caxis(clim)
xlabel('L [km]','fontsize',fs)
ylabel('H^z [km]','fontsize',fs)
title('\textbf{b) Constant total, $\Lambda = 3.6 \cdot 10^{-8}$}', 'fontsize',fs, 'interpreter','latex')
% mintick = 1E4*ceil(min(10.^strfun_max(:))/1E4);
% maxtick = 1E4*floor(max(10.^strfun_max(:))/1E4);
% ticks = linspace(log10(mintick),log10(maxtick),10);
% set(cb,'Ticks',ticks)
% set(cb,'TickLabels',1E4*round(10.^ticks/1E4)/1E6)
view(2)
axis tight
%hold all
%contour3(L/1000,H/1000,strfun_max,'k')

end