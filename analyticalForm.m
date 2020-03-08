function psim2s = analyticalForm (x, p, p_max, L, Hz, Q0, N2, r, eps)
% x in degrees (-180 to 180), p in hPa (via meshgrid), create 2D with meshgrid. L and H (actual geoetric height) in meters
format compact
x = x/360 * 2*pi*6371E3;
H = Hz * sqrt(N2/(r*eps));
H,L,Q0
rho0 = 0.67;
T0 = 250;
cp = 1005;
y = 100*(p_max - p)/rho0/9.81 * sqrt(N2/(r*eps));

ksi = x/L + 1i*y/H; % [unitless]
w = H*x/L + 1i*L*y/H; % [m]
z = x - 1i*y;         % [m]

E = zeros(size(x));
E_inside = Q0*ksi.*(1 - ksi.*(2*z + conj(w))/(3*(H+L)) )/(L+H)*L*H/N2;
inside = 1 - (x/L).^2 - (y/H).^2 > 0;
E(inside) = E_inside(inside);

s = z.*sqrt(1 - (L^2 - H^2)./(z.^2));
E_outside = 2*Q0./z .* (2*z ./ (z + s)).^2 .* ((z + 2*s)./(3*z))*L*H/N2/8;
E(~inside) = E_outside(~inside);

%display(['Analytic_sampled_from_gridpoints: ', num2str(max(real(E(:))))]);
Qtot = cp*T0*Q0*rho0*pi*6371E3*40/180*Hz*L*pi/(2*9.81);
analytic_exact = 4/3/pi * 9.81*Qtot/(cp*T0*N2)/sqrt(2*L^2 + 3*L*H + H^2)*H/Hz;%rho0*pi*6371E3*40/180*Q0*2*L*H/(3*N2)/sqrt(2*L^2 + 3*L*H + H^2);
analytic_m2s = max(max(real(E)));
display(['Analytic exact [10^9 kg/s]: ', num2str(analytic_exact/1E9)]);
display(['Analytic exact [m^2/s]: ', num2str(analytic_m2s)]);
display(['Qtot [PW]: ', num2str(Qtot/1E15)]);

heat = 3 * (1 - (x/L).^2 - (y/H).^2);
heat(~inside) = 0;

x = x/(2*pi*6371E3)*360;
psim2s = real(E);
psi = real(E)*rho0*pi*6371E3*40/180;
x_b = [-L:1E3:L,L];
y_b = H*sqrt(1-(x_b/L).^2);
x_b = [x_b, x_b];
y_b = [y_b, -y_b];
p_b = p(end,1)/2 - y_b / sqrt(N2/(r*eps)) *rho0*9.81/100;
x_b = x_b/(2*pi*6371E3)*360;

z = (1000 - p)*100/rho0/9.81/1E3;

figure;
contourf(x,z,psi/1E9,-30:2.5:30,'edgecolor','none')
view(2); 
axis tight; 
xlabel('Longitude','fontsize',15)
ylabel('z [km]','fontsize',15)
set(gca,'fontsize',15)
%set(gca,'ydir','reverse'); 
title('b) Analytic streamfunction [Tg/s]', 'fontsize',14)
colormap(jet(24))
colorbar('fontsize',15)
caxis([-30,30])
% ma = max(psi); mi  = min(psi);
% caxis([min([-ma,mi]), max([ma,-mi])])
hold all
ylim([0,15])
%contour3(x,p,psi,'k')
N = length(x_b)/2;
%plot3(x_b(1:N),p_b(1:N),max([ma,-mi])+zeros(1,N),'color','r','linewidth',2.0)
%plot3(x_b(N+1:end),p_b(N+1:end),max([ma,-mi])+zeros(1,N),'color','r','linewidth',2.0)

figure;
colormap(copper(6))
contourf(x,z,heat,0.5:0.5:3.0,'edgecolor','none')
view(2); 
axis tight; 
xlabel('Longitude','fontsize',15)
ylabel('z [km]','fontsize',15)
set(gca,'fontsize',15)
%set(gca,'ydir','reverse'); 
title('a) Heating [K/day]','fontsize',14)
colorbar('fontsize',15)
ylim([0,15])
caxis([0,3])
set(gca,'Color','k')

%
%figure;
%surf(x,y,abs(E),'edgecolor','none')
%view(2); 
%axis tight; 
%set(gca,'fontsize',15)
%colorbar('fontsize',15)
%
%figure;
%surf(x,y,arg(E),'edgecolor','none')
%view(2); 
%axis tight; 
%set(gca,'fontsize',15)
%colorbar('fontsize',15)

%
%figure;
%surf(x,y,imag(E),'edgecolor','none')
%view(2)
%colorbar
%title('Im(E)')
%
%figure;
%plot(x(1,:), real(E(100,:)))



end
