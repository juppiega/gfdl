function [] = testHeating ()

x = -180:4:180;
p = 0:25:1000;
y = -90:1:90;
[X,Y,P] = meshgrid(x, y, p);
pwidth = 350;
pcenter = 500;
xwidth = 40;
ywidth = 30;
amp = 6;

lon_fac = (X./xwidth).^2;
lat_fac = (Y./ywidth).^2;
p_fac = ((P-pcenter)/pwidth).^2;
c = 1 - lon_fac - p_fac;
Q = zeros(size(X));
ind = c - lat_fac > 0;
Q(ind) = 1.5*amp*(c(ind)-lat_fac(ind))./sqrt(c(ind));



lat_max = 30;
lat_ind = -lat_max <= y & y <= lat_max;
Q_aver = squeeze(mean(Q(lat_ind,:,:),1));
Q_exact = zeros(size(Q_aver))';


[X,P] = meshgrid(x,p);
lon_fac = (X./xwidth).^2;
p_fac = ((P-pcenter)/pwidth).^2;
ind = 1 - lon_fac - p_fac > 0;
Q_exact(ind) = amp*(1 - lon_fac(ind) - p_fac(ind));

figure;
surf(X,P,Q_aver','edgecolor','none');
view(2);
colorbar;
axis tight;
climits = caxis;
title('Numerical')

figure;
surf(X,P,Q_exact,'edgecolor','none');
view(2);
colorbar;
caxis(climits);
axis tight;
title('Analytic')

figure;
surf(X,P,Q_aver'-Q_exact,'edgecolor','none');
view(2);
colorbar;
axis tight;
title('Numerical - Analytic')


end
