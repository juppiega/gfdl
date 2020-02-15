function sst_frierson()

lon_width = 70; 
max_sst = 32+273.15;
lat_width = 40;
min_eq_sst = max_sst - 10;
filename = 'sst_walker_70.nc';

lat = ncread('atmos_average.fr_unmodified.nc','lat');
lon = ncread('atmos_average.fr_unmodified.nc','lon');

sst = ncread('atmos_average.fr_unmodified.nc','t_surf')-8;
sst_mean = mean(sst(:,:,end),1);
figure;
plot(lat, sst_mean)
ampl = max_sst - max(sst_mean);

sst = zeros(length(lon), length(lat));
for i = 1:length(lon)
   sst_pert = ampl*exp(-lat.^2/lat_width^2).*exp(-(lon(i)-180)^2/lon_width^2);
   sst(i,:) = sst_mean + sst_pert';
end

figure;
[X,Y] = meshgrid(lon,lat);
surf(X,Y,sst','edgecolor','none')
view(2);
colorbar
caxis([max_sst-15, max_sst])

sst_mean = mean(sst(:,:,end),1);
figure;
plot(lat, sst_mean)
ampl = max_sst - max(sst_mean);

copyfile('atmos_average.fr_unmodified.nc',filename);
ncwrite(filename, 'flux_oceanq', sst);

end