function q_flux_frierson(lat_width, zm_ampl, latb, lon_width, filename)
% All in degrees and W/m^2

if nargin == 0 % La Nina conditions
    lat_width = 22;
    zm_ampl = 60;
    latb = 18;
    lon_width = 100;
    filename = 'q_flux_walker_100.nc';
end

lat = ncread('atmos_average.fr_unmodified.nc','lat');
lon = ncread('atmos_average.fr_unmodified.nc','lon');

zonal_mean = zm_ampl*(1-2.*lat.^2/lat_width.^2) .* exp(-lat.^2/lat_width.^2) ./ cosd(lat);
figure;
plot(lat, zonal_mean)
ylabel('OHT Divergence [W/m^2]')
xlabel('Lat')

heat_avail = 360*zonal_mean;

oht = zeros(length(lon), length(lat));
oht_lon = (1-exp(-(lon-180).^2/lon_width^2));
A = trapz(lon, oht_lon);

for i = 1:length(lat)
    if -latb <= lat(i) && lat(i) <= latb
        multip = heat_avail(i)/A;
        oht(:,i) = multip * oht_lon;
    else
        oht(:,i) = zonal_mean(i);
    end
end

[X,Y] = meshgrid(lon, lat);
figure;
surf(X,Y,oht','edgecolor','none');
view(2)
axis tight
colorbar

zonal_mean_test = mean(oht,1);
figure;
plot(lat, zonal_mean_test);

copyfile('atmos_average.fr_unmodified.nc',filename);
oht_orig = ncread(filename,'flux_oceanq');
oht = repmat(oht,[1,1,size(oht_orig,3)]);
ncwrite(filename, 'flux_oceanq', oht);


% lat = -89:1:89; x = lat*pi/180;
% lat_width = lat_width * pi/180;
% 
% q = ampl*(1-2.*x.^2/lat_width.^2) .* exp(-x.^2/lat_width.^2) ./ cos(x);
% 
% figure;
% plot(lat, q)
% 
% A = 40E6.*cosd(lat);
% P = cumtrapz(lat/180*20E6,q.*A/1E15);
% 
% figure;
% plot(lat, P)


end