function plot_era()

fname = 'era_udiv_678.nc';
u_div = ncread(fname,'u_div');
p = ncread(fname,'lev');
lon = ncread(fname,'lon');
lat = ncread(fname,'lat');
aver_lat = 20;

ind = -aver_lat <= lat & lat <= aver_lat;
u_mean = squeeze(mean(u_div(:,ind,:,:),2));
streamfun = cumtrapz(p, u_mean, 2) * pi * 100 *6371E3 * 2 * aver_lat/180/9.81;
z = compute_standard_height(p);
%streamfun = streamfun(:,z<16);
%z = z(z<=16);

display(['ERA5 Max [10^9 kg/s]: ', num2str(max(streamfun(:))/1E9)])
display(['ERA5 Min [10^9 kg/s]: ', num2str(min(streamfun(:))/1E9)])

[LON,Z] = meshgrid(lon, z);
figure;
surf(LON,Z,streamfun','edgecolor','none')
%contour(LON,Z,streamfun')
view(2)
colorbar
m = max(abs(streamfun(:)));
caxis([-m,m])
axis tight
yl = get(gca, 'ylim');
ylim([yl(1),15])
%set(gca,'ydir','reverse')

end