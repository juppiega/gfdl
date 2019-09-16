function plot_era()

fname = 'era_udiv_1.nc';
u_div = ncread(fname,'u_div');
%div = ncread(fname,'div');
p = ncread(fname,'lev');
lon = ncread(fname,'lon');
lat = ncread(fname,'lat');
aver_lat = 15;

ind = -aver_lat <= lat & lat <= aver_lat;
coslat = zeros([1, sum(ind), 1]);
coslat(:) = cosd(lat(ind));
udiv_mean = squeeze(sum(u_div(:,ind,:,1).*coslat/sum(coslat),2));
%div_mean = squeeze(mean(div(:,ind,:,:),2));
streamfun = cumtrapz(p, udiv_mean, 2) * pi * 100 *6371E3 * 2 * aver_lat/180/9.81;
z = compute_standard_height(p);
%streamfun = streamfun(:,z<16);
%z = z(z<=16);

strfun_west_280 = streamfun(lon < 280, :);
display(['ERA5 Max west of 280 [10^9 kg/s]: ', num2str(max(strfun_west_280(:))/1E9)])
strfun_west_200 = streamfun(lon < 200, :);
display(['ERA5 Min west of 200 [10^9 kg/s]: ', num2str(min(strfun_west_200(:))/1E9)])

% [LON,Z] = meshgrid(lon, z);
% figure;
% colormap(jet)
% contourf(LON,Z,ustar_mean',20,'edgecolor','none')
% %contour(LON,Z,streamfun')
% view(2)
% colorbar
% axis tight
% yl = get(gca, 'ylim');
% ylim([yl(1),15])
% title('u_star')

[LON,Z] = meshgrid(lon, z);

% figure;
% colormap(jet)
% contourf(LON,Z,div_mean',20,'edgecolor','none')
% %contour(LON,Z,streamfun')
% view(2)
% colorbar
% axis tight
% m = 1E-5;
% caxis([-m,m])
% yl = get(gca, 'ylim');
% ylim([yl(1),15])
% title('divergence')

figure;
colormap(jet(16))
contourf(LON,Z,streamfun'/1E9,16,'edgecolor','none')
%contour(LON,Z,streamfun')
view(2)
colorbar('location','northoutside')
m = 20;
caxis([-m,m])
axis tight
yl = get(gca, 'ylim');
ylim([yl(1),15])
fs = 15;
title(['d) La Ni',char(0241),'a streamfunction [Tg/s]'],'fontsize',fs)
xlabel('Longitude','fontsize',fs)
ylabel('Altitude [km]','fontsize',fs)
set(gca,'fontsize',fs)
%set(gca,'ydir','reverse')

end

function deriv = Dx(X_in, lon)

dx = (lon(2)-lon(1))/360 * 2*pi*6371E3;
X = zeros(size(X_in,1)+2,size(X_in,2));
X(1,:) = X_in(end,:);
X(end,:) = X_in(1,:);
X(2:end-1,:) = X_in;
deriv = (X(3:end,:) - X(1:end-2,:))/(2*dx);

end

function deriv = Dy(X, lat)

dy = (lat(1)-lat(2))/180 * pi * 6731E3;
deriv = zeros(size(X));
deriv(:,2:end-1,:) = (X(:,1:end-2,:) - X(:,3:end,:)) / (2*dy);

end

