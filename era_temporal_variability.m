function era_temporal_variability(filename)

format compact
%pkg load netcdf
lon = ncread(filename,'longitude');
t_udiv = ncread(filename,'time');
months = ncread(filename,'months');
lat = ncread(filename,'latitude');
p = ncread(filename,'level');
prec_file = 'era_precip.nc';
P_tot = ncread(prec_file,'cp') + ncread(prec_file,'lsp');
t_prec = ncread(prec_file,'time');

info = ncinfo(filename);
num_years = info.Dimensions(4).Length;
%u_div = ncread(filename,'u_div');

aver_end_lat = 30;
averInd = -aver_end_lat <= lat & lat <= aver_end_lat;
psimax = zeros(num_years,1);
psimin = psimax;
streamfun_tot = zeros(length(lon), length(p));
coslat = cosd(lat(averInd));

for i = 1:num_years
  u_this = ncread(filename,'u_div',[1,1,1,i],[Inf,Inf,Inf,1]);
  mean_u = squeeze(mean(u_this(:,averInd,:),2));
  
 %figure;
  streamfun = cumtrapz(p,mean_u,2)* pi * 100 *6371E3 * 2 * aver_end_lat/180/9.81;
  streamfun_tot = streamfun_tot + streamfun;
%  [X,Z] = meshgrid(lon,p);
%  surf(X,Z,streamfun','edgecolor','none');
%  %contour(X,Z,streamfun','showtext','on');
%  view(2);
%  axis tight
%  colorbar
%  title('Streamfunction from udiv')
%  set(gca,'Ydir','reverse')
%  caxis([-50E9, 50E9])
  ind = 0 < lon & lon < 360;
  s = streamfun(ind,:);
  psimin(i) = min(s(:));
  psimax(i) = max(s(:));
  
  prec_ind = find(t_prec == t_udiv(i));
  Q = 2.5E6*P_tot(:,:,prec_ind)/86400*1E3;
  Q_mean = sum(Q(:,averInd).*coslat'/sum(coslat),2);
  resid_heat = Q_mean-mean(Q_mean); 
  region_lims = find(resid_heat < 0);
  [~,widest_region] = max(diff(region_lims));
  wb = region_lims(widest_region);
  eb = region_lims(widest_region+1);
  resid_heat(resid_heat < 0) = 0;
  [Q_tot(i), power(i), width(i)] = compute_domain_total_heat(resid_heat, -aver_end_lat, aver_end_lat, lon, lon(wb), lon(eb));
  Q_50_200(i) = mean(mean(Q(50<lon&lon<210,averInd)));
end

streamfun_tot = streamfun_tot / num_years;

streamfun_devsum = zeros(size(streamfun_tot));
for i = 1:num_years
  u_this = ncread(filename,'u_div',[1,1,1,i],[Inf,Inf,Inf,1]);
  mean_u = squeeze(mean(u_this(:,averInd,:),2));
  
 %figure;
  streamfun = cumtrapz(p,mean_u,2)* pi * 100 *6371E3 * 2 * aver_end_lat/180/9.81;
  streamfun_devsum = streamfun_devsum + (streamfun-streamfun_tot).^2/(num_years-1);
end
streamfun_std = sqrt(streamfun_devsum);


figure;
colormap(jet)
[X,Z] = meshgrid(lon,compute_standard_height(p));
surf(X,Z,streamfun_tot','edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
axis tight
colorbar
title('Streamfunction from udiv')
%set(gca,'Ydir','reverse')
caxis([-2.5E10, 2.5E10])
ylim([compute_standard_height(p(end)),15])
title(filename)
disp('streamfun_mean max and min')
max(streamfun_tot(:)),min(streamfun_tot(:))

figure;
colormap(jet(16))
[X,Z] = meshgrid(lon,compute_standard_height(p));
r = normcdf(abs(streamfun_tot./streamfun_std));
streamfun_tot(r<0.90) = nan;
%surf(X,Z,r,'edgecolor','none');
contourf(X,Z,streamfun_tot',16,'edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
axis tight
colorbar
title('JJA 40 years')
%set(gca,'Ydir','reverse')
caxis([-2.5E10, 2.5E10])
xlim([0,360])
ylim([compute_standard_height(p(end)),15])

v = Q_50_200';
%psimax = psimax./v; psimin = psimin./v;
figure;
plot(1:num_years, psimax, 1:num_years, -psimin)
legend('psimax','-psimin')
ylim([1E10,8E10])
title(filename)
% figure;plot((1:num_years)',[Q_tot'])
% figure;plot((1:num_years)',[power'])
% figure;plot((1:num_years)',[Q_50_200'])

figure;
plot(-psimin,psimax,'.')
corr(-psimin,psimax)

% figure;
% subplot(2,1,1)
% hist(psimax,20)
% title(['psimax', filename])
% subplot(2,1,2)
% hist(-psimin,20)
% title(['psimin', filename])

display('psimax mean std nsr:')
mean(psimax), std(psimax), std(psimax)/mean(psimax)

display('psimin mean std nsr:')
mean(psimin), std(psimin), std(psimin)/mean(psimin)

end

function [heat, power, width] = compute_domain_total_heat(Q_vert_int, slat, nlat, lon, lon1, lon2)

RE = 6371E3;
A = RE^2*(sind(nlat)-sind(slat))*(lon2-lon1)*pi/180;
heat = mean(Q_vert_int(lon1<lon&lon<lon2))*A;
power = heat/A;
width = lon2-lon1;
%display(['Heat from ', num2str(lon1), ' to ', num2str(lon2), ' (', num2str(lon2-lon1),') : ', num2str(heat/1E15),' PW (',num2str(heat/A),' W/m^2)'])
%max_Q_int = max(Q_vert_int(lon1<lon&lon<lon2));

end