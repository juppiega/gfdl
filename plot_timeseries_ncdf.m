function [psimax, psimin] = plot_timeseries_ncdf (experiment_name)

format compact
%pkg load netcdf
filename = ['atmos_average.',experiment_name,'.nc'];
ncread(filename,'pfull');
lon = ncread(filename,'lon');
lat = ncread(filename,'lat');
p = ncread(filename,'pfull');

u_div = ncread(filename,'u_div');

aver_end_lat = 30;
averInd = -aver_end_lat <= lat & lat <= aver_end_lat;
psimax = zeros(size(u_div,4),1);
psimin = psimax;
streamfun_tot = zeros(length(lon), length(p));

for i = 1:size(u_div,4)
  u_this = u_div(:,:,:,i);
  mean_u = squeeze(mean(u_this(:,averInd,:),2));
  
%  figure;
  streamfun = cumtrapz(p,mean_u,2)* pi * 100 *6371E3 * 2 * aver_end_lat/180/9.81;
  streamfun_tot = streamfun_tot + streamfun;
%  [X,Z] = meshgrid(lon,p);
%  surf(X,Z,streamfun','edgecolor','none');
%  %contour(X,Z,streamfun','showtext','on');
%  view(2);
%  axis tight
%  colorbar
%  title('Streamfunction from u*')
%  set(gca,'Ydir','reverse')
%  caxis([min(streamfun(:)), -min(streamfun(:))])
  ind = 0 < lon & lon < 360;
  s = streamfun(ind,:);
  psimin(i) = min(s(:));
  psimax(i) = max(s(:));
end

streamfun_tot = streamfun_tot / size(u_div,4);
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
title(experiment_name)

c = 5;
figure;
plot(c:size(u_div,4), psimax(c:end), c:size(u_div,4), -psimin(c:end))
legend('psimax','-psimin')
ylim([1E10,6E10])
title(experiment_name)

figure;
subplot(2,1,1)
hist(psimax,20)
title(['psimax', experiment_name])
subplot(2,1,2)
hist(-psimin,20)
title(['psimin', experiment_name])

display('psimax mean std:')
mean(psimax(c:end)), std(psimax(c:end))

display('psimin mean std:')
mean(psimin(c:end)), std(psimin(c:end))

end

