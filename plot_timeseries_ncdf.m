function [psimax, psimin] = plot_timeseries_ncdf (experiment_name)

format compact
%pkg load netcdf
for i = 1:length(experiment_name)
    filename{i} = ['atmos_average.',experiment_name{i},'.nc'];
end
ncread(filename{1},'pfull');
lon = ncread(filename{1},'lon');
lat = ncread(filename{1},'lat');
p = ncread(filename{1},'pfull');

u_div = ncread(filename{1},'u_div');
for i = 2:length(experiment_name)
    u_div = cat(4, u_div, ncread(filename{i},'u_div'));
end

c = 1;
N = size(u_div,4);
aver_end_lat = 30;
averInd = -aver_end_lat <= lat & lat <= aver_end_lat;
psimax = zeros(size(u_div,4)-c+1,1);
psimin = psimax;
strfun_mean_1 = zeros(length(lon), length(p));
strfun_mean_2 = strfun_mean_1;
%T_tot = zeros(1, length(p));
half = ceil(N/2);

for i = c:N
  u_this = u_div(:,:,:,i);
  mean_u = squeeze(mean(u_this(:,averInd,:),2));
  %T = ncread(filename,'temp',[1,1,1,i],[Inf,Inf,Inf,1]);
  %T = squeeze(mean(T(:,averInd,:),2));
  %T_tot = T_tot + mean(T,1);
%  figure;
  streamfun = cumtrapz(p,mean_u,2)* pi * 100 *6371E3 * 2 * aver_end_lat/180/9.81;
  if i <= half
      strfun_mean_1 = strfun_mean_1 + streamfun;
  else
      strfun_mean_2 = strfun_mean_2 + streamfun;
  end
%  [X,Z] = meshgrid(lon,p);
%  surf(X,Z,streamfun','edgecolor','none');
%  %contour(X,Z,streamfun','showtext','on');
%  view(2);
%  axis tight
%  colorbar
%  title('Streamfunction from udiv')
%  set(gca,'Ydir','reverse')
%  caxis([-30E9,30E9])
  ind = 0 < lon & lon < 360;
  s = streamfun(ind,:);
  psimin(i-c+1) = min(s(:));
  psimax(i-c+1) = max(s(:));
end

strfun_mean_1 = strfun_mean_1 / half/1E9;
strfun_mean_2 = strfun_mean_2 / half/1E9;
%T_tot = T_tot / size(u_div,4);

strfun_std_1 = zeros(size(strfun_mean_1));
strfun_std_2 = strfun_std_1;
for i = 1:N
  u_this = u_div(:,:,:,i);
  mean_u = squeeze(mean(u_this(:,averInd,:),2));
  streamfun = cumtrapz(p,mean_u,2)* pi * 100 *6371E3 * 2 * aver_end_lat/180/9.81;
  if i <= half
      strfun_std_1 = strfun_std_1 + (streamfun - strfun_mean_1*1E9).^2/(half-1);
  else
      strfun_std_2 = strfun_std_2 + (streamfun - strfun_mean_2*1E9).^2/(half-1);
  end
end
strfun_std_1 = sqrt(strfun_std_1);
strfun_std_2 = sqrt(strfun_std_2);
  
figure;
colormap(jet(16))
z = compute_standard_height(p);%compute_actual_z(p'*100,T_tot)'/1000;
[X,Z] = meshgrid(lon,z);
contourf(X,Z,strfun_mean_1',16,'edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
set(gca,'fontsize',15)
axis tight
colorbar
%set(gca,'Ydir','reverse')
caxis([-20, 20])
xlabel('Longitude','fontsize',15)
ylabel('Altitude [km]','fontsize',15)
ylim([compute_standard_height(p(end)),15])
title('Nonlinear streamfunction [Tg/s]','fontsize',15)
disp('streamfun_mean max and min')
max(strfun_mean_1(:)),min(strfun_mean_1(:))

figure;
colormap(jet(16))
z = compute_standard_height(p);%compute_actual_z(p'*100,T_tot)'/1000;
[X,Z] = meshgrid(lon,z);
s = sqrt((strfun_std_1.^2 + strfun_std_2.^2)/half);
t = 1E9*(strfun_mean_1-strfun_mean_2)./s;
alpha = 0.05;
r = tcdf(t,2*(half-1));
%strfun_mean_1(r<alpha/2 | r>(1-alpha/2)) = nan;
contourf(X,Z,strfun_mean_1',16,'edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
set(gca,'fontsize',15)
axis tight
colorbar
ylim([compute_standard_height(p(end)),15])
xlabel('Longitude','fontsize',15)
ylabel('Altitude [km]','fontsize',15)
ylim([compute_standard_height(p(end)),15])
title('strfun_1','fontsize',15)
xlim([0,360])
%set(gca,'Ydir','reverse')
caxis([-20, 20])

figure;
colormap(jet(16))
%strfun_mean_2(r<alpha/2 | r>(1-alpha/2)) = nan;
contourf(X,Z,strfun_mean_2',16,'edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
set(gca,'fontsize',15)
axis tight
colorbar
ylim([compute_standard_height(p(end)),15])
xlabel('Longitude','fontsize',15)
ylabel('Altitude [km]','fontsize',15)
ylim([compute_standard_height(p(end)),15])
title('strfun_2','fontsize',15)
xlim([0,360])
%set(gca,'Ydir','reverse')
caxis([-20, 20])

err = abs((strfun_mean_1./strfun_mean_2)-1);
alpha = 0.1;
strfun_mean_1(err>alpha) = nan;
ratio_non_nan = sum(~isnan(strfun_mean_1(:)))/numel(strfun_mean_1)
figure;
colormap(jet(16))
contourf(X,Z,strfun_mean_1',16,'edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
set(gca,'fontsize',15)
axis tight
colorbar
xlabel('Longitude','fontsize',15)
ylabel('Altitude [km]','fontsize',15)
ylim([compute_standard_height(p(end)),15])
title('error','fontsize',15)
xlim([0,360])
%set(gca,'Ydir','reverse')
caxis([-20, 20])

% figure;
% colormap(jet)
% z = compute_standard_height(p);%compute_actual_z(p'*100,T_tot)'/1000;
% [X,Z] = meshgrid(lon,z);
% surf(X,Z,r','edgecolor','none');
% %contour(X,Z,streamfun','showtext','on');
% view(2);
% set(gca,'fontsize',15)
% axis tight
% colorbar
% xlabel('Longitude','fontsize',15)
% ylabel('Altitude [km]','fontsize',15)
% ylim([compute_standard_height(p(end)),15])
% title('tcdf(t)','fontsize',15)
% xlim([0,360])
% %set(gca,'Ydir','reverse')
% caxis([0,1])
% 
% figure;
% colormap(jet)
% z = compute_standard_height(p);%compute_actual_z(p'*100,T_tot)'/1000;
% [X,Z] = meshgrid(lon,z);
% surf(X,Z,t','edgecolor','none');
% %contour(X,Z,streamfun','showtext','on');
% view(2);
% set(gca,'fontsize',15)
% axis tight
% colorbar
% xlabel('Longitude','fontsize',15)
% ylabel('Altitude [km]','fontsize',15)
% ylim([compute_standard_height(p(end)),15])
% title('t-stat','fontsize',15)
% xlim([0,360])
% %set(gca,'Ydir','reverse')
% %caxis([0,1])

% figure;
% plot(c:size(u_div,4), psimax, c:size(u_div,4), -psimin)
% legend('psimax','-psimin')
% ylim([1E10,8E10])
% title(experiment_name)
% 
% figure;
% subplot(2,1,1)
% hist(psimax,20)
% title(['psimax', experiment_name])
% subplot(2,1,2)
% hist(-psimin,20)
% title(['psimin', experiment_name])

% display('psimax mean std:')
% mean(psimax(c:end)), std(psimax(c:end))
% 
% display('psimin mean std:')
% mean(psimin(c:end)), std(psimin(c:end))

end

