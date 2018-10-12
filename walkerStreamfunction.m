function [] = walkerStreamfunction()

pkg load netcdf

lon = ncread('atmos_average.nc','lon');
lat = ncread('atmos_average.nc','lat');
load div_wind
u = read_file('ucomp');
v = read_file('vcomp');
omega = read_file('omega');
p_3d = zeros(1,1,size(u,3));
p_3d(1,1,:) = ncread('atmos_average.nc','pfull');
T = read_file('temp');
rho = p_3d./(287*T)*100;
w = -omega * 287 .* T ./ (100*p_3d*9.8);

figure;
[X,Y] = meshgrid(lon,lat);
usurf = u_div(:,:,end);
surf(X,Y,usurf','edgecolor','none');
view(2);
axis tight
colorbar
title('u div surf')
%set(gca,'Ydir','reverse')
caxis([-max(usurf(:)), max(usurf(:))])

figure;
[X,Y] = meshgrid(lon,lat);
vsurf = v_div(:,:,end);
surf(X,Y,vsurf','edgecolor','none');
view(2);
axis tight
colorbar
title('v div surf')
%set(gca,'Ydir','reverse')
caxis([-max(vsurf(:)), max(vsurf(:))])

figure;
[X,Y] = meshgrid(lon,lat); 
usurf = u(:,:,end); usurf = usurf - mean(usurf,1);
surf(X,Y,usurf','edgecolor','none');
view(2);
axis tight
colorbar
title('u surf')
%set(gca,'Ydir','reverse')
caxis([min(usurf(:)), -min(usurf(:))])

figure;
[X,Y] = meshgrid(lon,lat);
vsurf = v(:,:,end); vsurf = vsurf - mean(vsurf,1);
surf(X,Y,vsurf','edgecolor','none');
view(2);
axis tight
colorbar
title('v surf')
%set(gca,'Ydir','reverse')
caxis([-max(vsurf(:)), max(vsurf(:))])

figure;
[X,Y] = meshgrid(lon,lat); 
usurf = u(:,:,6); usurf = usurf - mean(usurf,1);
surf(X,Y,usurf','edgecolor','none');
view(2);
axis tight
colorbar
title('u near top')
%set(gca,'Ydir','reverse')
caxis([min(usurf(:)), -min(usurf(:))])

figure;
[X,Y] = meshgrid(lon,lat);
vsurf = v(:,:,6); vsurf = vsurf - mean(vsurf,1);
surf(X,Y,vsurf','edgecolor','none');
view(2);
axis tight
colorbar
title('v near top')
%set(gca,'Ydir','reverse')
caxis([-max(vsurf(:)), max(vsurf(:))])


%figure;
%[X,Y] = meshgrid(lon,lat);
%qsurf = exp(-((X-180)/40).^2 - (Y/20).^2);
%surf(X,Y,qsurf,'edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('v surf')
%set(gca,'Ydir','reverse')
%hold on;
%Z = zeros(size(X))+1;
%i = 1:4:size(X,1);
%j = 1:4:size(X,2);
%quiver3(X(i,j),Y(i,j),Z(i,j),usurf(j,i),vsurf(j,i),zeros(size(X(i,j))),'color','k')


zonal_mean_w = squeeze(mean(w,1));
p = ncread('atmos_average.nc','pfull');
figure;
subplot(2,2,1);
[X,Z] = meshgrid(lat,p);
surf(X,Z,zonal_mean_w','edgecolor','none');
view(2);
axis tight
colorbar
title('Zonal mean w')
set(gca,'Ydir','reverse')
caxis([-max(zonal_mean_w(:)), max(zonal_mean_w(:))])

zonal_mean_v = squeeze(mean(v,1));
p = ncread('atmos_average.nc','pfull');
subplot(2,2,2);
surf(X,Z,zonal_mean_v','edgecolor','none');
view(2);
axis tight
colorbar
title('Zonal mean v')
set(gca,'Ydir','reverse')
caxis([-max(zonal_mean_v(:)), max(zonal_mean_v(:))])

zonal_mean_u = squeeze(mean(u,1));
p = ncread('atmos_average.nc','pfull');
subplot(2,2,3);
surf(X,Z,zonal_mean_u','edgecolor','none');
view(2);
axis tight
colorbar
title('Zonal mean u')
set(gca,'Ydir','reverse')
caxis([-max(zonal_mean_u(:)), max(zonal_mean_u(:))])

pot_t = T.*(1E3./p_3d).^(287/1004); pot_t = pot_t - pot_t(1,1,:);
zonal_mean_T = squeeze(mean(pot_t,1));
p = ncread('atmos_average.nc','pfull');
[X,Z] = meshgrid(lat,p);
subplot(2,2,4);
surf(X,Z,zonal_mean_T','edgecolor','none');
view(2);
axis tight
colorbar
title('Zonal mean \theta - \theta_0')
set(gca,'Ydir','reverse')
%caxis([-max(zonal_mean_T(:)), max(zonal_mean_T(:))])



aver_end_lat = 12;
averInd = -aver_end_lat <= lat & lat <= aver_end_lat;
f = zeros(1,length(lat),1); f(1,:,1) = 2*7.29E-5*sind(lat);
fv = f.*v;
zonal_pert_u = u - mean(u,1);
mean_w = squeeze(mean(w(:,averInd,:),2));
mean_u = squeeze(mean(u(:,averInd,:),2));
mean_ustar = squeeze(mean(zonal_pert_u(:,averInd,:),2));
mean_theta = squeeze(mean(pot_t(:,averInd,:),2));
p_2d(1,:) = ncread('atmos_average.nc','pfull');
mean_Q = fliplr(max(-6*exp(-((lon-180)/40).^2)*((p_2d - 140).*(p_2d - 900).*(p_2d*100 - 2825.0/3.0) * 3.0/211250000.0),0)/100);
mean_udiv = squeeze(mean(u_div(:,averInd,:),2));
mean_fv = squeeze(mean(fv(:,averInd,:),2));
fv_dz = -diff(mean_fv,1,2) / (p(end)-p(end-1));

figure;
subplot(2,2,1)
[X,Z] = meshgrid(lon,p);
surf(X,Z,mean_Q','edgecolor','none');
view(2);
axis tight
colorbar
title('Q')
set(gca,'Ydir','reverse')
%caxis([min(mean_Q(:)), -min(mean_Q(:))])

subplot(2,2,2)
[X,Z] = meshgrid(lon,p);
surf(X,Z,mean_fv','edgecolor','none');
view(2);
axis tight
colorbar
title('Meridional mean fv')
set(gca,'Ydir','reverse')
caxis([min(mean_fv(:)), -min(mean_fv(:))])

subplot(2,2,4)
[X,Z] = meshgrid(lon,p);
surf(X,Z,mean_theta','edgecolor','none');
view(2);
axis tight
colorbar
title('Meridional mean \theta - \theta_0')
set(gca,'Ydir','reverse')
%caxis([-max(mean_u(:)), max(mean_u(:))])

subplot(2,2,3)
streamfun = cumtrapz(p,mean_udiv,2)*2*pi*6371E3/9.8;
[X,Z] = meshgrid(lon,p);
surf(X,Z,streamfun','edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
axis tight
colorbar
title('Streamfunction from div. wind')
set(gca,'Ydir','reverse')
caxis([min(streamfun(:)), -min(streamfun(:))])
%min(streamfun(:)),max(streamfun(:))

figure;
streamfun = cumtrapz(p,mean_u,2)*2*pi*6371E3/9.8;
[X,Z] = meshgrid(lon,p);
surf(X,Z,streamfun','edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
axis tight
colorbar
title('Streamfunction dx = 40 degrees')
set(gca,'Ydir','reverse')
caxis([min(streamfun(:)), -min(streamfun(:))])
disp('u_full')
min(streamfun(:)),max(streamfun(:))

figure;
streamfun = cumtrapz(p,mean_udiv,2)*2*pi*6371E3/9.8;
[X,Z] = meshgrid(lon,p);
surf(X,Z,streamfun','edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
axis tight
colorbar
title('Streamfunction from div. wind')
set(gca,'Ydir','reverse')
caxis([min(streamfun(:)), -min(streamfun(:))])
disp('udiv')
min(streamfun(:)),max(streamfun(:))

figure;
streamfun = cumtrapz(p,mean_ustar,2)*2*pi*6371E3/9.8;
[X,Z] = meshgrid(lon,p);
surf(X,Z,streamfun','edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
axis tight
colorbar
title('Streamfunction from u*')
set(gca,'Ydir','reverse')
caxis([min(streamfun(:)), -min(streamfun(:))])
disp('ustar')
min(streamfun(:)),max(streamfun(:))

figure;
z = 100*(p(end)-p)/0.5/9.81;
streamfun = cumtrapz(z,mean_ustar,2)*2*pi*6371E3/9.8;
[X,Z] = meshgrid(lon,z);
surf(X,Z,streamfun','edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
axis tight
colorbar
title('Streamfunction from u*')
%set(gca,'Ydir','reverse')
caxis([min(streamfun(:)), -min(streamfun(:))])
disp('ustar_z')
min(streamfun(:)),max(streamfun(:))

%figure;
%[X,Z] = meshgrid(lon,p);
%surf(X,Z,mean_udiv','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('Meridional mean u div')
%set(gca,'Ydir','reverse')
%caxis([-max(mean_udiv(:)), max(mean_udiv(:))])


%figure;
%[X,Z] = meshgrid(lon,(p(1:end-1)+p(2:end))/2);
%surf(X,Z,fv_dz','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('Meridional mean fv')
%set(gca,'Ydir','reverse')
%caxis([-max(fv_dz(:)), max(fv_dz(:))])




end

function var = read_file(varname)

var = ncread('atmos_average.nc',varname);
var = squeeze(mean(var(:,:,:,end-1:end),4));

end