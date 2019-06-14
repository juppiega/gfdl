function [] = walkerStreamfunction(experiment_name)

filename = ['atmos_average.',experiment_name,'.nc'];
div_wind_name = ['div_wind.',experiment_name,'.mat'];
if ~exist(filename, 'file')
    error('File not found!')
end
if ~exist(div_wind_name, 'file')
    display('Solving divergent wind')
    solveDivergentWind(experiment_name);
end

pkg load netcdf

lon = ncread(filename,'lon');
lat = ncread(filename,'lat');
load(div_wind_name)
u = read_file(filename,'ucomp');
v = read_file(filename,'vcomp');
omega = read_file(filename,'omega');
p_3d = zeros(1,1,size(u,3));
p_3d(1,1,:) = ncread(filename,'pfull');
T = read_file(filename,'temp');
rho = p_3d./(287*T)*100;
w = -omega * 287 .* T ./ (100*p_3d*9.8);
integrand = 287*T./(100*p_3d*9.8);
geopot = cumtrapz(p_3d(:)*100, -integrand, 3);
geopot = geopot - geopot(:,:,end);
z_mid = 0.5*(geopot(:,:,1:end-1) + geopot(:,:,2:end));
lapse_rate = -1e3*(T(:,:,1:end-1) - T(:,:,2:end)) ./ (geopot(:,:,1:end-1) - geopot(:,:,2:end));
unstable = lapse_rate > (9.8/1.004);

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
%
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
vsurf = v(:,:,end); %vsurf = vsurf - mean(vsurf,1);
surf(X,Y,vsurf','edgecolor','none');
view(2);
axis tight
colorbar
title('v surf')
%set(gca,'Ydir','reverse')
caxis([-max(vsurf(:)), max(vsurf(:))])
%
%figure;
%[X,Y] = meshgrid(lon,lat); 
%usurf = u(:,:,6); usurf = usurf - mean(usurf,1);
%surf(X,Y,usurf','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('u near top')
%%set(gca,'Ydir','reverse')
%caxis([min(usurf(:)), -min(usurf(:))])
%
%figure;
%[X,Y] = meshgrid(lon,lat);
%vsurf = v(:,:,6); vsurf = vsurf - mean(vsurf,1);
%surf(X,Y,vsurf','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('v near top')
%%set(gca,'Ydir','reverse')
%caxis([-max(vsurf(:)), max(vsurf(:))])


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
p = ncread(filename,'pfull');
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
p = ncread(filename,'pfull');
subplot(2,2,2);
surf(X,Z,zonal_mean_v','edgecolor','none');
view(2);
axis tight
colorbar
title('Zonal mean v')
set(gca,'Ydir','reverse')
caxis([-max(zonal_mean_v(:)), max(zonal_mean_v(:))])

zonal_mean_u = squeeze(mean(u,1));
p = ncread(filename,'pfull');
subplot(2,2,3);
surf(X,Z,zonal_mean_u','edgecolor','none');
view(2);
axis tight
colorbar
title('Zonal mean u')
set(gca,'Ydir','reverse')
caxis([-max(zonal_mean_u(:)), max(zonal_mean_u(:))])

pot_t = T.*(1E3./p_3d).^(287/1004); 
pot_t_0 = zeros(1,1,size(pot_t,3));
pot_t_0(:) = mean(mean(pot_t(:,:,:)));
pot_t = pot_t - pot_t_0;
zonal_mean_T = squeeze(mean(pot_t,1));
p = ncread(filename,'pfull');
[X,Z] = meshgrid(lat,p);
subplot(2,2,4);
surf(X,Z,zonal_mean_T','edgecolor','none');
view(2);
axis tight
colorbar
title('Zonal mean \theta - \theta_0')
set(gca,'Ydir','reverse')
%caxis([-max(zonal_mean_T(:)), max(zonal_mean_T(:))])



aver_end_lat = 20;
averInd = -aver_end_lat <= lat & lat <= aver_end_lat;
f = zeros(1,length(lat),1); f(1,:,1) = 2*7.29E-5*sind(lat);
fv = f.*v;
zonal_pert_u = u - mean(u,1);
mean_w = squeeze(mean(w(:,averInd,:),2));
mean_u = squeeze(mean(u(:,averInd,:),2));
mean_ustar = squeeze(mean(zonal_pert_u(:,averInd,:),2));
pot_t_zon = pot_t - mean(pot_t, 1);
mean_theta = squeeze(mean(pot_t_zon(:,averInd,:),2));
mean_lapse = squeeze(mean(lapse_rate(:,averInd,:),2));
mean_zmid = squeeze(mean(z_mid(:,averInd,:),2));
any_unstable = squeeze(any(unstable(:,averInd,:),2));
p_2d(1,:) = ncread(filename,'pfull');
mean_Q = fliplr(max(6*(1-((lon-180)/40).^2 - ((p_2d-500e2)/350e2).^2),0));
mean_udiv = squeeze(mean(u_div(:,averInd,:),2));
mean_fv = squeeze(mean(fv(:,averInd,:),2));
fv_dz = -diff(mean_fv,1,2) / (p(end)-p(end-1));

figure;
[X,Z] = meshgrid(lon,0.5*(p(1:end-1) + p(2:end)));
surf(X,Z,double(any_unstable)','edgecolor','none');
view(2);
axis tight
colorbar
title('Unstable regions present')
set(gca,'Ydir','reverse')

figure;
%subplot(2,2,2)
%[X,Z] = meshgrid(lon,0.5*(p(1:end-1) + p(2:end)));
%surf(X,Z,mean_lapse','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('Lapse rate')
%set(gca,'Ydir','reverse')
%caxis([0, 9.8/1.004])
[X,Z] = meshgrid(lon,p);
surf(X,Z,mean_fv','edgecolor','none');
view(2);
axis tight
h=colorbar;
set(h,'fontsize',15)
title('cov(fv)','fontsize',15)
set(gca,'Ydir','reverse')
set(gca,'fontsize',15)
xlabel('longitude','fontsize',15)
ylabel('p [hPa]','fontsize',15)
%caxis([min(mean_udiv(:)), -min(mean_udiv(:))])

figure
%subplot(2,2,1)
[X,Z] = meshgrid(lon,p);
surf(X,Z,mean_ustar','edgecolor','none');
view(2);
axis tight
h=colorbar;
set(h,'fontsize',15)
title('Meridional mean u*','fontsize',15)
set(gca,'fontsize',15)
xlabel('longitude','fontsize',15)
ylabel('p [hPa]','fontsize',15)
set(gca,'Ydir','reverse')
caxis([min(mean_ustar(:)), -min(mean_ustar(:))])

%subplot(2,2,4)
%[X,Z] = meshgrid(lon,p);
%surf(X,Z,mean_theta','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('Meridional mean \theta - \theta_0')
%set(gca,'Ydir','reverse')
%%caxis([-max(mean_u(:)), max(mean_u(:))])
%
%subplot(2,2,3)
%streamfun = cumtrapz(p,mean_udiv,2)*2*pi*6371E3/9.8;
%[X,Z] = meshgrid(lon,p);
%surf(X,Z,streamfun','edgecolor','none');
%%contour(X,Z,streamfun','showtext','on');
%view(2);
%axis tight
%colorbar
%title('Streamfunction from div. wind')
%set(gca,'Ydir','reverse')
%caxis([min(streamfun(:)), -min(streamfun(:))])
%%min(streamfun(:)),max(streamfun(:))

%figure;
%[X,Z] = meshgrid(lon,p);
%surf(X,Z,mean_fv','edgecolor','none');
%view(2);
%axis tight
%colorbar('fontsize',15)
%shading interp
%ylabel('p [hPa]','fontsize',15)
%xlabel('Longitude','fontsize',15)
%title('Meridional mean fv [m/s^2]','fontsize',15)
%set(gca,'Ydir','reverse','fontsize',15)
%caxis([min(mean_fv(:)), -min(mean_fv(:))])

%figure;
%streamfun = cumtrapz(p,mean_u,2)*2*pi*6371E3/9.8;
%[X,Z] = meshgrid(lon,p);
%surf(X,Z,streamfun','edgecolor','none');
%%contour(X,Z,streamfun','showtext','on');
%view(2);
%axis tight
%colorbar
%title('Streamfunction dx = 40 degrees')
%set(gca,'Ydir','reverse')
%caxis([min(streamfun(:)), -min(streamfun(:))])
%disp('u_full')
%min(streamfun(:)),max(streamfun(:))

figure;
streamfun = cumtrapz(p,mean_udiv,2)*2*pi*6371E3/9.8;
[X,Z] = meshgrid(lon,p);
surf(X,Z,streamfun','edgecolor','none');
%contour(X,Z,streamfun','showtext','on');
view(2);
axis tight
colorbar('fontsize',15)
title('Streamfunction from div. wind [kg/s]','fontsize',15)
set(gca,'Ydir','reverse')
set(gca,'fontsize',15)
caxis([min(streamfun(:)), -min(streamfun(:))])
xlabel('Longitude','fontsize',15)
disp('udiv')
min(streamfun(:)),max(streamfun(:))

%figure;
%streamfun = cumtrapz(p,mean_ustar,2)*2*pi*6371E3/9.8;
%[X,Z] = meshgrid(lon,p);
%surf(X,Z,streamfun','edgecolor','none');
%%contour(X,Z,streamfun','showtext','on');
%view(2);
%axis tight
%colorbar
%title('Streamfunction from u*')
%set(gca,'Ydir','reverse')
%caxis([min(streamfun(:)), -min(streamfun(:))])
%disp('ustar')
%min(streamfun(:)),max(streamfun(:))

%figure;
%z = 100*(p(end)-p)/0.5/9.81;
%streamfun = cumtrapz(z,mean_ustar,2);
%[X,Z] = meshgrid(lon,z);
%surf(X,Z,streamfun','edgecolor','none');
%%contour(X,Z,streamfun','showtext','on');
%view(2);
%axis tight
%colorbar
%title('Streamfunction from u* [m2/s]')
%%set(gca,'Ydir','reverse')
%caxis([min(streamfun(:)), -min(streamfun(:))])
%disp('ustar [m2/s]')
%min(streamfun(:)),max(streamfun(:))
%
%figure;
%streamfun = cumtrapz(z,mean_udiv,2);
%[X,Z] = meshgrid(lon,z);
%surf(X,Z,streamfun','edgecolor','none');
%%contour(X,Z,streamfun','showtext','on');
%view(2);
%axis tight
%colorbar
%title('Streamfunction from div. wind')
%set(gca,'Ydir','reverse')
%caxis([min(streamfun(:)), -min(streamfun(:))])
%disp('udiv [m2/s]')
%min(streamfun(:)),max(streamfun(:))

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

function var = read_file(filename,varname)

var = ncread(filename,varname);
var = squeeze(mean(var(:,:,:,end-3:end),4));

end