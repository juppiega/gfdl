function [] = walkerStreamfunction(experiment_name)

latb = 20;

filename = ['atmos_average.',experiment_name,'.nc'];
div_wind_name = ['div_wind.',experiment_name,'.mat'];
if ~exist(filename, 'file')
    error('File not found!')
end

try
  u_div = read_file(filename,'u_div');
catch
  if ~exist(div_wind_name, 'file')
    display('Solving divergent wind')
    solveDivergentWind(experiment_name);
  end
  load(div_wind_name)
end

pkg load netcdf

lon = ncread(filename,'lon');
lat = ncread(filename,'lat');
u = read_file(filename,'ucomp');
v = read_file(filename,'vcomp');
omega = read_file(filename,'omega');
p_3d = zeros(1,1,size(u,3));
p_3d(1,1,:) = ncread(filename,'pfull');
p = ncread(filename,'pfull');
T = read_file(filename,'temp');
newton_heat = compute_newton_heat(T, lat, lon, p);
Q_prescribed = compute_Q_prescribed(experiment_name, lat, lon, p);
try
sst_mean = compute_sst_mean(filename);
catch
end

try
sst = read_2d(filename,'t_surf');
lsp = read_2d(filename,'condensation_rain');
q_flux = read_2d(filename,'flux_oceanq');
cp = read_2d(filename,'convection_rain');
lh_conv = read_file(filename, 'dt_qg_convection')*2.5E6/1005*86400;
lh_ls = read_file(filename, 'dt_qg_condensation')*2.5E6/1005*86400;
rad_cool = read_file(filename, 'tdt_rad') * 86400;
catch
zer = zeros(size(u,1),size(u,2));
sst = zer; lsp = zer; q_flux = zer; cp = zer; lh_conv = zer; lh_ls=zer; rad_cool = zer; 
end
precip = 2.5E6*(lsp+cp);

try
prec_hovm = compute_hovmoller(filename, 'convection_rain', lat, latb)*2.5E6;
catch
end
%Teq = read_file(filename,'teq');
try
T_adv = read_file(filename,'t_hor_adv')*86400;
T_newton = read_file(filename,'tdt_ndamp')*86400;
catch
T_adv = zeros(size(T));
T_newton = zeros(size(T));
end

rho = p_3d./(287*T)*100;
w = -omega * 287 .* T ./ (100*p_3d*9.8);
integrand = 287*T./(100*p_3d*9.8);
geopot = cumtrapz(p_3d(:)*100, -integrand, 3);
geopot = geopot - geopot(:,:,end);
z_mid = 0.5*(geopot(:,:,1:end-1) + geopot(:,:,2:end));
lapse_rate = -1e3*(T(:,:,1:end-1) - T(:,:,2:end)) ./ (geopot(:,:,1:end-1) - geopot(:,:,2:end));
unstable = lapse_rate > (9.8/1.004);
Sp = 0.5*(integrand(:,:,1:end-1) + integrand(:,:,2:end)) .* (9.8/1004 - lapse_rate/1000);

try
Q_heat =  Q_prescribed;% + newton_heat;
Q_vert_int = trapz(p(p>150)*100/9.81, 1004*Q_heat(:,:,p>150),3);
%Q_heat = -1004*0.5*(omega(:,:,1:end-1) + omega(:,:,2:end)) .* Sp;
%phalf = 0.5*(p(1:end-1)+p(2:end));
%Q_vert_int = trapz(phalf(phalf>200)*100/9.81, Q_heat(:,:,phalf>200),3);
catch
Q_heat = zeros(size(T));
end

if exist('sst_mean','var')
figure;
plot(1:length(sst_mean), sst_mean-273.15);
ylabel('SST mean')
end

%figure;
%[X,Y] = meshgrid(lon,lat);
%usurf = u_div(:,:,end);
%surf(X,Y,usurf','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('u div surf')
%%set(gca,'Ydir','reverse')
%caxis([-max(usurf(:)), max(usurf(:))])

%figure;
%[X,Y] = meshgrid(lon,lat);
%vsurf = v_div(:,:,end);
%surf(X,Y,vsurf','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('v div surf')
%%set(gca,'Ydir','reverse')
%caxis([-max(vsurf(:)), max(vsurf(:))])
%%
%figure;
%[X,Y] = meshgrid(lon,lat); 
%usurf = u(:,:,end); %usurf = usurf - mean(usurf,1);
%surf(X,Y,usurf','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('u surf')
%%set(gca,'Ydir','reverse')
%caxis([min(usurf(:)), -min(usurf(:))])
%%
%figure;
%[X,Y] = meshgrid(lon,lat);
%vsurf = v(:,:,end); %vsurf = vsurf - mean(vsurf,1);
%surf(X,Y,vsurf','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('v surf')
%%set(gca,'Ydir','reverse')
%caxis([-max(vsurf(:)), max(vsurf(:))])
%%
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


figure;
ind = -latb<=lat&lat<=latb;
[X,Y] = meshgrid(lon, lat(ind));
surf(X,Y,sst(:,ind)','edgecolor','none');
view(2)
axis tight
colorbar
title('SST')

figure;
subplot(3,1,1)
[~,i] = min(abs(p-200));
surf(X,Y,u(:,ind,i)','edgecolor','none');
view(2)
axis tight
colorbar
title('u 200')
caxis([-10,40])

subplot(3,1,2)
[~,i] = min(abs(p-500));
surf(X,Y,u(:,ind,i)','edgecolor','none');
view(2)
axis tight
colorbar
title('u 500')
caxis([-10,30])

subplot(3,1,3)
[~,i] = min(abs(p-850));
surf(X,Y,u(:,ind,i)','edgecolor','none');
view(2)
axis tight
colorbar
title('u 850')
caxis([-10,10])

figure;
surf(X,Y,precip(:,ind)','edgecolor','none');
view(2)
axis tight
colorbar
title('Precipitation [W/m^2]')
caxis([0,400])

figure;
surf(X,Y,Q_vert_int(:,ind)','edgecolor','none');
view(2)
axis tight
colorbar
title('Q total [W/m^2]')
caxis([-300,300])

figure;
clat = cosd(lat)';
[X,Y] = meshgrid(lon(1:end-1), lat(ind));
u_x = (u_div(2:end,:,:)-u_div(1:end-1,:,:))/(lon(2)-lon(1))*180/pi/6371E3./clat;
omega_w = -cumtrapz(p, u_x,3);
[~,i] = min(abs(p-500));
surf(X,Y,omega_w(:,ind,i)','edgecolor','none');
view(2)
axis tight
colorbar
title('Omega Walker 500')
%caxis([0,400])

figure;
[X,Y] = meshgrid(lon, lat);
surf(X,Y,q_flux','edgecolor','none');
view(2)
axis tight
colorbar
title('q-flux')

mean_precip = mean(precip(:,ind),2);
figure;
plot(lon,mean_precip)
m = mean(mean_precip);
hold all
plot([lon(1), lon(end)],[m,m])
ylabel('Precip [W/m^2]')
ylim([0,400])

if mean(mean_precip)>0
display('Precipitation anomaly:')
resid_heat = mean_precip-m; resid_heat(resid_heat < 0) = 0;
bz = find((mean_precip(1:end-1)-m) .* (mean_precip(2:end) - m) < 0);
maxInt = 0;
md = 0;
for i = 1:size(bz)-1
  int_this = trapz(bz(i):bz(i+1), resid_heat(bz(i):bz(i+1))');
  if int_this > maxInt
    maxInt = int_this;
    md = i;
  end
end
wb = bz(md);
eb = bz(md+1);
%figure;
%plot(lon, resid_heat,'.')
[heat, power, width] = compute_domain_total_heat(resid_heat, -latb, latb, lon, lon(wb), lon(eb));

end

display('Q vertical integral:')
mean_Q_int = mean(Q_vert_int(:,-latb<=lat & lat<=latb),2);
m = mean(mean_Q_int);
resid_heat = mean_Q_int-m; resid_heat(resid_heat < 0) = 0;
bz = find((mean_Q_int(1:end-1)-m) .* (mean_Q_int(2:end) - m) < 0);
maxInt = 0;
md = 0;
for i = 1:size(bz)-1
  int_this = trapz(bz(i):bz(i+1), resid_heat(bz(i):bz(i+1))');
  if int_this > maxInt
    maxInt = int_this;
    md = i;
  end
end
wb = bz(md);
eb = bz(md+1);
%figure;
%plot(lon, resid_heat,'.')
[heat, power, width] = compute_domain_total_heat(resid_heat, -latb, latb, lon, lon(wb), lon(eb));

save(['mean_Q_int_',experiment_name,'.mat'],'mean_Q_int','lon','-mat7-binary');
figure;
plot(lon,mean_Q_int)
hold all
plot([lon(1), lon(end)],[m,m])
ylabel('Q vertical integral [W/m^2]')
ylim([-300,300])

if exist('prec_hovm','var')
figure;
[X,t] = meshgrid(lon, 1:size(prec_hovm,2));
surf(X, t, prec_hovm', 'edgecolor','none');
view(2)
axis tight
colorbar
title('Precipitation hovmoller')
end


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
%pot_t = pot_t - pot_t_0;
zonal_mean_T = squeeze(mean(pot_t,1));
p = ncread(filename,'pfull');
[X,Z] = meshgrid(lat,p);
subplot(2,2,4);
surf(X,Z,zonal_mean_T','edgecolor','none');
view(2);
axis tight
colorbar
title('Zonal mean \theta')
caxis([270,370])
set(gca,'Ydir','reverse')
%caxis([-max(zonal_mean_T(:)), max(zonal_mean_T(:))])



aver_end_lat = latb;
averInd = -aver_end_lat <= lat & lat <= aver_end_lat;
f = zeros(1,length(lat),1); f(1,:,1) = 2*7.29E-5*sind(lat);
fv = f.*v;
zonal_pert_u = u - mean(u,1);
mean_w = squeeze(mean(omega(:,averInd,:),2));
mean_Sp = squeeze(mean(Sp(:,averInd,:),2));
mean_T = squeeze(mean(T(:,averInd,:),2));
mean_Tstar = mean_T - mean(mean_T,1);
%mean_Teq = squeeze(mean(Teq(:,averInd,:),2));
mean_T_adv = squeeze(mean(T_adv(:,averInd,:),2));
mean_T_newton = squeeze(mean(T_newton(:,averInd,:),2));
mean_Q_heat = squeeze(mean(Q_heat(:,averInd,:),2));
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



%figure;
%[X,Z] = meshgrid(lon,0.5*(p(1:end-1) + p(2:end)));
%surf(X,Z,double(any_unstable)','edgecolor','none');
%view(2);
%axis tight
%colorbar
%title('Unstable regions present')
%set(gca,'Ydir','reverse')

T_surf = mean(T(:,:,end),1);
ind = -30<lat&lat<30;
figure
plot(lat(ind), T_surf(ind)-273.15)
title('T surface')

figure;
%subplot(2,2,2)
[X,Z] = meshgrid(lon,0.5*(p(1:end-1) + p(2:end)));
surf(X,Z,mean_Sp','edgecolor','none');
view(2);
axis tight
colorbar
title('S_p')
set(gca,'Ydir','reverse')
caxis([0, 1.5E-3])

Sp_zonal_pert = (mean_Sp ./ mean(mean_Sp,1))-1;
figure;
surf(X,Z,Sp_zonal_pert','edgecolor','none');
view(2);
axis tight
colorbar
title('S_p* [relative]')
set(gca,'Ydir','reverse')
caxis([-1, 1])

ind = 0 < lon & lon < 361;
S_prof = mean(mean_Sp(ind,:),1);
figure
plot(S_prof,0.5*(p(1:end-1) + p(2:end)))
title('S_p')
set(gca,'Ydir','reverse')
xlim([0, 1.5E-3])
ylim([150,1000])

%T_prof = mean(mean_T(ind,:),1);
%z = compute_actual_z(p'*100, T_prof);
%es = exp(77.3450 + 0.0057*T_prof - 7235 ./ T_prof) ./ T_prof.^8.2;
%q_sat = 0.622*es./(p'*100);
%mse = T_prof + 9.81*z/1004 + 2.5E6*q_sat/1004;
%figure
%plot(T_prof, p)
%title('T MC/WP')
%set(gca,'Ydir','reverse')
%ylim([150,1000])
%%xlim([300, 400])

%[X,Z] = meshgrid(lon,p);
%surf(X,Z,mean_fv','edgecolor','none');
%view(2);
%axis tight
%h=colorbar;
%set(h,'fontsize',15)
%title('cov(fv)','fontsize',15)
%set(gca,'Ydir','reverse')
%set(gca,'fontsize',15)
%xlabel('longitude','fontsize',15)
%ylabel('p [hPa]','fontsize',15)
%%caxis([min(mean_udiv(:)), -min(mean_udiv(:))])

%%[X,Z] = meshgrid(lon,p);
%%surf(X,Z,mean_fv','edgecolor','none');
%%view(2);
%%axis tight
%%h=colorbar;
%%set(h,'fontsize',15)
%%title('cov(fv)','fontsize',15)
%%set(gca,'Ydir','reverse')
%%set(gca,'fontsize',15)
%%xlabel('longitude','fontsize',15)
%%ylabel('p [hPa]','fontsize',15)
%%%caxis([min(mean_udiv(:)), -min(mean_udiv(:))])

figure
%subplot(2,2,1)
[X,Z] = meshgrid(lon,p);
surf(X,Z,mean_u','edgecolor','none');
view(2);
axis tight
h=colorbar;
set(h,'fontsize',15)
title('Meridional mean u','fontsize',15)
set(gca,'fontsize',15)
xlabel('longitude','fontsize',15)
ylabel('p [hPa]','fontsize',15)
set(gca,'Ydir','reverse')
caxis([-40,40])

%figure
%%subplot(2,2,1)
%[X,Z] = meshgrid(lon,p);
%surf(X,Z,mean_ustar','edgecolor','none');
%view(2);
%axis tight
%h=colorbar;
%set(h,'fontsize',15)
%title('Meridional mean u*','fontsize',15)
%set(gca,'fontsize',15)
%xlabel('longitude','fontsize',15)
%ylabel('p [hPa]','fontsize',15)
%set(gca,'Ydir','reverse')
%caxis([min(mean_ustar(:)), -min(mean_ustar(:))])

figure
%subplot(2,2,1)
[X,Z] = meshgrid(lon,p);
surf(X,Z,mean_Q_heat','edgecolor','none');
view(2);
axis tight
h=colorbar;
set(h,'fontsize',15)
title('Meridional mean Q [K/d]','fontsize',15)
set(gca,'fontsize',15)
xlabel('longitude','fontsize',15)
ylabel('p [hPa]','fontsize',15)
m = max(max(abs(mean_Q_heat(:,p>200))));
caxis([-m,m])
set(gca,'Ydir','reverse')



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
%surf(X,Z,mean_T'-mean_Teq','edgecolor','none');
%view(2);
%axis tight
%colorbar('fontsize',15)
%shading interp
%ylabel('p [hPa]','fontsize',15)
%xlabel('Longitude','fontsize',15)
%title('Meridional mean T-Teq [K]','fontsize',15)
%set(gca,'Ydir','reverse','fontsize',15)
%%caxis([min(mean_fv(:)), -min(mean_fv(:))])

%figure;
%subplot(1,2,1)
%[X,Z] = meshgrid(lon,p);
%surf(X,Z,mean_T_adv','edgecolor','none');
%view(2);
%axis tight
%colorbar('fontsize',15)
%shading interp
%ylabel('p [hPa]','fontsize',15)
%xlabel('Longitude','fontsize',15)
%title('Meridional mean dT_hor_adv [K/day]','fontsize',15)
%set(gca,'Ydir','reverse','fontsize',15)
%%caxis([min(mean_fv(:)), -min(mean_fv(:))])
%clim = [-0.25,0.25];
%caxis(clim)

%subplot(1,2,2)
%[X,Z] = meshgrid(lon,p);
%surf(X,Z,mean_Tstar','edgecolor','none');
%view(2);
%axis tight
%colorbar('fontsize',15)
%shading interp
%ylabel('p [hPa]','fontsize',15)
%xlabel('Longitude','fontsize',15)
%title('Meridional mean T* [K]','fontsize',15)
%set(gca,'Ydir','reverse','fontsize',15)
%caxis([min(mean_fv(:)), -min(mean_fv(:))])
%caxis(clim)

%ind = 200 < p & p < 700;
%mTstar = mean_Tstar(:,ind);
%mTadv = mean_T_adv(:,ind);
%figure;
%plot(mTstar(:), mTadv(:), '.')
%xlabel('Tstar')
%ylabel('Tadv')
%%corr(mTstar(:), mTadv(:))
%pol = polyfit(mTstar(:), mTadv(:), 1)
%x = get(gca,'xlim');
%y = pol(1)*x + pol(2);
%hold all
%plot(x,y)


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
streamfun = cumtrapz(p,mean_udiv,2) * pi * 100 *6371E3 * 2 * aver_end_lat/180/9.81;
[X,Z] = meshgrid(lon,p);
contourf(X,Z,streamfun',-2E10:2E9:2E10);
%contour(X,Z,streamfun','showtext','on');
view(2);
axis tight
colorbar('fontsize',15)
title('Streamfunction from div. wind [kg/s]','fontsize',15)
set(gca,'Ydir','reverse')
set(gca,'fontsize',15)
caxis([-2E10,2E10])
xlabel('Longitude','fontsize',15)
disp('udiv')
min(streamfun(:)),max(streamfun(:))

dlon = lon(2)-lon(1);
mean_w = -(streamfun(2:end,:) - streamfun(1:end-1,:))/(dlon/360*40E6)/(2*latb/180*20E6);
figure
[X,Z] = meshgrid(lon(1:end-1),p);
surf(X,Z,mean_w','edgecolor','none');
view(2);
axis tight
colorbar
title('Meridional mean omega')
set(gca,'Ydir','reverse')
caxis([-2E-3, 2E-3])

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
%ind = 100 < lon & lon < 260;
%s = streamfun(ind,:);
%min(s(:)),max(s(:))

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
if size(var,4) == 2
  var = squeeze(var(:,:,:,end));
  return
end

try
var = squeeze(mean(var(:,:,:,20:end),4));
catch
var = squeeze(mean(var(:,:,:,end-2:end),4));
end

end

function var = read_2d(filename, varname)

var = ncread(filename,varname);
var = squeeze(mean(var(:,:,20:end),3));

end

function sst_mean = compute_sst_mean(filename)

var = ncread(filename,'temp');
sst_mean = squeeze(mean(mean(var(:,:,end,:),1),2));

end

function hovm = compute_hovmoller(filename, varname, lat, latb)

var = ncread(filename,varname);
hovm = zeros(size(var,1), size(var,3));
num_t = size(var,3);
ind = -latb <= lat & lat <= latb;
for i = 1:num_t
  varmean = mean(var(:,ind,i),2);
  hovm(:,i) = varmean;
end

end

function [heat, power, width] = compute_domain_total_heat(Q_vert_int, slat, nlat, lon, lon1, lon2)

RE = 6371E3;
A = RE^2*(sind(nlat)-sind(slat))*(lon2-lon1)*pi/180;
heat = mean(Q_vert_int(lon1<lon&lon<lon2))*A;
power = heat/A;
width = lon2-lon1;
display(['Heat from ', num2str(lon1), ' to ', num2str(lon2), ' (', num2str(lon2-lon1),') : ', num2str(heat/1E15),' PW (',num2str(heat/A),' W/m^2)'])
max_Q_int = max(Q_vert_int(lon1<lon&lon<lon2));

end

function newton_heat = compute_newton_heat(T, lat, lon, p)

p = p';
Tstrat = 200; T0 = 315; dTy = 60; dThez = 10; ka = 1/(8*86400); ks = 1/(4*86400); sig_b = 0.7; p0 = 1000;
sig = p/p0;
k = ka + (ks-ka).*max(0,(sig-sig_b)./(1-sig_b)).*cosd(lat).^4;

Teq = zeros(length(lat), length(p));
ft = (T0 - dTy.*sind(lat).^2 -dThez*cosd(lat).^2.*log(sig)).*(sig).^(2.0/7);
Teq = max(Tstrat, ft);
Teq_3d = zeros(1,length(lat), length(p));
k_3d = Teq_3d;
Teq_3d(:) = Teq(:);
k_3d(:) = k(:);

newton_heat = -k_3d.*(T-Teq_3d);


end

function Q_prescribed = compute_Q_prescribed(experiment_name, lat, lon, p)

Q_prescribed = zeros(length(lon), length(lat), length(p));
lat = lat';
p_3d = zeros(1,1,length(p));
p_3d(:) = p;

if strfind(experiment_name,'ns')
  H = 250;
  L = 50;
  if strfind(experiment_name,'ca')
    Q0 = 3;
  else
    Q0 = 3;
  end
elseif strfind(experiment_name,'nt')
  H = 350;
  L = 50;
  if strfind(experiment_name,'ca')
    Q0 = 3;
  else
    Q0 = 3*(250/350);
  end
elseif strfind(experiment_name,'wt')
  H = 350;
  L = 80;
  if strfind(experiment_name,'ca')
    Q0 = 3.7;
  else
    Q0 = 3.7*(50/80)*(250/350);
  end
else
  return;
end
dy = 20;
pcent = 150 + H;

Q_prescribed = max(0, Q0.*(1 - (lat/dy).^2 - ((lon-180)/L).^2 - ((p_3d-pcent)/H).^2)/86400);

end