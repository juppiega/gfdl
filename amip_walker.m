function amip_walker(opt)

fname = ['amip_udiv_',num2str(opt),'.nc'];
analyzeStreamfun(fname, 30)
analyzeRain(opt)

end

function analyzeRain(opt)

filename = 'pr_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc';
prec_file = filename;
P_tot = ncread(prec_file,'pr')*1E-3*365*86400;
Nmon = size(P_tot,3);
if isnumeric(opt)
    month_ind = opt:12:Nmon;
    aver_lat = 30;
else
    t_data = ncread(prec_file,'time');
    timestamps = datenum('1979-01-01') + t_data;
    [el, la] = find_el_nino_months(timestamps);
    if strcmpi(opt,'elnino')
        month_ind = find(el);
    else
        month_ind = find(la);
    end
    aver_lat = 15;
end

P_tot = mean(P_tot(:,:,month_ind),3);

lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
ind = -aver_lat <= lat & lat <= aver_lat;

Q = 2.5E6*P_tot/365/86400*1E3;

coslat = cosd(lat(ind));
Q_mean = sum(Q(:,ind).*coslat'/sum(coslat),2);

figure;
plot(lon, Q_mean,'linewidth',1.5,'color','k');
fs = 15;
title(['P ',num2str(opt)],'fontsize',fs)
ylim([0 250])
xlim([0,360])
hold all
m = mean(Q_mean);
plot(get(gca,'xlim'), [m, m],'linewidth',1,'color','k','linestyle','-')
xlabel('Longitude','fontsize',fs)
ylabel('W/m^2','fontsize',fs)
set(gca,'fontsize',fs)

region_lims = find(Q_mean - m < 0);
[~,widest_region] = max(diff(region_lims));
wb = region_lims(widest_region);
eb = region_lims(widest_region+1);
%compute_domain_total_heat(Q_mean, -aver_lat, aver_lat, lon, lon(wb), lon(eb));
resid_heat = Q_mean-m; resid_heat(resid_heat < 0) = 0;
if strcmpi(opt,'elnino')
    [heat, power, width] =compute_domain_total_heat(Q_mean-m, -aver_lat, aver_lat, lon, 65, lon(eb));
else
    [heat, power, width] = compute_domain_total_heat(resid_heat, -aver_lat, aver_lat, lon, lon(wb), lon(eb));
end
%compute_domain_total_heat(resid_heat, -aver_lat, aver_lat, lon, -1, 361);
%compute_domain_total_heat(Q_mean, -aver_lat, aver_lat, lon, 67, lon(eb));

dim = [.3 .1 .1 .3];
str = ['$\mathrm{P: ',num2str(power,'%3.1f'),'~W/m^2}$ ',char(10),'$\mathrm{Q_{tot}: ',num2str(heat/1E15,'%4.2f'),'~PW}$',char(10),...
    '$\mathrm{\Delta \lambda: ',num2str(round(width)),'^\circ}$'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','fontsize',fs,'interpreter','latex');

[X,Y] = meshgrid(lon, lat);
figure;
surf(X,Y,P_tot','edgecolor','none')
view(2)
colorbar
axis tight
title(['total P [m/year]: ',num2str(opt)])
caxis([0 3.5])

end

function analyzeStreamfun(fname_div, latb)

u_div = ncread(fname_div,'u_div');
p = ncread(fname_div,'level');
lon = ncread(fname_div,'longitude');
lat = ncread(fname_div,'latitude');
aver_lat = latb;

ind = -aver_lat <= lat & lat <= aver_lat;
coslat = zeros([1, sum(ind), 1]);
coslat(:) = cosd(lat(ind));
udiv_mean = squeeze(sum(u_div(:,ind,:,1).*coslat/sum(coslat),2));
streamfun = cumtrapz(p, udiv_mean, 2) * pi *6371E3 * 2 * aver_lat/180/9.81;

strfun_west_280 = streamfun(lon < 280, :);
display(['AMIP Max west of 280 [10^9 kg/s]: ', num2str(max(strfun_west_280(:))/1E9)])
strfun_west_200 = streamfun(lon < 200, :);
display(['AMIP Min west of 200 [10^9 kg/s]: ', num2str(min(strfun_west_200(:))/1E9)])

z = compute_standard_height(p/100);
[LON,Z] = meshgrid(lon, z);
figure;
colormap(jet(16))
contourf(LON,Z,streamfun'/1E9,16,'edgecolor','none')
view(2)
colorbar('location','northoutside')
m = 40;
caxis([-m,m])
axis tight
yl = get(gca, 'ylim');
ylim([yl(1),15])
fs = 15;
title(['d) La Ni',char(0241),'a streamfunction [Tg/s]'],'fontsize',fs)
xlabel('Longitude','fontsize',fs)
ylabel('Altitude [km]','fontsize',fs)
set(gca,'fontsize',fs)


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