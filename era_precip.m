function era_precip(opt)

plotRainHeating(opt)

end

function plotRainHeating(opt)

filename = 'era_precip.nc';
prec_file = 'era_precip.nc';
cp = ncread(prec_file,'cp');
if isnumeric(opt)
    ind = 1:12:size(cp,3); ind(end) = [];
    month_ind = [];
    for i = 1:length(opt)
        month_ind = [month_ind, ind+opt(i)-1];
    end
    aver_lat = 30;
else
    t_data = ncread(prec_file,'time');
    timestamps = datenum('1900-01-01') + t_data/24;
    [el, la] = find_el_nino_months(timestamps);
    if strcmpi(opt,'elnino')
        month_ind = find(el);
    else
        month_ind = find(la);
    end
    aver_lat = 15;
end

lat = ncread(filename,'latitude');
lon = ncread(filename,'longitude');
ind = -aver_lat <= lat & lat <= aver_lat;
lsp = ncread(filename,'lsp'); 
P_tot = (cp(:,:,month_ind)+lsp(:,:,month_ind))*365;
P_tot_merid = squeeze(mean(P_tot(:,ind,:),2));
ind_lon = 50 <= lon & lon <= 200;
P_tot_merid = squeeze(mean(P_tot_merid(ind_lon,:),1));

figure
hist(P_tot_merid,10)
title('P MC/WP [m/yr]')
nsr = std(P_tot_merid) / mean(P_tot_merid)

cp = mean(cp(:,:,month_ind),3)*365;
lsp = mean(lsp(:,:,month_ind),3)*365;


P_tot = cp+lsp;

Q = 3*1E3*P_tot*9.81*2.5E6/(365*8*1005*350E2);

Q_watt = 2.5E6*P_tot/365/86400*1E3;
Q = Q_watt;

Q_mean = mean(Q(:,ind),2);

coslat = cosd(lat(ind));
Q_domain_mean = sum(Q(:,ind).*coslat'/sum(coslat),2);
Q_mean = Q_domain_mean;

% figure;
% histogram(P_tot, 50, 'normalization','pdf');
% xlabel('P')

dQ_dx = [diff(Q_mean); Q_mean(1)-Q_mean(end)];

figure;
plot(lon, Q_mean,'linewidth',1.5,'color','k');
fs = 15;
title(['a) El Ni',char(0241),'o Precipitation'],'fontsize',fs)
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
compute_domain_total_heat(Q_mean, -aver_lat, aver_lat, lon, lon(wb), lon(eb));
resid_heat = Q_mean-m; resid_heat(resid_heat < 0) = 0;
[heat, power, width] = compute_domain_total_heat(resid_heat, -aver_lat, aver_lat, lon, lon(wb), lon(eb));
compute_domain_total_heat(resid_heat, -aver_lat, aver_lat, lon, -1, 361);
[heat, power, width] =compute_domain_total_heat(Q_mean-m, -aver_lat, aver_lat, lon, 65, lon(eb));
%compute_domain_total_heat(Q_mean, -aver_lat, aver_lat, lon, 67, lon(eb));

dim = [.3 .1 .1 .3];
str = ['$\mathrm{P: ',num2str(power,'%3.1f'),'~W/m^2}$ ',char(10),'$\mathrm{Q_{tot}: ',num2str(heat/1E15,'%4.2f'),'~PW}$',char(10),...
    '$\mathrm{\Delta \lambda: ',num2str(round(width)),'^\circ}$'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','fontsize',fs,'interpreter','latex');

[X,Y] = meshgrid(lon, lat);
% figure;
% surf(X,Y,cp','edgecolor','none')
% view(2)
% colorbar
% axis tight
% title('cp')
% caxis([0 3.5])

figure;
surf(X,Y,P_tot','edgecolor','none')
view(2)
colorbar
axis tight
title(['total P [m/year]: ',num2str(opt)])
caxis([0 3.5])

% figure;
% surf(X,Y,Q','edgecolor','none')
% view(2)
% colorbar
% axis tight
% title('Q [K/d]')
% caxis([0 4])

% figure
% zonal_mean_pert = (cp+lsp) - mean(cp+lsp, 1);
% surf(X,Y,zonal_mean_pert','edgecolor','none')
% view(2)
% colorbar
% axis tight
% title('P - [P]')
% caxis([0 3.5])

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