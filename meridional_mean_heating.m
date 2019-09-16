function meridional_mean_heating(filename, slat, nlat)

compute_mean_heating(filename, slat, nlat)

end

function compute_mean_heating(Q_file, slat, nlat)

load(Q_file)
prec_file = 'era_precip.nc';
cp = ncread(prec_file,'cp');
if isnumeric(opt)
    ind = 1:12:size(cp,3); ind(end) = [];
    month_ind = [];
    for i = 1:length(opt)
        month_ind = [month_ind, ind+opt(i)-1];
    end
else
    t_data = ncread(prec_file,'time');
    timestamps = datenum('1900-01-01') + t_data/24;
    [el, la] = find_el_nino_months(timestamps);
    if strcmpi(opt,'elnino')
        month_ind = find(el);
    else
        month_ind = find(la);
    end
end

cp = mean(cp(:,:,month_ind),3);
lsp = ncread(prec_file,'lsp'); lsp = mean(lsp(:,:,month_ind),3);

lat = ncread(prec_file,'latitude');
lon = ncread(prec_file,'longitude');

z = compute_standard_height(p);
P_tot = lsp+cp;
median_P = median(P_tot(:));
prec_limit = 0;%100E-3*12/365;
Q = Q_vert_adv*86400;
for j = 1:length(lat)
    for i = 1:length(lon)
        if P_tot(i,j) < prec_limit
            Q(i,j,:) = 0;
        end
    end
end

% z = compute_standard_height(p)*1000;
% strat_heating = compute_heating(lsp, z, @strat_profile);
% conv_heating = compute_heating(cp, z, @conv_profile);
% Q = strat_heating + conv_heating; %Q = Q*86400;
% z = z / 1000;

ind = slat <= lat & lat <= nlat;
clat = cosd(lat(ind));
Q_mean = squeeze(sum(Q(:,ind,:).*clat'/sum(clat),2));
%Q_mean = max(Q_mean, 0);

cp = 1005;
g = 9.81;
Q_vert_int = trapz(p*100,Q_mean,2)/86400 * cp / g;
Q_prec = 2.5E6 * squeeze(sum(P_tot(:,ind).*clat'/sum(clat),2)) /86400*1E3;

%Q_vert_int = Q_prec;

figure;
plot(lon, smooth(Q_vert_int,5));
ylabel('Q [W/m^2]')
hold all
m = mean(Q_prec);
%plot(get(gca,'xlim'), [m,m])
plot(lon, smooth(Q_prec-m,5));
legend('Q_all','Q_latent*')
%ylim([0,250])
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 50, 200);
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 50, 230);
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 50, 240);
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 50, 250);
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 50, 260);
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 120, 270);

compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 71, 196);
compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 91, 223);
compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 106, 223);
compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 130, 230);
compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 66, 230);

% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 70, 190);
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 70, 220);
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 65, 210);
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 120, 240);
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, 60, 240);
% compute_domain_total_heat(Q_vert_int, slat, nlat, lon, -1, 361);

% z = 0:dz:16E3;
[X,Z] = meshgrid(double(lon), z);

figure;
%colormap(hot)
%contourf(X,Z,Q_mean',0:0.5:4,'edgecolor','none');
surf(X,Z,Q_mean','edgecolor','none');
view(2)
axis tight
colorbar
ylim([0,14])
title(Q_file)
caxis([-1,1])

if isnumeric(opt) 
    fname = ['Q_filtered_',sprintf('%d',opt),'_rain.mat'];
else
    fname = ['Q_filtered_',opt,'_rain.mat'];
end
    
Q = Q_mean;
save(fname,'Q','lon','z');

end

function [heat, max_Q_int] = compute_domain_total_heat(Q_vert_int, slat, nlat, lon, lon1, lon2)

RE = 6371E3;
A = RE^2*(sind(nlat)-sind(slat))*(lon2-lon1)*pi/180;
heat = mean(Q_vert_int(lon1<lon&lon<lon2))*A;
display(['Heat from ', num2str(lon1), ' to ', num2str(lon2), ': ', num2str(heat/1E15),' PW (',num2str(heat/A),' W/m^2)'])
max_Q_int = max(Q_vert_int(lon1<lon&lon<lon2));

end

function heating = compute_heating(rain, z, Q_fun)

Q_norm = Q_fun(z);
heating = zeros([size(rain), length(z)]);

for j = 1:size(heating,2)
    for i = 1:size(heating,1)
        heating(i,j,:) = Q_norm*rain(i,j)*100;
    end
end

end

function Q_i = strat_profile(z_i)

z = [0,2,3,4,6,7,8,9,11,13]*1E3;
Q = [0,-5,-4.5,-3.5,6,8.5,9.5,8.5,3.5,0];

Q_i = interp1(z, Q, z_i, 'pchip', 0);

end

function Q_i = conv_profile(z_i)

z = [0,2,3,4,5,6,7,8,9,10,14]*1E3;
Q = [0,6,7.5,7,6,3.5,2,1.5,1,0.5,0];

Q_i = interp1(z, Q, z_i, 'pchip', 0);

end