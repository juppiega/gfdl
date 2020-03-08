function write_Q_ncdf(filename)

copyfile('heat_template.nc', [filename,'.nc']);
load([filename,'.mat'])
lat_era = lat;
lat_hs = ncread('heat_template.nc','lat');
p_era = p;
p_hs = ncread('heat_template.nc','pfull');

filename = [filename,'.nc'];
Q = max(Q_tot, 0);
for i = 1:size(Q,1)
    for k = 1:size(Q,3)
        Q(i,:,k) = interp1(lat_era,squeeze(Q(i,:,k)),lat_hs);
    end
end

for i = 1:size(Q,1)
    for j = 1:size(Q,2)
        Q(i,j,:) = interp1(p_era,squeeze(Q(i,j,:)),p_hs,'linear',0);
    end
end

ind = -20 <= lat_hs & lat_hs <= 20;
mean_Q = squeeze(mean(Q(:,ind,:),2));
figure
[X,Z] = meshgrid(lon,p_hs);
surf(X,Z,mean_Q'*86400,'edgecolor','none');
view(2);
axis tight
colorbar
title('Q [K/day]')
set(gca,'Ydir','reverse')

ncwrite(filename, 'local_heating', Q, [1,1,1,1]);
ncwrite(filename, 'local_heating', Q, [1,1,1,2]);
ncwrite(filename, 'local_heating', Q, [1,1,1,3]);
ncwriteatt(filename, 'time',  'units', 'days since 0000-00-00 00:00:00');
%ncdisp(filename)

end