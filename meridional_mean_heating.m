function meridional_mean_heating()

compute_mean_heating([12,1,2])

end

function compute_mean_heating(months)

filename = 'era_precip.nc';
cp = ncread(filename,'cp');
ind = 1:12:size(cp,3); ind(end) = [];
month_ind = [];
for i = 1:length(months)
    month_ind = [month_ind, ind+months(i)-1];
end

cp = mean(cp(:,:,month_ind),3);
lsp = ncread(filename,'lsp'); lsp = mean(lsp(:,:,month_ind),3);

lat = ncread(filename,'latitude');
lon = ncread(filename,'longitude');

dz = 500; % m
strat_heating = compute_heating(lsp, dz, @strat_profile);
conv_heating = compute_heating(cp, dz, @conv_profile);
Q = strat_heating + conv_heating;

aver_lat = 20;
ind = -aver_lat <= lat & lat <= aver_lat;
Q_mean = squeeze(mean(Q(:,ind,:),2));

z = 0:dz:16E3;
[X,Z] = meshgrid(double(lon), z);

figure;
surf(X,Z,Q_mean', 'edgecolor','none');
view(2)
axis tight
colorbar
title(num2str(months,'%d%d%d'))
caxis([0,4])

fname = ['Q_',num2str(months,'%d%d%d'),'.mat'];
Q = Q_mean;
save(fname,'Q','lon','z');

end

function heating = compute_heating(rain, dz, Q_fun)

z = 0:dz:16E3;
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