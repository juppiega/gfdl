function era_precip()

plotRainHeating([6,7,8])
plotRainHeating([12,1,2])

end

function plotRainHeating(months)

filename = 'era_precip.nc';
cp = ncread(filename,'cp');
ind = 1:12:size(cp,3); ind(end) = [];
month_ind = [];
for i = 1:length(months)
    month_ind = [month_ind, ind+months(i)-1];
end

cp = mean(cp(:,:,month_ind),3)*365;
lsp = ncread(filename,'lsp'); lsp = mean(lsp(:,:,month_ind),3)*365;

lat = ncread(filename,'latitude');
lon = ncread(filename,'longitude');

P_tot = cp+lsp;

Q = 3*1E3*P_tot*9.81*2.5E6/(365*8*1005*350E2);

aver_lat = 20;
ind = -aver_lat <= lat & lat <= aver_lat;
Q_mean = mean(Q(:,ind),2);
P_tropics = mean(P_tot(:,ind),2);
P_sorted = sort(P_tropics(:));

% figure;
% histogram(P_tot, 50, 'normalization','pdf');
% xlabel('P')


figure;
plot(lon, Q_mean);
title(['Q0 merid. mean: ', num2str(months)])
ylim([0 2.5])
hold all
plot(get(gca,'xlim'), [mean(Q_mean), mean(Q_mean)])
ind = 50<lon & lon<230;
Q_tot = sum(Q_mean(ind));
display(['Q_tot 50 to 230 longitude: ', num2str(Q_tot)])

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
title(['total P: ',num2str(months)])
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