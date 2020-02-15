function [S_prof,p] = reanalysis_sp(filename)

% filename = 'ERA5.nc';
% 
% T = ncread(filename,'t',[1,1,1,1],[Inf,Inf,Inf,1]);
% sum_T = zeros(size(T));
% num_months = 0;
% info = ncinfo(filename);
% N = info.Dimensions(4).Length;
% 
% for k = 1:N
%     T_this = ncread(filename,'t',[1,1,1,k],[Inf,Inf,Inf,1]);
%     sum_T = sum_T + T_this;
%     num_months = num_months + 1;    
% end
% 
% T = squeeze(sum_T / num_months);

load(filename)

aver_lat = 20;
p = p*100;
p_3d = zeros(1,1,size(T,3));
p_3d(1,1,:) = p;
ind = -aver_lat <= lat & lat <= aver_lat;
T_pert = T - mean(T,1);
Q_tot = max(Q_tot, 0);
Q = Q_tot;% - mean(Q_tot,1);
mean_T_pert = squeeze(mean(T_pert(:,ind,:),2));
mean_T = squeeze(mean(T(:,ind,:),2));
mean_Q = squeeze(mean(Q(:,ind,:),2)); %mean_Q = mean_Q - mean(mean_Q,1);

h1 = p_3d(1,1,2:end-1) - p_3d(1,1,1:end-2); 
h2 = p_3d(1,1,3:end) - p_3d(1,1,2:end-1); 
T_p = -h2./(h1.*(h1+h2)).*T(:,:,1:end-2) + (h2-h1)./(h2.*h1).*T(:,:,2:end-1) + h1./(h2.*(h1+h2)).*T(:,:,3:end);
alpha = 287*T(:,:,2:end-1)./p_3d(1,1,2:end-1);
S_p = alpha/1004 - T_p;
mean_Sp = squeeze(mean(S_p(:,ind,:),2));
theta = T.*(1000E2./p_3d).^(287/1004);
theta_p = -h2./(h1.*(h1+h2)).*theta(:,:,1:end-2) + (h2-h1)./(h2.*h1).*theta(:,:,2:end-1) + h1./(h2.*(h1+h2)).*theta(:,:,3:end);
g = 9.81;
rho = 1./alpha;
N2 = -g^2 .* rho .* theta_p./theta(:,:,2:end-1);
mean_N2 = squeeze(mean(N2(:,ind,:),2));

figure;
colormap(jet)
[X,Z] = meshgrid(lon,p(2:end-1)/100);
surf(X,Z,mean_Sp','edgecolor','none');
view(2);
axis tight
colorbar
title('S_p')
set(gca,'Ydir','reverse')
caxis([0, 1E-3])

figure;
colormap(jet)
[X,Z] = meshgrid(lon,p(2:end-1)/100);
surf(X,Z,mean_N2','edgecolor','none');
view(2);
axis tight
colorbar
title('N2')
set(gca,'Ydir','reverse')
%caxis([0, 1E-3])

% figure;
% [X,Z] = meshgrid(lon,p/100);
% surf(X,Z,mean_T','edgecolor','none');
% view(2);
% axis tight
% colorbar
% title('T')
% set(gca,'Ydir','reverse')

% figure;
% [X,Z] = meshgrid(lon,p/100);
% surf(X,Z,mean_Q','edgecolor','none');
% view(2);
% axis tight
% colorbar
% title('Q')
% set(gca,'Ydir','reverse')

ind = 50 < lon & lon < 200;
T_prof = mean(mean_T(ind,:),1);
z = compute_actual_z(p', T_prof);
p_mid = 0.5*(p(1:end-1)+p(2:end));
lapse_rate = (T_prof(1:end-1) - T_prof(2:end)) ./ (z(1:end-1) - z(2:end))*1E3;

[T_ma,Sp_ma,p_ma]=moist_Sp(T_prof(end-1)-273.15);
z_ma = compute_actual_z(p_ma', T_ma')';
lr_ma = (T_ma(1:end-1) - T_ma(2:end)) ./ (z_ma(1:end-1) - z_ma(2:end))*1E3;

ind = 50 < lon & lon < 200;
S_prof = mean(mean_Sp(ind,:),1);
figure
plot(S_prof,p(2:end-1)/100, Sp_ma, p_ma/100)
title(['S_p ', filename(1:end-4)])
legend('Obs','Moist ad. of T_s = 24C')
set(gca,'Ydir','reverse')
xlim([0, 1E-3])
ylim([150,1000])

figure
plot(lapse_rate, p_mid/100, lr_ma, p_ma(1:end-1)/100)
title('Lapse rate MC/WP')
legend('Obs','Moist ad. of T_s = 24C')
%plot(T_prof,z)
%title('T Maritime cont./WP')
set(gca,'Ydir','reverse')
%ylabel('z')
%xlim([0, 2])
%ylim([0,15])
ylim([150,1000])

es = exp(77.3450 + 0.0057*T_prof - 7235 ./ T_prof) ./ T_prof.^8.2;
q_sat = 0.622*es./p';
mse = T_prof + 9.81*z/1004 + 2.5E6*q_sat/1004;

figure
plot(T_prof, p/100, T_ma, p_ma/100)
title('T MC/WP')
legend('Obs','Moist ad. of T_s = 24C')
%plot(T_prof,z)
%title('T Maritime cont./WP')
set(gca,'Ydir','reverse')
%ylabel('z')
%xlim([300, 400])
%ylim([0,15])
ylim([150,1000])

Q_prof = mean(mean_Q(ind,:),1)*86400;
figure
z = compute_standard_height(p/100);
plot(Q_prof,z)
title('Q Maritime cont./WP')
%set(gca,'Ydir','reverse')
ylabel('z')
xlim([0, 2])
ylim([0,15])
z_prof = max(Q_prof)*(1-((z - 6)/5).^2);
p_prof = max(Q_prof)*(1-((p/100 - 500)/350).^2);
hold all
plot(z_prof, z)
plot(p_prof, z)
legend('ERA5','Q(z)','Q(p)')

% Q_prof = mean(mean_Q(ind,:),1)*86400;
% figure
% plot(Q_prof,p/100)
% title('Q')
% set(gca,'Ydir','reverse')
% ylabel('p')
% xlim([0, 2])
% ylim([10,1000])
% hold all
% plot(p_prof, p/100)

end
