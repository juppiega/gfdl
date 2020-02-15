function solve_approx_H()

aver_lat = 30;
load streamfun_months.mat
psimax_obs = abs(psimax_e);
Hz = 5E3;
N2 = 1.5E-4;
Hmax = 105*86400*1.22E-2*Hz;
Hmin = 0.01*86400*1.22E-2*Hz;
H = Hmin:1E4:Hmax;
Nreps = H/Hz;
L = L/360 * 2*pi*6371E3;
psi_corr = zeros(size(H));
j = 1:12;
for i = 1:length(H)
    psimax_est = 4*9.81*1E15*Q(j)./sqrt(2*L(j).^2 + 3*L(j)*H(i) + H(i)^2)/(3*pi*250*N2*1004)*Nreps(i);
    psi_corr(i) = rms(psimax_est - psimax_obs(j)*1E9);%sqrt(mean((psimax_est./(psimax_obs(j)*1E9)-1).^2));%corr(psimax_est, psimax_obs(j));
end

fs = 15;
figure;
[~,i] = min(psi_corr);
psimax_est = 4*9.81*1E15*Q(j)./sqrt(2*L(j).^2 + 3*L(j)*H(i) + H(i)^2)/(3*pi*250*N2*1004)*Nreps(i);
%psimax_est = psimax_est * mean(psimax_e./psimax_est);
scatter(1E3./L(j),psimax_obs(j),'filled','o','markeredgecolor','r','markerfacecolor','r');
hold all
scatter(1E3./L(j),psimax_est/1E9,'o','markeredgecolor','r');
xlabel('1/L [1/km]','fontsize',fs);
ylabel('M^E [Tg/s]','fontsize',fs);
title('a) ERA5 and analytic strength','fontsize',fs);
legend('ERA5','Analytic','location','northwest')
set(gca,'fontsize',fs-2)

% figure;
% scatter(psimax_obs(j), psimax_est);

psimax_obs = abs(psimax_w);
psi_corr_w = zeros(size(H));
j = 4:12;
for i = 1:length(H)
    psimax_est = 4*9.81*1E15*Q(j)./sqrt(2*L(j).^2 + 3*L(j)*H(i) + H(i)^2)/(3*pi*250*N2*1004)*Nreps(i);
    psi_corr_w(i) = rms(psimax_est - psimax_obs(j)*1E9);%sqrt(mean((psimax_est./(psimax_obs(j)*1E9)-1).^2));%corr(psimax_est, psimax_obs(j));
end

t_damp = H/(sqrt(N2)*Hz)/86400;
figure;
semilogx(t_damp, psi_corr_w/1E9,'linewidth',2.0,'color','b')
hold all
semilogx(t_damp, psi_corr/1E9,'linewidth',2.0,'color','r')
[~,i] = min(psi_corr);
characteristic_damping_time_east = t_damp(i)
lambda = 1/Nreps(i).^2
[~,i] = min(psi_corr_w);
characteristic_damping_time_west = t_damp(i)
lambda = 1/Nreps(i).^2

t0 = characteristic_damping_time_east;
%hold all
%plot([t0,t0], get(gca,'ylim'))
xlabel('$\sqrt{r^{-1} \alpha^{-1}}$ [days]', 'interpreter','latex','fontsize',fs)
ylabel('RMSE [Tg/s]','fontsize',fs)
title('b) RMSE vs damping','fontsize',fs)
%title('r: analytical vs ERA5','fontsize',fs)
set(gca,'fontsize',fs-2)
legend('West','East','location','northwest')
%axis tight
xlim([min(t_damp), max(t_damp)])

end