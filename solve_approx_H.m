function solve_approx_H()

load streamfun_months.mat
H = 1E5:1E4:1E8;
L = L/360 * 2*pi*6371E3;
psi_corr = zeros(size(H));
j = 1:12;
for i = 1:length(H)
    psimax_est = Q(j)./sqrt(2*L(j).^2 + 3*L(j)*H(i) + H(i)^2);
    psi_corr(i) = corr(psimax_est, psimax_e(j));
end

t_damp = H/(1E-2*5E3)/86400;
figure;
semilogx(t_damp, psi_corr)
[~,i] = max(psi_corr);
characteristic_damping_time = t_damp(i)
t0 = characteristic_damping_time;
hold all
plot([t0,t0], get(gca,'ylim'))
xlabel('$\sqrt{T_N T_R}$', 'interpreter','latex')
ylabel('r: analytical vs obs.')

end