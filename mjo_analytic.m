function mjo_analytic()

timescale = 6371E3/(2.3E-11*(456E3)^2);
lengthscale = 6371E3;

alp = 1.5; 
G = 0.1; Da = 1.5; 
Ca = 0.8; 
X = 1.5; kappa = 2; gamma = 1; diff = 0.1;
f_ke = @(x,k) khairoutdinov2018(x, k, alp, G, Da, Ca, X, kappa, gamma, diff);

mu = 0.3;  
lam = 0.1; 
the = 0.8; 
cs = 0.8; 
%Ca = 0;
taus0 = 0.01; D0 = 0.3; 
ao = 1.5;
tauh = 0.1;
f_full = @(x,k) full_ocean(x, k, alp, G, Da, Ca, X, kappa, gamma, diff, mu, lam, the, cs, taus0, D0, ao, tauh);

% [Ca_mat, mu_mat] = meshgrid(0.2:0.2:1.4, 0:0.05:0.5);
% c_mat = zeros(size(mu_mat));
% grow_mat = c_mat;
% for j = 1:size(c_mat,2)
%     for i = 1:size(c_mat,1)
%         f = @(x,k) full_ocean(x, k, alp, G, Da, Ca_mat(i,j), X, kappa, gamma, diff, mu_mat(i,j), lam, the, Ca_mat(i,j), taus0, D0, ao, tauh);
%         [c_mat(i,j), grow_mat(i,j)] = solveWaveProperties(f, 1);
%     end
% end
% 
% %grow_mat(c_mat >= 1) = nan;
% %c_mat(c_mat >= 1) = nan;
% 
% figure;
% surf(Ca_mat, mu_mat, c_mat);view(2);colorbar;xlabel('Ca');ylabel('\mu');title('Phase speed');axis tight
% figure;
% surf(Ca_mat, mu_mat, grow_mat);view(2);colorbar;xlabel('Ca');ylabel('\mu');title('Growth rate');axis tight

k = 1:10;
[speed_ke18, growth_ke18, freq_ke18] = solveWaveProperties(f_ke, k);

[speed_fo, growth_fo, freq_fo] = solveWaveProperties(f_full, k);

% figure;
% scatter(k, freq_ke18, growth_ke18*256/growth_ke18(1),'filled','markerEdgeColor','r','markerFaceColor','r');
% hold all
% scatter(k, freq_fo, growth_fo*256/growth_ke18(1),'markerEdgeColor','b');
% xlabel('Zonal wavenumber')
% ylabel('f')
% legend('KE18 (no ocean, Ca = 0.8)','Full ocean')

figure;
plot(k, speed_ke18, k, speed_fo, 'linewidth',2.0)
title('Phase speed')
xlabel('Zonal wavenumber')
legend('KE18','Full ocean')

figure;
plot(k, growth_ke18, k, growth_fo, 'linewidth',2.0)
title('Growth rate')
xlabel('Zonal wavenumber')
legend('KE18','Full ocean')

figure;
plot(k, freq_ke18, k, freq_fo, 'linewidth',2.0)
title('Frequency')
xlabel('Zonal wavenumber')
legend('KE18','Full ocean')

end

function [speed, growth, freq] = solveWaveProperties(target_fun, k_in)
% Looks for the most rapidly growing eastward propagating mode

% [sig_re,sig_im] = meshgrid(-1:0.01:1,-2:0.01:2);
% sig = complex(sig_re, sig_im);
% k_this = 1;
% f_out = target_fun(sig, k_this);
% 
% figure;
% % subplot(3,1,1)
% % surf(sig_re, sig_im, real(f_out),'edgecolor','none'); view(2); colorbar;axis tight;title('re')
% % subplot(3,1,2)
% % surf(sig_re, sig_im, imag(f_out),'edgecolor','none'); view(2); colorbar;axis tight;title('im')
% % subplot(3,1,3)
% freq = -sig_im;
% surf(sig_re, freq, log10(abs(f_out)),'edgecolor','none'); view(2); colorbar;axis tight;
% title(['log(abs) k = ',num2str(k_this)]);xlabel('Growth rate'); ylabel('freq')

di = 0.001; dr = 0.001;
res = -5:dr:5; ims = -20:di:0;
res_est = zeros(size(ims));
for k = 1:length(k_in)
    k_this = k_in(k);
    for i = 1:length(ims)
        sig = complex(res,ims(i));
        d = diff(abs(target_fun(sig,k_this)));
        min_ind = find(d<0, 1, 'last');
        if ~isempty(min_ind)
            res_est(i) = res(min_ind);
        else
            res_est(i) = res(1);
        end
    end
    
    [growth(k), max_ind] = max(res_est);
    freq(k) = -ims(max_ind);
end

speed = freq ./ k_in;

end

function val = khairoutdinov2018(sig, k, alp, G, D, C, X, kappa, gamma, diff)

val = -(1+C)*(-alp*k*1i - G*k*k - D*sig) + ...
    (-alp*k*1i - k*k - sig.*sig - sig*X) .* (-diff*k*k + kappa*C - sig*gamma);

end

function val = full_ocean(s, k, a, G, Da, ca, x, kappa, y, diff, mu, lam, the, cs, taus0, D0, ao, tauh)

val = -1i*k*lam*mu*(1+ca+diff*k^2-ca*kappa+s*y) + ...
    (-s-tauh).*(-1i*k*(-(1+ao)*mu+a*(D0+s+taus0)).*(1+ca+diff*k^2-ca*kappa+s*y) + ... 
    k^2*(cs*(-1+G)*mu-(D0+s+taus0).*(G+ca*G+diff*k^2-ca*kappa+s*y)) + ...
    s.*(cs*mu*(Da-s-x) + mu*the*(1+ca+diff*k^2-ca*kappa+s*y) - ...
    (D0+s+taus0).*((1 + ca)*Da+(s+x).*(diff*k^2-ca*kappa+s*y))));


end