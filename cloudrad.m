function cloudrad()

%N = 80;
numcols = 200;

%cloudfrac = zeros(N,1);
%cloudfrac((1:N/2)+3) = 0.15*(1-(((1:N/2) - N/4)/(N/4)).^2);
%cloudfrac(3*N/4:end) = 1.0*(1-(((3*N/4:N) - 7*N/8)/(N/7.5)).^2);
cloudfrac = load('example_cloudfrac.dat');
<<<<<<< HEAD
=======
cloudfrac(58:end) = 0;
cloudfrac = cloudfrac*10;
>>>>>>> 6f14b6b46841960d3eb7859b829362777d428b4c
T = load('example_T.dat');
N = length(cloudfrac);

cols = compute_independent_columns(cloudfrac, numcols);
meanfrac = mean(cols,2);

figure;
plot(cloudfrac, 1:N, meanfrac, 1:N); set(gca,'ydir','reverse')
title('cloudfrac')

hl = 3E3; rho0 = 0.2E-3; reff = 10E-6;
z = flip(linspace(0,18E3,N+1));

cwp = rho0*hl*(exp(-z(2:end)/hl) - exp(-z(1:end-1)/hl));
tau = 1.5*cwp/reff/1E3;

figure;
plot(tau,1:N);set(gca,'ydir','reverse') 
title('tau')

lat = 0; albedo = 0.07;
compute_sw(cols, tau, lat, albedo)

<<<<<<< HEAD
=======
eps = 1 - exp(-cwp*140);
Ts = T(end)+0.5*(T(end)+T(end-1));
p_half = linspace(0,1000E2,N+1);
[tdt, tdt_cs] = compute_lw(cols, lat, eps, Ts, T, p_half);

figure;
subplot(1,2,1)
plot(cloudfrac, 1:N); set(gca,'ydir','reverse')
title('cloudfrac')
subplot(1,2,2)
plot((tdt-tdt_cs)*86400, 1:N); set(gca,'ydir','reverse')
title('Cloud heating [K/day]')

end

function [tdt, tdt_cs] = compute_lw(cols, lat, lw_emiss, Ts, T, p_half)

pstd_mks = 1E5;
grav = 9.81;
cp_air = 1005;
linear_tau = 0.1;
wv_exponent  = 4.0;
ir_tau_eq = 6.0;
ir_tau_pole = 1.5;
lw_tau_clear_0 = ir_tau_eq + (ir_tau_pole - ir_tau_eq)*sind(lat)^2;

N_half = length(p_half);
N = N_half - 1;
lw_tau_clear = zeros(length(p_half),1);
for k = 1:N_half
    lw_tau_clear(k) = lw_tau_clear_0 * ( linear_tau * p_half(k)/pstd_mks  ...
       + (1.0 - linear_tau) * (p_half(k)/pstd_mks)^wv_exponent );
end

lw_dtrans_clear = zeros(N,1);
for k = 1:N
    lw_dtrans_clear(k) = exp( -(lw_tau_clear(k+1) - lw_tau_clear(k)) );
end

cloudfree = zeros(N,1);
[lw_down_clear, lw_up_clear] = compute_lw_column(cloudfree, lw_dtrans_clear, lw_emiss, Ts, T);
olr_cs = lw_up_clear(1);
lw_surf_cs = lw_down_clear(end);
F_cs = lw_up_clear - lw_down_clear;
tdt_cs = grav/cp_air*(F_cs(2:end)-F_cs(1:end-1))./(p_half(2:end) - p_half(1:end-1));

numcols = size(cols,2);
olr = 0.0; lw_surf = 0.0; F = zeros(N+1,1);
for i = 1:numcols
    [lw_down, lw_up] = compute_lw_column(cols(:,i), lw_dtrans_clear, lw_emiss, Ts, T);
    olr = olr + lw_up(1)/numcols;
    lw_surf = lw_surf + lw_down(end)/numcols;
    F = F + (lw_up-lw_down)/numcols;
end
tdt = grav/cp_air*(F(2:end)-F(1:end-1))./(p_half(2:end) - p_half(1:end-1));

lw_cre = olr_cs - olr

end

function [lw_down, lw_up] = compute_lw_column(cf, lw_dtrans_clear, eps, Ts, T)

sig = 5.670374419E-8;
Tc = (1-eps'.*cf);
N = length(cf);
T_tot = Tc.*lw_dtrans_clear;

B = sig*T.^4;
lw_down = zeros(N+1,1);
lw_up = zeros(N+1,1);

lw_down(1) = 0;
for k = 1:N
    lw_down(k+1) = T_tot(k)*lw_down(k) + B(k)*(1-T_tot(k)); 
end

B_surf = sig*Ts^4;
lw_up(end) = B_surf;
for k = N:-1:1
    lw_up(k) = lw_up(k+1)*T_tot(k) + B(k)*(1-T_tot(k));
end

>>>>>>> 6f14b6b46841960d3eb7859b829362777d428b4c
end

function compute_sw(cols, tau, lat, albedo)

solar_constant = 1360;
del_sol = 1.4;
del_sw = 0;

g = 0.8; f = g^2; 
tau_sc = tau*(1-f);
<<<<<<< HEAD
g_sc = (g-f)/(1-f);
=======
>>>>>>> 6f14b6b46841960d3eb7859b829362777d428b4c
coszen = cosd(lat);

alp = 0.75*coszen;
ga = 0.5;
d = alp - ga;
p = alp + ga;
pm1 = p - 1;

<<<<<<< HEAD
ncols = size(cols,2);
=======
numcols = size(cols,2);
>>>>>>> 6f14b6b46841960d3eb7859b829362777d428b4c
N = size(cols,1);
R = zeros(N,1);
T = R;
T_dirs = R;
rdiffs = R;
tdiffs = R;
p2          = (1. - 3.*sind(lat).^2)/4.;
insolation  = 0.25 * solar_constant * (1.0 + del_sol * p2 + del_sw * sind(lat));
sw_toa_up = nan(N,1);
sw_surf_abs = nan(N,1);
sw_toa_up_cs = nan(N,1);


<<<<<<< HEAD
for i = 1:ncols
=======
for i = 1:numcols
>>>>>>> 6f14b6b46841960d3eb7859b829362777d428b4c
    %figure;
    %plot(cols(:,i),1:N); set(gca,'ydir','reverse')
    
    numclouds = 0;
    incloud = false;
    tau_sum = 0;
    for k = 1:N
        if cols(k,i) > 0
            if ~incloud
                incloud = true;
                numclouds = numclouds + 1;
                tau_sum = tau_sum + tau_sc(k);
            else
                tau_sum = tau_sum + tau_sc(k);
            end
        elseif cols(k,i) == 0 || k == N
            if incloud
                e = exp(-tau_sum/coszen);
                ts = (1-g)*tau_sum/(1-f);
                r = ts / (4.0/3.0 + ts);
                t = 1 - r;
                %R(numclouds) = d*t*e + d*r - d;
                T(numclouds) = p*t + d*r*e - pm1*e;
                R(numclouds) = 1 - T(numclouds);
                rdiffs(numclouds) = r;
                tdiffs(numclouds) = t;
                
                T_dirs(numclouds) = e;
                tau_sum = 0;
                incloud = false;
            end
        end
    end
    
    R(numclouds + 1) = albedo;
    T(numclouds + 1) = 1.0 - albedo;
    rdiffs(numclouds + 1) = albedo;
    tdiffs(numclouds + 1) = 1.0 - albedo;
    
    T_tot = T(1); T_dir = T_dirs(1); rdiff = rdiffs(1); tdiff = tdiffs(1);
    for n = 2:numclouds+1
        T_tot = T_dir*T(n) + (tdiffs(n)*(T_tot-T_dir + T_dir*R(n)*rdiff)) / (1-rdiff*rdiffs(n));
        T_dir = T_dir*T_dirs(n);
        rdiff_orig = rdiff;
        rdiff = rdiff + (tdiff*tdiff*rdiffs(n))/(1.0-rdiff*rdiffs(n));
        tdiff = tdiff*tdiffs(n)/(1-rdiff_orig*rdiffs(n));
    end
    
    sw_surf_abs(i) = T_tot*insolation;
    sw_toa_up(i) = insolation - sw_surf_abs(i);
    sw_toa_up_cs(i) = insolation * albedo;
end

format compact
insolation
albedo_tot = 1 - nanmean(sw_surf_abs)/insolation
sw_cre = nanmean(sw_toa_up_cs) - nanmean(sw_toa_up)

end


function cols = compute_independent_columns(cloudfrac, numcols)

N = length(cloudfrac);
cols = zeros(N, numcols);

for i = 1:numcols
    x = rand;
    if x > 1-cloudfrac(1)
        cols(1,i) = 1.0;
    end
    for k = 2:N
        if cols(k-1,i) == 0.0
            x = rand*(1-cloudfrac(k-1));
        end
        if x > 1-cloudfrac(k)
            cols(k,i) = 1.0;
        end
    end
end

end