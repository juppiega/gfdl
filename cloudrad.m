function cloudrad()

%N = 80;
numcols = 200;

%cloudfrac = zeros(N,1);
%cloudfrac((1:N/2)+3) = 0.15*(1-(((1:N/2) - N/4)/(N/4)).^2);
%cloudfrac(3*N/4:end) = 1.0*(1-(((3*N/4:N) - 7*N/8)/(N/7.5)).^2);
cloudfrac = load('example_cloudfrac.dat');
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

end

function compute_sw(cols, tau, lat, albedo)

solar_constant = 1360;
del_sol = 1.4;
del_sw = 0;

g = 0.8; f = g^2; 
tau_sc = tau*(1-f);
g_sc = (g-f)/(1-f);
coszen = cosd(lat);

alp = 0.75*coszen;
ga = 0.5;
d = alp - ga;
p = alp + ga;
pm1 = p - 1;

ncols = size(cols,2);
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


for i = 1:ncols
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