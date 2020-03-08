function [T,S_p,p]=moist_Sp(T)
% T in Celsius

T = T+273.15;
[p,T] = ode45(@moist_adiabat, 1000E2:-1E2:100E2, T);

figure
subplot(1,2,1)
plot(T-273.15,p/100)
set(gca,'ydir','reverse')
xlabel('T')
ylabel('p')

[T_p, rho] = moist_adiabat(p,T);
S_p = 1./(1004*rho) - T_p;
subplot(1,2,2)
plot(S_p,p/100)
set(gca,'ydir','reverse')
xlabel('S_p')
ylabel('p')

end

function [lapse_rate,rho] = moist_adiabat(p,T)

es = exp(77.3450 + 0.0057*T - 7235 ./ T) ./ T.^8.2;
r = 0.622*es./(p-es);
Rd = 297;
Rv = 461.5;
L = 2.5E6;
cp = 1004;
rho = p./(Rd*T);
g = 9.81;

T_z = g*(1 + L*r./(Rd*T)) ./ (cp + L^2*r./(Rv*T.^2));
lapse_rate = T_z./(rho*g);

end