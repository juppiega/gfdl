function [psi] = numericalStreamfunction_reanalysis (fname)     
% Longitude in degrees, p in hPa

r = 1/(15*86400);
eps = 1/(5*86400);
N2 = 1.5E-4;
rho0 = 0.67;
g = 9.81;
T0 = 300;
a = 6371E3;

vert_param = r*eps/N2;

load(fname)
z = flip(z * 1E3);
Q_mean = fliplr(Q);
%Q_mean(lon < 130,:) = 0;
x = double(lon*2*pi*a/360);
L = x(end)-x(1);
H = z(end)-z(1);
x_spacing = (x(2)-x(1))/L;
z = z / H;

psi = zeros(length(x)+2,length(z)+2); % Add +2 for the boundary points.
Q = zeros(size(psi));

Q(2:end-1,2:end-1) = Q_mean*g/T0/86400;
Q(1,:) = Q(end-1,:);
Q(end,:) = Q(2,:);

% Plot Q
figure;
[plot_x,plot_z] = meshgrid(lon,z*H/1000);
colormap(hot);
contourf(plot_x, plot_z, Q(2:end-1,2:end-1)'/g*T0*86400, 0:0.5:4, 'edgecolor','none'); 
view(2); colorbar; axis tight; title('Q [K/day]')
xlabel('Longitude','fontsize',15)
ylabel('z [km]','fontsize',15)
ylim([0,14])
set(gca,'fontsize',15)

% Q zonal mean
Q_zm = mean(Q(2:end-1, 2:end-1),1);
figure
plot(Q_zm/g*T0*86400, z*H/1000)
xlabel('Q [K/day]')
ylabel('z [km]')
xlim([0,1])

%%%%%%%%%
% SOLVER
%%%%%%%%%
lnz = log(z);
h1 = [lnz(1); diff(lnz)]';
h2 = [diff(lnz); lnz(end)-lnz(end-1)]';
% Compute the constant multiplier (A(p) in (2))
multiplier = -2/x_spacing^2 - bsxfun(@rdivide, 2*vert_param, (h1.*h2));
% Forcing
F = -Dx(Q, x_spacing)*L/N2;
% Solver error tolerance. Iteration stops, when the rms error between LHS and RHS in (1) is smaller than tol.
err = 1E30;
tol = 5E7;

% Main solver loop
while err > tol
    % Compute neighboring effects 
    x_neighbors = Dxx(psi, x_spacing) + 2*psi(2:end-1,2:end-1)/x_spacing^2;
    z_neighbors = Dzz(psi, h1, h2, z) + bsxfun(@rdivide, 2*psi(2:end-1,2:end-1), (h1.*h2));
    neighbor_points = x_neighbors + z_neighbors * vert_param * (L/H)^2;
    
    % Update psi
    psi(2:end-1,2:end-1) = bsxfun(@rdivide, (F - neighbor_points), multiplier);
    
    % Approximate F using updated psi and compute the rms error. 
    PDE_approx = Dxx(psi, x_spacing) + vert_param*(L/H)^2*Dzz(psi, h1, h2, z);
    residual = F - PDE_approx;
    err = sqrt(sum(residual(:).^2));

end

% Plot result
psi_p = -rho0*pi*6371E3*40/180*psi(2:end-1,2:end-1);
figure;
colormap(jet)
contourf(plot_x, plot_z, psi_p',-20E9:1E9:20E9,'edgecolor','none'); 
view(2); 
axis tight; 
xlabel('Longitude','fontsize',15)
ylabel('z [km]','fontsize',15)
set(gca,'fontsize',15) 
title('Streamfunction [kg/s]')
colorbar('fontsize',15)
yl = get(gca, 'ylim');
m = max(abs(psi_p(:)));
caxis([-m,m])
ylim([yl(1),15])

% u = -rho0*g* (psi(2:end-1,2:end)-psi(2:end-1,1:end-1))/(dp*100);
% u = u(:,2:end);
% figure;
% surf(plot_x, plot_p, u','edgecolor','none'); 
% view(2); 
% axis tight; 
% xlabel('Longitude','fontsize',15)
% ylabel('p [hPa]','fontsize',15)
% set(gca,'fontsize',15)
% set(gca,'ydir','reverse'); 
% title('u [m/s]')
% colorbar('fontsize',15)

display(['Numerical [kg/s]: ', num2str(max(rho0*pi*6371E3*40/180*psi(:))/1E9)])
%analyticalForm (plot_x, plot_p, L, H/sqrt(N2/(r*eps)), Q0, N2, r, eps);

end

% *********************
% Derivative functions 
% *********************
function deriv = Dx(X, dx)

deriv = (X(3:end,:) - X(1:end-2,:))/(2*dx);
deriv = deriv(:,2:end-1);

end

function deriv = Dxx(X, dx)

X(1,:) = X(end-1,:);
X(end,:) = X(2,:);
deriv = (X(1:end-2,:) + X(3:end,:) - 2*X(2:end-1,:))/dx^2;
deriv = deriv(:,2:end-1);

end

function deriv = Dzz(X, h1, h2, z)

h2X = bsxfun(@times, h2, X(:,1:end-2));
h1X = bsxfun(@times, h1, X(:,3:end));
h12X = bsxfun(@times, h1+h2, X(:,2:end-1));
deriv = bsxfun(@rdivide, 2*(h2X + h1X - h12X), (h1.*h2 .* (h1+h2)));
deriv = deriv(2:end-1,:);
deriv = bsxfun(@rdivide, deriv, z');

end
