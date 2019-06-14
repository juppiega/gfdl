function [psi] = numericalStreamfunction_reanalysis (fname)     
% Longitude in degrees, p in hPa

r = 1/(15*86400);
eps = 1/(5*86400);
N2 = 1.5E-4;
rho0 = 0.73;
g = 9.81;
T0 = 300;
a = 6371E3;

vert_param = r*eps/N2;

load(fname)
Q_mean = Q;
x = double(lon*2*pi*a/360);
L = x(end)-x(1);
H = z(end)-z(1);
x_spacing = (x(2)-x(1))/L; 
z_spacing = (z(2)-z(1))/H;

psi = zeros(length(x)+2,length(z)+2); % Add +2 for the boundary points.
Q = zeros(size(psi));

Q(2:end-1,2:end-1) = Q_mean*g/T0/86400;

% Plot Q
figure;
[plot_x,plot_z] = meshgrid(lon,z/1000);
surf(plot_x, plot_z, Q(2:end-1,2:end-1)'/g*T0*86400, 'edgecolor','none'); 
view(2); colorbar; axis tight; title('Q [K/day]')
xlabel('Longitude','fontsize',15)
ylabel('z [km]','fontsize',15)
set(gca,'fontsize',15)

%%%%%%%%%
% SOLVER
%%%%%%%%%
multiplier = -2/x_spacing^2 - 2*vert_param*(L/H)^2/z_spacing^2;
% Forcing
F = -Dx(Q, x_spacing)*L/N2;
% Solver error tolerance. Iteration stops, when the rms error between LHS and RHS in (1) is smaller than tol.
err = 1E30;
tol = 5E7;

% Main solver loop
while err > tol
    % Compute neighboring effects 
    x_neighbors = Dxx(psi, x_spacing) + 2*psi(2:end-1,2:end-1)/x_spacing^2;
    z_neighbors = Dzz(psi, z_spacing) + 2*psi(2:end-1,2:end-1)/z_spacing^2;
    neighbor_points = x_neighbors + z_neighbors * vert_param * (L/H)^2;
    
    % Update psi
    psi(2:end-1,2:end-1) = bsxfun(@rdivide, (F - neighbor_points), multiplier);
    
    % Approximate F using updated psi and compute the rms error. 
    PDE_approx = Dxx(psi, x_spacing) + vert_param*(L/H)^2*Dzz(psi, z_spacing);
    residual = F - PDE_approx;
    err = sqrt(sum(residual(:).^2));

end

% Plot result
psi_p = -rho0*pi*6371E3*40/180*psi(2:end-1,2:end-1);
figure;
surf(plot_x, plot_z, psi_p','edgecolor','none'); 
view(2); 
axis tight; 
xlabel('Longitude','fontsize',15)
ylabel('z [km]','fontsize',15)
set(gca,'fontsize',15) 
title('Streamfunction [kg/s]')
colorbar('fontsize',15)
yl = get(gca, 'ylim');
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
analyticalForm (plot_x, plot_p, L, H/sqrt(N2/(r*eps)), Q0, N2, r, eps);

end

% *********************
% Derivative functions 
% *********************
function deriv = Dx(X, dx)

deriv = (X(3:end,:) - X(1:end-2,:))/(2*dx);
deriv = deriv(:,2:end-1);

end

function deriv = Dxx(X, dx)

deriv = (X(1:end-2,:) + X(3:end,:) - 2*X(2:end-1,:))/dx^2;
deriv = deriv(:,2:end-1);

end

function deriv = Dzz(X, dz)

deriv = (X(:,1:end-2) + X(:,3:end) - 2*X(:,2:end-1))/dz^2;
deriv = deriv(2:end-1,:);

end
