function [psi] = numericalStreamfunction (sig_lon, p_max, sig_p, dlon, dp)     
% Longitude in degrees, p in hPa

r = 1/(2*86400);
eps = 1/(2*86400);
N2 = 1.5E-4;
g = 9.81;
T0 = 250;
rho0 = 0.67;
a = 6371E3;
Q0 = 3*g/T0/86400; % m/s^3
display(['Q0: ', num2str(Q0/g*T0*86400)]);
L = (sig_lon/360)*2*pi*a;
H = sig_p*100/(rho0*g)*sqrt(N2/(r*eps));
vert_param = r*eps*rho0^2*g^2/N2;

x = (-180+dlon:dlon:180-dlon)*2*pi*a/360;
x_spacing = 2*pi*a*dlon/360;
p = (dp:dp:1000)' * 100;
if p(end) < 1000E2
    p = [p;p(end)+dp*100];
end
[X,P] = meshgrid(x,p);
X = X'; P = P';

psi = zeros(length(x)+2,length(p)+2); % Add +2 for the boundary points.
Q = zeros(size(psi));

Q(2:end-1,2:end-1) = max(Q0 * (1 - (X/L).^2 - ((P/100 - p_max)/sig_p).^2), 0);

% Plot Q
figure;
[plot_x,plot_p] = meshgrid(x/(2*pi*a)*360,p/100);
surf(plot_x, plot_p, Q(2:end-1,2:end-1)', 'edgecolor','none'); 
view(2); colorbar; axis tight; set(gca,'ydir','reverse'); title('Q 12-month mean')
xlabel('Longitude','fontsize',15)
ylabel('p [hPa]','fontsize',15)
set(gca,'fontsize',15)

% ****************************************************************************************************************
% The finite difference method is of the form:  A(p)*psi + h(psi_neighbors) = -d_x(Q)/N^2 = F   (2)                
% where A is constant in x (varies with p) and for each iteration, and h is a function of neighboring (in x and p) psi values.
% For each iteration, h is recomputed, and psi is updated simply as: psi = (F - h) / A
% ****************************************************************************************************************
% Compute pressure grid spacings for vertical difference.
h1 = [p(1); diff(p)]';
h2 = [diff(p); p(end)-p(end-1)]';
% Compute the constant multiplier (A(p) in (2))
multiplier = -2/x_spacing^2 - bsxfun(@rdivide, 2*vert_param, (h1.*h2));
% Forcing
F = -Dx(Q, x_spacing)/N2;
% Solver error tolerance. Iteration stops, when the rms error between LHS and RHS in (1) is smaller than tol.
err = 1E30;
tol = 1E-11;

% Main solver loop
while err > tol
    % Compute neighboring effects (h in equation (2))
    x_neighbors = Dxx(psi, x_spacing) + 2*psi(2:end-1,2:end-1)/x_spacing^2;
    p_neighbors = Dpp(psi, h1, h2) + bsxfun(@rdivide, 2*psi(2:end-1,2:end-1), (h1.*h2));
    neighbor_points = x_neighbors + p_neighbors * vert_param;
    
    % Update psi
    psi(2:end-1,2:end-1) = bsxfun(@rdivide, (F - neighbor_points), multiplier);
    
    % Approximate F using updated psi and compute the rms error. 
    PDE_approx = Dxx(psi, x_spacing) + vert_param*Dpp(psi, h1, h2);
    residual = F - PDE_approx;
    err = sqrt(sum(residual(:).^2));

end

% Plot result
psi_p = -rho0*pi*6371E3*40/180*psi(2:end-1,2:end-1);
figure;
surf(plot_x, plot_p, psi_p','edgecolor','none'); 
view(2); 
axis tight; 
xlabel('Longitude','fontsize',15)
ylabel('p [hPa]','fontsize',15)
set(gca,'fontsize',15)
set(gca,'ydir','reverse'); 
title('Streamfunction [m^2/s]')
colorbar('fontsize',15)

u = -rho0*g* (psi(2:end-1,2:end)-psi(2:end-1,1:end-1))/(dp*100);
u = u(:,2:end);
figure;
surf(plot_x, plot_p, u','edgecolor','none'); 
view(2); 
axis tight; 
xlabel('Longitude','fontsize',15)
ylabel('p [hPa]','fontsize',15)
set(gca,'fontsize',15)
set(gca,'ydir','reverse'); 
title('u [m/s]')
colorbar('fontsize',15)

display(['Numerical [kg/s]: ', num2str(max(rho0*pi*6371E3*40/180*psi(:))/1E9)])
display(['Numerical [m^2/s]: ', num2str(max(psi(:)))])
psi_m2s = analyticalForm (plot_x, plot_p, p_max, L, H/sqrt(N2/(r*eps)), Q0, N2, r, eps);
psi_analytic = zeros(length(x)+2,length(p)+2);
psi_analytic(2:end-1,2:end-1) = -psi_m2s';
% 
% Q_numerical = -N2*cumtrapz(x(5:end-5)', PDE_approx(5:end-5,5:end-5), 1)*86400*T0/g;
% figure; surf(psi','edgecolor','none'); view(2); colorbar; cl = caxis; % caxis([0,5])
% analytic_approx = Dxx(psi_analytic, x_spacing) + vert_param*Dpp(psi_analytic, h1, h2);
% Q_analytic = -N2*cumtrapz(x(5:end-5)', analytic_approx(5:end-5,5:end-5), 1)*86400*T0/g;
% figure; surf(psi_analytic','edgecolor','none'); view(2); colorbar; caxis(cl)

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

function deriv = Dpp(X, h1, h2)

h2X = bsxfun(@times, h2, X(:,1:end-2));
h1X = bsxfun(@times, h1, X(:,3:end));
h12X = bsxfun(@times, h1+h2, X(:,2:end-1));
deriv = bsxfun(@rdivide, 2*(h2X + h1X - h12X), (h1.*h2 .* (h1+h2)));
deriv = deriv(2:end-1,:);

end

% ****************************************************
% Read in a field, and compute its meridional average.
% ****************************************************
function var = read_file(filename, varname)

var = ncread(filename,varname);
var = squeeze(mean(var(:,:,:,:),4)); % Average over time
var = squeeze(mean(var,2)); % Average meridionally

end
