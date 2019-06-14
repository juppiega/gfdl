function [psi] = solveLeviStreamfun (filename, x_spacing, vert_param, heating_fields)
% psi = solveLeviStreamfun (filename, x_spacing, vert_param, heating_fields)
%
% Purpose: Solves the Boussinesq streamfunction (psi), given diabatic heating (Q):
%          d_xx(psi) + vert_param * d_zz(psi) = -d_x(Q)/N^2                   (1)
%          Algorithm solves the Poisson equation using the Jacobi finite difference method.
%          Assumes periodic boundaries in x, psi(p=0) = 0 and d(psi(surface))/dp = 0;
%
% Arguments:
% filename       : Name of the channel model output file to be read.
% x_spacing      : Distance between two grid points in the x-direction [meters].
% vert_param     : (optional, Default: 1E-4) 
%                  Magnitude of coupling between vertical levels.
%                  ( = Newton_rate * Rayleigh_rate * reference_density^2 * gravity^2 / N^2)
%                  Typically of order 1E-4. Exactly zero for weak temperature gradient.
% heating_fields : (optional, Default: {'tdt_ls', 'tdt_conv'}) 
%                  Names of the fields to be used for heating (Q).
%
% Output:
% psi            : The computed Boussinesq streamfunction [m^2/s].            

% Input checks
if nargin == 0 || nargin == 1
    error('FILENAME and X_SPACING have to be specified. See "> help solveLeviStreamfun" for more information')
end
if nargin <= 2
    vert_param = 1E-4;
end
if nargin <= 3
    heating_fields = {'tdt_ls', 'tdt_conv'};
end

% Gravity and temperature for converting Q from K/s to m/s^3
g = 9.81;
T = read_file(filename,'temp'); % read_file also computes the meridional averages.

% Initialize fields psi and Q
psi = zeros(size(T,1)+2,size(T,2)+2); % Add +2 for the boundary points.
Q = zeros(size(psi));

% Read in the heating fields.
for i = 1:length(heating_fields)
    heating_field = read_file(filename, heating_fields{i});
    Q(2:end-1,2:end-1) = Q(2:end-1,2:end-1) + g*heating_field./T;
end

% Coordinates x and pressure
x = ncread(filename,'grid_xt')*x_spacing;
p = ncread(filename,'pfull')*100;

% Plot Q
figure;
[plot_x,plot_p] = meshgrid(x/1E3,p/100);
surf(plot_x, plot_p, Q(2:end-1,2:end-1)','edgecolor','none'); view(2); colorbar; axis tight; set(gca,'ydir','reverse'); title('Q 12-month mean')
xlabel('x [km]','fontsize',15)
ylabel('p [hPa]','fontsize',15)
set(gca,'fontsize',15)

u = read_file(filename,'ucomp');
streamfun = cumtrapz(p,u,2)*2*pi*6371E3/9.8;
figure;
surf(plot_x, plot_p, -streamfun','edgecolor','none'); view(2); colorbar; axis tight; set(gca,'ydir','reverse'); title('Streamfunction [kg/s]')
xlabel('x [km]','fontsize',15)
ylabel('p [hPa]','fontsize',15)
set(gca,'fontsize',15)

theta = T'.*(1000./plot_p).^(2/7);
rho = plot_p*100./(287*T');
theta_grad = (theta(2:end,:) - theta(1:end-1,:))./(plot_p(2:end,:) - plot_p(1:end-1,:))/100;
stability = -9.81^2 * rho(2:end,:) .* theta_grad ./ theta(2:end,:);
figure;
surf(plot_x(2:end,:), plot_p(2:end,:), stability,'edgecolor','none'); view(2); colorbar; axis tight; set(gca,'ydir','reverse'); title('N^2 [s^{-2}]')
xlabel('x [km]','fontsize',15)
ylabel('p [hPa]','fontsize',15)
caxis([0, 2E-4])
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
N2 = 1E-4;
F = -Dx(Q, x_spacing)/N2;
% Solver error tolerance. Iteration stops, when the rms error between LHS and RHS in (1) is smaller than tol.
err = 1E30;
tol = 1E-7;

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
psi = psi(2:end-1,2:end-1);
figure;
surf(plot_x, plot_p, psi','edgecolor','none'); 
view(2); 
axis tight; 
xlabel('x [km]','fontsize',15)
ylabel('p [hPa]','fontsize',15)
set(gca,'fontsize',15)
set(gca,'ydir','reverse'); 
title('Streamfunction [m^2/s]')
colorbar

end

% *********************
% Derivative functions 
% *********************
function deriv = Dx(X, dx)

X(1,:) = X(end-1,:);
X(end,:) = X(2,:);
deriv = (X(3:end,:) - X(1:end-2,:))/(2*dx);
deriv = deriv(:,2:end-1);

end

function deriv = Dxx(X, dx)

X(1,:) = X(end-1,:);
X(end,:) = X(2,:);
deriv = (X(1:end-2,:) + X(3:end,:) - 2*X(2:end-1,:))/dx^2;
deriv = deriv(:,2:end-1);

end

function deriv = Dpp(X, h1, h2)

X(:,end) = X(:,end-1);
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
