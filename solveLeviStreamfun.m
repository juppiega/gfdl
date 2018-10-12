function [] = solveLeviStreamfun (filename)

pkg load netcdf
tdt_conv = read_file(filename,'tdt_conv');
tdt_ls = read_file(filename,'tdt_ls');
tdt_lw = read_file(filename,'tdt_lw');
T = read_file(filename,'temp');
x_spacing = 25E3;
x = ncread(filename,'grid_xt')*x_spacing;
y = ncread(filename,'grid_yt');
p = ncread(filename,'pfull')*100;
h1 = [p(1); diff(p)]';
h2 = [diff(p); p(end)-p(end-1)]';

figure;
[plot_x,plot_p] = meshgrid(x,p);
surf(plot_x, plot_p, 9.81*(tdt_ls+tdt_conv)'./T'); view(2); colorbar; axis tight; set(gca,'ydir','reverse'); title('LS + CONV 12-month')


vert_param = 1E-4;

psi = zeros(size(T,1)+2,size(T,2)+2);
Q = zeros(size(psi));
Q(2:end-1,2:end-1) = 9.81*(tdt_ls + tdt_conv)./T;

err = 1E30;
tol = 1E-11;
multiplier = -2/x_spacing^2 - 2*vert_param./(h1.*h2);
F = -Dx(Q, x_spacing);
while err > tol

    x_neighbors = Dxx(psi, x_spacing) + 2*psi(2:end-1,2:end-1)/x_spacing^2;
    p_neighbors = Dpp(psi, h1, h2) + 2*psi(2:end-1,2:end-1)./(h1.*h2);
    neighbor_points = x_neighbors + p_neighbors * vert_param;
    psi(2:end-1,2:end-1) = (F - neighbor_points) ./ multiplier;
    
    PDE_approx = Dxx(psi, x_spacing) + vert_param*Dpp(psi, h1, h2);
    residual = F - PDE_approx;
    err = sqrt(sum(residual(:).^2));
end

psi = psi(2:end-1,2:end-1);
figure;
surf(plot_x, plot_p, psi'); view(2); colorbar; axis tight; set(gca,'ydir','reverse'); title('Streamfunction LS + CONV + Newt. cool.')

end

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
deriv = 2*(h2.*X(:,1:end-2) + h1.*X(:,3:end) - (h1+h2).*X(:,2:end-1)) ./ (h1.*h2 .* (h1+h2));
deriv = deriv(2:end-1,:);

end

function var = read_file(filename, varname)

var = ncread(filename,varname);
var = squeeze(mean(var(:,:,:,:),4)); % Average over time
var = squeeze(mean(var,2)); % Average meridionally

end
