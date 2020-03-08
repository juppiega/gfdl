function [] = solveDivergentWind (experiment_name, div, lat, lon, div_wind_name)


if ~isempty(experiment_name) && nargin == 1
    pkg load netcdf
    filename = ['atmos_average.',experiment_name,'.nc'];
    div_wind_name = ['div_wind.',experiment_name,'.mat'];
    div = read_file(filename, 'div');
    lon = ncread(filename,'lon');
    lat = ncread(filename,'lat'); 
end
if nargin < 5
    error('solveDivergentWind: too few arguments')
end
lon = lon * pi/180;
lat = lat * pi/180;
dlat = diff(lat); dlat = dlat(floor(length(dlat)/2));
dlon = lon(2) - lon(1);

u_div = zeros(size(div));
v_div = zeros(size(div));
velocity_pot = zeros(size(div));

lat_multip = zeros(1,length(lat)); lat_multip(1,:) = lat;

clat = repmat(cos(lat'),length(lon),1);
clat2 = clat.^2;
slat = repmat(sin(lat'),length(lon),1);

RE = 6371E3;

more off

for k = 1:size(div,3)
  
  psi = zeros(size(div,1)+2, size(div,2)+2);
  F_this = div(:,:,k) .* clat2 * (RE)^2;
  err = 1E30;
  tol = 1E-7;
  level = k
  
  while err > tol
    c1 = clat2 .* (D2lat(psi,dlat) + 2*psi(2:end-1,2:end-1)/dlat^2);
    sc = slat.*clat .* Dlat(psi,dlat);
    l1 = D2lon(psi,dlon) + 2*psi(2:end-1,2:end-1)/dlon^2;
    neighbour_points = c1 - sc + l1;
    multiplier = (-2*clat2/dlat^2 - 2/dlon^2);
    psi(2:end-1,2:end-1) = (F_this - neighbour_points) ./ multiplier;
    
    div_psi = (clat2.*D2lat(psi,dlat) - slat.*clat.*Dlat(psi,dlat) + D2lon(psi,dlon)) ./ (RE^2*clat2);
    residual = div(:,:,k) - div_psi;
    
    err = sqrt(sum(residual(clat>0.5).^2));
    %pause(0.1)
    
  end
  
  u_div(:,:,k) = Dlon(psi,dlon)./(RE*clat);
  v_div(:,:,k) = Dlat(psi,dlat)/RE;
  velocity_pot(:,:,k) = psi(2:end-1,2:end-1);
%  figure;
%  subplot(2,1,1)
%  surf(div(:,:,k),'edgecolor','none');
%  view(2);
%  colorbar;
%  axis tight
%  subplot(2,1,2)
%  surf(div_psi,'edgecolor','none');
%  view(2);
%  colorbar;
%  axis tight
  
end

save(div_wind_name,'u_div','v_div','velocity_pot');

end

function deriv = D2lat(X, dlat)
  deriv = (X(:,3:end) + X(:,1:end-2) - 2*X(:,2:end-1))/dlat^2;
  deriv = deriv(2:end-1,:);
end

function deriv = D2lon(X, dlon)
  X(1,:) = X(end-1,:); X(end,:) = X(2,:);
  deriv = (X(3:end,:) + X(1:end-2,:) - 2*X(2:end-1,:))/dlon^2;
  deriv = deriv(:,2:end-1);
end

function deriv = Dlat(X, dlat)
  deriv = (X(:,3:end) - X(:,1:end-2)) / (2*dlat);
  deriv = deriv(2:end-1,:);
end

function deriv = Dlon(X, dlon)
  X(1,:) = X(end-1,:); X(end,:) = X(2,:);
  deriv = (X(3:end,:) - X(1:end-2,:)) / (2*dlon);
  deriv = deriv(:,2:end-1);
end

function var = read_file(filename, varname)

var = ncread(filename,varname);
if size(var,4) == 2
  var = squeeze(var(:,:,:,end));
  return
end
var = squeeze(mean(var(:,:,:,end-5:end),4));

end