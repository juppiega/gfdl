function [] = debug_simulation (filename, varname, surface_vals, start_i, end_i)

pkg load netcdf
if strcmpi(varname,'Q_approx')
  w = ncread(filename,'omega');
  B = ncread(filename,'temp');
  var = w*1E-4 + B / (15*86400);
elseif strcmpi(varname,'ps')
  var = ncread(filename,varname);
  var = var - mean(var,1);
elseif strcmpi(varname,'p_pert')
  B = ncread(filename,'temp');
  p = ncread(filename,'pfull');
  var = -cumtrapz(p, B, 3)/9.81;
  %var = var - var(:,:,end,:);
  %var = var - mean(var,1);
else
  var = ncread(filename,varname);
  if surface_vals
    var = squeeze(var(:,:,end-1,:));
  end
  var = var - mean(var,1);
end



if nargin <= 3
    start_i = 1;
end
if nargin <= 4
    end_i = size(var,length(size(var)));
end

lon = ncread(filename,'lon');
lat = ncread(filename,'lat');
try
  p = ncread(filename,'pfull');
  [X_p,p] = meshgrid(lon,p);
catch

end
[X_lat,Y] = meshgrid(lon,lat);

eq = round(size(var,2)/2);

figure;

for k = start_i:end_i
    if strcmpi(varname,'ps')
      surf(X_lat,Y,squeeze(var(:,:,k))', 'edgecolor','none');
      %m = max(abs(var(:)));
      %caxis([-m, m]);
    end
    
    if surface_vals
      Z = var(:,:,k)';
      surf(X_lat, Y, Z, 'edgecolor','none');
    else
      surf(X_p, p, squeeze(var(:,eq,:,k))', 'edgecolor','none');
    end
    colorbar
    view(2)
    axis tight
    title(num2str(k));
    set(gca,'ydir','reverse')
    
    pause(0.5)
end

if strcmpi(varname,'ps')
  figure;
  plot(lon, var(:,eq,end),'linewidth',2.0)
  title('p_s equator')
  axis tight
  hold all
  plot([180 180],get(gca,'ylim'),'linestyle','--')
end

if surface_vals && strcmpi(varname,'ucomp')
  figure;
  plot(lon, var(:,eq,end),'linewidth',2.0)
  title('Surface u at equator')
  axis tight
  hold all
  plot([180 180],get(gca,'ylim'),'linestyle','--')
end

end
