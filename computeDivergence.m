function div = computeDivergence(u_in, v_in, lat, lon)

u = zeros(size(u_in,1)+2,size(u_in,2)+2,size(u_in,3));
v = u;
u(2:end-1,2:end-1,:) = u_in;
v(2:end-1,2:end-1,:) = v_in;
dlat = abs(lat(2) - lat(1)) * pi/180;
dlon = (lon(2) - lon(1)) * pi/180;
RE = 6371E3;
clat = cosd(lat)';

div = Dlon(u,dlon)./(RE*clat);
div = div + Dlat(v(:,2:end-1,:).*clat,dlat)./(RE*clat);

end

function deriv = Dlat(X, dlat)
  deriv = zeros(size(X));
  deriv(:,2:end-1,:) = (X(:,3:end,:) - X(:,1:end-2,:)) / (2*dlat);
  deriv(:,end,:) = (X(:,end,:) - X(:,end-1,:)) / (2*dlat);
  deriv(:,1,:) = (X(:,2,:) - X(:,1,:)) / (2*dlat);
  deriv = deriv(2:end-1,:,:);
end

function deriv = Dlon(X, dlon)
  X(1,:,:) = X(end-1,:,:); X(end,:,:) = X(2,:,:);
  deriv = (X(3:end,:,:) - X(1:end-2,:,:)) / (2*dlon);
  deriv = deriv(:,2:end-1,:);
end