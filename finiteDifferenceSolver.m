function [] = finiteDifferenceSolver()

accuracy = 1E-5;
maxIter = 5000;
global b; % Num of boundary points
b = 3;

Lx = 1;
Ly = Lx;
ztop = 2;
H = ztop/10;%2E3;
z0 = ztop/2;
Q0 = 1E-5;
eps = 1;
r = 0;
beta = 1;%2.3E-11;
Nsq = 1E-4;
hx = Lx/10;%100E3;
hz = ztop/40;
x = -7*Lx:hx:7*Lx;
y = -5*Ly:hx:5*Ly;
z = 0:hz:ztop;
epsBetaYsq = nan(length(y)-6,1,1);
epsBetaYsq(:) = (eps/beta)^2 + y(4:end-3).^2;

Q = zeros(length(y),length(x),length(z));
for i = 1:length(y)
    for j = 1:length(x)
        for k = 1:length(z)
            Q(i,j,k) = x(j)^4 + y(i)^4;%Q0 * exp(-(x(j)/Lx)^2 - (y(i)/Ly)^2 - ((z(k)-z0)/H)^2);
        end
    end
end

%F = (bsxfun(@times,Lap_Dzz(Q,hx,hz),epsBetaYsq) + 3*(eps/beta)*Dxzz(Q,hx,hz) + ...
%    5*yDyzz(Q,hx,hz,y) + 4*Dzz(Q,hz)) / Nsq;
%F = eps/Nsq * Dzz(Q,hz);
F = ones(size(Biharmonic(Q,hx)));
Zplot = permute(F(50,:,:),[3,2,1]);
figure; surf(Zplot(:,:,1),'edgecolor','none');view(2);axis tight;colorbar

Zplot = permute(Q(50,:,:),[3,2,1]);
figure; surf(Zplot(:,:,1),'edgecolor','none');view(2);axis tight;colorbar;



B = zeros(length(y),length(x),length(z))/10;
B(b+1:size(B,1)-b,b+1:size(B,2)-b,2:size(B,3)-1) = 1; ind_internal = find(B);
B(:) = 0;%randn(length(y),length(x),length(z))/100;
residual = 1E30;
B_updated = B;
divider = (640*eps/(108*beta^2 * hx^4) + r/Nsq * 8*epsBetaYsq / (hx^2 * hz^2) - r/Nsq * 8/hz^2);
%divider = -4/hx^2 - 2*r*eps/Nsq/hz^2;
numIter = 0;
while residual > accuracy && numIter < maxIter
     coeff_bi = Biharmonic(B,hx) - 640*B(b+1:size(B,1)-b,b+1:size(B,2)-b,2:size(B,3)-1) / (108*hx^4);
     coeff_lap = Lap_Dzz(B,hx,hz) - 8*B(b+1:size(B,1)-b,b+1:size(B,2)-b,2:size(B,3)-1) / (hx^2 * hz^2);
     coeff_vertical = Dzz(B,hz) + 2*B(b+1:size(B,1)-b,b+1:size(B,2)-b,2:size(B,3)-1) / hz^2;
%     
     neighbor_points = eps*coeff_bi/beta^2;% + ... %2*Lap_Dx(B,hx)/beta + ...
         %r/Nsq * (bsxfun(@times,coeff_lap,epsBetaYsq) + 3*(eps/beta)*Dxzz(B,hx,hz) + ...
         %5*yDyzz(B,hx,hz,y) + 4*coeff_vertical);

%    coeff_lap = Laplacian(B,hx) + 4*B(3:size(B,1)-2,3:size(B,2)-2,2:size(B,3)-1)/hx^2;
%    coeff_vertical = Dzz(B,hz) + 2*B(3:size(B,1)-2,3:size(B,2)-2,2:size(B,3)-1)/hz^2;
%    neighbor_points = coeff_lap + r*eps/Nsq * coeff_vertical;
    
    B_updated = bsxfun(@rdivide,(F - neighbor_points),divider);
    
    %err = eps/beta^2 * Biharmonic(B,hx) - F;
    
    residual = mean(abs(B(ind_internal)-B_updated(:)))
    B(ind_internal) = B_updated;
    
    numIter = numIter + 1
    if residual > 10
        break;
    end
end
numIter

err = Biharmonic(B,hx) - F;

Zplot = permute(err(50,:,:),[3,2,1]);
figure; surf(Zplot(:,:,1),'edgecolor','none');view(2);axis tight;colorbar

Zplot = permute(B(50,:,:),[3,2,1]);
figure; surf(Zplot(:,:,1),'edgecolor','none');view(2);axis tight;colorbar

end

function [derivative] = Dx(X,hx)

derivative = (X(:,3:end,:) - X(:,1:end-2,:)) / (2*hx);

end 

function [derivative] = Dy(X,hx)

derivative = (X(3:end,:,:) - X(1:end-2,:,:)) / (2*hx);

end

function [derivative] = Dzz_internal(X,hz)

derivative = (X(:,:,3:end) - 2*X(:,:,2:end-1) + X(:,:,1:end-2)) / (hz*hz);

end

function [lap] = Laplacian_internal(X,hx)

lap = (X(2:end-1,3:end,:) + X(2:end-1,1:end-2,:) + X(1:end-2,2:end-1,:) + ...
    X(3:end,2:end-1,:) - 4*X(2:end-1,2:end-1,:)) / (hx*hx);

end

function lap = Laplacian(X,hx)

global b
lap = Laplacian_internal(X,hx);
lap = lap(b:end-b+1,b:end-b+1,2:end-1);

end

function [derivative] = Dzz(X,hz)

global b
derivative = Dzz_internal(X,hz);
derivative = derivative(b+1:end-b,b+1:end-b,:);

end

function [derivative] = Dxzz(X,hx,hz)

global b
derivative = Dx(Dzz_internal(X,hz),hx);
derivative = derivative(b+1:end-b,b:end-b+1,:);

end

function [derivative] = yDyzz(X,hx,hz,y)

global b
derivative = Dy(Dzz_internal(X,hz),hx);
derivative = derivative(b:end-b+1,b+1:end-b,:);
Y_cube = nan(length(y)-6,1,1);
Y_cube(:) = y(b+1:end-b);
derivative = bsxfun(@times,derivative,Y_cube);

end

function bi = Biharmonic(X,hx)

bi = (640 * X(4:end-3,4:end-3,:) ...
     -177 * (X(3:end-4,4:end-3,:) + X(5:end-2,4:end-3,:) + X(4:end-3,3:end-4,:) + X(4:end-3,5:end-2,:)) +...
       11 * (X(1:end-6,4:end-3,:) + X(7:end,4:end-3,:) + X(4:end-3,1:end-6,:) + X(4:end-3,7:end,:)) + ...
        3 * (X(2:end-5,7:end,:) + X(6:end-1,7:end,:) + X(2:end-5,1:end-6,:) + X(6:end-1,1:end-6,:) +...
             X(1:end-6,2:end-5,:) + X(1:end-6,6:end-1,:) + X(7:end,2:end-5,:) + X(7:end,6:end-1,:))) / (108*hx^4);
bi = bi(:,:,2:end-1);

end

function derivative = Lap_Dzz(X,hx,hz)

global b
derivative = Laplacian_internal(Dzz_internal(X,hz),hx);
derivative = derivative(b:end-b+1,b:end-b+1,:);

end

function derivative = Lap_Dx(X,hx)

global b
derivative = Dx(Laplacian_internal(X,hx),hx);
derivative = derivative(b:end-b+1,b-1:end-b+2,2:end-1);

end