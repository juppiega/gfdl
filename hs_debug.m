function [] = hs_debug (filename)
pkg load netcdf
T = ncread(filename,'temp');
u = ncread(filename,'ucomp');
w = ncread(filename,'omega');
h = ncread(filename,'height');
p = ncread(filename,'pres_full');

endi = 10;
%figure;for i=1:endi;var=squeeze(T(:,32,:,i));surf(flipud(var'));view(2);colorbar;title(num2str(i));pause(0.5);end
figure;for i=1:endi;var=squeeze(u(:,32,:,i));surf(flipud(var'));view(2);colorbar;title(num2str(i));pause(0.5);end
%figure;for i=1:endi;var=squeeze(w(:,32,:,i));surf(flipud(var'));view(2);colorbar;title(num2str(i));pause(0.5);end
figure;for i=1:endi;var=squeeze(h(:,32,:,i));var=var-mean(var,1);surf(flipud(var'));view(2);colorbar;title(num2str(i));pause(0.5);end
gp_grad = -9.80*(var(3:end,:) - var(1:end-2,:));
figure;surf(flipud(gp_grad'));view(2);colorbar;title('gp grad')

figure;for i=1:endi;var=squeeze(p(:,32,:,i));var=var-mean(var,1);surf(flipud(var'));view(2);colorbar;title(num2str(i));pause(0.5);end
p_grad = -(var(3:end,:) - var(1:end-2,:))/0.5;
figure;surf(flipud(p_grad'));view(2);colorbar;title('p grad')

figure;surf(flipud(p_grad'+gp_grad'));view(2);colorbar;title('total grad')

endfunction
