function walker_era(opt)

if isnumeric(opt) 
    result_name = ['era_div_7918_',sprintf('%d',opt),'.nc'];
else
    result_name = ['era_div_7918_',opt,'.nc'];
end
delete(result_name);

eraFiles = dir('adaptor.mars*');
div = ncread(eraFiles(1).name,'d',[1,1,1,1],[Inf,Inf,Inf,1]);
sum_div = zeros(size(div));
num_months = 0;
for i = 1:length(eraFiles)
    filename = eraFiles(i).name;
    if isnumeric(opt)
        info = ncinfo(filename);
        N = info.Dimensions(4).Length;
        ind = 1:12:N; %ind(end) = [];
        month_ind = [];
        for j = 1:length(opt)
            month_ind = [month_ind, ind+opt(j)-1];
        end
    else
        t_data = ncread(filename,'time');
        timestamps = double(datenum('1900-01-01') + t_data/24);
        [el, la] = find_el_nino_months(timestamps);
        if strcmpi(opt,'elnino')
            month_ind = find(el);
        else
            month_ind = find(la);
        end
    end
    for k = 1:length(month_ind)
        d = ncread(filename,'d',[1,1,1,month_ind(k)],[Inf,Inf,Inf,1]);
        sum_div = sum_div + d;
        num_months = num_months + 1;    
    end
end

mean_div = sum_div / num_months;

aver_lat = 20;
lat = double(ncread(filename,'latitude'));
lon = double(ncread(filename,'longitude'));
lev = double(ncread(filename,'level'));
ind = -aver_lat <= lat & lat <= aver_lat;
div = squeeze(mean(mean_div(:,ind,:,:),2));
z = compute_standard_height(lev);
[LON,Z] = meshgrid(lon, z);
figure;
colormap(jet)
contourf(LON,Z,div',20,'edgecolor','none')
%contour(LON,Z,streamfun')
view(2)
colorbar
axis tight
m = 1E-5;
caxis([-m,m])
yl = get(gca, 'ylim');
ylim([yl(1),15])
title('divergence_direct')

nccreate(result_name,'mean_div',...
           'Dimensions',{'longitude',1440,'latitude',721,'level',37,'time',1})
ncwrite(result_name, 'mean_div', mean_div);
nccreate(result_name,'latitude','Dimensions',{'latitude', length(lat)});
nccreate(result_name,'longitude','Dimensions',{'longitude', length(lon)});
nccreate(result_name,'level','Dimensions',{'level', length(lev)});
ncwrite(result_name,'latitude',lat);
ncwrite(result_name,'longitude',lon);
ncwrite(result_name,'level',lev);

end

function [data,done] = load_variable(filename)
persistent previous_month

if isempty(previous_month)
    previous_month = 0;
end

data = ncread(filename, 'd', [1,1,1,previous_month+1], [Inf,Inf,Inf,1]);
info = ncinfo(filename);
done = previous_month+1 == info.Dimensions(4).Length;
previous_month = previous_month + 1;

end

