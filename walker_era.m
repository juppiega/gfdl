function walker_era(months)

result_name = ['era_div_7918_',num2str(months,'%d%d%d'),'.nc'];
delete(result_name);

eraFiles = dir('adaptor.mars*');
div = ncread(eraFiles(1).name,'d',[1,1,1,1],[Inf,Inf,Inf,1]);
sum_div = zeros(size(div));
num_months = 0;
for i = 1:length(eraFiles)
    filename = eraFiles(i).name;
    info = ncinfo(filename);
    time_len = info.Dimensions(4).Length;
    ind = 1:12:time_len; %ind(end) = [];
    month_ind = [];
    for k = 1:length(months)
        month_ind = [month_ind, ind+months(k)-1];
    end
    for k = 1:length(month_ind)
        d = ncread(filename,'d',[1,1,1,month_ind(k)],[Inf,Inf,Inf,1]);
        sum_div = sum_div + d;
        num_months = num_months + 1;    
    end
end

mean_div = sum_div / num_months;
nccreate(result_name,'mean_div',...
           'Dimensions',{'longitude',1440,'latitude',721,'level',37,'time',1})
ncwrite(result_name, 'mean_div', mean_div);
lat = ncread(filename,'latitude');
lon = ncread(filename,'longitude');
lev = ncread(filename,'level');
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

