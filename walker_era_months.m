function walker_era_months(opt)

if isnumeric(opt) 
    result_name = ['era_div_7918_',sprintf('%d',opt),'.nc'];
else
    error('Only numeric month input, e.g. opt=[6,7,8]')
end
delete(result_name);

eraFiles = dir('adaptor.mars*');
div = ncread(eraFiles(1).name,'d',[1,1,1,1],[Inf,Inf,Inf,1]);
mean_div = zeros(size(div));
month_ind_cell = cell(length(eraFiles),1);
num_months = 0;
for i = 1:length(eraFiles)
    filename = eraFiles(i).name;
    time_this = datenum('1900-01-01')+double(ncread(filename,'time'))/24;
    [datestr(time_this(1)),' to ', datestr(time_this(end))]
    info = ncinfo(filename);
    N = info.Dimensions(4).Length;
    ind = 1:12:N; %ind(end) = [];
    month_ind = [];
    for j = 1:length(opt)
        month_ind = [month_ind, ind+opt(j)-1];
    end
    num_months = num_months + length(month_ind);
    month_ind_cell{i} = month_ind;
end

ntime = num_months/length(opt);
nccreate(result_name,'mean_div',...
       'Dimensions',{'longitude',1440,'latitude',721,'level',37,'time',ntime})
nccreate(result_name,'time','Dimensions',{'time', ntime}); 
nccreate(result_name,'months','Dimensions',{'months', length(opt)});
ncwrite(result_name,'months',opt);

year_ind = 1;
for i = 1:length(eraFiles)
    filename = eraFiles(i).name;
    month_ind = sort(month_ind_cell{i});
    for k = 1:length(opt):length(month_ind)
        for kk = 0:length(opt)-1
            d = ncread(filename,'d',[1,1,1,month_ind(k)+kk],[Inf,Inf,Inf,1]);
            mean_div = mean_div + d/length(opt);
        end
        ncwrite(result_name, 'mean_div', mean_div, [1,1,1,year_ind]);
        t = ncread(filename,'time',month_ind(k),1);
        ncwrite(result_name, 'time', t, year_ind);
        mean_div(:) = 0;
        year_ind
        year_ind = year_ind + 1;
    end
end
% figure;
% plot(1:length(mean_div_vec), mean_div_vec);

lat = double(ncread(filename,'latitude'));
lon = double(ncread(filename,'longitude'));
lev = double(ncread(filename,'level'));

nccreate(result_name,'latitude','Dimensions',{'latitude', length(lat)});
nccreate(result_name,'longitude','Dimensions',{'longitude', length(lon)});
nccreate(result_name,'level','Dimensions',{'level', length(lev)});
ncwrite(result_name,'latitude',lat);
ncwrite(result_name,'longitude',lon);
ncwrite(result_name,'level',lev);

end