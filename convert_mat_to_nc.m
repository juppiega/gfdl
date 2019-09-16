function convert_mat_to_nc()

result_name = 'era_uv_elnino.nc';
load('ERA5_Q_elnino.mat')
delete(result_name)
nccreate(result_name,'u',...
           'Dimensions',{'longitude',1440,'latitude',721,'level',length(p),'time',1})
ncwrite(result_name, 'u', u);
nccreate(result_name,'v',...
           'Dimensions',{'longitude',1440,'latitude',721,'level',length(p),'time',1})
ncwrite(result_name, 'v', v);
nccreate(result_name,'latitude','Dimensions',{'latitude', length(lat)});
nccreate(result_name,'longitude','Dimensions',{'longitude', length(lon)});
nccreate(result_name,'level','Dimensions',{'level', length(p)});
ncwrite(result_name,'latitude',lat);
ncwrite(result_name,'longitude',lon);
ncwrite(result_name,'level',p);

result_name = 'era_uv_lanina.nc';
load('ERA5_Q_lanina.mat')
delete(result_name)
nccreate(result_name,'u',...
           'Dimensions',{'longitude',1440,'latitude',721,'level',length(p),'time',1})
ncwrite(result_name, 'u', u);
nccreate(result_name,'v',...
           'Dimensions',{'longitude',1440,'latitude',721,'level',length(p),'time',1})
ncwrite(result_name, 'v', v);
nccreate(result_name,'latitude','Dimensions',{'latitude', length(lat)});
nccreate(result_name,'longitude','Dimensions',{'longitude', length(lon)});
nccreate(result_name,'level','Dimensions',{'level', length(p)});
ncwrite(result_name,'latitude',lat);
ncwrite(result_name,'longitude',lon);
ncwrite(result_name,'level',p);

end