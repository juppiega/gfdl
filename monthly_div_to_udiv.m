function monthly_div_to_udiv(filename)

info = ncinfo(filename);
num_months = info.Dimensions(4).Length;

outfile = [filename(1:end-3),'_udiv.nc'];
delete(outfile);
nccreate(outfile,'u_div',...
       'Dimensions',{'longitude',1440,'latitude',721,'level',37,'time',num_months})
lat = double(ncread(filename,'latitude'));
lon = double(ncread(filename,'longitude'));
lev = double(ncread(filename,'level'));

nccreate(outfile,'latitude','Dimensions',{'latitude', length(lat)});
nccreate(outfile,'longitude','Dimensions',{'longitude', length(lon)});
nccreate(outfile,'level','Dimensions',{'level', length(lev)}); 
t = ncread(filename,'time'); mo = ncread(filename,'months');
nccreate(outfile,'time','Dimensions',{'time', length(t)});
nccreate(outfile,'months','Dimensions',{'months', length(mo)});
ncwrite(outfile,'latitude',lat);
ncwrite(outfile,'longitude',lon);
ncwrite(outfile,'level',lev);
ncwrite(outfile,'time',t);
ncwrite(outfile,'months',mo);

for i = 1:num_months
    d = ncread(filename,'mean_div',[1,1,1,i],[Inf,Inf,Inf,1]);
    ncwrite('era_div_in.nc', 'mean_div', d, [1,1,1,1]);
    i
    stat = system('unset LD_LIBRARY_PATH && ncl era_streamfunction.ncl');
    if stat
        error('Something wrong executing ncl');
    end
    u_div = ncread('era_udiv_out.nc','u_div',[1,1,1,1],[Inf,Inf,Inf,1]);
    ncwrite(outfile,'u_div',u_div,[1,1,1,i])
end

end