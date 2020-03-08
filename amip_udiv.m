function amip_udiv()

uFile = 'ua_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc';
vFile = 'va_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc';
u = ncread(uFile,'ua');
v = ncread(vFile,'va');
time = ncread(uFile,'time') + datenum('1979-01-01');
lat = ncread(uFile,'lat');
lon = ncread(uFile,'lon');
p = ncread(uFile,'plev');
Nmon = size(u,4);

for i = 1:12
    disp(['Solving month ',num2str(i)])
    ind = i:12:Nmon;
    outfile = ['amip_udiv_',num2str(i),'.nc'];
    if exist(outfile,'file')
        continue
    end
    umean = mean(u(:,:,:,ind),4);
    vmean = mean(v(:,:,:,ind),4);
    
    div = computeDivergence(umean,vmean,lat,lon);
    div(isnan(div)) = 0;
    solveDivWindNcl(div, lat, lon, p, outfile)
end

[el, la] = find_el_nino_months(time);

umean = mean(u(:,:,:,el),4);
vmean = mean(v(:,:,:,el),4);
div = computeDivergence(umean,vmean,lat,lon);
div(isnan(div)) = 0;
outfile = ['amip_udiv_elnino.nc'];
solveDivWindNcl(div, lat, lon, p, outfile)

umean = mean(u(:,:,:,la),4);
vmean = mean(v(:,:,:,la),4);
div = computeDivergence(umean,vmean,lat,lon);
div(isnan(div)) = 0;
outfile = ['amip_udiv_lanina.nc'];
solveDivWindNcl(div, lat, lon, p, outfile)


end

function solveDivWindNcl(div, lat, lon, lev, outfile)

delete(outfile);
nccreate(outfile,'u_div',...
       'Dimensions',{'longitude',length(lon),'latitude',length(lat),'level',length(lev),'time',1})
nccreate(outfile,'mean_div',...
       'Dimensions',{'longitude',length(lon),'latitude',length(lat),'level',length(lev),'time',1})

nccreate(outfile,'latitude','Dimensions',{'latitude', length(lat)});
nccreate(outfile,'longitude','Dimensions',{'longitude', length(lon)});
nccreate(outfile,'level','Dimensions',{'level', length(lev)}); 
nccreate(outfile,'time','Dimensions',{'time', 1}); 
ncwrite(outfile,'latitude',lat);
ncwrite(outfile,'longitude',lon);
ncwrite(outfile,'level',lev);
ncwrite(outfile,'time',1);

delete('era_div_in.nc')
copyfile(outfile, 'era_div_in.nc');

ncwrite('era_div_in.nc', 'mean_div', div);
stat = system('unset LD_LIBRARY_PATH && ncl era_streamfunction.ncl');
if stat
    error('Something wrong executing ncl');
end
u_div = ncread('era_udiv_out.nc','u_div');
ncwrite(outfile,'u_div',u_div)

end