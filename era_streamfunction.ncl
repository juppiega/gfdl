
function compute_streamfunction(u, p_in, south_lat, north_lat, lat, lon)
local num_levels, u_div_meridional_mean, p, streamfun, dims, i, RE, GRAV, PI
begin
  u(:, :, ind(lat.lt.south_lat .or. lat.gt.north_lat), :) = 0.0
  u_div_meridional_mean = transpose(rm_single_dims(dim_avg_n(u, 2)))

  p = p_in(::-1)
  u_div_meridional_mean = u_div_meridional_mean(:,::-1)
  streamfun = new(dimsizes(u_div_meridional_mean), typeof(p), u@_FillValue)
  dims = dimsizes(streamfun)
  num_levels = dims(1)
  streamfun(:,0) = 0
  do i = 1, num_levels-1
    streamfun(:,i) = streamfun(:,i-1) + (p(i)-p(i-1)) * (u_div_meridional_mean(:,i) + u_div_meridional_mean(:,i-1)) / 2
  end do

  RE = 6371D3
  GRAV = 9.81
  PI = 3.14159265358979
  streamfun := RE * (north_lat - south_lat) * transpose(streamfun) * 180.0/PI / GRAV

  streamfun!0 = "lev"
  streamfun!1 = "lon"
  streamfun&lev = p
  streamfun&lon = lon

  system("/bin/rm -f era_streamfun.nc")
  ncdf = addfile("era_streamfun.nc" ,"c")
  ncdf->streamfun  = streamfun

  return(streamfun)
  
end


function compute_divergent_streamfunction(south_lat, north_lat, filename)
local f, u, v, p, lat, lon, div, u_div, v_div, dims, wks, res, streamfun
begin

  f    = addfile (filename, "r")
  div  = f->mean_div
  p    = f->level
  lat  = f->latitude
  lon  = f->longitude               
  
  u_div    = new ( dimsizes(div), typeof(div) )
  v_div    = new ( dimsizes(div), typeof(div) )

  ;print(dimsizes(div))

  dv2uvf(div, u_div, v_div)

  u_div!1 = "lev"
  u_div!2 = "lat"
  u_div!3 = "lon"
  u_div&lev = p
  u_div&lat = lat
  u_div&lon = lon

  system("/bin/rm -f era_udiv.nc")
  ncdf = addfile("era_udiv_1212.nc" ,"c")
  ncdf->u_div  = u_div

  streamfun = compute_streamfunction(u_div, p, south_lat, north_lat, lat, lon)

  dims = dimsizes(streamfun)
  num_levels = dims(0)

  wks   = gsn_open_wks ("png", str_concat((/filename, "_divergent_streamfuntion"/)) ) 
  res                 = True                        
  res@tiMainString    = "Divergent Walker Streamfunction" 
  res@cnFillOn             = True                  ; turn on color fill
  res@cnFillPalette        = "BlWhRe"              ; set color map 
  res@trYMinF = 70
         
  ;plot  = gsn_csm_pres_hgt(wks, streamfun, res ) 

  return(streamfun)

end


begin
  filename = "era_div_7918_1212.nc"
  
  streamfun_div = compute_divergent_streamfunction(-20, 20, filename)
end



