
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
  PI = get_pi("double")
  streamfun := RE * (north_lat - south_lat) * transpose(streamfun) * 180.0/PI / GRAV

  streamfun!0 = "lev"
  streamfun!1 = "lon"
  streamfun&lev = p
  streamfun&lon = lon

  return(streamfun)
  
end


function compute_divergent_streamfunction(south_lat, north_lat, filename)
local f, u, v, p, lat, lon, div, u_div, v_div, dims, wks, res, streamfun
begin

  f    = addfile (filename, "r")
  u    = f->ua_unmsk
  v    = f->va_unmsk
  p    = f->plev19
  lat  = f->lat
  lon  = f->lon

  div = uv2dvF_Wrap(u,v)                
  
  u_div    = new ( dimsizes(u), typeof(u), u@_FillValue )
  v_div    = new ( dimsizes(v), typeof(v), v@_FillValue )

  dv2uvf(div, u_div, v_div)

  streamfun = compute_streamfunction(u_div, p, south_lat, north_lat, lat, lon)

  dims = dimsizes(streamfun)
  num_levels = dims(0)

  wks   = gsn_open_wks ("png", str_concat((/filename, "_divergent_streamfuntion"/)) ) 
  res                 = True                        
  res@tiMainString    = "Divergent Walker Streamfunction" 
  res@cnFillOn             = True                  ; turn on color fill
  res@cnFillPalette        = "BlWhRe"              ; set color map 
  res@trYMinF = 70
         
  plot  = gsn_csm_pres_hgt(wks, streamfun, res ) 

  return(streamfun)

end

function compute_zonal_pert_streamfunction(south_lat, north_lat, filename)
local f, u, v, p, lat, lon, div, u_div, v_div, dims, wks, res, streamfun
begin

  f    = addfile (filename, "r")
  u    = f->ua_unmsk
  v    = f->va_unmsk
  p    = f->plev19
  lat  = f->lat
  lon  = f->lon            
  
  u_star = new ( dimsizes(u), typeof(u), u@_FillValue )
  
  dims = dimsizes(u)
  u_star = u - conform_dims(dimsizes(u), dim_avg_Wrap(u), (/0,1,2/))

  streamfun = compute_streamfunction(u_star, p, south_lat, north_lat, lat, lon)

  num_levels = dims(0)

  wks   = gsn_open_wks ("png", str_concat((/filename, "_zonal_pert_streamfuntion"/)) ) 
  res                 = True                        
  res@tiMainString    = "Walker Streamfunction from u*" 
  res@cnFillOn             = True                  ; turn on color fill
  res@cnFillPalette        = "BlWhRe"              ; set color map   
  res@trYMinF = 70
           
  plot  = gsn_csm_pres_hgt(wks, streamfun, res ) 

  return(streamfun)

end



begin
  filename = "atmos_cmip.1979-2014.djf.tmn.nc"
  
  streamfun_zon = compute_zonal_pert_streamfunction(-15, 15, filename)
  streamfun_div = compute_divergent_streamfunction(-15, 15, filename)
end



