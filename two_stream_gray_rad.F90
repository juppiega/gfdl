!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cloud_mod

   use sat_vapor_pres_mod,    only:  compute_qs

   use fms_mod,               only: open_file, check_nml_error, &
                                    mpp_pe, close_file, error_mesg, FATAL

   use constants_mod,         only: stefan, cp_air, grav, pstd_mks

   use    diag_manager_mod,   only: register_diag_field, send_data

   use    time_manager_mod,   only: time_type, &
                                    operator(+), operator(-), operator(/=)

implicit none
private

public :: cloudfrac, cloud_init, cloud_end, compute_clouds

character(len=10), parameter :: mod_name = 'cloud'
logical :: initialized     = .false.

real :: crit_rh_aloft = 0.7
real :: crit_rh_surf = 0.9
real :: crit_rh_shape = 4.0
real :: liq_dens_surf = 0.21E-3 ! Cloud water at the surface [kg / m^3]
real :: ref_liquid_scale_height = 700 ! Cloud water reference scale height [m]
real :: reff_liq = 10E-6 ! Liquid drop effective radius [m]
real :: g_liq = 0.8 ! Asymmetry parameter
real :: k_liq = 90 ! Absorption coefficient [m^2/kg(absorber)]

namelist/cloudrad_nml/ crit_rh_aloft, crit_rh_surf, crit_rh_shape, &
                       liq_dens_surf, ref_liquid_scale_height, reff_liq

!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_cloudfrac, id_lwp, id_tau, id_lw_emiss
real, allocatable, dimension(:,:,:) :: cloudfrac, lwp, tau, reff, lw_emiss

real :: missing_value = -999.

contains

!============================================================================================
subroutine cloud_init(is, ie, js, je, num_levels, axes, Time)
implicit none
!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
!-------------------------------------------------------------------------------------

allocate (cloudfrac          (ie-is+1, je-js+1, num_levels))
allocate (lwp          (ie-is+1, je-js+1, num_levels))
allocate (tau          (ie-is+1, je-js+1, num_levels))
allocate (reff          (ie-is+1, je-js+1, num_levels))
allocate (lw_emiss          (ie-is+1, je-js+1, num_levels))

id_cloudfrac = &
        register_diag_field ( mod_name, 'cloudfrac', axes(1:3), Time, &
               'Cloud fraction', &
               '[0-1]', missing_value=missing_value               )

id_lwp = &
        register_diag_field ( mod_name, 'lwp', axes(1:3), Time, &
               'Liquid water path', &
               'kg/m^2', missing_value=missing_value               )

id_tau = &
        register_diag_field ( mod_name, 'tau', axes(1:3), Time, &
               'Cloud extinction optical depth', &
               'unitless', missing_value=missing_value               )

id_lw_emiss = &
        register_diag_field ( mod_name, 'lw_emiss', axes(1:3), Time, &
               'Longwave emissivity', &
               'unitless', missing_value=missing_value               )

initialized = .true.

end subroutine

!============================================================================================
subroutine compute_clouds(tg, sphum, pfull, phalf, z_full, z_half, Time_diag)
implicit none
real, intent(in), dimension(:,:,:) :: tg, sphum, pfull, phalf, z_full, z_half
type(time_type), intent(in)         :: Time_diag

call compute_cloudfrac(tg, sphum, pfull, z_full, Time_diag)
call compute_cloud_water_path(sphum, phalf, z_half, Time_diag)
reff = reff_liq ! TODO: ICE MODIFY
call compute_tau(Time_diag)
call compute_emissivity(Time_diag)

end subroutine

!============================================================================================
subroutine compute_emissivity(Time_diag)
implicit none
real, parameter :: D = 1.66 ! Diffusivity factor
logical :: used
type(time_type), intent(in)         :: Time_diag

lw_emiss = 1 - exp(-D*k_liq*lwp) ! TODO: ICE MODIFY

if ( id_lw_emiss > 0 ) then
   used = send_data ( id_lw_emiss, lw_emiss, Time_diag)
endif

end subroutine

!============================================================================================
subroutine compute_tau(Time_diag)
implicit none
type(time_type), intent(in)         :: Time_diag

real, parameter :: rho_water = 1E3
logical :: used

tau = 1.5 * lwp / (reff*rho_water) ! TODO: ICE MODIFY

if ( id_tau > 0 ) then
   used = send_data ( id_tau, tau, Time_diag)
endif

end subroutine

!============================================================================================
subroutine compute_cloud_water_path(sphum, phalf, z_half, Time_diag)
implicit none
real, intent(in), dimension(:,:,:) :: sphum, phalf, z_half
type(time_type), intent(in)         :: Time_diag

integer :: i,j,k,Nlev
real, dimension(size(sphum,3)) :: q_col, cwp_col
real, dimension(size(phalf,3)) :: phalf_col, zhalf_col
real :: water_integral, hl
logical :: used

Nlev = size(sphum,3)
do j = 1,size(sphum,2)
  do i = 1,size(sphum,1)

    phalf_col = phalf(i,j,:)
    zhalf_col = z_half(i,j,:)
    water_integral = 0.0
    q_col = sphum(i,j,:)
    do k = 1,Nlev
      water_integral = water_integral + q_col(k)*(phalf_col(k+1)-phalf_col(k))/GRAV
    end do

    hl = ref_liquid_scale_height * log(1.0 + water_integral)
    cwp_col = hl*liq_dens_surf * (exp(-zhalf_col(2:Nlev+1)/hl) - exp(-zhalf_col(1:Nlev)/hl))
    lwp(i,j,:) = cwp_col

  end do
end do

if ( id_lwp > 0 ) then
   used = send_data ( id_lwp, lwp, Time_diag)
endif

end subroutine

!============================================================================================
subroutine compute_cloudfrac(tg, sphum, pfull, z_full, Time_diag)
implicit none
real, intent(in), dimension(:,:,:) :: tg, sphum, pfull, z_full
type(time_type), intent(in)         :: Time_diag
! Locals
integer :: i,j
logical :: used
real, dimension(size(tg,3)) :: tg_col, rh_col, rh_crit_col, q_col, qsat_col, pfull_col

do j = 1,size(tg,2)
  do i = 1,size(tg,1)

    tg_col = tg(i,j,:)
    q_col = sphum(i,j,:)
    pfull_col = pfull(i,j,:)

    call compute_qs (tg_col, pfull_col, qsat_col)
    rh_col = q_col / qsat_col

    rh_crit_col = crit_rh_aloft + (crit_rh_surf - crit_rh_aloft) * exp(1.0 - (pstd_mks/pfull_col)**crit_rh_shape)

    cloudfrac(i,j,:) = max(0.0, 1.0 - sqrt(1.0 - min(1.0, (rh_col-rh_crit_col)/(1.0-rh_crit_col))))

  end do
end do

if ( id_cloudfrac > 0 ) then
   used = send_data ( id_cloudfrac, cloudfrac, Time_diag)
endif

end

!============================================================================================
subroutine cloud_end()
implicit none

  deallocate(cloudfrac)

end subroutine

!============================================================================================
end module

module two_stream_gray_rad_mod

! ==================================================================================
! ==================================================================================

   use fms_mod,               only: open_file, check_nml_error, &
                                    mpp_pe, close_file

   use constants_mod,         only: stefan, cp_air, grav, pstd_mks

   use    diag_manager_mod,   only: register_diag_field, send_data

   use    time_manager_mod,   only: time_type, &
                                    operator(+), operator(-), operator(/=)
   use cloud_mod
 
!==================================================================================
implicit none
private
!==================================================================================

! version information 
 
character(len=128) :: version='$Id: two_stream_gray_rad.F90,v 1.1.2.1 2013/01/24 14:44:49 pjp Exp $'
character(len=128) :: tag='Two-stream gray atmosphere'

!==================================================================================

! public interfaces

public :: two_stream_gray_rad_init, two_stream_gray_rad_down, two_stream_gray_rad_up, two_stream_gray_rad_end, &
          compute_radiation_with_clouds
!==================================================================================


! module variables
logical :: initialized     = .false.
real    :: solar_constant  = 1360.0
real    :: del_sol         = 1.4
real    :: del_sw          = 0.0
real    :: ir_tau_eq       = 6.0
real    :: ir_tau_pole     = 1.5
real    :: atm_abs         = 0.0
real    :: sw_diff         = 0.0
real    :: linear_tau      = 0.1
real    :: albedo_value    = 0.06
real    :: wv_exponent     = 4.0 
real    :: solar_exponent  = 4.0 

real, allocatable, dimension(:,:)   :: insolation, p2, albedo, lw_tau_0, sw_tau_0
real, allocatable, dimension(:,:)   :: b_surf
real, allocatable, dimension(:,:,:) :: b, tdt_rad, tdt_solar
real, allocatable, dimension(:,:,:) :: lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux 
real, allocatable, dimension(:,:,:) :: lw_tau, sw_tau, lw_dtrans
real, allocatable, dimension(:,:)   :: olr, net_lw_surf, toa_sw_in

real, save :: pi, deg_to_rad , rad_to_deg

namelist/two_stream_gray_rad_nml/ solar_constant, del_sol, &
           ir_tau_eq, ir_tau_pole, atm_abs, sw_diff, &
           linear_tau, del_sw, albedo_value, wv_exponent, &
           solar_exponent
        
!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_olr, id_swdn_sfc, id_swdn_toa, id_net_lw_surf, id_lwdn_sfc, id_lwup_sfc, &
           id_tdt_rad, id_tdt_solar, id_flux_rad, id_flux_lw, id_flux_sw, id_swdn_sfc_cs, &
           id_swup_toa, id_swup_toa_cs

character(len=10), parameter :: mod_name = 'two_stream'

real :: missing_value = -999.


contains


! ==================================================================================
! ==================================================================================


subroutine two_stream_gray_rad_init(is, ie, js, je, num_levels, axes, Time)

!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/)
integer :: ierr, io, unit
!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

unit = open_file ('input.nml', action='read')
ierr=1
do while (ierr /= 0)
   read  (unit, nml=two_stream_gray_rad_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'two_stream_gray_rad_nml')
enddo
10 call close_file (unit)

unit = open_file ('logfile.out', action='append')
if ( mpp_pe() == 0 ) then
  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
  write (unit, nml=two_stream_gray_rad_nml)
endif
call close_file (unit)

pi         = 4. * atan(1.)
deg_to_rad = pi/180.
rad_to_deg = 180./pi

call cloud_init(is, ie, js, je, num_levels, axes, Time)

initialized = .true.

allocate (b                (ie-is+1, je-js+1, num_levels))
allocate (tdt_rad          (ie-is+1, je-js+1, num_levels))
allocate (tdt_solar        (ie-is+1, je-js+1, num_levels))

allocate (lw_dtrans        (ie-is+1, je-js+1, num_levels))
allocate (lw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (sw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (lw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (lw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (lw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (sw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (rad_flux         (ie-is+1, je-js+1, num_levels+1))

allocate (b_surf           (ie-is+1, je-js+1))
allocate (lw_tau_0         (ie-is+1, je-js+1))
allocate (sw_tau_0         (ie-is+1, je-js+1))
allocate (olr              (ie-is+1, je-js+1))
allocate (net_lw_surf      (ie-is+1, je-js+1))
allocate (toa_sw_in        (ie-is+1, je-js+1))

allocate (insolation       (ie-is+1, je-js+1))
allocate (albedo           (ie-is+1, je-js+1))
allocate (p2               (ie-is+1, je-js+1))


!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------

    id_olr = &
    register_diag_field ( mod_name, 'olr', axes(1:2), Time, &
               'outgoing longwave radiation', &
               'W/m2', missing_value=missing_value               )
    id_swdn_sfc_cs = &
    register_diag_field ( mod_name, 'swdn_sfc_cs', axes(1:2), Time, &
               'Absorbed SW at surface (clear sky)', &
               'W/m2', missing_value=missing_value               )
    id_swdn_sfc = &
    register_diag_field ( mod_name, 'swdn_sfc', axes(1:2), Time, &
               'Absorbed SW at surface', &
               'W/m2', missing_value=missing_value               )
    id_swdn_toa = &
    register_diag_field ( mod_name, 'swdn_toa', axes(1:2), Time, &
               'SW flux down at TOA', &
               'W/m2', missing_value=missing_value               )
    id_swup_toa = &
    register_diag_field ( mod_name, 'swup_toa', axes(1:2), Time, &
               'SW flux up at TOA', &
               'W/m2', missing_value=missing_value               )
    id_swup_toa_cs = &
    register_diag_field ( mod_name, 'swup_toa_cs', axes(1:2), Time, &
               'SW flux up at TOA (clear sky)', &
               'W/m2', missing_value=missing_value               )
    id_lwup_sfc = &
    register_diag_field ( mod_name, 'lwup_sfc', axes(1:2), Time, &
               'LW flux up at surface', &
               'W/m2', missing_value=missing_value               )

    id_lwdn_sfc = &
    register_diag_field ( mod_name, 'lwdn_sfc', axes(1:2), Time, &
               'LW flux down at surface', &
               'W/m2', missing_value=missing_value               )

    id_net_lw_surf = &
    register_diag_field ( mod_name, 'net_lw_surf', axes(1:2), Time, &
               'Net upward LW flux at surface', &
               'W/m2', missing_value=missing_value               )

    id_tdt_rad = &
        register_diag_field ( mod_name, 'tdt_rad', axes(1:3), Time, &
               'Temperature tendency due to radiation', &
               'K/s', missing_value=missing_value               )

    id_tdt_solar = &
        register_diag_field ( mod_name, 'tdt_solar', axes(1:3), Time, &
               'Temperature tendency due to solar radiation', &
               'K/s', missing_value=missing_value               )

    id_flux_rad = &
        register_diag_field ( mod_name, 'flux_rad', axes(half), Time, &
               'Total radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_lw = &
        register_diag_field ( mod_name, 'flux_lw', axes(half), Time, &
               'Net longwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw = &
        register_diag_field ( mod_name, 'flux_sw', axes(half), Time, &
               'Net shortwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )

return
end subroutine two_stream_gray_rad_init

! ==================================================================================

subroutine two_stream_gray_rad_down (is, js, Time_diag, lat, p_half, t,         &
                           net_surf_sw_down, surf_lw_down)

! Begin the radiation calculation by computing downward fluxes.
! This part of the calculation does not depend on the surface temperature.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in), dimension(:,:)    :: lat
real, intent(out), dimension(:,:)   :: net_surf_sw_down
real, intent(out), dimension(:,:)   :: surf_lw_down
real, intent(in), dimension(:,:,:)  :: t, p_half

integer :: i, j, k, n

logical :: used

n = size(t,3)

! insolation at TOA
p2          = (1. - 3.*sin(lat)**2)/4.
insolation  = 0.25 * solar_constant * (1.0 + del_sol * p2 + del_sw * sin(lat))

! LW optical thickness
lw_tau_0    = ir_tau_eq + (ir_tau_pole - ir_tau_eq)*sin(lat)**2

! SW optical thickness
sw_tau_0    = (1.0 - sw_diff*sin(lat)**2)*atm_abs

! constant albedo 
albedo(:,:) = albedo_value

! compute optical depths for each model level
do k = 1, n+1
 
  lw_tau(:,:,k) = lw_tau_0 * ( linear_tau * p_half(:,:,k)/pstd_mks     &
       + (1.0 - linear_tau) * (p_half(:,:,k)/pstd_mks)**wv_exponent )

  sw_tau(:,:,k) = sw_tau_0 * (p_half(:,:,k)/pstd_mks)**solar_exponent

end do

! longwave source function

b = stefan*t**4

! longwave differential transmissivity
do k = 1, n
   lw_dtrans(:,:,k) = exp( -(lw_tau(:,:,k+1) - lw_tau(:,:,k)) )
end do

! compute downward longwave flux by integrating downward
lw_down(:,:,1)      = 0.
do k = 1, n
   lw_down(:,:,k+1) = lw_down(:,:,k)*lw_dtrans(:,:,k) + b(:,:,k)*(1. - lw_dtrans(:,:,k))
end do

! compute downward shortwave flux
do k = 1, n+1
   sw_down(:,:,k)   = insolation(:,:) * exp(-sw_tau(:,:,k)) 
end do

surf_lw_down     = lw_down(:, :, n+1)
toa_sw_in        = sw_down(:, :, 1)
net_surf_sw_down = sw_down(:, :, n+1) * (1. - albedo)

!------- downward lw flux surface -------
if ( id_lwdn_sfc > 0 ) then
   used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag)
endif
!------- incoming sw flux toa -------
if ( id_swdn_toa > 0 ) then
   used = send_data ( id_swdn_toa, toa_sw_in, Time_diag)
endif
!------- downward sw flux surface -------
if ( id_swdn_sfc > 0 ) then
   used = send_data ( id_swdn_sfc, net_surf_sw_down, Time_diag)
endif

return
end subroutine two_stream_gray_rad_down

! ==================================================================================

subroutine two_stream_gray_rad_up (is, js, Time_diag, lat, p_half, t_surf, t, tdt)

! Now complete the radiation calculation by computing the upward and net fluxes.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat
real, intent(in) , dimension(:,:)   :: t_surf
real, intent(in) , dimension(:,:,:) :: t, p_half
real, intent(inout), dimension(:,:,:) :: tdt


integer :: i, j, k, n

logical :: used

n = size(t,3)
b_surf            = stefan*t_surf**4

! compute upward longwave flux by integrating upward
lw_up(:,:,n+1)    = b_surf
do k = n, 1, -1
   lw_up(:,:,k)   = lw_up(:,:,k+1)*lw_dtrans(:,:,k) + b(:,:,k)*(1.0 - lw_dtrans(:,:,k))
end do

! compute upward shortwave flux (here taken to be constant)
do k = 1, n+1
   sw_up(:,:,k)   = albedo(:,:) * sw_down(:,:,n+1)
end do

! net fluxes (positive up)
lw_flux  = lw_up - lw_down
sw_flux  = sw_up - sw_down
rad_flux = lw_flux + sw_flux

do k = 1, n
   tdt_rad(:,:,k)   = ( rad_flux(:,:,k+1) - rad_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

   tdt_solar(:,:,k) = ( sw_flux(:,:,k+1) - sw_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

   tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
end do

olr         = lw_up(:,:,1)
net_lw_surf = lw_flux(:, :, n+1)    

!------- outgoing lw flux toa (olr) -------
if ( id_olr > 0 ) then
   used = send_data ( id_olr, olr, Time_diag)
endif
!------- upward lw flux surface -------
if ( id_lwup_sfc > 0 ) then
   used = send_data ( id_lwup_sfc, b_surf, Time_diag)
endif
!------- net upward lw flux surface -------
if ( id_net_lw_surf > 0 ) then
   used = send_data ( id_net_lw_surf, net_lw_surf, Time_diag)
endif
!------- temperature tendency due to radiation ------------
if ( id_tdt_rad > 0 ) then
   used = send_data ( id_tdt_rad, tdt_rad, Time_diag)
endif
!------- temperature tendency due to solar radiation ------------
if ( id_tdt_solar > 0 ) then
   used = send_data ( id_tdt_solar, tdt_solar, Time_diag)
endif
!------- total radiative flux (at half levels) -----------
if ( id_flux_rad > 0 ) then
   used = send_data ( id_flux_rad, rad_flux, Time_diag)
endif
!------- longwave radiative flux (at half levels) --------
if ( id_flux_lw > 0 ) then 
   used = send_data ( id_flux_lw, lw_flux, Time_diag)
endif
if ( id_flux_sw > 0 ) then
   used = send_data ( id_flux_sw, sw_flux, Time_diag)
endif

return
end subroutine two_stream_gray_rad_up

! ==================================================================================

                                                                      
subroutine two_stream_gray_rad_end
                                                                                                      
deallocate (b, tdt_rad, tdt_solar) 
deallocate (lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux)
deallocate (b_surf, olr, net_lw_surf, toa_sw_in, lw_tau_0, sw_tau_0)
deallocate (lw_dtrans, lw_tau, sw_tau)
deallocate (insolation, p2, albedo)

call cloud_end()

end subroutine two_stream_gray_rad_end

! ==================================================================================

end module two_stream_gray_rad_mod




