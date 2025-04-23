! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! Note:  This version has a namelist item for the max number of
! gps observations that can be read in, but it is currently commented out.
! Search below for the 'NAMELIST' string to find where to comment the
! code back in.  There are 2 places.  This maximum must be the total of
! all gps obs in all input files, which if you are reading multiple obs_seq
! files (e.g. for the obs_diag program) might be a larger number than 100K.

!>@todo we should have a local vs nonlocal forward operator for GPS RO,
!>so we don't have to add the metadata for the local operator.  big space
!>and time savings.  also, we should add GPSRO_BENDING_ANGLE if someone
!>can contribute a forward operator for it.

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! TEMPERATURE,             QTY_TEMPERATURE,        COMMON_CODE
! SPECIFIC_HUMIDITY,       QTY_SPECIFIC_HUMIDITY,  COMMON_CODE
! PRESSURE,                QTY_PRESSURE,           COMMON_CODE
! ADSB_HEIGHT,             QTY_ADSB_HEIGHT
! END DART PREPROCESS TYPE DEFINITIONS


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_adsb_mod, only : get_expected_adsb_height, interactive_adsb_data, &
!                              read_adsb_data, write_adsb_data
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(ADSB_HEIGHT)
!            call get_expected_adsb_height(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(ADSB_HEIGHT)
!            call read_adsb_data(obs_def%key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(ADSB_HEIGHT)
!            call write_adsb_data(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(ADSB_HEIGHT)
!            call interactive_adsb_data(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_adsb_mod

use        types_mod, only : r8, missing_r8, RAD2DEG, DEG2RAD, PI
use    utilities_mod, only : register_module, error_handler, E_ERR, &
                             nmlfileunit, check_namelist_read,      &
                             find_namelist_in_file, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, set_location, get_location, &
                             is_vertical, &
                             VERTISHEIGHT
use  assim_model_mod, only : interpolate

use     obs_kind_mod, only : QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY, &
                             QTY_PRESSURE, QTY_U_WIND_COMPONENT, QTY_ADSB_HEIGHT

use  ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: set_adsb_data, get_adsb_data, write_adsb_data, read_adsb_data, &
          get_expected_adsb_height, interactive_adsb_data

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

! Storage for the special information required for ADS-B observations
!
integer :: max_adsb_obs = 100000

type adsb_type
   private
   real(r8)            :: lon_rec
   real(r8)            :: lat_rec
   real(r8)            :: alt_rec
   real(r8)            :: aircraft_arc_ang       ! angle between receiver radius vector and 
                                                 ! aircraft radius vector   
   real(r8)            :: ray_direction(3)       ! z-component is AoA (beta)
   real(r8)            :: rfict
   real(r8)            :: step_size
end type adsb_type

type(adsb_type), allocatable :: adsb_data(:)

namelist /obs_def_adsb_nml/ max_adsb_obs

character(len=129) :: string1, string2
integer  :: ii
integer  :: keycount 

contains

!------------------------------------------------------------------------------
subroutine initialize_module
!------------------------------------------------------------------------------
!
! initialize global gps private key number and allocate space for obs data
integer :: rc, iunit

call register_module(source, revision, revdate)
module_initialized = .true.

! global count of all gps observations from any input file
keycount = 0

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_def_adsb_nml", iunit)
read(iunit, nml = obs_def_adsb_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_adsb_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_def_adsb_nml)
if (do_nml_term()) write(     *     , nml=obs_def_adsb_nml)

! find max number of gps obs which can be stored, and initialize type
allocate(adsb_data(max_adsb_obs), stat = rc)
if (rc /= 0) then
   write(string1, *) 'initial allocation failed for adsb observation data,', &
                       'itemcount = ', max_adsb_obs
   call error_handler(E_ERR,'initialize_module', string1, &
                      source, revision, revdate)
endif

end subroutine initialize_module



subroutine set_adsb_data(adsbkey, nx, ny, nz, rfict0, ds, lon_rec, lat_rec, alt_rec, aircraft_arc_ang)
!------------------------------------------------------------------------------
!
! increment key and set all private data for this observation

integer,          intent(out) :: adsbkey
real(r8),         intent(in)  :: nx, ny, nz, rfict0, ds, lon_rec, lat_rec, alt_rec, aircraft_arc_ang

if ( .not. module_initialized ) call initialize_module

keycount = keycount + 1
adsbkey = keycount

if(adsbkey > max_adsb_obs) then
   write(string1, *) 'key (',adsbkey,') exceeds max_adsb_obs (',max_adsb_obs,')'
   string2 = 'Increase max_adsb_obs in input.nml &obs_def_adsb_nml namelist.'
   call error_handler(E_ERR,'read_adsb_data', string1, &
                      source, revision, revdate, text2=string2)
endif

adsb_data(adsbkey)%ray_direction(1) = nx
adsb_data(adsbkey)%ray_direction(2) = ny
adsb_data(adsbkey)%ray_direction(3) = nz
adsb_data(adsbkey)%lon_rec = lon_rec
adsb_data(adsbkey)%lat_rec = lat_rec
adsb_data(adsbkey)%alt_rec = alt_rec
adsb_data(adsbkey)%rfict     = rfict0
adsb_data(adsbkey)%step_size = ds
adsb_data(adsbkey)%aircraft_arc_ang   = aircraft_arc_ang

end subroutine set_adsb_data

subroutine get_adsb_data(adsbkey, nx, ny, nz, rfict0, ds, lon_rec, lat_rec, alt_rec, aircraft_arc_ang)
!------------------------------------------------------------------------------
!
! return all private data for this observation

integer,          intent(in)  :: adsbkey
real(r8),         intent(out) :: nx, ny, nz, rfict0, ds, lon_rec, lat_rec, alt_rec, aircraft_arc_ang

if ( .not. module_initialized ) call initialize_module

if (adsbkey < 1 .or. adsbkey > keycount) then
   write(string1, *) 'key (',adsbkey,') out of valid range (1<=key<=',keycount,')'
   call error_handler(E_ERR,'get_adsb_data', string1, &
                      source, revision, revdate)
endif

nx = adsb_data(adsbkey)%ray_direction(1)
ny = adsb_data(adsbkey)%ray_direction(2)
nz = adsb_data(adsbkey)%ray_direction(3)
lon_rec = adsb_data(adsbkey)%lon_rec
lat_rec = adsb_data(adsbkey)%lat_rec
alt_rec = adsb_data(adsbkey)%alt_rec
rfict0 = adsb_data(adsbkey)%rfict     
ds = adsb_data(adsbkey)%step_size
aircraft_arc_ang = adsb_data(adsbkey)%aircraft_arc_ang
 
end subroutine get_adsb_data



subroutine write_adsb_data(adsbkey, ifile, fform)
!------------------------------------------------------------------------------
!

integer,          intent(in)           :: adsbkey, ifile
character(len=*), intent(in), optional :: fform


if ( .not. module_initialized ) call initialize_module

! Write the 5 character identifier for verbose formatted output
! Write out the obs_def key for this observation
if (ascii_file_format(fform)) then
   write(ifile,11) adsbkey
   write(ifile, *) (adsb_data(adsbkey)%ray_direction(ii), ii=1, 3), &
                  adsb_data(adsbkey)%lon_rec, adsb_data(adsbkey)%lat_rec, &
                  adsb_data(adsbkey)%alt_rec, adsb_data(adsbkey)%rfict, &
                   adsb_data(adsbkey)%step_size, adsb_data(adsbkey)%aircraft_arc_ang
                  
11  format('adsb', i8)
else
   write(ifile) adsbkey
   write(ifile, *) (adsb_data(adsbkey)%ray_direction(ii), ii=1, 3), &
                  adsb_data(adsbkey)%lon_rec, adsb_data(adsbkey)%lat_rec, &
                  adsb_data(adsbkey)%alt_rec, adsb_data(adsbkey)%rfict, &
                   adsb_data(adsbkey)%step_size, adsb_data(adsbkey)%aircraft_arc_ang
endif

end subroutine write_adsb_data


 subroutine read_adsb_data(adsbkey, ifile, fform)
!------------------------------------------------------------------------------
!
! Every ADS-B observation has its own (metadata) adsbkey.
! When you read multiple adsb observation sequence files, it is necessary 
! to track the total number of metadata adsbkeys read, not just the number 
! in the current file.
! 

integer,          intent(out)          :: adsbkey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

integer :: keyin    ! the metadata key in the current obs sequence

real(r8) :: nx, ny, nz, rfict0, ds, lon_rec, lat_rec, alt_rec, aircraft_arc_ang
character(len=8) :: header

if (ascii_file_format(fform)) then
   read(ifile, FMT='(a8, i8)') header, keyin    ! throw away keyin
   if(header /= 'adsb') then
       call error_handler(E_ERR,'read_adsb_data', &
       'Expected header "adsb" in input file', source, revision, revdate)
   endif
   read(ifile, *) nx, ny, nz, lon_rec, lat_rec, alt_rec, rfict0, ds,  aircraft_arc_ang
else
   read(ifile) keyin          ! read and throw away
   read(ifile) nx, ny, nz, lon_rec, lat_rec, alt_rec, rfict0, ds, aircraft_arc_ang
endif


! increment key and set all private data for this observation
call set_adsb_data(adsbkey, nx, ny, nz, rfict0, ds, lon_rec, lat_rec, alt_rec, aircraft_arc_ang)

end subroutine read_adsb_data


subroutine interactive_adsb_data(adsbkey)
!----------------------------------------------------------------------
!
! Interactively prompt for the info needed to create a ADS-B refractivity 
! observation.  Increments the key number and returns it.

integer, intent(out) :: adsbkey

real(r8) :: nx, ny, nz, rfict0, ds, aircraft_arc_ang, lat_rec, lon_rec, alt_rec

if ( .not. module_initialized ) call initialize_module

write(*, *)
write(*, *) 'Beginning to inquire information on ADS-B type.'
write(*, *)

100 continue



    ! FIXME:  i have no idea what valid values are for any
   !  of the following items, so i cannot add any error checking or
   !  guidance for the user.

   write(*, *)
   write(*, *) 'Enter X, Y, Z value for ray direction'
   write(*, *)
   read(*,*) nx, ny, nz

   write(*, *)
   write(*, *) 'Enter local curvature radius'
   write(*, *)
   read(*,*) rfict0

   write(*, *)
   write(*, *) 'Enter receiver lat, lon and height'
   write(*, *)
   read(*,*) lat_rec, lon_rec, alt_rec


   write(*, *)
   write(*, *) 'Enter step size'
   write(*, *)
   read(*,*) ds

   write(*, *)
   write(*, *) 'Enter ray top'
   write(*, *)
   read(*,*) aircraft_arc_ang

   nx = 0.0_r8
   ny = 0.0_r8
   nz = 0.0_r8
   rfict0 = 0.0_r8
   ds = 0.0_r8
   lat_rec = 0.0_r8
   lon_rec = 0.0_r8
   alt_rec = 0.0_r8
   aircraft_arc_ang = 0.0


! increment key and set all private data for this observation
call set_adsb_data(adsbkey, nx, ny, nz, rfict0, ds, lon_rec, lat_rec, alt_rec, aircraft_arc_ang)

write(*, *)
write(*, *) 'End of specialized section for gps observation data.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_adsb_data

!> Distributed version of get_expected_adsb_height
 subroutine get_expected_adsb_height(state_handle, ens_size,  location, adsbkey, ray_height, istatus)
!------------------------------------------------------------------------------
!
! Purpose: Calculate the ray altitude after integration
!
!------------------------------------------------------------------------------
!
! inputs:
!    state_vector:    DART state vector
!
! output parameters:
!    ray_height: altitude of ray after integration
!    istatus:  =0 normal; =1 outside of domain.
!------------------------------------------------------------------------------
!  Author: Ollie Lewis 
!  Version 1.1: Sep 17 2024
!
!------------------------------------------------------------------------------

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: adsbkey
real(r8),            intent(out) :: ray_height(ens_size)
integer,             intent(out) :: istatus(ens_size)

! local variables.  first set is per ensemble member.
! second set is independent of a particular ensemble member.

real(r8) :: ref_ray(ens_size), ref_up(ens_size), ref_low(ens_size), ref_grad(ens_size)
real(r8) :: ref00(ens_size)
real(r8) :: xx(ens_size), yy(ens_size), zz(ens_size)   ! cartesian ray positions   
real(r8) :: nxx(ens_size), nyy(ens_size), nzz(ens_size)  ! ray direction in cartesian coords
real(r8) :: height1(ens_size)   ! height of ray is different for each ensemble member
real(r8) :: height0(ens_size)   
integer  :: this_istatus(ens_size)

real(r8) :: nx, ny, nz          ! direction of ray at receiver              
real(r8) :: xo(ens_size), yo(ens_size), zo(ens_size)
real(r8) :: lat1(ens_size), lon1(ens_size)
real(r8) :: arc(ens_size)
real(r8) :: aircraft_arc_ang
real(r8) :: lon(ens_size), lat(ens_size), height(ens_size)
integer  :: iter
logical  :: return_now

if ( .not. module_initialized ) call initialize_module



lon      = adsb_data(adsbkey)%lon_rec                       ! degree: 0 to 360
lat      = adsb_data(adsbkey)%lat_rec                       ! degree: -90 to 90
height0   = adsb_data(adsbkey)%alt_rec                    ! (m)

! to use track_status() start out with istatus all success.
istatus = 0

nx = adsb_data(adsbkey)%ray_direction(1)
ny = adsb_data(adsbkey)%ray_direction(2)
nz = adsb_data(adsbkey)%ray_direction(3)

aircraft_arc_ang = adsb_data(adsbkey)%aircraft_arc_ang

! convert location of the perigee from geodetic to Cartesian coordinate

call geo2carte (height0, lat, lon, ens_size, xo, yo, zo, adsb_data(adsbkey)%rfict )

iter = 0

xx = xo
yy = yo
zz = zo
nzz = nz
nxx = nx
nyy = ny

INTEGRATE: do 

   iter = iter + 1

   !  integrate to one direction of the ray for one step
   ! HK These are now different for each ensemble member
   xx = xx + adsb_data(adsbkey)%step_size * nxx
   yy = yy + adsb_data(adsbkey)%step_size * nyy
   zz = zz + adsb_data(adsbkey)%step_size * nzz

   ! convert the location of the point to geodetic coordinates 
   ! height(m), lat, lon(deg)

  call carte2geo(xx, yy, zz, ens_size, height1, lat1, lon1, adsb_data(adsbkey)%rfict )  
   
   ! compute the angle subtended by the ray at the centre
   ! of the Earth (this defines the end point of the ray)

   call haversine(lon, lat, lon1, lat1, ens_size, arc)

   if (any(arc >= adsb_data(adsbkey)%aircraft_arc_ang)) exit INTEGRATE
   

   ! get the refractivity at this ray point(ref00)
   call ref_local(state_handle, ens_size, lat1, lon1, height1, ref00, this_istatus)
   call track_status(ens_size, this_istatus, ray_height, istatus, return_now)
   if (return_now) return

   ref_ray = ref00

   ! get the refractivity above ray at this ray point
   call ref_local(state_handle, ens_size, lat1, lon1, height1+1.0_r8, ref00, this_istatus)
   call track_status(ens_size, this_istatus, ray_height, istatus, return_now)
   if (return_now) return

   ref_up = ref00

   ! get the refractivity below ray at this ray point
   call ref_local(state_handle, ens_size, lat1, lon1, height1-1.0_r8, ref00, this_istatus)
   call track_status(ens_size, this_istatus, ray_height, istatus, return_now)
   if (return_now) return

   ref_low = ref00

   ! compute refractivity gradient (assume only vertical gradient currently)
   ref_grad = ref_ray * (ref_up - ref_low) / 2.0_r8

   ! update vertical component of ray direction
   nxx = nxx + (ref_grad * cos(lat1*DEG2RAD) * cos(lon1*DEG2RAD)) * adsb_data(adsbkey)%step_size
   nyy = nyy + (ref_grad * cos(lat1*DEG2RAD) * sin(lon1*DEG2RAD)) * adsb_data(adsbkey)%step_size
   nzz = nzz + (ref_grad * sin(lat1*DEG2RAD)) * adsb_data(adsbkey)%step_size

end do INTEGRATE

! finish the integration of the height of the ray

where(istatus == 0) ray_height = height1    ! in m


! make sure return is missing_r8 if failure.
!>@todo is the first line necessary?  i believe the second one
!> is a necessary test.
where (istatus /= 0) ray_height = missing_r8
where (istatus == 0 .and. ray_height < 0.0_r8)
   istatus = 5
   ray_height = missing_r8
endwhere

end subroutine get_expected_adsb_height

!> Distributed version.  i removed the unused 'location' argument.
! i changed this to be lat, lon, height - which is a more logical layout of
! the location information.  it used to be height, lat, lon - be careful if
! anyone else calls this routine.  it isn't public so no code outside this
! module can be calling it.   nsc 30 oct 2015

subroutine ref_local(state_handle, ens_size, lats, lons, height, ref00, istatus0)
!------------------------------------------------------------------------------
!
! Calculate local refractivity at any ADS-B ray point (lat, lon, height)
!
! inputs:
!    lat, lon, height:  ADS-B ray observation location (units: degrees, degrees, meters)
!
! output:
!    ref00: modeled local refractivity at ray point(unit: N-1, ~1.0e-4 to e-6)
!
!------------------------------------------------------------------------------

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
real(r8), intent(in) :: lons(ens_size), lats(ens_size)
real(r8), intent(in) :: height(ens_size)
real(r8),            intent(out) :: ref00(ens_size)
integer,             intent(out) :: istatus0(ens_size)

real(r8), parameter::  rd = 287.05_r8, rv = 461.51_r8, c1 = 77.6d-6 , &
                       c2 = 3.73d-1,  rdorv = rd/rv

real(r8) :: lon2(ens_size)
real(r8) :: t(ens_size), q(ens_size), p(ens_size), ew(ens_size)
real(r8) :: t_dummy(ens_size), q_dummy(ens_size), p_dummy(ens_size)
integer  :: this_istatus(ens_size)
logical  :: return_now
integer  :: imem
type(location_type) :: location

if ( .not. module_initialized ) call initialize_module

! for sequential sets of calls to interpolate, start with a 0 istatus
! and each call to track_status() will set new failures.
istatus0 = 0
ref00 = missing_r8

lon2 = lons
where (lons > 360.0_r8)
    lon2 = lons - 360.0_r8
end where

where (lons < 0.0_r8)
    lon2 = lons + 360.0_r8
end where


!  required variable units for calculation of refractivity
!   t :  Kelvin, from top to bottom
!   q :  kg/kg, from top to bottom
!   p :  mb

do imem = 1, ens_size
     

      location = set_location( lon2(imem), lats(imem), height(imem), VERTISHEIGHT )

      call interpolate(state_handle, ens_size, location,  QTY_TEMPERATURE, t_dummy, this_istatus)
      call track_status(ens_size, this_istatus, ref00, istatus0, return_now)
      if (return_now) return

      call interpolate(state_handle, ens_size, location, QTY_SPECIFIC_HUMIDITY, q_dummy, this_istatus)
      call track_status(ens_size, this_istatus, ref00, istatus0, return_now)
      if (return_now) return

      call interpolate(state_handle, ens_size, location,  QTY_PRESSURE, p_dummy, this_istatus)
      call track_status(ens_size, this_istatus, ref00, istatus0, return_now)
      if (return_now) return
      

      t(imem) = t_dummy(imem)
      q(imem) = q_dummy(imem)
      p(imem) = p_dummy(imem)

enddo


where (istatus0 == 0) 
p     = p * 0.01_r8      ! to mb

ew    = q * p/(rdorv + (1.0_r8-rdorv)*q )
ref00 = 1.0_r8 + c1*p/t + c2*ew/(t**2)              ! 1 + N/1e6
endwhere

end subroutine ref_local


 subroutine geo2carte (s1, s2, s3, ens_size, x1, x2, x3, rfict0) 
!------------------------------------------------------------------------------
!
!  Converts geodetical coordinates to cartesian with a reference sphere
!------------------------------------------------------------------------------
!  input parameters:
!   s - geodetical coordinates
!        (height (m), latitude (degree), longitude (degree))
!                     -90 to 90           0 to 360
!  output parameters:
!   x - cartesian coordinates (m) connected with the earth(x, y, z-coordinate)
!------------------------------------------------------------------------------
implicit none
integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: rfict0    ! units: m
real(r8), intent(in)  :: s1(ens_size), s2(ens_size), s3(ens_size)
real(r8), intent(out) :: x1(ens_size), x2(ens_size), x3(ens_size)
real(r8) :: g3, g4
integer  :: imem

if ( .not. module_initialized ) call initialize_module

do imem = 1, ens_size
   g3 = s1(imem) + rfict0
   g4 = g3 * cos(s2(imem)*DEG2RAD) 
   x1(imem) = g4 * cos(s3(imem)*DEG2RAD)
   x2(imem) = g4 * sin(s3(imem)*DEG2RAD)
   x3(imem) = g3 * sin(s2(imem)*DEG2RAD)
end do

end subroutine geo2carte


 subroutine carte2geo (x1, x2, x3, ens_size, s1, s2, s3, rfict0)
!------------------------------------------------------------------------------
!
!  Converts cartesian coordinates to geodetical.
!
!   input parameters:
!        x - cartesian coordinates (x, y, z-coordinate, unit: m)
!
!   output parameters:
!        s - geodetical coordinates
!            (height (m), latitude (deg), longitude (deg))
!                          -90 to 90         0 to 360
!------------------------------------------------------------------------------
implicit none
integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: rfict0
real(r8), intent(in)  :: x1(ens_size), x2(ens_size), x3(ens_size)
real(r8), intent(out) :: s1(ens_size), s2(ens_size), s3(ens_size)

real(r8), parameter :: crcl  = 2.0_r8 * PI, &
                       crcl2 = 4.0_r8 * PI

real(r8) :: rho, sphi, azmth
integer  :: imem

if ( .not. module_initialized ) call initialize_module

do imem = 1, ens_size
   rho   = sqrt (x1(imem)**2 + x2(imem)**2 + x3(imem)**2 ) 
   sphi  = x3(imem)/rho
   s1(imem)    = rho - rfict0
   s2(imem)    = asin (sphi) 
   azmth = atan2 (x2(imem), x1(imem))
   s3(imem)    = mod((azmth + crcl2), crcl)

   s2(imem)    = s2(imem) * RAD2DEG
   s3(imem)    = s3(imem) * RAD2DEG
end do

end  subroutine carte2geo

subroutine haversine (lon1, lat1, lon2, lat2, ens_size, arc)
!------------------------------------------------------------------------------
!
!  Computes the angle subtended by the ray at the centre of the Earth
!
!   input parameters:
!        lon1, lat1 - longitude and latitude of observer
!        lon2, lat2 - longitude and latitude of ray point
!
!   output parameters:
!        arc - angle subtended by the ray at the centre of the Earth
!------------------------------------------------------------------------------
implicit none
integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: lon1(ens_size), lat1(ens_size), lon2(ens_size), lat2(ens_size)
real(r8), intent(out) :: arc(ens_size)
integer  :: imem

real(r8) :: a, delta_lat, delta_lon

if ( .not. module_initialized ) call initialize_module

do imem = 1, ens_size
   delta_lon = (lon2(imem) - lon1(imem)) * DEG2RAD
   delta_lat = (lat2(imem) - lat1(imem)) * DEG2RAD

   a = sin(delta_lat/2.0_r8)**2 + cos(lat1(imem) * DEG2RAD)*cos(lat2(imem)*DEG2RAD)*(sin(delta_lon/2.0_r8)**2)

   arc = 2 * atan2(sqrt(a), sqrt(1.0_r8 - a))
end do

end  subroutine haversine


end module obs_def_adsb_mod

! END DART PREPROCESS MODULE CODE

