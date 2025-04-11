! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! Forward operator to compute total precipitable water in a column,
! in centimeters, over the ocean.   Can be used as an example of a
! forward operator that loops over either fixed pressure levels or
! over model levels.

! This code is only correct for TPW over the ocean where the surface
! is at 0m elevation.  For TWP over land you must be able to compute
! the surface elevation for any given lat/lon location.  This code
! assumes the model can return the surface pressure for any given
! lat/lon location, and the specific humidity for a given location
! where the vertical is either pressure or model levels.

! keep in mind that fortran allows only 31 characters in parameter
! definitions (which is what this string is going to be used for).
! if the platform name gets longer than 5 chars, consider going
! to something like xxx_TOTAL_PRECIP_WATER to give you room to
! put in more descriptive platform names.

! BEGIN DART PREPROCESS TYPE DEFINITIONS
!  CTT,                 QTY_CTT,
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_ctt_mod, only : get_expected_ctt
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(CTT)
!            call get_expected_ctt(state_handle, ens_size, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(CTT)
!           continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(CTT)
!           continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(CTT)
!           continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_ctt_mod

use        types_mod, only : r8, missing_r8, RAD2DEG, DEG2RAD, PI
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, do_nml_file, do_nml_term, &
                             check_namelist_read, find_namelist_in_file
use     location_mod, only : location_type, set_location, get_location, &
                             write_location, read_location, &
                             VERTISLEVEL, VERTISPRESSURE, VERTISSURFACE
use time_manager_mod, only : time_type, read_time, write_time, &
                             set_time, set_time_missing
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : QTY_SURFACE_PRESSURE, QTY_VAPOR_MIXING_RATIO, &
                             QTY_CLOUDWATER_MIXING_RATIO, QTY_RAINWATER_MIXING_RATIO, &
                             QTY_ICE_MIXING_RATIO, QTY_GRAUPEL_MIXING_RATIO, &
                             QTY_TEMPERATURE, QTY_PRESSURE
use ensemble_manager_mod,  only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public ::  get_expected_ctt

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

character(len=129) :: msgstring

integer :: max_pressure_intervals = 1000   ! increase as needed

! default samples the atmosphere between the surface and 200 hPa 
! at the model level numbers.  if model_levels is set false,
! then the default samples at 40 heights, evenly divided in
! linear steps in pressure between the surface and top.

logical  :: model_levels = .true.        ! if true, use model levels, ignores num_pres_int
real(r8) :: pressure_top = 5000.0       ! top pressure in pascals
logical  :: separate_surface_level = .true.  ! false: level 1 of 3d grid is sfc
                                             ! true: sfc is separate from 3d grid
integer  :: num_pressure_intervals = 40  ! number of intervals if model_levels is F


namelist /obs_def_ctt_nml/ model_levels, pressure_top,  &
                           separate_surface_level, num_pressure_intervals

contains

!------------------------------------------------------------------------------
subroutine initialize_module()

! should be called once by the other routines in this file to be sure
! the namelist has been read and any initialization code is run.

integer :: iunit, rc

if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_def_ctt_nml", iunit)
read(iunit, nml = obs_def_ctt_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_ctt_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_def_ctt_nml)
if (do_nml_term()) write(     *     , nml=obs_def_ctt_nml)

end subroutine initialize_module


!------------------------------------------------------------------------------
subroutine get_expected_ctt(state_handle, ens_size, location, ctt, istatus)

!------------------------------------------------------------------------------
! Purpose:  To calculate the cloud top temperature at a particular location.
! inputs:
!    state_vector:    DART state vector
!    location:        Observation location
!
! output parameters:
!    ctt:     cloud top temperature
!    istatus: 0 if ok, a positive value for error
!------------------------------------------------------------------------------
!  Author: Ollie Lewis ,  Version 1.1: Apr 10, 2025 
!  
!------------------------------------------------------------------------------

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: ctt(ens_size)
integer,             intent(out) :: istatus(ens_size)

! local variables
real(r8) :: lon, lat, height, obsloc(3)
real(r8) :: pressure(ens_size)
real(r8) :: qv(ens_size), qc(ens_size), qr(ens_size), qi(ens_size), qg(ens_size), &
            qtot(ens_size), qtot_prev(ens_size)
real(r8) :: psfc(ens_size)
real(r8) :: t(ens_size)
type(location_type) :: location2
integer  :: which_vert, k, lastk, first_non_surface_level

integer  :: iter
integer  :: this_istatus(ens_size)
integer  :: imem
logical  :: return_now

if ( .not. module_initialized ) call initialize_module

ctt = missing_r8
istatus = 0

! location is the lat/lon where we need to compute the column quantity
obsloc   = get_location(location)
lon      = obsloc(1)                       ! degree: 0 to 360
lat      = obsloc(2)                       ! degree: -90 to 90

qtot = 0.0_r8

which_vert = VERTISSURFACE
height = 0.0
location2 = set_location(lon, lat, height,  which_vert)

! interpolate the surface pressure and specific humidity at the desired location
! assumes the values returned from the interpolation will be in these units:
!   surface pressure :  Pa
!   moisture         :  kg/kg
call interpolate(state_handle, ens_size, location2, QTY_SURFACE_PRESSURE, pressure, this_istatus)
call track_status(ens_size, this_istatus, ctt, istatus, return_now)
if (return_now) return

! save this for use below
psfc = pressure

! there are two options for constructing the column of values.  if 'model_levels'
! is true, we query the model by vertical level number.  the 'separate_surface_level'
! flag should be set to indicate if the lowest level of the 3d grid is the
! surface or if the surface values are a separate quantity below the 3d grid.

if (model_levels) then
   
   ! some models have a 3d grid of values and the lowest level contains
   ! the surface quantities.  others have a separate field for the
   ! surface values and the 3d grid starts at some given elevation.
   ! if the namelist value 'separate_surface_level'  is true, we will
   ! ask to interpolate a surface pressure first and then work up the
   ! 3d column starting at level 1.  if it is false, we assume level 1
   ! was the surface pressure and we start here at level 2.

   if (separate_surface_level) then
      first_non_surface_level = 1
   else
      first_non_surface_level = 2
   endif

   ! construct a pressure column on model levels

   ! call the model until the interpolation call fails (above the top level)
   ! (this is not a fatal error unless the first call fails).
   ! also exit the loop if the pressure is above the namelist-specified pressure top

   lastk = 2

   do imem = 1, ens_size

      LEVELS: do k=first_non_surface_level, 10000   ! something unreasonably large
      
         ! call the model_mod to get the pressure and specific humidity at each level 
         ! from the model and fill out the pressure and qv arrays.  the model must
         ! support a vertical type of level number.

         which_vert = VERTISLEVEL
         location2 = set_location(lon, lat, real(k, r8),  which_vert)

         call interpolate(state_handle, ens_size, location2, QTY_PRESSURE, pressure, this_istatus)
         call track_status(ens_size, this_istatus, pressure, istatus, return_now)
         if (pressure(imem) < pressure_top) exit LEVELS
         if (return_now) return

         ! Computing mixing ratios to compute total mixing ratio for threshold calculation

         qtot_prev = qtot

         call interpolate(state_handle, ens_size, location2, QTY_VAPOR_MIXING_RATIO, qv, this_istatus)
         call track_status(ens_size, this_istatus, qv, istatus, return_now)

         call interpolate(state_handle, ens_size, location2, QTY_CLOUDWATER_MIXING_RATIO, qc, this_istatus)
         call track_status(ens_size, this_istatus, qc, istatus, return_now)

         call interpolate(state_handle, ens_size, location2, QTY_RAINWATER_MIXING_RATIO, qr, this_istatus)
         call track_status(ens_size, this_istatus, qr, istatus, return_now)

!         call interpolate(state_handle, ens_size, location2, QTY_ICE_MIXING_RATIO, qi, this_istatus)
!         call track_status(ens_size, this_istatus, qi, istatus, return_now)

!         call interpolate(state_handle, ens_size, location2, QTY_GRAUPEL_MIXING_RATIO, qg, this_istatus)
!         call track_status(ens_size, this_istatus, qg, istatus, return_now)


         qtot = qv(imem) + qc(imem) + qr(imem)

         write(*,*) qtot, imem, k, ctt(imem), pressure(imem)*1.0d-2, pressure_top


         if (qtot(imem) < 1.0d-6 .and. qtot_prev(imem) >= 1.0d-6) then
            call interpolate(state_handle, ens_size, location2, QTY_TEMPERATURE, t, this_istatus)
            call track_status(ens_size, this_istatus, t, istatus, return_now)

            ctt(imem) = t(imem)

            write(*,*) qtot, imem, k, ctt(imem)

         end if


   
         if (return_now) return
      
         lastk = lastk + 1
      enddo LEVELS

      lastk = lastk - 1

      if (lastk == 1) then 

         location2 = set_location(lon, lat, real(first_non_surface_level, r8),  which_vert)
         call interpolate(state_handle, ens_size, location2, QTY_TEMPERATURE, t, this_istatus)
         call track_status(ens_size, this_istatus, ctt, istatus, return_now)

         ctt(imem) = t(imem)

         if (return_now) return
      endif

   end do

endif

end subroutine get_expected_ctt

end module obs_def_ctt_mod

! END DART PREPROCESS MODULE CODE

