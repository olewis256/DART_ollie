! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_adsb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read 
use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, &
                                  nc_get_global_attribute
use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              increment_time, get_time, get_date, operator(-), GREGORIAN
use      location_mod, only : VERTISHEIGHT, set_location
use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,         &
                               static_init_obs_sequence, init_obs, destroy_obs,   &
                               write_obs_seq, init_obs_sequence,                  &
                               destroy_obs_sequence, set_obs_values, set_obs_def, &
                               set_copy_meta_data, set_qc, set_qc_meta_data, get_num_obs
use      obs_kind_mod, only : ADSB_HEIGHT
use   obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                               set_obs_def_error_variance, set_obs_def_location,        &
                               set_obs_def_key
use obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen, getvar_real_2d, &
                              getvar_int_2d, query_varname, set_missing_name
use    obs_def_adsb_mod, only : set_adsb_data

implicit none

character(len=129),  parameter :: adsb_netcdf_file = '../data/adsb_input.nc'
character(len=129), parameter :: adsb_out_file    = 'obs_seq.in'
character(len=129)            :: adsb_outfile

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer  :: iunit, io
integer  :: ncid, nstn, n, i, oday, osec, nused, index, ntime, it, obs_num
logical  :: file_exist, first_obs
real(r8) :: qc
real(r8) :: pwv_miss = -999.
!real(r8) :: pwverr_miss = -999.

character(len=129), allocatable :: sitenames(:)
real(r8), allocatable :: lat(:), lon(:), elev(:), toff(:)
real(r8), allocatable :: aircraft_arc_ang(:,:), &
                         alt_air(:,:), alt_err(:,:), lat_air(:,:), lon_air(:,:), roc(:), &
                         azim(:,:), aoa(:,:)
real(r8)              :: xdir, ydir, zdir                         
real(r8)              :: obs_val(1), qc_val(1)
real(r8)              :: obs_window

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

integer, parameter :: MAX_NAME = 256
character(len=MAX_NAME) :: varname(5)

! For the data resource, check Readme in data/ directory.
! Sumoinet data provides observation times only as offset [sec] from 00Z 
! at each gregorian day. Thus, the actual observation date is read from 
! global attributes available in the file.
! The Suominet data has two different data regions - either over CONUS or the
! whole globe. We add the region name to the output file name (as well as
! the observation time).
! In this converter, we can process either hourly conus data or daily global 
! data depending on the global_data flag.
! There are other surface variables available at the stations if one wants
! to read from this dataset. For now, we only get pwv.

integer  :: max_num_obs                 = 1000000
real(r8) :: obs_window_hr               = 6          ! read data only every 6 hr from 00Z.
logical  :: global_data                 = .true.

namelist /convert_adsb_nml/ global_data, max_num_obs, obs_window_hr

character(len=6) :: ftail 
character(len=19):: sdate
character(len= 8):: ymd            ! YYYYMMDD
character(len=10):: ymdh           ! YYYYMMDDHH
integer          :: iyear, iday, ihour
integer          :: iyr, imo, idy, ihr, imn, isc

!------------
! start of executable code
!------------

call initialize_utilities('convert_adsb')

call find_namelist_in_file("input.nml", "convert_adsb_nml", iunit)
read(iunit, nml = convert_adsb_nml, iostat = io)
call check_namelist_read(iunit, io, "mpas_obs_preproc_nml")


first_obs = .true.

ncid = nc_open_file_readonly(adsb_netcdf_file, 'convert_adsb')
call nc_get_global_attribute(ncid, 'start_date', sdate)
read(sdate,'(i4,1x,i3,1x,i2)')  iyear, iday, ihour

call set_calendar_type(GREGORIAN)
comp_day0 = set_date(iyear,1,1,0,0,0)
call get_time(comp_day0, osec, oday)
!print *,'DATE for ',iyear,' 01-01_00:00 => (oday, osec): ',oday, osec
oday = oday + iday - 1
comp_day0 = set_time(ihour*3600,oday)
call get_date(comp_day0, iyr, imo, idy, ihr, imn, isc)

! Final output file name
if( global_data ) then
   ftail = '.globe'
   write(ymd,'(I4,2I2.2)') iyr, imo, idy
   adsb_outfile = trim(adsb_out_file) // ftail // '.' // trim(ymd)
   print *,'OBS_DATE: ',oday,' gregorian days which is ',ymd
else
   ftail = '.conus'
   write(ymdh,'(I4,3I2.2)') iyr, imo, idy, ihr
   adsb_outfile = trim(adsb_out_file) // ftail // '.' // trim(ymdh)
   print *,'OBS_DATE: ',oday,' gregorian days which is ',ymdh
endif
print *,'Output file: ',trim(adsb_outfile)


call getdimlen(ncid, "nstn",      nstn)
call getdimlen(ncid, "ntime", ntime)
call set_missing_name("missing_value")

obs_window = obs_window_hr * 3600.0_r8        ! obs frequency (from hours to seconds)

allocate( lat(nstn))         
allocate( lon(nstn))
allocate(elev(nstn))         
allocate(toff(ntime))
allocate( azim(ntime,nstn))
allocate( aoa(ntime,nstn))
allocate( aircraft_arc_ang(ntime,nstn))
allocate( alt_air(ntime,nstn))
allocate( alt_err(ntime,nstn))
allocate( lon_air(ntime,nstn))
allocate( lat_air(ntime,nstn))
allocate( roc(nstn))

! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.
qc = 0.0_r8

! read in the data arrays

! we have gpspw data files which have different names for the 
! lat/lon/elev/obs arrays in the netcdf file.  there doesn't seem
! to be a global attr to say which one is in use, so for now try
! both options.  

varname(1) = 'lat_receiver'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real(ncid, varname(index),  lat            ) ! station latitude

varname(1) = 'lon_receiver'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real(ncid, varname(index),  lon            ) ! station longitude

varname(1) = 'alt_receiver'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real(ncid, varname(index),  elev           ) ! station elevation

varname(1) = 'azim'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real_2d(ncid, varname(index),  azim       )   ! ray x-direction

varname(1) = 'aoa'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real_2d(ncid, varname(index),  aoa       )   ! ray z-direction

varname(1) = 'aircraft_arc_ang'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real_2d(ncid, varname(index),  aircraft_arc_ang       )   ! aircraft arc angle

varname(1) = 'alt_air'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real_2d(ncid, varname(index),  alt_air       )   

varname(1) = 'alt_err'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real_2d(ncid, varname(index),  alt_err       )   

varname(1) = 'lon_air'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real_2d(ncid, varname(index),  lon_air      )   

varname(1) = 'lat_air'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real_2d(ncid, varname(index),  lat_air       )   

varname(1) = 'roc'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real(ncid, varname(index),  roc       )   ! radius of curvature of Earth at receiver

varname(1) = 'time_offset'    ! "PWV window midpoint time delta from start_time"
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real(ncid, varname(index),  toff        ) ! obs time offset in seconds

! if user says to use them, read in QCs if present

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=adsb_outfile, exist=file_exist)

if ( file_exist ) then

  ! existing file found, append to it
   call read_obs_seq(adsb_outfile, 0, 0, nstn, obs_seq)

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, max_num_obs)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'ADS-B refraction observation')
  end do
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Data QC')
  end do

endif

obs_num = 1

where (lon > 360.0_r8)
    lon = lon - 360.0_r8
end where

where (lon < 0.0_r8)
    lon = lon + 360.0_r8
end where

where (lon_air > 360.0_r8)
    lon_air = lon_air - 360.0_r8
end where

where (lon_air < 0.0_r8)
    lon_air = lon_air + 360.0_r8
end where

timloop: do it = 1, ntime

  ! compute time of observation
  time_obs = increment_time(comp_day0, nint(toff(it)))

  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)
  print*,'time_obs for it = ',it,oday, osec, ntime

  

  nused = 0
  stnloop: do n = 1, nstn

    call tanvec01(lon(n), lat(n), azim(it, n), aoa(it, n), xdir,ydir, zdir)

    call set_adsb_data(obs_num, xdir, ydir, zdir, roc(n), 100.0_r8, lon(n),lat(n),elev(n), aircraft_arc_ang(it,n))
    call set_obs_def_location(obs_def,set_location(lon_air(it,n),lat_air(it,n),alt_air(it,n),VERTISHEIGHT))
    call set_obs_def_type_of_obs(obs_def, ADSB_HEIGHT)
    call set_obs_def_time(obs_def, set_time(osec, oday))
    call set_obs_def_error_variance(obs_def, alt_err(it,n) * alt_err(it,n))
    call set_obs_def_key(obs_def, obs_num)
    call set_obs_def(obs, obs_def)
 
    obs_val(1) = alt_air(it,n)
    call set_obs_values(obs, obs_val)

    qc_val(1)  = 0
    call set_qc(obs, qc_val)
    call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    nused = nused + 1
    obs_num = obs_num + 1
    


  end do stnloop
  print*, 'nused =',nused

end do timloop

deallocate( lat )         
deallocate( lon )
deallocate( elev )         
deallocate( toff )
deallocate( azim )
deallocate( aoa )

deallocate( aircraft_arc_ang )
deallocate( alt_air )
deallocate( alt_err )
deallocate( lon_air )
deallocate( lat_air )
deallocate( roc )

call nc_close_file(ncid, 'convert_adsb')

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, adsb_outfile)


! end of main program
call finalize_utilities()

contains

subroutine tanvec01(lon0, lat0, azim0, elev0, ux, uy, uz)

use   types_mod, only : r8, deg2rad

implicit none

real(r8), intent(in)  :: lon0, lat0, azim0, elev0
real(r8), intent(out) :: ux, uy, uz

real(r8) :: zz0(3), lon, lat, azim, rtp(3), rnm(3), rno(3), uon(3), vlen0, rn(3), ufin(3), elev

zz0(1) = 0.0_r8  ;  zz0(2) = 0.0_r8  ;  zz0(3) = 1.0_r8
lon = lon0 * deg2rad  ;  lat = lat0 * deg2rad  ;  azim = azim0 * deg2rad ; elev = elev0 * deg2rad

rtp(1) = cos(lat) * cos(lon)
rtp(2) = cos(lat) * sin(lon)
rtp(3) = sin(lat)

!  compute unit vector normal to merdion plane through tangent point
call vprod(rtp, zz0, rnm)
vlen0 = sqrt(rnm(1)*rnm(1) + rnm(2)*rnm(2) + rnm(3)*rnm(3))
rnm(:) = rnm(:) / vlen0

!  compute unit vector toward north from perigee point
call vprod(rnm, rtp, rno)
vlen0 = sqrt(rno(1)*rno(1) + rno(2)*rno(2) + rno(3)*rno(3))
rno(:) = rno(:) / vlen0

!  rotate the vector rno around rtp for a single azim to get tangent vector
call spin(rno, rtp, azim, uon)
vlen0 = sqrt(uon(1)*uon(1) + uon(2)*uon(2) + uon(3)*uon(3))
uon(:) = uon(:) / vlen0  

!  compute unit vector normal to tangent ray and rtp
call vprod(uon, rtp, rn)
vlen0 = sqrt(rn(1)*rn(1) + rn(2)*rn(2) + rn(3)*rn(3))
rn(:) = rn(:) / vlen0

!  rotate the vector rno around rn for a single elevation to get ray direction
call spin(uon, rn, elev, ufin)
vlen0 = sqrt(ufin(1)*ufin(1) + ufin(2)*ufin(2) + ufin(3)*ufin(3))
ux = ufin(1) / vlen0  ;  uy = ufin(2) / vlen0  ;  uz = ufin(3) / vlen0


return
end subroutine tanvec01

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   spin - subroutine that rotates vector v1 around vs clockwise by 
!          by a specified angle.
!
!    v1 - vector to rotate
!    vs - vector to rotate about
!     a - angle to rotate v1 around
!    v2 - output vector after rotation
!
!     created Hui Liu NCAR/MMM
!     modified June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spin(v1, vs, a, v2)

use   types_mod, only : r8, deg2rad

implicit none

real(r8), intent(in)  :: v1(3), vs(3), a
real(r8), intent(out) :: v2(3)

real(r8) :: vsabs, vsn(3), a1, a2, a3, s(3,3) 

! Calculation of the unit vector for the rotation
vsabs  = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
vsn(:) = vs(:) / vsabs

! Calculation the rotation matrix
a1 = cos(a)  ;  a2 = 1.0_r8 - a1  ;  a3 = sin(a)
s(1,1) = a2 * vsn(1) * vsn(1) + a1
s(1,2) = a2 * vsn(1) * vsn(2) - a3 * vsn(3)
s(1,3) = a2 * vsn(1) * vsn(3) + a3 * vsn(2)
s(2,1) = a2 * vsn(2) * vsn(1) + a3 * vsn(3)
s(2,2) = a2 * vsn(2) * vsn(2) + a1
s(2,3) = a2 * vsn(2) * vsn(3) - a3 * vsn(1)
s(3,1) = a2 * vsn(3) * vsn(1) - a3 * vsn(2)
s(3,2) = a2 * vsn(3) * vsn(2) + a3 * vsn(1)
s(3,3) = a2 * vsn(3) * vsn(3) + a1

!  Compute the rotated vector
v2(:) = s(:,1) * v1(1) + s(:,2) * v1(2) + s(:,3) * v1(3)

return
end subroutine spin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   vprod - subroutine that computes the vector product of two vectors
!
!    x - first vector to take product of
!    y - second vector to take product of
!    z - vector product
!
!     created Hui Liu NCAR/MMM
!     modified June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vprod(x, y, z)
 
use   types_mod, only : r8

implicit none

real(r8), intent(in)  :: x(3), y(3)
real(r8), intent(out) :: z(3)

z(1) = x(2)*y(3) - x(3)*y(2)
z(2) = x(3)*y(1) - x(1)*y(3)
z(3) = x(1)*y(2) - x(2)*y(1)

return
end subroutine vprod

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
