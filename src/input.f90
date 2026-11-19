module input
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env, only: error_unit, iostat_end
    use hdf5
    use parameters

    private get_dimsize_2d, get_dimsize_3d, &
            read_4flt_1d, read_8flt_1d, &
            read_1int_2d, read_4int_2d, read_4flt_2d, &
            read_4flt_3d, &
            read_4flt_4d, abort_hdf

contains

subroutine analyze_file_name( &
    cgranule, year, month, cmon, ym, imon, &
    obsstr, obsnext, avgstr_y, avgstr, avgend, &
    fname, &
    ooutyr, ooutmon, &
    istat, ofirst, ofinal)

    implicit none

    character(len =   6), intent(out) :: cgranule
    character(len =   4), intent(out) :: year
    character(len =   2), intent(out) :: month
    character(len =   3), intent(out) :: cmon
    character(len =   7), intent(out) :: ym
    character(len =  16), intent(out) :: obsstr, obsnext
    character(len =  16), intent(out) :: avgstr_y, avgstr, avgend
    character(len = 256), intent(out) :: fname
    integer, intent(out) :: imon
    logical, intent(out) :: ooutyr, ooutmon

    integer, intent(inout) :: istat
    logical, intent(inout) :: ofirst, ofinal

    integer :: slash_ix, parse_status, iyear, imon_next, iyear_next
    character(len = 256), save :: fname_next
    character(len =   4) :: year_next, hour, hour_next
    character(len =   2) :: month_next, day, day_next
    character(len =   7) :: ym_next
    character(len =  16) :: obsend

    ooutyr  = .false.
    ooutmon = .false.

    if ( ofirst ) then
       read( input_file_idx, '(A)', iostat=istat ) fname
       if ( istat .eq. iostat_end ) then
          write(error_unit,*) "### No input HDF5 file. [Abort]"
          stop 1
       elseif ( istat .ne. 0 ) then
          write(error_unit,*) "### Filelist read error. IOSTAT =", istat
          stop 1
       endif
       read( input_file_idx, '(A)', iostat=istat ) fname_next
       if ( istat .eq. iostat_end ) then
          ooutyr  = .true.
          ooutmon = .true.
          ofinal  = .true.
       elseif ( istat .ne. 0 ) then
          write(error_unit,*) "### Filelist read error. IOSTAT =", istat
          stop 1
       endif
    else
       fname = trim( fname_next )  ! substitute t-1 older filename
       read( input_file_idx, '(A)', iostat=istat ) fname_next
       if ( istat .eq. iostat_end ) then
          ooutyr  = .true.
          ooutmon = .true.
          ofinal  = .true.
       elseif ( istat .ne. 0 ) then
          write(error_unit,*) "### Filelist read error. IOSTAT =", istat
          stop 1
       endif
    endif
    fname = trim( adjustl(fname) )
    if ( verbose_logging ) write( *,* ) "*** FILENAME (IN) : ", trim( fname )
    slash_ix = scan( fname, "/", .TRUE. )
    cgranule = fname( slash_ix + 47 : slash_ix + 52 )
    year     = "20" // fname( slash_ix + 21 : slash_ix + 22 )
    month    = fname( slash_ix + 23 : slash_ix + 24 )
    ym       = year // month
    day      = fname( slash_ix + 25 : slash_ix + 26 )
    hour     = fname( slash_ix + 28 : slash_ix + 31 )
    read ( year,*,iostat=parse_status ) iyear
    if ( parse_status .ne. 0 ) then
       write(error_unit,*) "### Invalid year in input filename: ", trim(fname)
       stop 1
    endif
    read ( month,*,iostat=parse_status ) imon        ! CHARACTER => INTEGER
    if ( parse_status .ne. 0 ) then
       write(error_unit,*) "### Invalid month in input filename: ", trim(fname)
       stop 1
    endif
    if ( imon .lt. 1 .or. imon .gt. 12 ) then
       write(error_unit,*) "### Invalid month in input filename: ", trim(fname)
       stop 1
    endif
    write( cmon,'(i2.2)' ) imon  ! INTEGER => CHARACTER
    if ( month .eq. '01' ) cmon = "jan"
    if ( month .eq. '02' ) cmon = "feb"
    if ( month .eq. '03' ) cmon = "mar"
    if ( month .eq. '04' ) cmon = "apr"
    if ( month .eq. '05' ) cmon = "may"
    if ( month .eq. '06' ) cmon = "jun"
    if ( month .eq. '07' ) cmon = "jul"
    if ( month .eq. '08' ) cmon = "aug"
    if ( month .eq. '09' ) cmon = "sep"
    if ( month .eq. '10' ) cmon = "oct"
    if ( month .eq. '11' ) cmon = "nov"
    if ( month .eq. '12' ) cmon = "dec"
    if ( ofirst ) then
       obsstr = trim( year ) // trim( month ) // trim( day ) // " " // trim( hour ) // "00"
       obsend = trim( obsstr )
       avgstr = trim( obsstr )
       avgstr_y = trim( obsstr )
       avgend = trim( obsend )
       ofirst = .false.
    else
       obsend = trim( year ) // trim( month ) // trim( day ) // " " // trim( hour ) // "00"
       avgend = trim( obsend )
    endif
    if ( verbose_logging ) write( *,* ) "GRANULE = ", cgranule
    if ( verbose_logging ) write( *,* ) "YYYY = ", year
    if ( verbose_logging ) write( *,* ) "MM =   ", month, " ", imon, " ", cmon
    if ( verbose_logging ) write( *,* ) "YM =   ", ym
    if ( verbose_logging ) write( *,* ) "ISTAT = ", istat
    if ( verbose_logging ) write( *,* ) "OBSSTR = ", obsstr
    if ( verbose_logging ) write( *,* ) "AVGSTR_Y = ", avgstr_y
    if ( ofinal ) then
       if ( verbose_logging ) write(*,*) "OFINAL = .TRUE."
    else
       if ( verbose_logging ) write(*,*) "OFINAL = .FALSE."
       ! next file
       fname_next = trim( adjustl(fname_next) )
       slash_ix   = scan( fname_next, "/", .TRUE. )
       year_next  = "20" // fname_next( slash_ix + 21 : slash_ix + 22 )
       month_next = fname_next( slash_ix + 23 : slash_ix + 24 )
       ym_next    = year_next // month_next
       read ( year_next,*,iostat=parse_status ) iyear_next
       if ( parse_status .ne. 0 ) then
          write(error_unit,*) "### Invalid year in input filename: ", trim(fname_next)
          stop 1
       endif
       read ( month_next,*,iostat=parse_status ) imon_next
       if ( parse_status .ne. 0 ) then
          write(error_unit,*) "### Invalid month in input filename: ", trim(fname_next)
          stop 1
       endif
       if ( imon_next .lt. 1 .or. imon_next .gt. 12 ) then
          write(error_unit,*) "### Invalid month in input filename: ", trim(fname_next)
          stop 1
       endif
       day_next   = fname_next( slash_ix + 25 : slash_ix + 26 )
       hour_next  = fname_next( slash_ix + 28 : slash_ix + 31 )
       obsnext = trim( year_next ) // trim( month_next ) // trim( day_next ) // " " // trim( hour_next ) // "00"
       if ( verbose_logging ) write( *,* ) "*** FILENAME (IN NEXT) : ", trim( obsnext )
       if ( obsnext .le. obsend ) then
          write(error_unit,*) "### Input filelist is not in strictly chronological order."
          write(error_unit,*) "### Current: ", trim(fname)
          write(error_unit,*) "### Next:    ", trim(fname_next)
          stop 1
       endif
    endif
    ! Output flag
    if ( ofinal ) then
       ooutyr  = .true.
       ooutmon = .true.
       if ( verbose_logging ) write(*,*) "OOUTYR = .TRUE."
       if ( verbose_logging ) write(*,*) "OOUTMON = .TRUE."
    else
       if ( year .ne. year_next ) then
          ooutyr = .true.
          if ( verbose_logging ) write(*,*) "OOUTYR = .TRUE."
       else
          if ( verbose_logging ) write(*,*) "OOUTYR = .FALSE."
       endif
       if ( ym .ne. ym_next ) then
          ooutmon = .true.
          if ( verbose_logging ) write(*,*) "OOUTMON = .TRUE."
       else
          if ( verbose_logging ) write(*,*) "OOUTMON = .FALSE."
       endif
    endif

end subroutine analyze_file_name

subroutine read_HDF5_files( &
        varname, fname, &
        nbin, nray, &
        lat, lon, lsflg, height, temp, pres, shum, rhum, &
        cot, cwp, ciwp, dbz, cmask, &
        cldavel, clddvel, cldtvel1, cldtvel2, cldtvelt, &
        cwc1, cwc2, ciwc, cer1, cer2, cerliq, cerice, &
        numc1, numc2, numliq, numice, &
        ccat_cpratl, ctyp_cpratl, cphase1, cphase2 )
    implicit none
    character(*), dimension(ndset), intent(in) :: varname
    character(*), intent(in) :: fname

    integer, intent(out) :: nbin, nray
    real(dp), dimension(:), allocatable, intent(out) :: lat
    real(dp), dimension(:), allocatable, intent(out) :: lon
    integer(kind=4), dimension(:), allocatable, intent(out) :: lsflg  ! 0=water; 1=land
    real(sp), dimension(:,:), allocatable, intent(out) :: height
    real(sp), dimension(:,:), allocatable, intent(out) :: temp
    real(sp), dimension(:,:), allocatable, intent(out) :: pres
    real(sp), dimension(:,:), allocatable, intent(out) :: shum
    real(sp), dimension(:,:), allocatable, intent(out) :: rhum
    real(sp), dimension(:),   allocatable, intent(out) :: cot
    real(sp), dimension(:),   allocatable, intent(out) :: cwp
    real(sp), dimension(:),   allocatable, intent(out) :: ciwp
    real(sp), dimension(:,:), allocatable, intent(out) :: dbz
!    integer(kind=1), dimension(:,:), allocatable, intent(out) :: cmask
    integer(kind=4), dimension(:,:), allocatable, intent(out) :: cmask
    real(sp), dimension(:,:), allocatable, intent(out) :: cldavel
    real(sp), dimension(:,:), allocatable, intent(out) :: clddvel
    real(sp), dimension(:,:), allocatable, intent(out) :: cldtvel1
    real(sp), dimension(:,:), allocatable, intent(out) :: cldtvel2
    real(sp), dimension(:,:), allocatable, intent(out) :: cldtvelt
    real(sp), dimension(:,:), allocatable, intent(out) :: cwc1
    real(sp), dimension(:,:), allocatable, intent(out) :: cwc2
    real(sp), dimension(:,:), allocatable, intent(out) :: ciwc
    real(sp), dimension(:,:), allocatable, intent(out) :: cer1
    real(sp), dimension(:,:), allocatable, intent(out) :: cer2
    real(sp), dimension(:,:), allocatable, intent(out) :: cerliq
    real(sp), dimension(:,:), allocatable, intent(out) :: cerice
    real(sp), dimension(:,:), allocatable, intent(out) :: numc1
    real(sp), dimension(:,:), allocatable, intent(out) :: numc2
    real(sp), dimension(:,:), allocatable, intent(out) :: numliq
    real(sp), dimension(:,:), allocatable, intent(out) :: numice
    integer(kind=4), dimension(:,:), allocatable, intent(out) :: ccat_cpratl
!    integer(kind=1), dimension(:,:), allocatable, intent(out) :: ctyp_cpratl
    integer(kind=4), dimension(:,:), allocatable, intent(out) :: ctyp_cpratl
!    integer(kind=1), dimension(:,:), allocatable, intent(out) :: cphase1
    integer(kind=4), dimension(:,:), allocatable, intent(out) :: cphase1
!    integer(kind=1), dimension(:,:), allocatable, intent(out) :: cphase2
    integer(kind=4), dimension(:,:), allocatable, intent(out) :: cphase2

    integer(HSIZE_T), dimension(1) :: dim1d   ! dim1d(1)=nray
    integer(HSIZE_T), dimension(2) :: dim2d   ! dim2d(1)=nbin, dim2d(2)=nray
    integer(HID_T) :: file_id                 ! shared file identifier for one granule
    integer :: hdferr                         ! HDF5 status
    integer :: astat                          ! status ID from data assignment

    !--- Open HDF5 once for this granule
    call h5open_f( hdferr )
    if ( hdferr .lt. 0 ) then
       write(error_unit,*) "### read_HDF5_files: H5OPEN FAILED", hdferr, trim(fname)
       error stop 1
    endif
    call h5fopen_f( fname, H5F_ACC_RDONLY_F, file_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5FOPEN", hdferr, fname, "" )
    endif

    !--- Retrieve the Dimensions of a Datafield in the HDF5
    call get_dimsize_2d( file_id, fname, varname(4), & ! in
                        & dim2d )             ! out
    nbin = dim2d(1)
    nray = dim2d(2)
    dim1d(1) = int(nray, kind=HSIZE_T)
    !if ( nbin .gt. nray ) then
    if ( nbin .ne. 206 ) then
       write(error_unit,*) "### MISMATCHED TO RETRIEVE DIMSIZES:"
       write(error_unit,*) "    NBIN, NRAY = ", nbin, nray
       write(error_unit,*) "    DIM2D(1) and DIM2D(2) = ", dim2d(1), dim2d(2)
       write(error_unit,*) "### ABORT"
       stop
    else
       if ( verbose_logging ) write( *,* ) "*** NBIN, NRAY = ", nbin, nray
       if ( verbose_logging ) write( *,* ) "    DIM2D(1) and DIM2D(2) = ", dim2d(1), dim2d(2)
    endif

    !--- Get Geolocation Data (/ScienceData/Geo/)
    ! latitude
    allocate( lat( nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE lat", astat
       stop
    endif
    call read_8flt_1d( file_id, fname, trim(varname(1)), dim1d, &  ! in
                     & lat )                              ! out

    ! longitude
    allocate( lon( nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE lon", astat
       stop
    endif
    call read_8flt_1d( file_id, fname, trim(varname(2)), dim1d, &  ! in
                     & lon )                              ! out

    ! land_water_flag
    allocate( lsflg( nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE lsflg", astat
       stop
    endif
    call read_4int_1d( file_id, fname, trim(varname(3)), dim1d, &  ! in
                     & lsflg )                            ! out

    ! height
    allocate( height( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE height", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(4)), dim2d, &  ! in
                     & height )                           ! out
    height(:,:) = height(:,:) * 1e-3_dp  ! Unit: m => km
    where (.not. ieee_is_finite(height))
       height = undef 
    end where

    !--- Retrieve the Specific Data (/ScienceData/Data/)
    ! "GRID_temperature_1km"
    allocate( temp( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE temp", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(5)), dim2d, &  ! in
                     & temp )                             ! out

    ! "GRID_pressure_1km"
    allocate( pres( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE pres", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(6)), dim2d, &  ! in
                     & pres )                             ! out
    !pres(:,:) = pres(:,:) * 1e-2_dp  ! Unit: Pa => hPa
                                      ! Note: Unit is already "hPa"

    ! "GRID_specificHumidity_1km"
    allocate( shum( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE shum", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(7)), dim2d, &  ! in
                     & shum )                             ! out

    ! "GRID_relativeHumidity_1km"
    allocate( rhum( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE rhum", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(8)), dim2d, &  ! in
                     & rhum )                             ! out

    ! "optical_thickness_1km" or "optical_thickness_10km"
    allocate( cot( nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cot", astat
       stop
    endif
    call read_4flt_1d( file_id, fname, trim(varname(9)), dim1d, &  ! in
                     & cot )                              ! out

    ! "cloud_water_path_1km" or "cloud_water_path_10km"
    allocate( cwp( nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cwp", astat
       stop
    endif
    call read_4flt_1d( file_id, fname, trim(varname(10)), dim1d, &  ! in
                     & cwp )                              ! out

    ! "cloud_ice_water_path_1km" or "cloud_ice_water_path_10km"
    allocate( ciwp( nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE ciwp", astat
       stop
    endif
    call read_4flt_1d( file_id, fname, trim(varname(11)), dim1d, &  ! in
                     & ciwp )                              ! out

    ! "cloud_radar_reflectivity_1km" or "cloud_radar_reflectivity_10km"
    allocate( dbz( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE dbz", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(12)), dim2d, &  ! in
                     & dbz )                               ! out
    where (.not. ieee_is_finite(dbz))
       dbz = undef 
    end where

    ! "cloud_mask_cpr_atlid_1km" or "cloud_mask_cpr_atlid_10km"
    allocate( cmask( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cmask", astat
       stop
    endif
!    call read_1int_2d( file_id, fname, trim(varname(13)), dim2d, &  ! in
!                     & cmask )                             ! out
    call read_4int_2d( file_id, fname, trim(varname(13)), dim2d, &  ! in
                     & cmask )                             ! out

    ! "cloud_air_velocity_1km" or "cloud_air_velocity_10km"
    allocate( cldavel( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cldavel", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(14)), dim2d, &  ! in
                     & cldavel )                           ! out
    where ( .not. ieee_is_finite(cldavel) )
       cldavel = undef
    end where
    where ( cldavel < -100._sp .or. cldavel > 100._sp )
       cldavel = undef
    end where

    ! "cloud_doppler_velocity_1km" or "cloud_doppler_velocity_10km"
    allocate( clddvel( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE clddvel", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(15)), dim2d, &  ! in
                     & clddvel )                           ! out
    where ( .not. ieee_is_finite(clddvel) )
       clddvel = undef
    end where
    where ( clddvel < -100._sp .or. clddvel > 100._sp )
       clddvel = undef
    end where

    ! "cloud_terminal_velocity1_1km" or "cloud_terminal_velocity1_10km"
    allocate( cldtvel1( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cldtvel1", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(16)), dim2d, &  ! in
                     & cldtvel1 )                          ! out
    where ( .not. ieee_is_finite(cldtvel1) )
       cldtvel1 = undef
    end where
    where ( cldtvel1 < -100._sp .or. cldtvel1 > 100._sp )
       cldtvel1 = undef
    end where

    ! "cloud_terminal_velocity2_1km" or "cloud_terminal_velocity2_10km"
    allocate( cldtvel2( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cldtvel2", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(17)), dim2d, &  ! in
                     & cldtvel2 )                          ! out
    where ( .not. ieee_is_finite(cldtvel2) )
       cldtvel2 = undef
    end where
    where ( cldtvel2 < -100._sp .or. cldtvel2 > 100._sp )
       cldtvel2 = undef
    end where

    ! "total_cloud_terminal_velocity_1km" or "total_cloud_terminal_velocity_10km"
    allocate( cldtvelt( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cldtvelt", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(18)), dim2d, &  ! in
                     & cldtvelt )                          ! out
    where ( .not. ieee_is_finite(cldtvelt) )
       cldtvelt = undef
    end where
    where ( cldtvelt < -100._sp .or. cldtvelt > 100._sp )
       cldtvelt = undef
    end where

    ! "cloud_water_content1_1km" or "cloud_water_content1_10km"
    allocate( cwc1( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cwc1", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(19)), dim2d, &  ! in
                     & cwc1 )                              ! out

    ! "cloud_water_content2_1km" or "cloud_water_content2_10km"
    allocate( cwc2( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cwc2", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(20)), dim2d, &  ! in
                     & cwc2 )                              ! out

    ! "cloud_ice_content_1km" or "cloud_ice_content_10km"
    allocate( ciwc( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE ciwc", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(21)), dim2d, &  ! in
                     & ciwc )                              ! out

    ! "cloud_effective_radius1_1km" or "cloud_effective_radius1_10km"
    allocate( cer1( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cer1", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(22)), dim2d, &  ! in
                     & cer1 )                              ! out

    ! "cloud_effective_radius2_1km" or "cloud_effective_radius2_10km"
    allocate( cer2( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cer2", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(23)), dim2d, &  ! in
                     & cer2 )                              ! out

    ! "cloud_water_effective_radius_1km" or "cloud_water_effective_radius_10km"
    allocate( cerliq( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cerliq", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(24)), dim2d, &  ! in
                     & cerliq )                            ! out

    ! "cloud_ice_effective_radius_1km" or "cloud_ice_effective_radius_10km"
    allocate( cerice( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cerice", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(25)), dim2d, &  ! in
                     & cerice )                            ! out

    ! "cloud_number_concentration1_1km" or "cloud_number_concentration1_10km"
    allocate( numc1( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE numc1", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(26)), dim2d, &  ! in
                     & numc1 )                             ! out
    numc1(:,:) = numc1(:,:) * 1e-6_sp  ! Unit: # m-3 => # cm-3

    ! "cloud_number_concentration2_1km" or "cloud_number_concentration2_10km"
    allocate( numc2( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE numc2", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(27)), dim2d, &  ! in
                     & numc2 )                             ! out
    numc2(:,:) = numc2(:,:) * 1e-6_sp  ! Unit: # m-3 => # cm-3

    ! "total_cloud_water_number_concentration_1km" or "total_cloud_water_number_concentration_10km"
    allocate( numliq( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE numliq", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(28)), dim2d, &  ! in
                     & numliq )                            ! out
    numliq(:,:) = numliq(:,:) * 1e-6_sp  ! Unit: # m-3 => # cm-3

    ! "total_cloud_ice_number_concentration_1km" or "total_cloud_ice_number_concentration_10km"
    allocate( numice( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE numice", astat
       stop
    endif
    call read_4flt_2d( file_id, fname, trim(varname(29)), dim2d, &  ! in
                     & numice )                            ! out
    numice(:,:) = numice(:,:) * 1e-6_sp  ! Unit: # m-3 => # cm-3

    ! "cloud_particle_category_cpr_atlid_1km" or "cloud_particle_category_cpr_atlid_10km"
    allocate( ccat_cpratl( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE ccat_cpratl", astat
       stop
    endif
    call read_4int_2d( file_id, fname, trim(varname(30)), dim2d, &  ! in
                     & ccat_cpratl )                       ! out
    where ( ccat_cpratl .eq. -2147483647 .or. ccat_cpratl .lt. 0 .or. ccat_cpratl .gt. 12 )
       ccat_cpratl = undef_i4
    end where

    ! "cloud_particle_type_cpr_atlid_1km" or "cloud_particle_type_cpr_atlid_10km"
    allocate( ctyp_cpratl( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE ctyp_cpratl", astat
       stop
    endif
!    call read_1int_2d( file_id, fname, trim(varname(31)), dim2d, &  ! in
!                     & ctyp_cpratl )                       ! out
    call read_4int_2d( file_id, fname, trim(varname(31)), dim2d, &  ! in
                     & ctyp_cpratl )                       ! out

    ! "cloud_phase1_1km" or "cloud_phase1_10km"
    allocate( cphase1( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cphase1", astat
       stop
    endif
!    call read_1int_2d( file_id, fname, trim(varname(32)), dim2d, &  ! in
!                     & cphase1 )                           ! out
    call read_4int_2d( file_id, fname, trim(varname(32)), dim2d, &  ! in
                     & cphase1 )                           ! out

    ! "cloud_phase2_1km" or "cloud_phase2_10km"
    allocate( cphase2( nbin,nray ), stat=astat )
    if ( astat /= 0 ) then
       write(error_unit,*) "### CANNOT ALLOCATE cphase2", astat
       stop
    endif
!    call read_1int_2d( file_id, fname, trim(varname(33)), dim2d, &  ! in
!                     & cphase2 )                           ! out
    call read_4int_2d( file_id, fname, trim(varname(33)), dim2d, &  ! in
                     & cphase2 )                           ! out

    !--- Close the shared granule file and HDF5 interface
    call h5fclose_f( file_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5FCLOSE", hdferr, fname, "", file_id=file_id )
    endif
    call h5close_f( hdferr )
    if ( hdferr .lt. 0 ) then
       write(error_unit,*) "### read_HDF5_files: H5CLOSE FAILED", hdferr, trim(fname)
       error stop 1
    endif

end subroutine read_HDF5_files

subroutine get_dimsize_2d( file_id, filename, dataname, & ! in
                         & dims )                ! out
    implicit none

    integer hdferr
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dataname
    integer(hid_t), intent(in) :: file_id ! Shared file identifier
    integer(hid_t) :: data_space_id ! data space identifier
    integer(hid_t) :: dset_id       ! Dataset identifier

    integer(hsize_t), dimension(2)              :: maxdims
    integer(hsize_t), dimension(2), intent(out) :: dims


    ! Access to Dataset
    call h5dopen_f( file_id, dataname, dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DOPEN", hdferr, filename, dataname, file_id=file_id )
    endif

    ! Get Data Space
    call h5dget_space_f( dset_id, data_space_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DGET_SPACE", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Get the Size of Dataset
    call h5sget_simple_extent_dims_f( data_space_id, dims, maxdims, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5SGET_SIMPLE_EXTENT_DIMS", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id, space_id=data_space_id )
    endif
    if ( verbose_logging ) write( *,* ) "*** DIMENSION SIZE (2D): ", dims(1:2)
    if ( verbose_logging ) write( *,* ) "                MAXDIMS: ", maxdims(1:2)

    ! Close Data Space
    call h5sclose_f( data_space_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5SCLOSE", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Close the Datasets
    call h5dclose_f( dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DCLOSE", hdferr, filename, dataname, file_id=file_id )
    endif

    return
end subroutine get_dimsize_2d


subroutine get_dimsize_3d( file_id, filename, dataname, & ! in
                         & dims )                ! out
    implicit none

    integer hdferr
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dataname
    integer(hid_t), intent(in) :: file_id ! Shared file identifier
    integer(hid_t) :: data_space_id ! data space identifier
    integer(hid_t) :: dset_id       ! Dataset identifier

    integer(hsize_t), dimension(3)              :: maxdims
    integer(hsize_t), dimension(3), intent(out) :: dims


    ! Access to Dataset
    call h5dopen_f( file_id, dataname, dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DOPEN", hdferr, filename, dataname, file_id=file_id )
    endif

    ! Get Data Space
    call h5dget_space_f( dset_id, data_space_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DGET_SPACE", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Get the Size of Dataset
    call h5sget_simple_extent_dims_f( data_space_id, dims, maxdims, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5SGET_SIMPLE_EXTENT_DIMS", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id, space_id=data_space_id )
    endif
    if ( verbose_logging ) write( *,* ) "*** DIMENSION SIZE (3D): ", dims(1:3)
    if ( verbose_logging ) write( *,* ) "                MAXDIMS: ", maxdims(1:3)

    ! Close Data Space
    call h5sclose_f( data_space_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5SCLOSE", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Close the Datasets
    call h5dclose_f( dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DCLOSE", hdferr, filename, dataname, file_id=file_id )
    endif

    return
end subroutine get_dimsize_3d

subroutine read_4int_1d( file_id, filename, dataname, dims, & ! in
                       & data )                      ! out
    implicit none

    integer hdferr
    integer(hid_t), intent(in) :: file_id ! Shared file identifier
    integer(hid_t) :: dset_id  ! Dataset identifier

    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: dataname
    integer(hsize_t), intent(in)  :: dims(1)
    integer(kind=4),  intent(inout) :: data(:) ! Data buffers

    ! Access to Dataset
    call h5dopen_f( file_id, dataname, dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DOPEN", hdferr, filename, dataname, file_id=file_id )
    endif
    if ( verbose_logging ) write( *,* ) "*** ACCESS TO ... ", dataname

    ! Read Data
    call h5dread_f( dset_id, H5T_NATIVE_INTEGER, data, dims, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DREAD", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Close the Datasets
    call h5dclose_f( dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DCLOSE", hdferr, filename, dataname, file_id=file_id )
    endif

end subroutine read_4int_1d

subroutine read_4flt_1d( file_id, filename, dataname, dims, & ! in
                       & data )                      ! out
    implicit none

    integer hdferr
    integer(hid_t), intent(in) :: file_id ! Shared file identifier
    integer(hid_t) :: dset_id  ! Dataset identifier

    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: dataname
    integer(hsize_t), intent(in)  :: dims(1)
    real(sp),         intent(inout) :: data(:) ! Data buffers

    ! Access to Dataset
    call h5dopen_f( file_id, dataname, dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DOPEN", hdferr, filename, dataname, file_id=file_id )
    endif
    if ( verbose_logging ) write( *,* ) "*** ACCESS TO ... ", dataname

    ! Read Data
    call h5dread_f( dset_id, H5T_NATIVE_REAL, data, dims, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DREAD", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Close the Datasets
    call h5dclose_f( dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DCLOSE", hdferr, filename, dataname, file_id=file_id )
    endif

end subroutine read_4flt_1d

subroutine read_8flt_1d( file_id, filename, dataname, dims, & ! in
                       & data )                      ! out
    implicit none

    integer hdferr
    integer(hid_t), intent(in) :: file_id ! Shared file identifier
    integer(hid_t) :: dset_id  ! Dataset identifier

    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: dataname
    integer(hsize_t), intent(in)  :: dims(1)
    real(dp),         intent(inout) :: data(:) ! Data buffers

    ! Access to Dataset
    call h5dopen_f( file_id, dataname, dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DOPEN", hdferr, filename, dataname, file_id=file_id )
    endif
    if ( verbose_logging ) write( *,* ) "*** ACCESS TO ... ", dataname

    ! Read Data
    call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DREAD", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Close the Datasets
    call h5dclose_f( dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DCLOSE", hdferr, filename, dataname, file_id=file_id )
    endif

end subroutine read_8flt_1d

subroutine read_1int_2d( file_id, filename, dataname, dims, & ! in
                       & data )                      ! out
    implicit none

    integer hdferr
    integer(hid_t), intent(in) :: file_id ! Shared file identifier
    integer(hid_t) :: dset_id  ! Dataset identifier

    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: dataname
    integer(hsize_t), intent(in)  :: dims(2)
    integer(kind=1),  intent(out) :: data(:,:) ! Data buffers

    ! Access to Dataset
    call h5dopen_f( file_id, dataname, dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DOPEN", hdferr, filename, dataname, file_id=file_id )
    endif
    if ( verbose_logging ) write( *,* ) "*** ACCESS TO ... ", dataname

    ! Read Data
    call h5dread_f( dset_id, H5T_NATIVE_INTEGER, data, dims, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DREAD", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Close the Datasets
    call h5dclose_f( dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DCLOSE", hdferr, filename, dataname, file_id=file_id )
    endif

end subroutine read_1int_2d

subroutine read_4int_2d( file_id, filename, dataname, dims, & ! in
                       & data )                      ! out
    implicit none

    integer hdferr
    integer(hid_t), intent(in) :: file_id ! Shared file identifier
    integer(hid_t) :: dset_id  ! Dataset identifier

    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: dataname
    integer(hsize_t), intent(in)  :: dims(2)
    integer(kind=4),  intent(out) :: data(:,:) ! Data buffers

    ! Access to Dataset
    call h5dopen_f( file_id, dataname, dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DOPEN", hdferr, filename, dataname, file_id=file_id )
    endif
    if ( verbose_logging ) write( *,* ) "*** ACCESS TO ... ", dataname

    ! Read Data
    call h5dread_f( dset_id, H5T_NATIVE_INTEGER, data, dims, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DREAD", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Close the Datasets
    call h5dclose_f( dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DCLOSE", hdferr, filename, dataname, file_id=file_id )
    endif

end subroutine read_4int_2d

subroutine read_4flt_2d( file_id, filename, dataname, dims, & ! in
                       & data )                      ! out
    implicit none

    integer hdferr
    integer(hid_t), intent(in) :: file_id ! Shared file identifier
    integer(hid_t) :: dset_id  ! Dataset identifier

    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: dataname
    integer(hsize_t), intent(in)  :: dims(2)
    real(sp),         intent(inout) :: data(:,:) ! Data buffers

    ! Access to Dataset
    call h5dopen_f( file_id, dataname, dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DOPEN", hdferr, filename, dataname, file_id=file_id )
    endif
    if ( verbose_logging ) write( *,* ) "*** ACCESS TO ... ", dataname

    ! Read Data
    call h5dread_f( dset_id, H5T_NATIVE_REAL, data, dims, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DREAD", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Close the Datasets
    call h5dclose_f( dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DCLOSE", hdferr, filename, dataname, file_id=file_id )
    endif

end subroutine read_4flt_2d

subroutine read_4flt_3d( file_id, filename, dataname, dims, & ! in
                       & data )                      ! out
    implicit none

    integer hdferr
    integer(hid_t), intent(in) :: file_id ! Shared file identifier
    integer(hid_t) :: dset_id  ! Dataset identifier

    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: dataname
    integer(hsize_t), intent(in)  :: dims(3)
    real(sp),         intent(inout) :: data(:,:,:) ! Data buffers

    ! Access to Dataset
    call h5dopen_f( file_id, dataname, dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DOPEN", hdferr, filename, dataname, file_id=file_id )
    endif
    if ( verbose_logging ) write( *,* ) "*** ACCESS TO ... ", dataname

    ! Read Data
    call h5dread_f( dset_id, H5T_NATIVE_REAL, data, dims, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DREAD", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Close the Datasets
    call h5dclose_f( dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DCLOSE", hdferr, filename, dataname, file_id=file_id )
    endif

end subroutine read_4flt_3d

subroutine read_4flt_4d( file_id, filename, dataname, dims, & ! in
                       & data )                      ! out
    implicit none

    integer hdferr
    integer(hid_t), intent(in) :: file_id ! Shared file identifier
    integer(hid_t) :: dset_id  ! Dataset identifier

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dataname
    integer(hsize_t), intent(in) :: dims(4)
    real(sp),         intent(inout) :: data(:,:,:,:) ! Data buffers

    ! Access to Dataset
    call h5dopen_f( file_id, dataname, dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DOPEN", hdferr, filename, dataname, file_id=file_id )
    endif
    if ( verbose_logging ) write( *,* ) "*** ACCESS TO ... ", dataname

    ! Read Data
    call h5dread_f( dset_id, H5T_NATIVE_REAL, data, dims, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DREAD", hdferr, filename, dataname, &
                       file_id=file_id, dset_id=dset_id )
    endif

    ! Close the Datasets
    call h5dclose_f( dset_id, hdferr )
    if ( hdferr .lt. 0 ) then
       call abort_hdf( "H5DCLOSE", hdferr, filename, dataname, file_id=file_id )
    endif

end subroutine read_4flt_4d

subroutine abort_hdf( operation, hdferr, filename, dataname, &
                      file_id, dset_id, space_id )
    implicit none
    character(len=*), intent(in) :: operation, filename, dataname
    integer, intent(in) :: hdferr
    integer(HID_T), intent(in), optional :: file_id, dset_id, space_id
    integer :: cleanup_err

    write(error_unit,*) "### HDF5 ERROR"
    write(error_unit,*) "    operation = ", trim(operation)
    write(error_unit,*) "    error     = ", hdferr
    write(error_unit,*) "    file      = ", trim(filename)
    if ( len_trim(dataname) .gt. 0 ) then
       write(error_unit,*) "    dataset   = ", trim(dataname)
    endif

    if ( present(space_id) ) then
       if ( space_id .ge. 0 ) call h5sclose_f( space_id, cleanup_err )
    endif
    if ( present(dset_id) ) then
       if ( dset_id .ge. 0 ) call h5dclose_f( dset_id, cleanup_err )
    endif
    if ( present(file_id) ) then
       if ( file_id .ge. 0 ) call h5fclose_f( file_id, cleanup_err )
    endif
    call h5close_f( cleanup_err )
    error stop 1
end subroutine abort_hdf

end module input
