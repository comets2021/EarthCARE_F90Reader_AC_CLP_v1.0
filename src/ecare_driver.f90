!=======================================================================!
  !--- Purpose:
  !      * Create Level3 product from the EarthCARE Level2 granule files:
  !        - Product Name: CPR+ATLID/2B.AC_CLP
  !
  !--- Software Required:
  !      * HDF5 Library
  !        (https://portal.hdfgroup.org/display/support)
  !      * NetCDF4 Library
  !        (https://www.unidata.ucar.edu/software/netcdf/)
  !        - netCDF:         Version 4.9.2 or later
  !        - netCDF-Fortran: Version 4.6.1 or later
  !
  !--- History:
  !      * Dec 15, 2025 - Takuro Michibata: Test Run
  !      * Jan 16, 2026 - Takuro Michibata: Modify P-Level Output
  !      * Jan 31, 2026 - Takuro Michibata: Add ISCCP-like Histograms
  !      * Feb 20, 2026 - Takuro Michibata: Add Cloud Fraction
  !      * Jul  6, 2026 - Takuro Michibata: Update for vCb product
  !      * Nov 19, 2026 - Takuro Michibata: First Release
  !
  !--- Contact:
  !    michibata@riam.kyushu-u.ac.jp
  !
!=======================================================================!
program ecare_main

   use, intrinsic :: iso_fortran_env, only: error_unit
   use parameters
   use numbers
   use input
   use output

   implicit none

   integer :: imon
   integer :: n_arg, istat, arg_status, filelist_length
   integer, parameter :: filelist_max_length = 4096
   integer :: env_status, env_length, files_processed
   integer(kind=8) :: tick_start, tick_now, tick_rate
   real(dp) :: elapsed_seconds
   character(len = 32) :: log_environment
   character(len = filelist_max_length) :: filelist
   character(len =   3) :: version
   character(len = 128), dimension(ndset) :: varname
   character(len =   6) :: cgranule
   character(len = 256) :: fname
   character(len =   4) :: year
   character(len =   2) :: month
   character(len =   3) :: cmon
   character(len =   7) :: ym
   character(len =  16) :: obsstr, obsnext
   character(len =  16) :: avgstr_y, avgstr, avgend
   character(len =   7) :: processing_start_month, annual_start_month

   logical :: ofirst, ooutyr, ooutmon, ofinal

   integer :: nray                                    ! # of profiles
   integer :: nbin                                    ! # of vertical bins

   real(dp), dimension(:), allocatable :: lat         ! Latitude
   real(dp), dimension(:), allocatable :: lon         ! Longitude
   integer(kind=4), dimension(:), allocatable :: lsflg ! Land/Sea Flag (0=water; 1=land)
   real(sp), dimension(:,:), allocatable :: height    ! Height [m]
   real(sp), dimension(:,:), allocatable :: temp      ! Temperature [K]
   real(sp), dimension(:,:), allocatable :: pres      ! Pressure [Pa]
   real(sp), dimension(:,:), allocatable :: shum      ! Specific Humidity [kg kg-1]
   real(sp), dimension(:,:), allocatable :: rhum      ! Relative Humidity [%]
   real(sp), dimension(:),   allocatable :: cot       ! Optical Thickness for Water+Ice Clouds
   real(sp), dimension(:),   allocatable :: cwp       ! Liquid Water Path [g m-2]
   real(sp), dimension(:),   allocatable :: ciwp      ! Ice Water Path [g m-2]
   real(sp), dimension(:,:), allocatable :: dbz       ! Attenuated Radar Reflectivity [dBZ]
!   integer(kind=1),  dimension(:,:), allocatable :: cmask     ! C4 Cloud Mask (CPR or ATLID) 
   integer(kind=4),  dimension(:,:), allocatable :: cmask     ! C4 Cloud Mask (CPR or ATLID) 
   real(sp), dimension(:,:), allocatable :: cldavel   ! Air Vertical Velocity [m s-1]
   real(sp), dimension(:,:), allocatable :: clddvel   ! Cloud Doppler Velocity [m s-1]
   real(sp), dimension(:,:), allocatable :: cldtvel1  ! Cloud Terminal Velocity 1 [m s-1]
   real(sp), dimension(:,:), allocatable :: cldtvel2  ! Cloud Terminal Velocity 2 [m s-1]
   real(sp), dimension(:,:), allocatable :: cldtvelt  ! Total Cloud Terminal Velocity [m s-1]
                                                        ! => Sedimentation velocity 
   real(sp), dimension(:,:), allocatable :: cwc1      ! Cloud Water Content 1 [g m-3]
   real(sp), dimension(:,:), allocatable :: cwc2      ! Cloud Water Content 2 [g m-3]
   real(sp), dimension(:,:), allocatable :: ciwc      ! Ice Water Content [g m-3]
   real(sp), dimension(:,:), allocatable :: cer1      ! Cloud Effective Radius 1 [micron]
   real(sp), dimension(:,:), allocatable :: cer2      ! Cloud Effective Radius 2 [micron]
   real(sp), dimension(:,:), allocatable :: cerliq    ! Cloud Effective Radius Liq [micron]
   real(sp), dimension(:,:), allocatable :: cerice    ! Cloud Effective Radius Ice [micron]
   real(sp), dimension(:,:), allocatable :: numc1     ! Cloud Number Concentration 1 [# m-3]
   real(sp), dimension(:,:), allocatable :: numc2     ! Cloud Number Concentration 2 [# m-3]
   real(sp), dimension(:,:), allocatable :: numliq    ! Cloud Number Concentration Liq [# m-3]
   real(sp), dimension(:,:), allocatable :: numice    ! Cloud Number Concentration Ice [# m-3]
   integer(kind=4), dimension(:,:), allocatable :: ccat_cpratl ! Cloud Particle Category CPR+ATLID
!   integer(kind=1), dimension(:,:), allocatable :: ctyp_cpratl ! Cloud Particle Type CPR+ATLID
   integer(kind=4), dimension(:,:), allocatable :: ctyp_cpratl ! Cloud Particle Type CPR+ATLID
!   integer(kind=1), dimension(:,:), allocatable :: cphase1   ! Cloud Phase 1
   integer(kind=4), dimension(:,:), allocatable :: cphase1   ! Cloud Phase 1
!   integer(kind=1), dimension(:,:), allocatable :: cphase2   ! Cloud Phase 2
   integer(kind=4), dimension(:,:), allocatable :: cphase2   ! Cloud Phase 2

   !-----------------------------------------------------------------!
   !  INLINE GENERATORS
   !-----------------------------------------------------------------!
   !--- Main Routine
   real               :: glon( ilon_max + 1 )         ! grid value (longitude)
   real               :: glat( ilat_max + 1 )         ! grid value (latitude)
   real               :: gplev( dimz + 1 )            ! grid value (P-levels)
   real               :: gdbz( idbz_max + 1 )         ! grid value (CFAD dBZ)
   real               :: glev( ilev_max + 1 )         ! grid value (CFAD Lev)
   integer            :: ilon, ilat, iplev            ! grid boundary number for lon/lat/plev
   integer            :: idbz, ilev                   ! dimension number for dbz/lev for CFAD
   real(sp), dimension( dimx,dimy ) :: regmap

   !************************************************************************************************
   !--- I/O Settings (Arguments, Filename, Path, etc.)
   ! Consistensy Check for the Arguments
   n_arg = command_argument_count()
   if( n_arg /= 2 ) then
      write( error_unit,* ) "### Wrong Arguments!"
      write( error_unit,* ) "Usage: ./ecare_driver [FILELIST] [VERSION]"
      stop 1
   endif

   ! Normal logging is the default.
   log_environment = ""
   call get_environment_variable( "ECARE_VERBOSE", log_environment, &
                                  length=env_length, status=env_status )
   if ( env_status .eq. 0 .and. env_length .gt. 0 ) then
      select case ( trim(adjustl(log_environment)) )
      case ( "1", "true", "TRUE", "yes", "YES", "on", "ON", "verbose", "VERBOSE" )
         verbose_logging = .true.
      case default
         verbose_logging = .false.
      end select
   endif

   ! Get Arguments for Reading HDF5 Filelist
   call get_command_argument( 1, length=filelist_length, status=arg_status )
   if ( arg_status .ne. 0 ) then
      write(error_unit,*) "### Cannot determine filelist argument length. STATUS =", arg_status
      stop 1
   endif
   if ( filelist_length .le. 0 ) then
      write(error_unit,*) "### Empty filelist argument."
      stop 1
   endif
   if ( filelist_length .gt. len(filelist) ) then
      write(error_unit,*) "### Filelist argument is too long."
      write(error_unit,*) "### Length / maximum =", filelist_length, len(filelist)
      stop 1
   endif
   call get_command_argument( 1, value=filelist, status=arg_status )
   if ( arg_status .ne. 0 ) then
      write(error_unit,*) "### Cannot read filelist argument. STATUS =", arg_status
      stop 1
   endif
   call get_command_argument( 2, value=version, status=arg_status )
   if ( arg_status .ne. 0 ) then
      write(error_unit,*) "### Cannot read version argument. STATUS =", arg_status
      stop 1
   endif
   open(input_file_idx, file = trim(filelist), status = "old", action = "read", iostat = istat )
   if ( istat .ne. 0 ) then
      write(error_unit,*) "### Cannot open input filelist. IOSTAT =", istat
      stop 1
   endif
   files_processed = 0
   call system_clock( count=tick_start, count_rate=tick_rate )
   if ( .not. verbose_logging ) then
      write(*,'(A)') "EarthCARE F90Reader: normal logging (set ECARE_VERBOSE=1 for details)"
   endif

   ! Initialized (Monthly Output)
   n_obs2d          = 0.0_dp;  n_obs3d          = 0.0_dp;  n_cld3d          = 0.0_dp
   n_c4_valid       = 0;       n_c4_cloud       = 0;       n_c4_mask_invalid = 0
   n_c4_type        = 0;       n_c4_category    = 0;       n_c4_type_qc = 0
   n_c4_category_invalid = 0
   n_clr            = 0;       n_cld            = 0;       n_drz            = 0;       n_prc = 0
   n_lwp            = 0;       sum_lwp          = 0.0_dp
   n_iwp            = 0;       sum_iwp          = 0.0_dp
   n_cot            = 0;       sum_cot          = 0.0_dp
   n_cth            = 0;       sum_cth          = 0.0_dp
   n_cth_c4         = 0;       sum_cth_c4       = 0.0_dp
   n_ctp_c4         = 0;       sum_ctp_c4       = 0.0_dp
   n_pctgr          = 0;       n_ptype          = 0
   sum_temp         = 0.0_dp;  sum_shum         = 0.0_dp;  sum_rhum         = 0.0_dp
   n_cldavel        = 0;       sum_cldavel      = 0.0_dp
   n_clddvel        = 0;       sum_clddvel      = 0.0_dp
   n_cldtvel1       = 0;       sum_cldtvel1     = 0.0_dp
   n_cldtvel2       = 0;       sum_cldtvel2     = 0.0_dp
   n_cldtvelt       = 0;       sum_cldtvelt     = 0.0_dp
   n_cwc1           = 0;       sum_cwc1         = 0.0_dp
   n_cwc2           = 0;       sum_cwc2         = 0.0_dp
   n_ciwc           = 0;       sum_ciwc         = 0.0_dp
   n_cer1           = 0;       sum_cer1         = 0.0_dp
   n_cer2           = 0;       sum_cer2         = 0.0_dp
   n_cerliq         = 0;       sum_cerliq       = 0.0_dp
   n_cerice         = 0;       sum_cerice       = 0.0_dp
   n_numc1          = 0;       sum_numc1        = 0.0_dp
   n_numc2          = 0;       sum_numc2        = 0.0_dp
   n_numliq         = 0;       sum_numliq       = 0.0_dp
   n_numice         = 0;       sum_numice       = 0.0_dp
   n_cfad_obs       = 0;       n_cfad_bin       = 0
   n_isccp          = 0;       sum_isccpdv      = 0
   n_isccp_c4       = 0;       n_isccpdv_c4     = 0;       sum_isccpdv_c4 = 0
   ! Initialized (All-Term Output)
   n_obs2d_all      = 0.0_dp;  n_obs3d_all      = 0.0_dp;  n_cld3d_all      = 0.0_dp
   n_c4_valid_all   = 0;       n_c4_cloud_all   = 0;       n_c4_mask_invalid_all = 0
   n_c4_type_all    = 0;       n_c4_category_all = 0;      n_c4_type_qc_all = 0
   n_c4_category_invalid_all = 0
   n_clr_all        = 0.0_dp;  n_cld_all        = 0.0_dp;  n_drz_all        = 0.0_dp;  n_prc_all = 0.0_dp
   n_lwp_all        = 0;       sum_lwp_all      = 0.0_dp
   n_iwp_all        = 0;       sum_iwp_all      = 0.0_dp
   n_cot_all        = 0;       sum_cot_all      = 0.0_dp
   n_cth_all        = 0;       sum_cth_all      = 0.0_dp
   n_cth_c4_all     = 0;       sum_cth_c4_all   = 0.0_dp
   n_ctp_c4_all     = 0;       sum_ctp_c4_all   = 0.0_dp
   n_pctgr_all      = 0;       n_ptype_all      = 0
   sum_temp_all     = 0.0_dp;  sum_shum_all     = 0.0_dp;  sum_rhum_all     = 0.0_dp
   n_cldavel_all    = 0;       sum_cldavel_all  = 0.0_dp
   n_clddvel_all    = 0;       sum_clddvel_all  = 0.0_dp
   n_cldtvel1_all   = 0;       sum_cldtvel1_all = 0.0_dp
   n_cldtvel2_all   = 0;       sum_cldtvel2_all = 0.0_dp
   n_cldtvelt_all   = 0;       sum_cldtvelt_all = 0.0_dp
   n_cwc1_all       = 0;       sum_cwc1_all     = 0.0_dp
   n_cwc2_all       = 0;       sum_cwc2_all     = 0.0_dp
   n_ciwc_all       = 0;       sum_ciwc_all     = 0.0_dp
   n_cer1_all       = 0;       sum_cer1_all     = 0.0_dp
   n_cer2_all       = 0;       sum_cer2_all     = 0.0_dp
   n_cerliq_all     = 0;       sum_cerliq_all   = 0.0_dp
   n_cerice_all     = 0;       sum_cerice_all   = 0.0_dp
   n_numc1_all      = 0;       sum_numc1_all    = 0.0_dp
   n_numc2_all      = 0;       sum_numc2_all    = 0.0_dp
   n_numliq_all     = 0;       sum_numliq_all   = 0.0_dp
   n_numice_all     = 0;       sum_numice_all   = 0.0_dp
   n_cfad_obs_all   = 0;       n_cfad_bin_all   = 0
   n_isccp_all      = 0;       sum_isccpdv_all  = 0
   n_isccp_c4_all   = 0;       n_isccpdv_c4_all = 0;       sum_isccpdv_c4_all = 0
   regmap           = 0

   ! Grid Settings
   do ilon = 1, ilon_max + 1
      glon(ilon) = (ilon-1)*gintv - 180.0_sp
   enddo
   do ilat = 1, ilat_max + 1
      glat(ilat) = (ilat-1)*gintv - 90.0_sp
   enddo
   if ( dimz .eq. 23 .or. dimz .eq. 40 ) then
      gplev(1)  = 0.0_sp
      gplev(2:) = plev(:)
   else
      write(error_unit,*) "### Vertical Levels NOT Supported.", dimz
      stop
   endif

   do idbz = 1, idbz_max + 1
      gdbz(idbz) = (idbz-1)*gintvdbz + cfad_dbzmin
   enddo
   do ilev = 1, ilev_max + 1
      glev(ilev) = (ilev-1)*gintvlev + cfad_levmin
   enddo

   ! Populate Variable Name in Original HDF5 Files
   if ( i_scan .eq. 1 ) then
      if ( version == "vCb" ) then
      varname( 1) = "/ScienceData/Geo/latitude"
      varname( 2) = "/ScienceData/Geo/longitude"
      varname( 3) = "/ScienceData/Geo/land_water_flag"
      varname( 4) = "/ScienceData/Geo/height"
      varname( 5) = "/ScienceData/Data/GRID_temperature_1km"
      varname( 6) = "/ScienceData/Data/GRID_pressure_1km"  ! Pa => hPa
      varname( 7) = "/ScienceData/Data/GRID_specificHumidity_1km"
      varname( 8) = "/ScienceData/Data/GRID_relativeHumidity_1km"
      varname( 9) = "/ScienceData/Data/optical_thickness_1km"
      varname(10) = "/ScienceData/Data/liquid_water_path_1km"
      varname(11) = "/ScienceData/Data/ice_water_path_1km"
      varname(12) = "/ScienceData/Data/cloud_radar_reflectivity_1km"
      varname(13) = "/ScienceData/Data/cloud_mask_cpr_atlid_1km"
      varname(14) = "/ScienceData/Data/cloud_air_velocity_1km"
      varname(15) = "/ScienceData/Data/cloud_doppler_velocity_1km"
      varname(16) = "/ScienceData/Data/cloud_terminal_velocity1_1km"
      varname(17) = "/ScienceData/Data/cloud_terminal_velocity2_1km"
      varname(18) = "/ScienceData/Data/total_cloud_terminal_velocity_1km"
      varname(19) = "/ScienceData/Data/cloud_water_content1_1km"
      varname(20) = "/ScienceData/Data/cloud_water_content2_1km"
      varname(21) = "/ScienceData/Data/ice_water_content_1km"
      varname(22) = "/ScienceData/Data/cloud_effective_radius1_1km"
      varname(23) = "/ScienceData/Data/cloud_effective_radius2_1km"
      varname(24) = "/ScienceData/Data/water_effective_radius_1km"
      varname(25) = "/ScienceData/Data/ice_effective_radius_1km"
      varname(26) = "/ScienceData/Data/cloud_number_concentration1_1km"
      varname(27) = "/ScienceData/Data/cloud_number_concentration2_1km"
      varname(28) = "/ScienceData/Data/total_cloud_water_number_concentration_1km"
      varname(29) = "/ScienceData/Data/total_cloud_ice_number_concentration_1km"
      varname(30) = "/ScienceData/Data/cloud_particle_category_cpr_atlid_1km"
      varname(31) = "/ScienceData/Data/cloud_particle_type_cpr_atlid_1km"
      varname(32) = "/ScienceData/Data/cloud_phase1_1km"
      varname(33) = "/ScienceData/Data/cloud_phase2_1km"
      else
      varname( 1) = "/ScienceData/Geo/latitude"
      varname( 2) = "/ScienceData/Geo/longitude"
      varname( 3) = "/ScienceData/Geo/land_water_flag"
      varname( 4) = "/ScienceData/Geo/height"
      varname( 5) = "/ScienceData/Data/GRID_temperature_1km"
      varname( 6) = "/ScienceData/Data/GRID_pressure_1km"
      varname( 7) = "/ScienceData/Data/GRID_specificHumidity_1km"
      varname( 8) = "/ScienceData/Data/GRID_relativeHumidity_1km"
      varname( 9) = "/ScienceData/Data/optical_thickness_1km"
      varname(10) = "/ScienceData/Data/cloud_water_path_1km"
      varname(11) = "/ScienceData/Data/cloud_ice_water_path_1km"
      varname(12) = "/ScienceData/Data/cloud_radar_reflectivity_1km"
      varname(13) = "/ScienceData/Data/cloud_mask_cpr_atlid_1km"
      varname(14) = "/ScienceData/Data/cloud_air_velocity_1km"
      varname(15) = "/ScienceData/Data/cloud_doppler_velocity_1km"
      varname(16) = "/ScienceData/Data/cloud_terminal_velocity1_1km"
      varname(17) = "/ScienceData/Data/cloud_terminal_velocity2_1km"
      varname(18) = "/ScienceData/Data/total_cloud_terminal_velocity_1km"
      varname(19) = "/ScienceData/Data/cloud_water_content1_1km"
      varname(20) = "/ScienceData/Data/cloud_water_content2_1km"
      varname(21) = "/ScienceData/Data/cloud_ice_content_1km"
      varname(22) = "/ScienceData/Data/cloud_effective_radius1_1km"
      varname(23) = "/ScienceData/Data/cloud_effective_radius2_1km"
      varname(24) = "/ScienceData/Data/cloud_water_effective_radius_1km"
      varname(25) = "/ScienceData/Data/cloud_ice_effective_radius_1km"
      varname(26) = "/ScienceData/Data/cloud_number_concentration1_1km"
      varname(27) = "/ScienceData/Data/cloud_number_concentration2_1km"
      varname(28) = "/ScienceData/Data/total_cloud_water_number_concentration_1km"
      varname(29) = "/ScienceData/Data/total_cloud_ice_number_concentration_1km"
      varname(30) = "/ScienceData/Data/cloud_particle_category_cpr_atlid_1km"
      varname(31) = "/ScienceData/Data/cloud_particle_type_cpr_atlid_1km"
      varname(32) = "/ScienceData/Data/cloud_phase1_1km"
      varname(33) = "/ScienceData/Data/cloud_phase2_1km"
      endif
   elseif ( i_scan .eq. 2 ) then
      if ( version == "vCb" ) then
      varname( 1) = "/ScienceData/Geo/latitude"
      varname( 2) = "/ScienceData/Geo/longitude"
      varname( 3) = "/ScienceData/Geo/land_water_flag"
      varname( 4) = "/ScienceData/Geo/height"
      varname( 5) = "/ScienceData/Data/GRID_temperature_1km"
      varname( 6) = "/ScienceData/Data/GRID_pressure_1km"  ! Pa => hPa
      varname( 7) = "/ScienceData/Data/GRID_specificHumidity_1km"
      varname( 8) = "/ScienceData/Data/GRID_relativeHumidity_1km"
      varname( 9) = "/ScienceData/Data/optical_thickness_10km"
      varname(10) = "/ScienceData/Data/liquid_water_path_10km"
      varname(11) = "/ScienceData/Data/ice_water_path_10km"
      varname(12) = "/ScienceData/Data/cloud_radar_reflectivity_10km"
      varname(13) = "/ScienceData/Data/cloud_mask_cpr_atlid_10km"
      varname(14) = "/ScienceData/Data/cloud_air_velocity_10km"
      varname(15) = "/ScienceData/Data/cloud_doppler_velocity_10km"
      varname(16) = "/ScienceData/Data/cloud_terminal_velocity1_10km"
      varname(17) = "/ScienceData/Data/cloud_terminal_velocity2_10km"
      varname(18) = "/ScienceData/Data/total_cloud_terminal_velocity_10km"
      varname(19) = "/ScienceData/Data/cloud_water_content1_10km"
      varname(20) = "/ScienceData/Data/cloud_water_content2_10km"
      varname(21) = "/ScienceData/Data/ice_water_content_10km"
      varname(22) = "/ScienceData/Data/cloud_effective_radius1_10km"
      varname(23) = "/ScienceData/Data/cloud_effective_radius2_10km"
      varname(24) = "/ScienceData/Data/water_effective_radius_10km"
      varname(25) = "/ScienceData/Data/ice_effective_radius_10km"
      varname(26) = "/ScienceData/Data/cloud_number_concentration1_10km"
      varname(27) = "/ScienceData/Data/cloud_number_concentration2_10km"
      varname(28) = "/ScienceData/Data/total_cloud_water_number_concentration_10km"
      varname(29) = "/ScienceData/Data/total_cloud_ice_number_concentration_10km"
      varname(30) = "/ScienceData/Data/cloud_particle_category_cpr_atlid_10km"
      varname(31) = "/ScienceData/Data/cloud_particle_type_cpr_atlid_10km"
      varname(32) = "/ScienceData/Data/cloud_phase1_10km"
      varname(33) = "/ScienceData/Data/cloud_phase2_10km"
      else
      varname( 1) = "/ScienceData/Geo/latitude"
      varname( 2) = "/ScienceData/Geo/longitude"
      varname( 3) = "/ScienceData/Geo/land_water_flag"
      varname( 4) = "/ScienceData/Geo/height"
      varname( 5) = "/ScienceData/Data/GRID_temperature_1km"
      varname( 6) = "/ScienceData/Data/GRID_pressure_1km"
      varname( 7) = "/ScienceData/Data/GRID_specificHumidity_1km"
      varname( 8) = "/ScienceData/Data/GRID_relativeHumidity_1km"
      varname( 9) = "/ScienceData/Data/optical_thickness_10km"
      varname(10) = "/ScienceData/Data/cloud_water_path_10km"
      varname(11) = "/ScienceData/Data/cloud_ice_water_path_10km"
      varname(12) = "/ScienceData/Data/cloud_radar_reflectivity_10km"
      varname(13) = "/ScienceData/Data/cloud_mask_cpr_atlid_10km"
      varname(14) = "/ScienceData/Data/cloud_air_velocity_10km"
      varname(15) = "/ScienceData/Data/cloud_doppler_velocity_10km"
      varname(16) = "/ScienceData/Data/cloud_terminal_velocity1_10km"
      varname(17) = "/ScienceData/Data/cloud_terminal_velocity2_10km"
      varname(18) = "/ScienceData/Data/total_cloud_terminal_velocity_10km"
      varname(19) = "/ScienceData/Data/cloud_water_content1_10km"
      varname(20) = "/ScienceData/Data/cloud_water_content2_10km"
      varname(21) = "/ScienceData/Data/cloud_ice_content_10km"
      varname(22) = "/ScienceData/Data/cloud_effective_radius1_10km"
      varname(23) = "/ScienceData/Data/cloud_effective_radius2_10km"
      varname(24) = "/ScienceData/Data/cloud_water_effective_radius_10km"
      varname(25) = "/ScienceData/Data/cloud_ice_effective_radius_10km"
      varname(26) = "/ScienceData/Data/cloud_number_concentration1_10km"
      varname(27) = "/ScienceData/Data/cloud_number_concentration2_10km"
      varname(28) = "/ScienceData/Data/total_cloud_water_number_concentration_10km"
      varname(29) = "/ScienceData/Data/total_cloud_ice_number_concentration_10km"
      varname(30) = "/ScienceData/Data/cloud_particle_category_cpr_atlid_10km"
      varname(31) = "/ScienceData/Data/cloud_particle_type_cpr_atlid_10km"
      varname(32) = "/ScienceData/Data/cloud_phase1_10km"
      varname(33) = "/ScienceData/Data/cloud_phase2_10km"
      endif
   endif

   !--- Main Loop:
   istat   = 0
   ofirst  = .true.
   ooutyr  = .false.
   ooutmon = .false.
   ofinal  = .false.
   processing_start_month = ""
   annual_start_month = ""

   do
      if ( ofinal ) then
         close(input_file_idx)
      else
         ! Analyze Input Filename
         if ( verbose_logging ) write(*,*) "=== Call analyze_file_name ==="
         call analyze_file_name( &
            cgranule, year, month, cmon, ym, imon, &
            obsstr, obsnext, avgstr_y, avgstr, avgend, &
            fname, &
            ooutyr, ooutmon, &
            istat, ofirst, ofinal)

         ! Read HDF5 Data
         if ( verbose_logging ) write(*,*) "=== Call read_HDF5_files ==="
         call read_HDF5_files( &
            varname, fname, &
            nbin, nray, &
            lat, lon, lsflg, height, temp, pres, shum, rhum, &
            cot, cwp, ciwp, dbz, cmask, &
            cldavel, clddvel, cldtvel1, cldtvel2, cldtvelt, &
            cwc1, cwc2, ciwc, cer1, cer2, cerliq, cerice, numc1, numc2, numliq, numice, &
            ccat_cpratl, ctyp_cpratl, cphase1, cphase2 )

         ! Accumulate Statistics
         if ( verbose_logging ) write(*,*) "=== Call output_level2 ==="
         call output_level2( &
            year, month, cgranule, imon, &
            nbin, nray, &
            lat, lon, lsflg, height, temp, pres, shum, rhum, &
            cot, cwp, ciwp, dbz, cmask, &
            cldavel, clddvel, cldtvel1, cldtvel2, cldtvelt, &
            cwc1, cwc2, ciwc, cer1, cer2, cerliq, cerice, numc1, numc2, numliq, numice, &
            ccat_cpratl, ctyp_cpratl, cphase1, cphase2, &
            regmap )

         if ( len_trim(processing_start_month) .eq. 0 ) then
            processing_start_month = year // "-" // month
         endif
         if ( len_trim(annual_start_month) .eq. 0 ) then
            annual_start_month = year // "-" // month
         endif

         ! Release the Array
         if ( verbose_logging ) write(*,*) "=== Release the Array ==="
         deallocate(lat,   stat=istat); call report_deallocation('dealloc lat  :', istat)
         deallocate(lon,   stat=istat); call report_deallocation('dealloc lon  :', istat)
         deallocate(lsflg, stat=istat); call report_deallocation('dealloc lsflg:', istat)
         deallocate(height,stat=istat); call report_deallocation('dealloc height:', istat)
         deallocate(temp,  stat=istat); call report_deallocation('dealloc temp :', istat)
         deallocate(pres,  stat=istat); call report_deallocation('dealloc pres :', istat)
         deallocate(shum,  stat=istat); call report_deallocation('dealloc shum :', istat)
         deallocate(rhum,  stat=istat); call report_deallocation('dealloc rhum :', istat)
         deallocate(cot,   stat=istat); call report_deallocation('dealloc cot  :', istat)
         deallocate(cwp,   stat=istat); call report_deallocation('dealloc cwp  :', istat)
         deallocate(ciwp,  stat=istat); call report_deallocation('dealloc ciwp :', istat)
         deallocate(dbz,   stat=istat); call report_deallocation('dealloc dbz  :', istat)
         deallocate(cmask, stat=istat); call report_deallocation('dealloc cmask:', istat)
         deallocate(cldavel, stat=istat); call report_deallocation('dealloc cldavel:', istat)
         deallocate(clddvel, stat=istat); call report_deallocation('dealloc clddvel:', istat)
         deallocate(cldtvel1,stat=istat); call report_deallocation('dealloc cldtvel1:', istat)
         deallocate(cldtvel2,stat=istat); call report_deallocation('dealloc cldtvel2:', istat)
         deallocate(cldtvelt,stat=istat); call report_deallocation('dealloc cldtvelt:', istat)
         deallocate(cwc1,   stat=istat); call report_deallocation('dealloc cwc1 :', istat)
         deallocate(cwc2,   stat=istat); call report_deallocation('dealloc cwc2 :', istat)
         deallocate(ciwc,   stat=istat); call report_deallocation('dealloc ciwc :', istat)
         deallocate(cer1,   stat=istat); call report_deallocation('dealloc cer1 :', istat)
         deallocate(cer2,   stat=istat); call report_deallocation('dealloc cer2 :', istat)
         deallocate(cerliq, stat=istat); call report_deallocation('dealloc cerliq:', istat)
         deallocate(cerice, stat=istat); call report_deallocation('dealloc cerice:', istat)
         deallocate(numc1,  stat=istat); call report_deallocation('dealloc numc1:', istat)
         deallocate(numc2,  stat=istat); call report_deallocation('dealloc numc2:', istat)
         deallocate(numliq, stat=istat); call report_deallocation('dealloc numliq:', istat)
         deallocate(numice, stat=istat); call report_deallocation('dealloc numice:', istat)
         deallocate(ccat_cpratl, stat=istat); call report_deallocation('dealloc ccat :', istat)
         deallocate(ctyp_cpratl, stat=istat); call report_deallocation('dealloc ctyp :', istat)
         deallocate(cphase1, stat=istat); call report_deallocation('dealloc cphase1:', istat)
         deallocate(cphase2, stat=istat); call report_deallocation('dealloc cphase2:', istat)

         files_processed = files_processed + 1
         if ( .not. verbose_logging .and. &
              ( mod(files_processed,100) .eq. 0 .or. ofinal ) ) then
            call system_clock( count=tick_now )
            elapsed_seconds = real(tick_now-tick_start,dp) / real(tick_rate,dp)
            write(*,'(A,I0,A,F10.2,A)') "Progress: ", files_processed, &
                                       " files, elapsed ", elapsed_seconds, " s"
         endif

         ! Check Output
         if ( ooutyr .or. ooutmon ) then
            continue
         else
            cycle
         endif

      endif ! if ( ofinal )

      !--- Output Level3 Statistics
      if ( verbose_logging ) write(*,*) "=== Call output_level3 ==="
      call output_level3( &
         ooutyr, ooutmon, ofinal, ym, year, imon, &
         annual_start_month, processing_start_month, &
         glon, glat, gplev, gdbz, glev )

      if ( .not. ofinal ) then
         if ( ooutyr ) then
            ooutmon = .false.  ! reset flag
            ooutyr  = .false.  ! reset flag
            avgstr_y = trim( obsnext )
            avgstr   = trim( obsnext )
            ! reset for next year
            n_obs2d      = 0.0_dp;  n_obs3d      = 0.0_dp;  n_cld3d      = 0.0_dp
            n_c4_valid   = 0;       n_c4_cloud   = 0;       n_c4_mask_invalid = 0
            n_c4_type    = 0;       n_c4_category = 0;      n_c4_type_qc = 0
            n_c4_category_invalid = 0
            n_clr        = 0;       n_cld        = 0;       n_drz        = 0;       n_prc = 0
            n_lwp        = 0;       sum_lwp      = 0.0_dp
            n_iwp        = 0;       sum_iwp      = 0.0_dp
            n_cot        = 0;       sum_cot      = 0.0_dp
            n_cth        = 0;       sum_cth      = 0.0_dp
            n_cth_c4     = 0;       sum_cth_c4   = 0.0_dp
            n_ctp_c4     = 0;       sum_ctp_c4   = 0.0_dp
            n_pctgr      = 0;       n_ptype      = 0
            sum_temp     = 0.0_dp;  sum_shum     = 0.0_dp;  sum_rhum     = 0.0_dp
            n_cldavel    = 0;       sum_cldavel  = 0.0_dp
            n_clddvel    = 0;       sum_clddvel  = 0.0_dp
            n_cldtvel1   = 0;       sum_cldtvel1 = 0.0_dp
            n_cldtvel2   = 0;       sum_cldtvel2 = 0.0_dp
            n_cldtvelt   = 0;       sum_cldtvelt = 0.0_dp
            n_cwc1       = 0;       sum_cwc1     = 0.0_dp
            n_cwc2       = 0;       sum_cwc2     = 0.0_dp
            n_ciwc       = 0;       sum_ciwc     = 0.0_dp
            n_cer1       = 0;       sum_cer1     = 0.0_dp
            n_cer2       = 0;       sum_cer2     = 0.0_dp
            n_cerliq     = 0;       sum_cerliq   = 0.0_dp
            n_cerice     = 0;       sum_cerice   = 0.0_dp
            n_numc1      = 0;       sum_numc1    = 0.0_dp
            n_numc2      = 0;       sum_numc2    = 0.0_dp
            n_numliq     = 0;       sum_numliq   = 0.0_dp
            n_numice     = 0;       sum_numice   = 0.0_dp
            n_cfad_obs   = 0;       n_cfad_bin   = 0
            n_isccp      = 0;       sum_isccpdv  = 0
            n_isccp_c4   = 0;       n_isccpdv_c4 = 0;       sum_isccpdv_c4 = 0
            annual_start_month = ""
            cycle ! back to reading HDF file
         else
            ooutmon = .false.
            avgstr = trim( obsnext )
            cycle ! back to reading HDF file
         endif
      else
         close(input_file_idx)
         exit
      endif

   end do ! back to reading HDF file

   if ( .not. verbose_logging ) then
      call system_clock( count=tick_now )
      elapsed_seconds = real(tick_now-tick_start,dp) / real(tick_rate,dp)
      write(*,'(A,I0,A,F10.2,A)') "Completed: ", files_processed, &
                                  " files, elapsed ", elapsed_seconds, " s"
   endif

   stop

contains

   subroutine report_deallocation(label, status)
      character(len=*), intent(in) :: label
      integer, intent(in) :: status

      if ( verbose_logging ) write(*,*) label, status
      if ( status .ne. 0 ) write(error_unit,*) "### ", trim(label), " failed:", status
   end subroutine report_deallocation

end program ecare_main
