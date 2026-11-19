module parameters

  integer, parameter :: sp = kind(0.0)                      ! single precision
  integer, parameter :: dp = kind(0.0d0)                    ! double precision

  integer, parameter :: ndset = 33                          ! # of read dataset
  !integer, parameter :: i_scan = 1                          ! 1 =  1 km scan
  integer, parameter :: i_scan = 2                          ! 2 = 10 km scan
  integer, parameter :: nvar2d =  9                         ! total output 2D vars
  integer, parameter :: nvar3d = 20                         ! total output 3D vars
  integer, parameter :: nvar3dcld = 17                      ! total output 3D-CF vars

  integer, parameter :: nmon = 12
  integer, parameter :: nctgr = 20                          ! dimsize of particle category
  !integer, parameter :: ntype = 20                          ! dimsize of particle type
  integer, parameter :: ntype = 14                          ! dimsize of particle type
  integer, parameter :: ntype_c4 = 14                       ! C4 particle type 1--14
  integer, parameter :: ncategory_c4 = 12                   ! C4 particle category 1--12
  integer, parameter :: ntype_qc_c4 = 4                     ! Type 15, 16, 17, other unsupported

  integer, parameter :: input_file_idx = 11

  ! Runtime logging mode.  Normal mode is the default; set
  ! ECARE_VERBOSE=1 to retain the legacy diagnostic output.
  logical :: verbose_logging
  data verbose_logging / .false. /

  !--- I/O option for users
  ! output to ASCII files?
  logical :: outascii
  !data outascii / .true. /
  data outascii / .false. /

  ! read/write 3dim/4dim variables?
  logical :: odump3d
  data odump3d / .true. /
  !data odump3d / .false. /

  ! output directory
  character(len = 128) :: outpath_lev2, outpath_lev3
  data                    outpath_lev2  / "../output/Level2" /
  data                    outpath_lev3  / "../output/Level3" /

  !--- User Defined Configs ---
  real, parameter    :: dbzmin = -35.0_sp                          ! minimum detection dBZ sensitivity
  integer, parameter :: cmaskmin = 20                              ! minimum threshold of cloud mask
  real, parameter    :: reso = 2.50_sp                             ! horizontal resolution
  !real, parameter    :: reso = 2.00_sp                             ! horizontalresolution [degr]
  !real, parameter    :: reso = 1.50_sp                             ! horizontalresolution [degr]
  real, parameter    :: cfad_ddbz = 1.0_sp                         ! resolution for CFAD dBZ
  !real, parameter    :: cfad_dlev = 0.25_sp                        ! resolution for CFAD Height (Unit = km)
  real, parameter    :: cfad_dlev = 0.5_sp                        ! resolution for CFAD Height (Unit = km)
  real, parameter    :: cfad_dbzmin = -50.0_sp                      ! minimum value of CFAD dbZ bin
  real, parameter    :: cfad_dbzmax = 20.0_sp                      ! maximum value of CFAD dbZ bin
  real, parameter    :: cfad_levmin =  0.0_sp                      ! minimum value of CFAD Lev bin (Unit = km)
  real, parameter    :: cfad_levmax = 20.0_sp                      ! maximum value of CFAD Lev bin (Unit = km)
  real, parameter    :: gintv  = reso/2.0_sp                       ! grid interval for lon/lat
  real, parameter    :: gintvdbz = cfad_ddbz/2.0_sp                ! grid interval for radar (Unit = dBZ)
  real, parameter    :: gintvlev = cfad_dlev/2.0_sp                ! grid interval for altitude (Unit = km)
  integer, parameter :: ilon_max = int( 360/(reso/2.0_sp) )        ! total # of grid for longitude
  integer, parameter :: ilat_max = int( 180/(reso/2.0_sp) )        ! total # of grid for latitude
  integer, parameter :: idbz_max = int( ( cfad_dbzmax - cfad_dbzmin ) &
                                 & / ( cfad_ddbz / 2.0_sp ) )      ! total # of grid for CFAD dBZ
  integer, parameter :: ilev_max = int( ( cfad_levmax - cfad_levmin ) &
                                 & / ( cfad_dlev / 2.0_sp ) )      ! total # of grid for CFAD Lev
  integer, parameter :: dimx = int( 360 / reso )                   ! dimension size for longitude
  integer, parameter :: dimy = int( 180 / reso )                   ! dimension size for latitude
  !integer, parameter :: dimz = 23                                  ! dimension size for P-level
  integer, parameter :: dimz = 40                                  ! dimension size for P-level
                                                                   ! dimz=23 or dimz=40
  !real(sp), parameter :: plev(dimz) = (/ &
  !    70._sp, 100._sp, 125._sp, 150._sp, 200._sp, 250._sp, 300._sp, 350._sp, &
  !   400._sp, 450._sp, 500._sp, 550._sp, 600._sp, 650._sp, 700._sp, 750._sp, &
  !   775._sp, 800._sp, 850._sp, 900._sp, 925._sp, 950._sp, 1000._sp /)
  real(sp), parameter :: plev(dimz) = (/ &
     10._sp,  70._sp, 100._sp, 125._sp, 150._sp, 175._sp, 200._sp, 250._sp, &
    275._sp, 300._sp, 325._sp, 350._sp, 400._sp, 450._sp, 500._sp, 525._sp, &
    550._sp, 575._sp, 600._sp, 625._sp, 650._sp, 675._sp, 700._sp, 725._sp, &
    750._sp, 775._sp, 800._sp, 825._sp, 850._sp, 875._sp, 900._sp, 910._sp, &
    925._sp, 940._sp, 950._sp, 960._sp, 970._sp, 980._sp, 990._sp, 1000._sp /)
  integer, parameter :: dimdbz = idbz_max / 2                      ! dimension size for CFAD dBZ
  integer, parameter :: dimlev = ilev_max / 2                      ! dimension size for CFAD Lev
  integer, parameter :: dimcot = 6                                 ! dimension size for ISCCP COT
  integer, parameter :: dimctp = 7                                 ! dimension size for ISCCP CTP
  real(sp), parameter :: gcot(dimcot+1) = (/ &
     0._sp, 1.3_sp, 3.6_sp, 9.4_sp, 23._sp, 60._sp, 1000._sp /)
  real(sp), parameter :: gctp(dimctp+1) = (/ &
    10._sp, 180._sp, 310._sp, 440._sp, 560._sp, 680._sp, 800._sp, 1100._sp /)
  !integer, parameter :: nreg = 3                                  ! 1=Tropics; 2=Midlatitudes; 3=Global
  integer, parameter :: nreg = 16                                  ! Based on Michibata et al. (2019, GMD)
                                                                   !-  1. Tropical Warm Pool (5S-15N,70E-135E)
                                                                   !-  2. ITCZ (5N-15N,140E-140W)
                                                                   !-  3. SPCZ (15S-5S,150E-130W)
                                                                   !-  4. North East Pacific (25N-50N,160W-135W)
                                                                   !-  5. California Stratocumulus (15N-35N,130W-110W)
                                                                   !-  6. Peruvian (10S-30S,90W-70W)
                                                                   !-  7. North Atlantic (30N-60N,45W-10W)
                                                                   !-  8. Namibian (30S-0S,25W-15E)
                                                                   !-  9. Australian (15S-35S,80E-115E)
                                                                   !- 10. Japan (25N-50N,125E-150E)
                                                                   !- 11. Eqt. cold tongue (5S-5N,130W-85W)
                                                                   !- 12. Eastern Asia (20N-40N,105E-120E)
                                                                   !- 13. Southern Ocean (40S-60S)
                                                                   !- 14. Arctic (66N-90N)
                                                                   !- 15. Antarctic (66S-90S)
                                                                   !- 16. Global
  real, parameter    :: rair = 287.04         ! gas constant of air: 287.04 [J/kg/K]
  real, parameter    :: undef = -999.0000_sp  ! missing undef value
  integer(kind=4), parameter :: undef_i4 = -999_4  ! missing undef value

end module parameters
