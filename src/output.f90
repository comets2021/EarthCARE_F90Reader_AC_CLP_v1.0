module output
   use, intrinsic :: iso_fortran_env, only: error_unit
   use netcdf
   use parameters
   use numbers
   use pdfs

   private get_regid, calc_pdf, calc_ratio_1d, calc_ratio_2d, increment_c4_count

contains

subroutine output_level2( &
   year, month, cgranule, imon, &
   nbin, nray, &
   lat, lon, lsflg, height, temp, pres, shum, rhum, &
   cot, cwp, ciwp, dbz, cmask, &
   cldavel, clddvel, cldtvel1, cldtvel2, cldtvelt, &
   cwc1, cwc2, ciwc, cer1, cer2, cerliq, cerice, numc1, numc2, numliq, numice, &
   ccat_cpratl, ctyp_cpratl, cphase1, cphase2, &
   regmap )

   use, intrinsic :: ieee_arithmetic
   implicit none

   character(*), intent(in) :: year, month, cgranule
   integer, intent(in) :: nbin, nray, imon
   real(dp), dimension(:), allocatable, intent(inout) :: lat
   real(dp), dimension(:), allocatable, intent(inout) :: lon
   integer(kind=4), dimension(:), allocatable, intent(in) :: lsflg
   real(sp), dimension(:,:), allocatable, intent(in) :: height ! Note: Unit is "km"
   real(sp), dimension(:,:), allocatable, intent(in) :: temp
   real(sp), dimension(:,:), allocatable, intent(in) :: pres   ! Note: Unit is "hPa"
   real(sp), dimension(:,:), allocatable, intent(in) :: shum
   real(sp), dimension(:,:), allocatable, intent(in) :: rhum
   real(sp), dimension(:),   allocatable, intent(in) :: cot
   real(sp), dimension(:),   allocatable, intent(in) :: cwp
   real(sp), dimension(:),   allocatable, intent(in) :: ciwp
   real(sp), dimension(:,:), allocatable, intent(in) :: dbz
!   integer(kind=1), dimension(:,:), allocatable, intent(in) :: cmask
   integer(kind=4), dimension(:,:), allocatable, intent(in) :: cmask
   real(sp), dimension(:,:), allocatable, intent(in) :: cldavel
   real(sp), dimension(:,:), allocatable, intent(in) :: clddvel
   real(sp), dimension(:,:), allocatable, intent(in) :: cldtvel1
   real(sp), dimension(:,:), allocatable, intent(in) :: cldtvel2
   real(sp), dimension(:,:), allocatable, intent(in) :: cldtvelt
   real(sp), dimension(:,:), allocatable, intent(in) :: cwc1
   real(sp), dimension(:,:), allocatable, intent(in) :: cwc2
   real(sp), dimension(:,:), allocatable, intent(in) :: ciwc
   real(sp), dimension(:,:), allocatable, intent(in) :: cer1
   real(sp), dimension(:,:), allocatable, intent(in) :: cer2
   real(sp), dimension(:,:), allocatable, intent(in) :: cerliq
   real(sp), dimension(:,:), allocatable, intent(in) :: cerice
   real(sp), dimension(:,:), allocatable, intent(in) :: numc1   ! Note: Unit is "# cm-3"
   real(sp), dimension(:,:), allocatable, intent(in) :: numc2   ! Note: Unit is "# cm-3"
   real(sp), dimension(:,:), allocatable, intent(in) :: numliq  ! Note: Unit is "# cm-3"
   real(sp), dimension(:,:), allocatable, intent(in) :: numice  ! Note: Unit is "# cm-3"
   integer(kind=4), dimension(:,:), allocatable, intent(in) :: ccat_cpratl
!   integer(kind=1), dimension(:,:), allocatable, intent(in) :: ctyp_cpratl
   integer(kind=4), dimension(:,:), allocatable, intent(in) :: ctyp_cpratl
!   integer(kind=1), dimension(:,:), allocatable, intent(in) :: cphase1
   integer(kind=4), dimension(:,:), allocatable, intent(in) :: cphase1
!   integer(kind=1), dimension(:,:), allocatable, intent(in) :: cphase2   
   integer(kind=4), dimension(:,:), allocatable, intent(in) :: cphase2   
   real(sp), dimension(:,:), intent(inout) :: regmap

   !--- Local
   integer :: i, k                      ! i=lon/lat loop; k=level loop
   integer :: ireg                      ! dimension number for region for CFAD
   integer :: itype                     ! dimension number for cloud praticle type
   integer :: ix, iy, iz                ! dimension number for lon/lat/alt
   integer :: idbz, ilev                ! dimension number for dbz/lev for CFAD
   integer :: icot, ictp                ! dimension number for cot/ctp for ISCCP Histogram
   integer :: icot_c4, ictp_c4          ! dimension number for C4 cot/ctp for ISCCP Histogram
   integer :: n_dum                     ! dummy
   integer :: n_excluded_type15         ! excluded non-cloud echo bins (Type 15)
   integer :: n_excluded_type16         ! excluded fully attenuated bins (Type 16)
   integer :: n_excluded_type17         ! excluded non-cloud echo bins (Type 17)
   integer :: n_excluded_unknown        ! excluded unsupported particle-type bins
   integer :: n_invalid_geolocation     ! excluded profiles with invalid geolocation
   real    :: rho                       ! air density [kg m-3]
   real    :: dbzmax                    ! maximum dBZ within all C4 cloud layers
   real    :: cth                       ! cloud top height [m]
   real    :: ctp                       ! cloud top pressure [hPa]
   real    :: cth_c4                    ! C4 cloud top height [m]
   real    :: ctp_c4                    ! C4 cloud top pressure [hPa]
   real    :: mean_dum                  ! dummy
   real    :: mean_dum_c4               ! profile-mean Doppler velocity for C4 ISCCP
   real(dp), parameter :: geolocation_fill = 9.96921e36_dp
   character(len = 128) :: outfname
   logical :: octop
   logical :: bad_granule
   logical :: has_cloud_layer
   logical :: has_valid_cloud_mask
   logical :: has_valid_incloud_dbz

   bad_granule = .false.
   n_excluded_type15 = 0
   n_excluded_type16 = 0
   n_excluded_type17 = 0
   n_excluded_unknown = 0
   n_invalid_geolocation = 0

   !--- Inline Generators (Cumulating Level2)
   profile_loop: do i = 1, nray
      octop = .false.
      cth = undef
      ctp = undef
      cth_c4 = undef
      ctp_c4 = undef
      dbzmax = undef
      has_cloud_layer = .false.
      has_valid_cloud_mask = .false.
      has_valid_incloud_dbz = .false.
      icot = 0; ictp = 0
      icot_c4 = 0; ictp_c4 = 0
      n_dum = 0
      mean_dum = 0.0_sp
      mean_dum_c4 = undef

      ! Search Gridbox
      if ( .not. ieee_is_finite(lon(i)) .or. &
         & .not. ieee_is_finite(lat(i)) ) then
         n_invalid_geolocation = n_invalid_geolocation + 1
         bad_granule = .true.
         cycle profile_loop
      endif
      if ( lon(i) .eq. geolocation_fill .or. lat(i) .eq. geolocation_fill ) then
         n_invalid_geolocation = n_invalid_geolocation + 1
         bad_granule = .true.
         cycle profile_loop
      endif
      if ( lon(i) < -180.0_sp .or. lon(i) > 180.0_sp &
      .or. lat(i) <  -90.0_sp .or. lat(i) >  90.0_sp ) then
         n_invalid_geolocation = n_invalid_geolocation + 1
         bad_granule = .true.
         cycle profile_loop
      endif
      if ( lon(i) .eq. 180.00000_sp ) lon(i) = lon(i) - 1e-2_sp
      if ( lat(i) .eq.  90.00000_sp ) lat(i) = lat(i) - 1e-2_sp
      ix = int( ( lon(i) + 180.0_sp ) / reso ) + 1
      iy = int( ( lat(i) +  90.0_sp ) / reso ) + 1
      if ( ix .gt. dimx  .or.  iy .gt. dimy  .or.  ix .le. 0 .or.  iy .le. 0 ) then
         n_invalid_geolocation = n_invalid_geolocation + 1
         bad_granule = .true.
         cycle profile_loop
      endif

      ireg = 0
      if ( nreg .eq. 3 ) then
         ! Define Region ID (Simple)
         if ( lat(i) .ge. -35.0_sp  .and.  lat(i) .le.  35.0_sp ) ireg = 1  ! Tropics
         if ( lat(i) .gt.  35.0_sp  .and.  lat(i) .le.  90.0_sp ) ireg = 2  ! Mid-Latitudes
         if ( lat(i) .lt. -35.0_sp  .and.  lat(i) .ge. -90.0_sp ) ireg = 2  ! Mid-Latitudes
         if ( ireg .eq. 0 ) then
            write(error_unit,*) "### ERROR IN REGION ID: ", ireg, lat(i)
            stop
         endif
      elseif ( nreg .eq. 16 ) then
         ! Define Region ID (Updated Michibata et al., 2019, GMD)
         call get_regid( lon(i),  &   ! [IN]    longitude
                       & lat(i),  &   ! [IN]    latitude
                       & ireg       ) ! [INOUT] region ID
         regmap( ix,iy ) = real( ireg )
      else
         write(error_unit,*) "### ERROR IN REGION DEFINITION: ", nreg
         stop
      endif

      !--- Statistics for 3D Variables
      do k = 1, nbin

         ! Diagnose the cloud layer independently of CPR reflectivity.
         ! AC_CLP vCb C4 cloud mask: 0=clear, 1=cloud.
         if ( ieee_is_finite(height(k,i)) .and. &
            & ieee_is_finite(pres(k,i))   .and. &
            & pres(k,i) .gt. 0.0_sp .and. pres(k,i) .lt. 1100.0_sp ) then
            if ( pres(k,i) .gt. 0.5_sp*(plev(dimz-1)+plev(dimz)) ) then
               iz = dimz
            else
               do iz = 1, dimz-1
                  if ( pres(k,i) .lt. 0.5_sp*(plev(iz)+plev(iz+1)) ) exit
               end do
            endif
            if ( iz .le. 0 .or. iz .gt. dimz ) then
               write(error_unit,*) "### ERROR: C4 iz diagnosis ", iz, pres(k,i)
               stop 1
            endif

            ! C4 3D diagnostics use native vertical bins and the same pressure
            ! layer assignment as legacy cld3d, without a CPR dBZ threshold.
            if ( cmask(k,i) .eq. 0 .or. cmask(k,i) .eq. 1 ) then
               call increment_c4_count( n_c4_valid(ix,iy,iz,imon), "n_c4_valid" )
               if ( cmask(k,i) .eq. 1 ) then
                  call increment_c4_count( n_c4_cloud(ix,iy,iz,imon), "n_c4_cloud" )
                  itype = ctyp_cpratl(k,i)
                  if ( itype .ge. 1 .and. itype .le. ntype_c4 ) then
                     call increment_c4_count( n_c4_type(ix,iy,iz,itype,imon), "n_c4_type" )
                  endif
                  if ( ccat_cpratl(k,i) .ge. 1 .and. &
                     & ccat_cpratl(k,i) .le. ncategory_c4 ) then
                     call increment_c4_count( &
                        n_c4_category(ix,iy,iz,ccat_cpratl(k,i),imon), "n_c4_category" )
                  endif
               endif
            else
               call increment_c4_count( n_c4_mask_invalid(ix,iy,iz,imon), &
                                        "n_c4_mask_invalid" )
            endif

            ! Core QC counts use all analysis-valid bins, independent of the
            ! C4 mask.  Classes 1--3 are Types 15--17; class 4 is any other
            ! unsupported raw Type code.  Official Types 0--14 are not QC.
            itype = ctyp_cpratl(k,i)
            if ( itype .ge. 15 .and. itype .le. 17 ) then
               call increment_c4_count( n_c4_type_qc(ix,iy,iz,itype-14,imon), &
                                        "n_c4_type_qc" )
            elseif ( itype .lt. 0 .or. itype .gt. 17 ) then
               call increment_c4_count( n_c4_type_qc(ix,iy,iz,4,imon), &
                                        "n_c4_type_qc" )
            endif
            if ( ccat_cpratl(k,i) .lt. 0 .or. &
               & ccat_cpratl(k,i) .gt. ncategory_c4 ) then
               call increment_c4_count( n_c4_category_invalid(ix,iy,iz,imon), &
                                        "n_c4_category_invalid" )
            endif

            if ( cmask(k,i) .eq. 0 .or. cmask(k,i) .eq. 1 ) then
               has_valid_cloud_mask = .true.
            endif
            if ( cmask(k,i) .eq. 1 ) then
               if ( .not. has_cloud_layer .or. height(k,i)*1.d3 .gt. cth_c4 ) then
                  cth_c4 = height(k,i)*1.d3
                  ctp_c4 = pres(k,i)
               endif
               has_cloud_layer = .true.
            endif
         endif

         ! Exclude "NaN/Inf" values
         if ( .not. ieee_is_finite(height(k,i)) .or. &
              .not. ieee_is_finite(dbz(k,i)) ) cycle

         if ( pres(k,i) .gt. 0.0_sp .and. pres(k,i) .lt. 1100.0_sp ) then
            if ( pres(k,i) .gt. 0.5_sp*(plev(dimz-1)+plev(dimz)) ) then
               iz = dimz
            else
               do iz = 1, dimz-1
                  if ( pres(k,i) .lt. 0.5_sp*(plev(iz)+plev(iz+1)) ) exit
               end do
            endif
            if ( iz .le. 0 .or. iz .gt. dimz ) then
               write(error_unit,*) "### ERROR: iz diagnosis ", iz, pres(k,i)
               stop
            endif

            rho = pres(k,i)*1.d2 / ( rair * temp(k,i) )
            ! Cloud-Top Height [m]
            if ( .not. octop  .and. &
               & dbz(k,i) .ge. dbzmin ) then
               octop = .true.
               cth = height(k,i)*1.d3
               ctp = pres(k,i)
            endif
            ! Maximum valid CPR reflectivity within all C4 cloud layers.
            ! The HDF5 valid_range and the cloud-layer condition are
            ! deliberately independent of the -35 dBZ classification limit.
            if ( cmask(k,i) .eq. 1      .and. &
               & dbz(k,i) .ge. -60.0_sp .and. &
               & dbz(k,i) .le.  40.0_sp       ) then
               if ( .not. has_valid_incloud_dbz ) then
                  dbzmax = dbz(k,i)
                  has_valid_incloud_dbz = .true.
               else
                  dbzmax = max( dbzmax, dbz(k,i) )
               endif
            endif

            ! # of observation
            n_obs3d(ix,iy,iz,imon) = n_obs3d(ix,iy,iz,imon) + 1

            ! # of clouds
            itype = ctyp_cpratl(k,i)
            if ( dbz(k,i) .ge. dbzmin  .and. &
               & itype .ne. 0                ) then
               if ( itype .ge. 1 .and. itype .le. 14 ) then
                  n_cld3d(ix,iy,iz,0,imon) = n_cld3d(ix,iy,iz,0,imon) + 1   ! NOTE: cloudy-sky here.
                  n_cld3d(ix,iy,iz,itype,imon) = n_cld3d(ix,iy,iz,itype,imon) + 1
               elseif ( itype .eq. 15 ) then
                  n_excluded_type15 = n_excluded_type15 + 1
               elseif ( itype .eq. 16 ) then
                  n_excluded_type16 = n_excluded_type16 + 1
               elseif ( itype .eq. 17 ) then
                  n_excluded_type17 = n_excluded_type17 + 1
               else
                  n_excluded_unknown = n_excluded_unknown + 1
               endif
            endif

            ! Temperature [K]
            sum_temp(ix,iy,iz,imon) = sum_temp(ix,iy,iz,imon) + temp(k,i)

            ! Spcific Humidity [kg kg-1]
            sum_shum(ix,iy,iz,imon) = sum_shum(ix,iy,iz,imon) + shum(k,i)

            ! Relative Humidity [%]
            sum_rhum(ix,iy,iz,imon) = sum_rhum(ix,iy,iz,imon) + rhum(k,i)

            ! Air Vertical Velocity [m s-1]
            if ( cldavel(k,i) .gt. -20.0_sp .and. cldavel(k,i) .lt. 20.0_sp ) then
               n_cldavel(ix,iy,iz,imon) = n_cldavel(ix,iy,iz,imon) + 1
               sum_cldavel(ix,iy,iz,imon) = sum_cldavel(ix,iy,iz,imon) + cldavel(k,i)
            endif

            ! Cloud Doppler Velocity [m s-1]
            if ( clddvel(k,i) .gt. -20.0_sp .and. clddvel(k,i) .lt. 20.0_sp ) then
               n_clddvel(ix,iy,iz,imon) = n_clddvel(ix,iy,iz,imon) + 1
               sum_clddvel(ix,iy,iz,imon) = sum_clddvel(ix,iy,iz,imon) + clddvel(k,i)
               ! for ISCCP-like Histogram
               n_dum = n_dum + 1
               mean_dum = mean_dum + clddvel(k,i)
            endif

            ! Cloud Terminal Velocity 1 [m s-1]
            if ( cldtvel1(k,i) .gt. -20.0_sp .and. cldtvel1(k,i) .lt. 20.0_sp ) then
               n_cldtvel1(ix,iy,iz,imon) = n_cldtvel1(ix,iy,iz,imon) + 1
               sum_cldtvel1(ix,iy,iz,imon) = sum_cldtvel1(ix,iy,iz,imon) + cldtvel1(k,i)
            endif

            ! Cloud Terminal Velocity 2 [m s-1]
            if ( cldtvel2(k,i) .gt. -20.0_sp .and. cldtvel2(k,i) .lt. 20.0_sp ) then
               n_cldtvel2(ix,iy,iz,imon) = n_cldtvel2(ix,iy,iz,imon) + 1
               sum_cldtvel2(ix,iy,iz,imon) = sum_cldtvel2(ix,iy,iz,imon) + cldtvel2(k,i)
            endif

            ! Total Cloud Terminal Velocity (Sedimentation) [m s-1]
            if ( cldtvelt(k,i) .gt. -20.0_sp .and. cldtvelt(k,i) .lt. 0.0_sp ) then
               n_cldtvelt(ix,iy,iz,imon) = n_cldtvelt(ix,iy,iz,imon) + 1
               sum_cldtvelt(ix,iy,iz,imon) = sum_cldtvelt(ix,iy,iz,imon) + cldtvelt(k,i)
            endif

            ! Cloud Water Content 1 [g m-3]
            if ( cwc1(k,i) .gt. 0.0_sp .and. cwc1(k,i) .lt. 10.0_sp ) then
               n_cwc1(ix,iy,iz,imon) = n_cwc1(ix,iy,iz,imon) + 1
               sum_cwc1(ix,iy,iz,imon) = sum_cwc1(ix,iy,iz,imon) + cwc1(k,i)
            endif

            ! Cloud Water Content 2 [g m-3]
            if ( cwc2(k,i) .gt. 0.0_sp .and. cwc2(k,i) .lt. 10.0_sp ) then
               n_cwc2(ix,iy,iz,imon) = n_cwc2(ix,iy,iz,imon) + 1
               sum_cwc2(ix,iy,iz,imon) = sum_cwc2(ix,iy,iz,imon) + cwc2(k,i)
            endif

            ! Ice Water Content [g m-3]
            if ( ciwc(k,i) .gt. 0.0_sp .and. ciwc(k,i) .lt. 10.0_sp ) then
               n_ciwc(ix,iy,iz,imon) = n_ciwc(ix,iy,iz,imon) + 1
               sum_ciwc(ix,iy,iz,imon) = sum_ciwc(ix,iy,iz,imon) + ciwc(k,i)
            endif

            ! Cloud Effective Radius 1 [micron]
            if ( cer1(k,i) .gt. 0.0_sp .and. cer1(k,i) .lt. 10000.0_sp ) then
               n_cer1(ix,iy,iz,imon) = n_cer1(ix,iy,iz,imon) + 1
               sum_cer1(ix,iy,iz,imon) = sum_cer1(ix,iy,iz,imon) + cer1(k,i)
            endif

            ! Cloud Effective Radius 2 [micron]
            if ( cer2(k,i) .gt. 0.0_sp .and. cer2(k,i) .lt. 10000.0_sp ) then
               n_cer2(ix,iy,iz,imon) = n_cer2(ix,iy,iz,imon) + 1
               sum_cer2(ix,iy,iz,imon) = sum_cer2(ix,iy,iz,imon) + cer2(k,i)
            endif

            ! Cloud Effective Radius Liq [micron]
            if ( cerliq(k,i) .gt. 0.0_sp .and. cerliq(k,i) .lt. 10000.0_sp ) then
               n_cerliq(ix,iy,iz,imon) = n_cerliq(ix,iy,iz,imon) + 1
               sum_cerliq(ix,iy,iz,imon) = sum_cerliq(ix,iy,iz,imon) + cerliq(k,i)
            endif

            ! Cloud Effective Radius Ice [micron]
            if ( cerice(k,i) .gt. 0.0_sp .and. cerice(k,i) .lt. 10000.0_sp ) then
               n_cerice(ix,iy,iz,imon) = n_cerice(ix,iy,iz,imon) + 1
               sum_cerice(ix,iy,iz,imon) = sum_cerice(ix,iy,iz,imon) + cerice(k,i)
            endif

            ! Cloud Number Concentration 1 [# cm-3]
            if ( numc1(k,i) .gt. 0.0_sp .and. numc1(k,i) .lt. 1.e21 ) then
               n_numc1(ix,iy,iz,imon) = n_numc1(ix,iy,iz,imon) + 1
               sum_numc1(ix,iy,iz,imon) = sum_numc1(ix,iy,iz,imon) + numc1(k,i)
            endif

            ! Cloud Number Concentration 2 [# cm-3]
            if ( numc2(k,i) .gt. 0.0_sp .and. numc2(k,i) .lt. 1.e21 ) then
               n_numc2(ix,iy,iz,imon) = n_numc2(ix,iy,iz,imon) + 1
               sum_numc2(ix,iy,iz,imon) = sum_numc2(ix,iy,iz,imon) + numc2(k,i)
            endif

            ! Cloud Number Concentration Liq [# cm-3]
            if ( numliq(k,i) .gt. 0.0_sp .and. numliq(k,i) .lt. 1.e21 ) then
               n_numliq(ix,iy,iz,imon) = n_numliq(ix,iy,iz,imon) + 1
               sum_numliq(ix,iy,iz,imon) = sum_numliq(ix,iy,iz,imon) + numliq(k,i)
            endif

            ! Cloud Number Concentration Ice [# cm-3]
            if ( numice(k,i) .gt. 0.0_sp .and. numice(k,i) .lt. 1.e21 ) then
               n_numice(ix,iy,iz,imon) = n_numice(ix,iy,iz,imon) + 1
               sum_numice(ix,iy,iz,imon) = sum_numice(ix,iy,iz,imon) + numice(k,i)
            endif

         endif

         ! CFAD
         if ( height(k,i) .ge. 1.0_dp       .and. &  ! remove ground clutter
            & height(k,i) .lt. cfad_levmax  .and. &
            & dbz(k,i)    .ge. cfad_dbzmin  .and. &
            & dbz(k,i)    .lt. cfad_dbzmax        ) then
            idbz = int( ( dbz(k,i) - cfad_dbzmin ) / cfad_ddbz ) + 1
            ilev = int( ( height(k,i) - cfad_levmin ) / cfad_dlev ) + 1
            ! Regional CFAD
            if ( ireg .ne. 0 ) then
               n_cfad_obs( ireg,imon ) = n_cfad_obs( ireg,imon ) + 1
               n_cfad_bin( idbz,ilev,ireg,imon ) = n_cfad_bin( idbz,ilev,ireg,imon ) + 1
            endif
            ! Global CFAD
            n_cfad_obs( nreg,imon ) = n_cfad_obs( nreg,imon ) + 1
            n_cfad_bin( idbz,ilev,nreg,imon ) = n_cfad_bin( idbz,ilev,nreg,imon ) + 1
         endif

      enddo

      !--- Statistics for 2D Variables
      ! # of geolocated profiles (kept independent of classification validity)
      n_obs2d( ix,iy,imon ) = n_obs2d( ix,iy,imon ) + 1

      if ( .not. has_cloud_layer .and. has_valid_cloud_mask ) then
         ! Explicit clear-sky detection by the C4 cloud mask.  Profiles with
         ! only undefined mask codes are not counted as clear.
         n_clr( ix,iy,imon ) = n_clr( ix,iy,imon ) + 1
      elseif ( has_cloud_layer .and. has_valid_incloud_dbz ) then
         if ( dbzmax .ge. -35.0_sp .and. dbzmax .lt. -15.0_sp ) then
            n_cld( ix,iy,imon ) = n_cld( ix,iy,imon ) + 1
         elseif ( dbzmax .ge. -15.0_sp .and. dbzmax .lt. 0.0_sp ) then
            n_drz( ix,iy,imon ) = n_drz( ix,iy,imon ) + 1
         elseif ( dbzmax .ge. 0.0_sp ) then
            n_prc( ix,iy,imon ) = n_prc( ix,iy,imon ) + 1
         endif
      endif

      ! LWP [g m-2]
      if ( cwp(i) .gt. 0.0_sp .and. cwp(i) .lt. 20000.0_sp ) then
         n_lwp( ix,iy,imon ) = n_lwp( ix,iy,imon ) + 1
         sum_lwp( ix,iy,imon ) = sum_lwp( ix,iy,imon ) + cwp(i)
      endif

      ! IWP [g m-2]
      if ( ciwp(i) .gt. 0.0_sp .and. ciwp(i) .lt. 20000.0_sp ) then
         n_iwp( ix,iy,imon ) = n_iwp( ix,iy,imon ) + 1
         sum_iwp( ix,iy,imon ) = sum_iwp( ix,iy,imon ) + ciwp(i)
      endif

      ! COT
      if ( cot(i) .gt. 0.0_sp .and. cot(i) .lt. 1000.0_sp ) then
         n_cot( ix,iy,imon ) = n_cot( ix,iy,imon ) + 1
         sum_cot( ix,iy,imon ) = sum_cot( ix,iy,imon ) + cot(i)
      endif

      ! CTH [m]
      if ( cth .gt. 0.0_sp  .and. cth .lt. 20000.0_sp ) then
         n_cth( ix,iy,imon ) = n_cth( ix,iy,imon ) + 1
         sum_cth( ix,iy,imon ) = sum_cth( ix,iy,imon ) + cth
      endif

      ! C4 Cloud Top Height [m]
      if ( has_cloud_layer ) then
         n_cth_c4( ix,iy,imon ) = n_cth_c4( ix,iy,imon ) + 1
         sum_cth_c4( ix,iy,imon ) = sum_cth_c4( ix,iy,imon ) + cth_c4
      endif

      ! C4 Cloud Top Pressure [hPa]
      if ( has_cloud_layer ) then
         n_ctp_c4( ix,iy,imon ) = n_ctp_c4( ix,iy,imon ) + 1
         sum_ctp_c4( ix,iy,imon ) = sum_ctp_c4( ix,iy,imon ) + ctp_c4
      endif

      if ( n_dum .gt. 0 ) mean_dum_c4 = mean_dum / real(n_dum)

      !--- Statistics for ISCCP-like Histogram
      if ( cot(i) .gt. gcot(1) .and. cot(i) .lt. gcot(dimcot+1) .and. &
         & ctp    .ge. gctp(1) .and. ctp    .lt. gctp(dimctp+1) .and. &
         & n_dum  .gt. 0 ) then
         mean_dum = mean_dum / real(n_dum)  ! mean Doppler velocity
         ! Diagnose COT bin (icot)
         do icot = 1, dimcot
            if ( cot(i) .ge. gcot(icot) .and. cot(i) .lt. gcot(icot+1) ) then
               exit
            endif
         enddo
         ! Diagnose CTP bin (ictp)
         do ictp = 1, dimctp
            if ( ctp .ge. gctp(ictp) .and. ctp .lt. gctp(ictp+1) ) then
               exit
            endif
         enddo

         ! Regional ISCCP Histogram
         if ( ireg .ne. 0 ) then
            n_isccp( icot,ictp,ireg,imon ) = n_isccp( icot,ictp,ireg,imon ) + 1
            sum_isccpdv( icot,ictp,ireg,imon ) = sum_isccpdv( icot,ictp,ireg,imon ) + mean_dum
         endif
         ! Global ISCCP Histogram
         n_isccp( icot,ictp,nreg,imon ) = n_isccp( icot,ictp,nreg,imon ) + 1
         sum_isccpdv( icot,ictp,nreg,imon ) = sum_isccpdv( icot,ictp,nreg,imon ) + mean_dum
      endif

      !--- Statistics for C4 ISCCP-like Histogram
      ! Frequency does not require a valid Doppler velocity.
      if ( cot(i) .gt. gcot(1) .and. cot(i) .lt. gcot(dimcot+1) .and. &
         & ctp_c4 .ge. gctp(1) .and. ctp_c4 .lt. gctp(dimctp+1) ) then
         ! Diagnose COT bin (icot_c4)
         do icot_c4 = 1, dimcot
            if ( cot(i) .ge. gcot(icot_c4) .and. cot(i) .lt. gcot(icot_c4+1) ) then
               exit
            endif
         enddo
         ! Diagnose C4 CTP bin (ictp_c4)
         do ictp_c4 = 1, dimctp
            if ( ctp_c4 .ge. gctp(ictp_c4) .and. ctp_c4 .lt. gctp(ictp_c4+1) ) then
               exit
            endif
         enddo

         ! Regional C4 ISCCP Histogram
         if ( ireg .ne. 0 ) then
            n_isccp_c4( icot_c4,ictp_c4,ireg,imon ) = &
               n_isccp_c4( icot_c4,ictp_c4,ireg,imon ) + 1
            if ( n_dum .gt. 0 ) then
               n_isccpdv_c4( icot_c4,ictp_c4,ireg,imon ) = &
                  n_isccpdv_c4( icot_c4,ictp_c4,ireg,imon ) + 1
               sum_isccpdv_c4( icot_c4,ictp_c4,ireg,imon ) = &
                  sum_isccpdv_c4( icot_c4,ictp_c4,ireg,imon ) + mean_dum_c4
            endif
         endif
         ! Global C4 ISCCP Histogram
         n_isccp_c4( icot_c4,ictp_c4,nreg,imon ) = &
            n_isccp_c4( icot_c4,ictp_c4,nreg,imon ) + 1
         if ( n_dum .gt. 0 ) then
            n_isccpdv_c4( icot_c4,ictp_c4,nreg,imon ) = &
               n_isccpdv_c4( icot_c4,ictp_c4,nreg,imon ) + 1
            sum_isccpdv_c4( icot_c4,ictp_c4,nreg,imon ) = &
               sum_isccpdv_c4( icot_c4,ictp_c4,nreg,imon ) + mean_dum_c4
         endif
      endif

   enddo profile_loop

   if ( n_invalid_geolocation .gt. 0 ) then
      write(error_unit,'(A,A,A,A,A,A,A,I0)') &
         & "### WARNING: profiles with invalid geolocation were excluded; year/month/granule=", &
         & trim(year), "/", trim(month), "/", trim(cgranule), &
         & "; excluded_profiles=", n_invalid_geolocation
   endif

   if ( n_excluded_type15 + n_excluded_type16 + n_excluded_type17 + n_excluded_unknown .gt. 0 ) then
      write(error_unit,'(A,A,A,A,A,A,A,I0,A,I0,A,I0,A,I0)') &
         & "### WARNING: bins satisfying the existing cld3d pressure/dBZ/type gate were excluded; year/month/granule=", &
         & trim(year), "/", trim(month), "/", trim(cgranule), "; Type15(non-cloud)=", n_excluded_type15, &
         & "; Type17(non-cloud)=", n_excluded_type17, "; Type16(fully-attenuated)=", n_excluded_type16, &
         & "; unknown-code=", n_excluded_unknown
   endif

   !--- Output to Ascii File
   if ( outascii ) then
      outfname = trim( outpath_lev2 ) // "/" // trim(year) // "/" // trim(month) // "/" // "ecare_" // cgranule // ".txt"
      open( 21, file = trim( outfname ), status = "unknown", action = "write" )
      if ( verbose_logging ) write( *,* ) "*** FILENAME (OUT): ", trim( outfname )
      if ( verbose_logging ) write( *,* ) "*** Output to text file..."
      do i = 1, nray
         write( 21,* ) i, lon(i), lat(i), cwp(i), ciwp(i), cot(i), cth, dbzmax
         do k = 1, nbin
            write( 21,* ) k, height(k,i), temp(k,i), pres(k,i), shum(k,i), rhum(k,i), dbz(k,i), cmask(k,i), &
                        & cldavel(k,i), clddvel(k,i), cldtvel1(k,i), cldtvel2(k,i), cldtvelt(k,i), &
                        & cwc1(k,i), cwc2(k,i), ciwc(k,i), cer1(k,i), cer2(k,i), cerliq(k,i), cerice(k,i), &
                        & numc1(k,i), numc2(k,i), numliq(k,i), numice(k,i), &
                        & ccat_cpratl(k,i), ctyp_cpratl(k,i), cphase1(k,i), cphase2(k,i)
         enddo
      enddo
      close( 21 )
   endif

   !--- Log/Warning Message
   if ( bad_granule ) then
      outfname = "./log/" // trim(year) // trim(month) // "_" // cgranule // ".log"
      open( 22, file = trim( outfname ), status = "unknown", action = "write" )
      write( 22,* ) "### NBIN, NRAY ==> ", nbin, nray
      write( 22,* ) "### Invalid geolocation profiles excluded ==> ", n_invalid_geolocation
      write( 22,* ) "### Criteria: non-finite, HDF5 FillValue, out-of-range, or invalid grid index"
      close( 22 )
   endif

end subroutine output_level2

subroutine get_regid( r8lon, r8lat, id )
   implicit none

   real(dp), intent(in)    :: r8lon
   real(dp), intent(in)    :: r8lat
   integer, intent(inout) :: id
   real :: lon, lat

   lon = r8lon
   lat = r8lat
   if( lon < 0 ) lon = lon + 360.0_dp  !! -180 to 180 --> 0 to 360
   !--- Region definition following Michibata et al. (2019,GMD)
   ! 1. Tropical Warm Pool (5S-15N,70E-135E)
   if ( lon >=  70.0_sp .and. lon <= 135.0_sp .and. lat >=  -5.0_sp .and. lat <=  15.0_sp ) id =  1
   ! 2. ITCZ (5N-15N,140E-140W)
   if ( lon >= 140.0_sp .and. lon <= 220.0_sp .and. lat >=   5.0_sp .and. lat <=  15.0_sp ) id =  2
   ! 3. SPCZ (15S-5S,150E-130W)
   if ( lon >= 150.0_sp .and. lon <= 230.0_sp .and. lat >= -15.0_sp .and. lat <=  -5.0_sp ) id =  3
   ! 4. North East Pacific (25N-50N,160W-135W)
   if ( lon >= 200.0_sp .and. lon <= 225.0_sp .and. lat >=  25.0_sp .and. lat <=  50.0_sp ) id =  4
   ! 5. California StCu (15N-35N,130W-110W)
   if ( lon >= 230.0_sp .and. lon <= 250.0_sp .and. lat >=  15.0_sp .and. lat <=  35.0_sp ) id =  5
   ! 6. Peruvian (10S-30S,90W-70W)
   if ( lon >= 270.0_sp .and. lon <= 290.0_sp .and. lat >= -30.0_sp .and. lat <= -10.0_sp ) id =  6
   ! 7. North Atlantic (30N-60N,45W-10W)
   if ( lon >= 315.0_sp .and. lon <= 350.0_sp .and. lat >=  30.0_sp .and. lat <=  60.0_sp ) id =  7
   ! 8. Namibian (30S-0S,25W-15E)
   if ( lon >= 335.0_sp .and. lon <= 360.0_sp .and. lat >= -30.0_sp .and. lat <=   0.0_sp ) id =  8
   if ( lon >=   0.0_sp .and. lon <=  15.0_sp .and. lat >= -30.0_sp .and. lat <=   0.0_sp ) id =  8
   ! 9. Australian (15S-35S,80E-115E)
   if ( lon >=  80.0_sp .and. lon <= 115.0_sp .and. lat >= -35.0_sp .and. lat <= -15.0_sp ) id =  9
   ! 10. Japan (25N-50N,125E-150E)
   if ( lon >= 125.0_sp .and. lon <= 150.0_sp .and. lat >=  25.0_sp .and. lat <=  50.0_sp ) id = 10
   ! 11. Eqt. cold tongue (5S-5N,130W-85W)
   if ( lon >= 230.0_sp .and. lon <= 275.0_sp .and. lat >=  -5.0_sp .and. lat <=   5.0_sp ) id = 11
   ! 12. Eastern Asia (20N-40N,105E-120E)
   if ( lon >= 105.0_sp .and. lon <= 120.0_sp .and. lat >=  20.0_sp .and. lat <=  40.0_sp ) id = 12
   ! 13. Southern Ocean (40S-60S)
   if ( lon >=   0.0_sp .and. lon <= 360.0_sp .and. lat >= -60.0_sp .and. lat <= -40.0_sp ) id = 13
   ! 14. Arctic (66N-90N)
   if ( lon >=   0.0_sp .and. lon <= 360.0_sp .and. lat >= 66.0_sp .and. lat <= 90.0_sp ) id = 14
   ! 15. Antarctic (66S-90S)
   if ( lon >=   0.0_sp .and. lon <= 360.0_sp .and. lat >= -90.0_sp .and. lat <= -66.0_sp ) id = 15

   return
end subroutine get_regid

subroutine output_level3( &
   ooutyr, ooutmon, ofinal, ym, year, imon, &
   annual_start_month, processing_start_month, &
   glon, glat, gplev, gdbz, glev )

   use summary_processor, only: output_summary

   implicit none
   logical, intent(in) :: ooutyr, ooutmon, ofinal
   character(*), intent(in) :: ym, year
   character(*), intent(in) :: annual_start_month, processing_start_month
   integer, intent(in) :: imon
   real(sp), dimension(:), intent(in) :: glon                 ! grid value (longitude)
   real(sp), dimension(:), intent(in) :: glat                 ! grid value (latitude)
   real(sp), dimension(:), intent(in) :: gplev                ! grid value (P-levels)
   real(sp), dimension(:), intent(in) :: gdbz( idbz_max + 1 ) ! grid value (CFAD dBZ)
   real(sp), dimension(:), intent(in) :: glev( ilev_max + 1 ) ! grid value (CFAD Lev)

   call calc_pdf(ooutyr, ooutmon, ofinal, imon)

   call output_summary( &
      ooutyr, ooutmon, ofinal, &
      ym, year, imon, &
      annual_start_month, processing_start_month, &
      glon, glat, gplev, gdbz, glev, &
      n_obs2d, n_obs2d_all, &
      n_obs3d, n_obs3d_all, &
      n_cld3d, n_cld3d_all  )

end subroutine output_level3

subroutine calc_pdf(ooutyr, ooutmon, ofinal, imon)
   implicit none
   logical, intent(in) :: ooutyr, ooutmon, ofinal
   integer, intent(in) :: imon

   integer :: iz, im
   integer :: idbz, ilev  ! dimension number for dbz/lev for CFAD
   integer :: icot, ictp  ! dimension number for cot/ctp for ISCCP
   real(sp) :: smpl_isccp( size(n_isccp,3) )  ! expected to have dimension NREG
   real(sp) :: smpl_isccp_c4( size(n_isccp_c4,3) )
   real(sp) :: smpl_pdf( size(n_obs2d,1), size(n_obs2d,2) )

   ! Initialize Variables
   if ( imon .eq. 1 ) then
      mean_lwp = 0;      mean_iwp = 0;      mean_cot = 0;      mean_cth = 0
      mean_cth_c4 = 0;   mean_ctp_c4 = 0
      pdf_clr = 0;       pdf_cld = 0;       pdf_drz = 0;       pdf_prc = 0
      mean_temp = 0;     mean_shum = 0;     mean_rhum = 0
      mean_cldavel = 0;  mean_clddvel = 0;  mean_cldtvel1 = 0; mean_cldtvel2 = 0
      mean_cldtvelt = 0; mean_cwc1 = 0;     mean_cwc2 = 0;     mean_ciwc = 0
      mean_cer1 = 0;     mean_cer2 = 0;     mean_cerliq = 0;   mean_cerice = 0
      mean_numc1 = 0;    mean_numc2 = 0;    mean_numliq = 0;   mean_numice = 0
   endif

   !--- Monthly Statistics: 2D Variables
   ! Explicit clear/cloud classification population.  Cloudy profiles with no
   ! valid in-cloud dBZ, Zmax<-35, and undetermined mask profiles are excluded.
   smpl_pdf(:,:) = real( n_clr(:,:,imon) + n_cld(:,:,imon) &
                       & + n_drz(:,:,imon) + n_prc(:,:,imon), kind=sp )
   where ( smpl_pdf(:,:) .eq. 0.0_sp )
      pdf_clr(:,:,imon) = undef
      pdf_cld(:,:,imon) = undef
      pdf_drz(:,:,imon) = undef
      pdf_prc(:,:,imon) = undef
   else where 
      pdf_clr(:,:,imon) = 100.0_sp * real( n_clr(:,:,imon), kind=sp ) / smpl_pdf(:,:)
      pdf_cld(:,:,imon) = &
         100.0_sp * real( n_cld(:,:,imon), kind=sp ) / smpl_pdf(:,:)
      pdf_drz(:,:,imon) = &
         100.0_sp * real( n_drz(:,:,imon), kind=sp ) / smpl_pdf(:,:)
      pdf_prc(:,:,imon) = &
         100.0_sp * real( n_prc(:,:,imon), kind=sp ) / smpl_pdf(:,:)
   end where

   ! LWP
   mean_lwp(:,:,imon) = &
      calc_ratio_2d( real( sum_lwp(:,:,imon) ), &
                   & real( n_lwp(:,:,imon) ) )

   ! IWP
   mean_iwp(:,:,imon) = &
      calc_ratio_2d( real( sum_iwp(:,:,imon) ), &
                   & real( n_iwp(:,:,imon) ) )

   ! COT
   mean_cot(:,:,imon) = &
      calc_ratio_2d( real( sum_cot(:,:,imon) ), &
                   & real( n_cot(:,:,imon) ) )

   ! CTH
   mean_cth(:,:,imon) = &
      calc_ratio_2d( real( sum_cth(:,:,imon) ), &
                   & real( n_cth(:,:,imon) ) )

   ! C4 CTH
   mean_cth_c4(:,:,imon) = &
      calc_ratio_2d( real( sum_cth_c4(:,:,imon) ), &
                   & real( n_cth_c4(:,:,imon) ) )

   ! C4 CTP
   mean_ctp_c4(:,:,imon) = &
      calc_ratio_2d( real( sum_ctp_c4(:,:,imon) ), &
                   & real( n_ctp_c4(:,:,imon) ) )

   !--- Monthly Statistics: 3D Variables
   do iz = 1, dimz

      ! Temperature [K]
      mean_temp(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_temp(:,:,iz,imon) ), &
                      & real( n_obs3d(:,:,iz,imon) ) )

      ! Specific Humidity [kg kg-1]
      mean_shum(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_shum(:,:,iz,imon) ), &
                      & real( n_obs3d(:,:,iz,imon) ) )

      ! Relative Humidity [%]
      mean_rhum(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_rhum(:,:,iz,imon) ), &
                      & real( n_obs3d(:,:,iz,imon) ) )

      ! Air Vertical Velocity [m s-1]
      mean_cldavel(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_cldavel(:,:,iz,imon) ), &
                      & real( n_cldavel(:,:,iz,imon) ) )

      ! Cloud Doppler Velocity [m s-1]
      mean_clddvel(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_clddvel(:,:,iz,imon) ), &
                      & real( n_clddvel(:,:,iz,imon) ) )

      ! Cloud Terminal Velocity 1 [m s-1]
      mean_cldtvel1(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_cldtvel1(:,:,iz,imon) ), &
                      & real( n_cldtvel1(:,:,iz,imon) ) )

      ! Cloud Terminal Velocity 2 [m s-1]
      mean_cldtvel2(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_cldtvel2(:,:,iz,imon) ), &
                      & real( n_cldtvel2(:,:,iz,imon) ) )

      ! Total Cloud Terminal Velocity [m s-1]
      mean_cldtvelt(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_cldtvelt(:,:,iz,imon) ), &
                      & real( n_cldtvelt(:,:,iz,imon) ) )

      ! Cloud Water Content 1 [g m-3]
      mean_cwc1(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_cwc1(:,:,iz,imon) ), &
                      & real( n_cwc1(:,:,iz,imon) ) )

      ! Cloud Water Content 2 [g m-3]
      mean_cwc2(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_cwc2(:,:,iz,imon) ), &
                      & real( n_cwc2(:,:,iz,imon) ) )

      ! Ice Water Content [g m-3]
      mean_ciwc(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_ciwc(:,:,iz,imon) ), &
                      & real( n_ciwc(:,:,iz,imon) ) )

      ! Cloud Effective Radius 1 [micron]
      mean_cer1(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_cer1(:,:,iz,imon) ), &
                      & real( n_cer1(:,:,iz,imon) ) )

      ! Cloud Effective Radius 2 [micron]
      mean_cer2(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_cer2(:,:,iz,imon) ), &
                      & real( n_cer2(:,:,iz,imon) ) )

      ! Cloud Effective Radius Liq [micron]
      mean_cerliq(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_cerliq(:,:,iz,imon) ), &
                      & real( n_cerliq(:,:,iz,imon) ) )

      ! Cloud Effective Radius Ice [micron]
      mean_cerice(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_cerice(:,:,iz,imon) ), &
                      & real( n_cerice(:,:,iz,imon) ) )

      ! Cloud Number Concentration 1 [# cm-3]
      mean_numc1(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_numc1(:,:,iz,imon) ), &
                      & real( n_numc1(:,:,iz,imon) ) )

      ! Cloud Number Concentration 2 [# cm-3]
      mean_numc2(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_numc2(:,:,iz,imon) ), &
                      & real( n_numc2(:,:,iz,imon) ) )

      ! Cloud Number Concentration Liq [# cm-3]
      mean_numliq(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_numliq(:,:,iz,imon) ), &
                      & real( n_numliq(:,:,iz,imon) ) )

      ! Cloud Number Concentration Ice [# cm-3]
      mean_numice(:,:,iz,imon) = &
         calc_ratio_2d( real( sum_numice(:,:,iz,imon) ), &
                      & real( n_numice(:,:,iz,imon) ) )

   enddo

   !--- Annual Statistics: 2D Variables
   if ( ooutyr ) then
      ! Averaging
      n_c4_valid(:,:,:,0) = 0
      n_c4_cloud(:,:,:,0) = 0
      n_c4_mask_invalid(:,:,:,0) = 0
      n_c4_type(:,:,:,:,0) = 0
      n_c4_category(:,:,:,:,0) = 0
      n_c4_type_qc(:,:,:,:,0) = 0
      n_c4_category_invalid(:,:,:,0) = 0
      do im = 1, 12
         if ( any(n_c4_valid(:,:,:,im) > huge(0_4)-n_c4_valid(:,:,:,0)) .or. &
            & any(n_c4_cloud(:,:,:,im) > huge(0_4)-n_c4_cloud(:,:,:,0)) .or. &
            & any(n_c4_mask_invalid(:,:,:,im) > &
            &     huge(0_4)-n_c4_mask_invalid(:,:,:,0)) .or. &
            & any(n_c4_type(:,:,:,:,im) > huge(0_4)-n_c4_type(:,:,:,:,0)) .or. &
            & any(n_c4_category(:,:,:,:,im) > huge(0_4)-n_c4_category(:,:,:,:,0)) .or. &
            & any(n_c4_type_qc(:,:,:,:,im) > huge(0_4)-n_c4_type_qc(:,:,:,:,0)) .or. &
            & any(n_c4_category_invalid(:,:,:,im) > &
            &     huge(0_4)-n_c4_category_invalid(:,:,:,0)) ) then
            write(error_unit,*) "### C4 annual int32 count overflow risk"
            stop 1
         endif
         n_c4_valid(:,:,:,0) = n_c4_valid(:,:,:,0) + n_c4_valid(:,:,:,im)
         n_c4_cloud(:,:,:,0) = n_c4_cloud(:,:,:,0) + n_c4_cloud(:,:,:,im)
         n_c4_mask_invalid(:,:,:,0) = n_c4_mask_invalid(:,:,:,0) + &
                                      n_c4_mask_invalid(:,:,:,im)
         n_c4_type(:,:,:,:,0) = n_c4_type(:,:,:,:,0) + n_c4_type(:,:,:,:,im)
         n_c4_category(:,:,:,:,0) = n_c4_category(:,:,:,:,0) + &
                                    n_c4_category(:,:,:,:,im)
         n_c4_type_qc(:,:,:,:,0) = n_c4_type_qc(:,:,:,:,0) + n_c4_type_qc(:,:,:,:,im)
         n_c4_category_invalid(:,:,:,0) = n_c4_category_invalid(:,:,:,0) + &
                                           n_c4_category_invalid(:,:,:,im)
      enddo
      n_obs2d(:,:,0)   = sum( n_obs2d(:,:,1:12),   DIM=3 )
      n_obs3d(:,:,:,0) = sum( n_obs3d(:,:,:,1:12), DIM=4 )
      n_cld3d(:,:,:,:,0) = sum( n_cld3d(:,:,:,:,1:12), DIM=5 )
      n_clr(:,:,0) = sum( n_clr(:,:,1:12), DIM=3 )

      ! Explicit clear/cloud classification population
      smpl_pdf(:,:) = real( n_clr(:,:,0) &
                          & + sum(n_cld(:,:,1:12), DIM=3) &
                          & + sum(n_drz(:,:,1:12), DIM=3) &
                          & + sum(n_prc(:,:,1:12), DIM=3), kind=sp )
      where ( smpl_pdf(:,:) .eq. 0.0_sp )
         pdf_clr(:,:,0) = undef
         pdf_cld(:,:,0) = undef
         pdf_drz(:,:,0) = undef
         pdf_prc(:,:,0) = undef
      else where
         pdf_clr(:,:,0) = 100.0_sp * real( n_clr(:,:,0), kind=sp ) / smpl_pdf(:,:)
         pdf_cld(:,:,0) = &
            100.0_sp * real( sum(n_cld(:,:,1:12), DIM=3), kind=sp ) / smpl_pdf(:,:)
         pdf_drz(:,:,0) = &
            100.0_sp * real( sum(n_drz(:,:,1:12), DIM=3), kind=sp ) / smpl_pdf(:,:)
         pdf_prc(:,:,0) = &
            100.0_sp * real( sum(n_prc(:,:,1:12), DIM=3), kind=sp ) / smpl_pdf(:,:)
      end where

      ! LWP
      mean_lwp(:,:,0) = &
         calc_ratio_2d( real( sum(sum_lwp(:,:,1:12), DIM=3) ), &
                      & real( sum(n_lwp(:,:,1:12), DIM=3) ) )

      ! IWP
      mean_iwp(:,:,0) = &
         calc_ratio_2d( real( sum(sum_iwp(:,:,1:12), DIM=3) ), &
                      & real( sum(n_iwp(:,:,1:12), DIM=3) ) )

      ! COT
      mean_cot(:,:,0) = &
         calc_ratio_2d( real( sum(sum_cot(:,:,1:12), DIM=3) ), &
                      & real( sum(n_cot(:,:,1:12), DIM=3) ) )

      ! CTH
      mean_cth(:,:,0) = &
         calc_ratio_2d( real( sum(sum_cth(:,:,1:12), DIM=3) ), &
                      & real( sum(n_cth(:,:,1:12), DIM=3) ) )

      ! C4 CTH
      mean_cth_c4(:,:,0) = &
         calc_ratio_2d( real( sum(sum_cth_c4(:,:,1:12), DIM=3) ), &
                      & real( sum(n_cth_c4(:,:,1:12), DIM=3) ) )

      ! C4 CTP
      mean_ctp_c4(:,:,0) = &
         calc_ratio_2d( real( sum(sum_ctp_c4(:,:,1:12), DIM=3) ), &
                      & real( sum(n_ctp_c4(:,:,1:12), DIM=3) ) )

   !--- Annual Statistics: 3D Variables
      do iz = 1, dimz

      ! Temperature [K]
      mean_temp(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_temp(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_obs3d(:,:,iz,1:12), DIM=3) ) )

      ! Specific Humidity [kg kg-1]
      mean_shum(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_shum(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_obs3d(:,:,iz,1:12), DIM=3) ) )

      ! Relative Humidity [%]
      mean_rhum(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_rhum(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_obs3d(:,:,iz,1:12), DIM=3) ) )

      ! Air Vertical Velocity [m s-1]
      mean_cldavel(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_cldavel(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_cldavel(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Doppler Velocity [m s-1]
      mean_clddvel(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_clddvel(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_clddvel(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Terminal Velocity 1 [m s-1]
      mean_cldtvel1(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_cldtvel1(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_cldtvel1(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Terminal Velocity 2 [m s-1]
      mean_cldtvel2(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_cldtvel2(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_cldtvel2(:,:,iz,1:12), DIM=3) ) )

      ! Total Cloud Terminal Velocity [m s-1]
      mean_cldtvelt(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_cldtvelt(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_cldtvelt(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Water Content 1 [g m-3]
      mean_cwc1(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_cwc1(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_cwc1(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Water Content 2 [g m-3]
      mean_cwc2(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_cwc2(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_cwc2(:,:,iz,1:12), DIM=3) ) )

      ! Ice Water Content [g m-3]
      mean_ciwc(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_ciwc(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_ciwc(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Effective Radius 1 [micron]
      mean_cer1(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_cer1(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_cer1(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Effective Radius 2 [micron]
      mean_cer2(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_cer2(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_cer2(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Effective Radius Liq [micron]
      mean_cerliq(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_cerliq(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_cerliq(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Effective Radius Ice [micron]
      mean_cerice(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_cerice(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_cerice(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Number Concentration 1 [# cm-3]
      mean_numc1(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_numc1(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_numc1(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Number Concentration 2 [# cm-3]
      mean_numc2(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_numc2(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_numc2(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Number Concentration Liq [# cm-3]
      mean_numliq(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_numliq(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_numliq(:,:,iz,1:12), DIM=3) ) )

      ! Cloud Number Concentration Ice [# cm-3]
      mean_numice(:,:,iz,0) = &
         calc_ratio_2d( real( sum(sum_numice(:,:,iz,1:12), DIM=3) ), &
                      & real( sum(n_numice(:,:,iz,1:12), DIM=3) ) )

      enddo
   endif ! ooutyr

   !--- Annual Output
   if ( ooutyr ) then
      ! update for accumulating whole analysis period
      if ( any(n_c4_valid(:,:,:,0) > huge(0_4)-n_c4_valid_all) .or. &
         & any(n_c4_cloud(:,:,:,0) > huge(0_4)-n_c4_cloud_all) .or. &
         & any(n_c4_mask_invalid(:,:,:,0) > huge(0_4)-n_c4_mask_invalid_all) .or. &
         & any(n_c4_type(:,:,:,:,0) > huge(0_4)-n_c4_type_all) .or. &
         & any(n_c4_category(:,:,:,:,0) > huge(0_4)-n_c4_category_all) .or. &
         & any(n_c4_type_qc(:,:,:,:,0) > huge(0_4)-n_c4_type_qc_all) .or. &
         & any(n_c4_category_invalid(:,:,:,0) > &
         &     huge(0_4)-n_c4_category_invalid_all) ) then
         write(error_unit,*) "### C4 all-period int32 count overflow risk"
         stop 1
      endif
      n_c4_valid_all = n_c4_valid_all + n_c4_valid(:,:,:,0)
      n_c4_cloud_all = n_c4_cloud_all + n_c4_cloud(:,:,:,0)
      n_c4_mask_invalid_all = n_c4_mask_invalid_all + n_c4_mask_invalid(:,:,:,0)
      n_c4_type_all = n_c4_type_all + n_c4_type(:,:,:,:,0)
      n_c4_category_all = n_c4_category_all + n_c4_category(:,:,:,:,0)
      n_c4_type_qc_all = n_c4_type_qc_all + n_c4_type_qc(:,:,:,:,0)
      n_c4_category_invalid_all = n_c4_category_invalid_all + &
                                  n_c4_category_invalid(:,:,:,0)
      where ( n_obs2d(:,:,0) .ne. undef )
         n_obs2d_all(:,:) = n_obs2d_all(:,:) + n_obs2d(:,:,0)
      end where
      where ( n_obs3d(:,:,:,0) .ne. undef )
         n_obs3d_all(:,:,:) = n_obs3d_all(:,:,:) + n_obs3d(:,:,:,0)
      end where
      where ( n_cld3d(:,:,:,:,0) .ne. undef )
         n_cld3d_all(:,:,:,:) = n_cld3d_all(:,:,:,:) + n_cld3d(:,:,:,:,0)
      end where
      ! 2D variables
      n_cld_all(:,:) = n_cld_all(:,:) + sum( n_cld(:,:,1:12), DIM=3 )
      n_drz_all(:,:) = n_drz_all(:,:) + sum( n_drz(:,:,1:12), DIM=3 )
      n_prc_all(:,:) = n_prc_all(:,:) + sum( n_prc(:,:,1:12), DIM=3 )
      n_clr_all(:,:) = n_clr_all(:,:) + n_clr(:,:,0)
      n_lwp_all(:,:) = n_lwp_all(:,:) + sum( n_lwp(:,:,1:12), DIM=3 )
      sum_lwp_all(:,:) = sum_lwp_all(:,:) + sum( sum_lwp(:,:,1:12), DIM=3 )
      n_iwp_all(:,:) = n_iwp_all(:,:) + sum( n_iwp(:,:,1:12), DIM=3 )
      sum_iwp_all(:,:) = sum_iwp_all(:,:) + sum( sum_iwp(:,:,1:12), DIM=3 )
      n_cot_all(:,:) = n_cot_all(:,:) + sum( n_cot(:,:,1:12), DIM=3 )
      sum_cot_all(:,:) = sum_cot_all(:,:) + sum( sum_cot(:,:,1:12), DIM=3 )
      n_cth_all(:,:) = n_cth_all(:,:) + sum( n_cth(:,:,1:12), DIM=3 )
      sum_cth_all(:,:) = sum_cth_all(:,:) + sum( sum_cth(:,:,1:12), DIM=3 )
      n_cth_c4_all(:,:) = n_cth_c4_all(:,:) + sum( n_cth_c4(:,:,1:12), DIM=3 )
      sum_cth_c4_all(:,:) = sum_cth_c4_all(:,:) + sum( sum_cth_c4(:,:,1:12), DIM=3 )
      n_ctp_c4_all(:,:) = n_ctp_c4_all(:,:) + sum( n_ctp_c4(:,:,1:12), DIM=3 )
      sum_ctp_c4_all(:,:) = sum_ctp_c4_all(:,:) + sum( sum_ctp_c4(:,:,1:12), DIM=3 )
      ! 3D variables
      sum_temp_all(:,:,:) = sum_temp_all(:,:,:) + sum( sum_temp(:,:,:,1:12), DIM=4 )
      sum_shum_all(:,:,:) = sum_shum_all(:,:,:) + sum( sum_shum(:,:,:,1:12), DIM=4 )
      sum_rhum_all(:,:,:) = sum_rhum_all(:,:,:) + sum( sum_rhum(:,:,:,1:12), DIM=4 )
      n_cldavel_all(:,:,:) = n_cldavel_all(:,:,:) + sum( n_cldavel(:,:,:,1:12), DIM=4 )
      sum_cldavel_all(:,:,:) = sum_cldavel_all(:,:,:) + sum( sum_cldavel(:,:,:,1:12), DIM=4 )
      n_clddvel_all(:,:,:) = n_clddvel_all(:,:,:) + sum( n_clddvel(:,:,:,1:12), DIM=4 )
      sum_clddvel_all(:,:,:) = sum_clddvel_all(:,:,:) + sum( sum_clddvel(:,:,:,1:12), DIM=4 )
      n_cldtvel1_all(:,:,:) = n_cldtvel1_all(:,:,:) + sum( n_cldtvel1(:,:,:,1:12), DIM=4 )
      sum_cldtvel1_all(:,:,:) = sum_cldtvel1_all(:,:,:) + sum( sum_cldtvel1(:,:,:,1:12), DIM=4 )
      n_cldtvel2_all(:,:,:) = n_cldtvel2_all(:,:,:) + sum( n_cldtvel2(:,:,:,1:12), DIM=4 )
      sum_cldtvel2_all(:,:,:) = sum_cldtvel2_all(:,:,:) + sum( sum_cldtvel2(:,:,:,1:12), DIM=4 )
      n_cldtvelt_all(:,:,:) = n_cldtvelt_all(:,:,:) + sum( n_cldtvelt(:,:,:,1:12), DIM=4 )
      sum_cldtvelt_all(:,:,:) = sum_cldtvelt_all(:,:,:) + sum( sum_cldtvelt(:,:,:,1:12), DIM=4 )
      n_cwc1_all(:,:,:) = n_cwc1_all(:,:,:) + sum( n_cwc1(:,:,:,1:12), DIM=4 )
      sum_cwc1_all(:,:,:) = sum_cwc1_all(:,:,:) + sum( sum_cwc1(:,:,:,1:12), DIM=4 )
      n_cwc2_all(:,:,:) = n_cwc2_all(:,:,:) + sum( n_cwc2(:,:,:,1:12), DIM=4 )
      sum_cwc2_all(:,:,:) = sum_cwc2_all(:,:,:) + sum( sum_cwc2(:,:,:,1:12), DIM=4 )
      n_ciwc_all(:,:,:) = n_ciwc_all(:,:,:) + sum( n_ciwc(:,:,:,1:12), DIM=4 )
      sum_ciwc_all(:,:,:) = sum_ciwc_all(:,:,:) + sum( sum_ciwc(:,:,:,1:12), DIM=4 )
      n_cer1_all(:,:,:) = n_cer1_all(:,:,:) + sum( n_cer1(:,:,:,1:12), DIM=4 )
      sum_cer1_all(:,:,:) = sum_cer1_all(:,:,:) + sum( sum_cer1(:,:,:,1:12), DIM=4 )
      n_cer2_all(:,:,:) = n_cer2_all(:,:,:) + sum( n_cer2(:,:,:,1:12), DIM=4 )
      sum_cer2_all(:,:,:) = sum_cer2_all(:,:,:) + sum( sum_cer2(:,:,:,1:12), DIM=4 )
      n_cerliq_all(:,:,:) = n_cerliq_all(:,:,:) + sum( n_cerliq(:,:,:,1:12), DIM=4 )
      sum_cerliq_all(:,:,:) = sum_cerliq_all(:,:,:) + sum( sum_cerliq(:,:,:,1:12), DIM=4 )
      n_cerice_all(:,:,:) = n_cerice_all(:,:,:) + sum( n_cerice(:,:,:,1:12), DIM=4 )
      sum_cerice_all(:,:,:) = sum_cerice_all(:,:,:) + sum( sum_cerice(:,:,:,1:12), DIM=4 )
      n_numc1_all(:,:,:) = n_numc1_all(:,:,:) + sum( n_numc1(:,:,:,1:12), DIM=4 )
      sum_numc1_all(:,:,:) = sum_numc1_all(:,:,:) + sum( sum_numc1(:,:,:,1:12), DIM=4 )
      n_numc2_all(:,:,:) = n_numc2_all(:,:,:) + sum( n_numc2(:,:,:,1:12), DIM=4 )
      sum_numc2_all(:,:,:) = sum_numc2_all(:,:,:) + sum( sum_numc2(:,:,:,1:12), DIM=4 )
      n_numliq_all(:,:,:) = n_numliq_all(:,:,:) + sum( n_numliq(:,:,:,1:12), DIM=4 )
      sum_numliq_all(:,:,:) = sum_numliq_all(:,:,:) + sum( sum_numliq(:,:,:,1:12), DIM=4 )
      n_numice_all(:,:,:) = n_numice_all(:,:,:) + sum( n_numice(:,:,:,1:12), DIM=4 )
      sum_numice_all(:,:,:) = sum_numice_all(:,:,:) + sum( sum_numice(:,:,:,1:12), DIM=4 )
   endif ! ooutyr

   !--- All-Term Output: 2D Variables
   if ( ofinal ) then
      ! Explicit clear/cloud classification population
      smpl_pdf(:,:) = real( n_clr_all(:,:) + n_cld_all(:,:) &
                          & + n_drz_all(:,:) + n_prc_all(:,:), kind=sp )
      where ( smpl_pdf(:,:) .eq. 0.0_sp )
         pdf_clr_all(:,:) = undef
         pdf_cld_all(:,:) = undef
         pdf_drz_all(:,:) = undef
         pdf_prc_all(:,:) = undef
      else where
         pdf_clr_all(:,:) = 100.0_sp * real( n_clr_all(:,:), kind=sp ) / smpl_pdf(:,:)
         pdf_cld_all(:,:) = 100.0_sp * real( n_cld_all(:,:), kind=sp ) / smpl_pdf(:,:)
         pdf_drz_all(:,:) = 100.0_sp * real( n_drz_all(:,:), kind=sp ) / smpl_pdf(:,:)
         pdf_prc_all(:,:) = 100.0_sp * real( n_prc_all(:,:), kind=sp ) / smpl_pdf(:,:)
      end where

      ! LWP
      mean_lwp_all(:,:) = &
         calc_ratio_2d( real( sum_lwp_all(:,:) ), &
                      & real( n_lwp_all(:,:) ) )

      ! IWP
      mean_iwp_all(:,:) = &
         calc_ratio_2d( real( sum_iwp_all(:,:) ), &
                      & real( n_iwp_all(:,:) ) )

      ! COT
      mean_cot_all(:,:) = &
         calc_ratio_2d( real( sum_cot_all(:,:) ), &
                      & real( n_cot_all(:,:) ) )

      ! CTH
      mean_cth_all(:,:) = &
         calc_ratio_2d( real( sum_cth_all(:,:) ), &
                      & real( n_cth_all(:,:) ) )

      ! C4 CTH
      mean_cth_c4_all(:,:) = &
         calc_ratio_2d( real( sum_cth_c4_all(:,:) ), &
                      & real( n_cth_c4_all(:,:) ) )

      ! C4 CTP
      mean_ctp_c4_all(:,:) = &
         calc_ratio_2d( real( sum_ctp_c4_all(:,:) ), &
                      & real( n_ctp_c4_all(:,:) ) )

      !--- All-Term Output: 3D Variables
      do iz = 1, dimz

      ! Temperature [K]
      mean_temp_all(:,:,iz) = &
         calc_ratio_2d( real( sum_temp_all(:,:,iz) ), &
                      & real( n_obs3d_all(:,:,iz) ) )

      ! Specific Humidity [kg kg-1]
      mean_shum_all(:,:,iz) = &
         calc_ratio_2d( real( sum_shum_all(:,:,iz) ), &
                      & real( n_obs3d_all(:,:,iz) ) )

      ! Relative Humidity [%]
      mean_rhum_all(:,:,iz) = &
         calc_ratio_2d( real( sum_rhum_all(:,:,iz) ), &
                      & real( n_obs3d_all(:,:,iz) ) )

      ! Air Vertical Velocity [m s-1]
      mean_cldavel_all(:,:,iz) = &
         calc_ratio_2d( real( sum_cldavel_all(:,:,iz) ), &
                      & real( n_cldavel_all(:,:,iz) ) )

      ! Cloud Doppler Velocity [m s-1]
      mean_clddvel_all(:,:,iz) = &
         calc_ratio_2d( real( sum_clddvel_all(:,:,iz) ), &
                      & real( n_clddvel_all(:,:,iz) ) )

      ! Cloud Terminal Velocity 1 [m s-1]
      mean_cldtvel1_all(:,:,iz) = &
         calc_ratio_2d( real( sum_cldtvel1_all(:,:,iz) ), &
                      & real( n_cldtvel1_all(:,:,iz) ) )

      ! Cloud Terminal Velocity 2 [m s-1]
      mean_cldtvel2_all(:,:,iz) = &
         calc_ratio_2d( real( sum_cldtvel2_all(:,:,iz) ), &
                      & real( n_cldtvel2_all(:,:,iz) ) )

      ! Total Cloud Terminal Velocity [m s-1]
      mean_cldtvelt_all(:,:,iz) = &
         calc_ratio_2d( real( sum_cldtvelt_all(:,:,iz) ), &
                      & real( n_cldtvelt_all(:,:,iz) ) )

      ! Cloud Water Content 1 [g m-3]
      mean_cwc1_all(:,:,iz) = &
         calc_ratio_2d( real( sum_cwc1_all(:,:,iz) ), &
                      & real( n_cwc1_all(:,:,iz) ) )

      ! Cloud Water Content 2 [g m-3]
      mean_cwc2_all(:,:,iz) = &
         calc_ratio_2d( real( sum_cwc2_all(:,:,iz) ), &
                      & real( n_cwc2_all(:,:,iz) ) )

      ! Ice Water Content [g m-3]
      mean_ciwc_all(:,:,iz) = &
         calc_ratio_2d( real( sum_ciwc_all(:,:,iz) ), &
                      & real( n_ciwc_all(:,:,iz) ) )

      ! Cloud Effective Radius 1 [micron]
      mean_cer1_all(:,:,iz) = &
         calc_ratio_2d( real( sum_cer1_all(:,:,iz) ), &
                      & real( n_cer1_all(:,:,iz) ) )

      ! Cloud Effective Radius 2 [micron]
      mean_cer2_all(:,:,iz) = &
         calc_ratio_2d( real( sum_cer2_all(:,:,iz) ), &
                      & real( n_cer2_all(:,:,iz) ) )

      ! Cloud Effective Radius Liq [micron]
      mean_cerliq_all(:,:,iz) = &
         calc_ratio_2d( real( sum_cerliq_all(:,:,iz) ), &
                      & real( n_cerliq_all(:,:,iz) ) )

      ! Cloud Effective Radius Ice [micron]
      mean_cerice_all(:,:,iz) = &
         calc_ratio_2d( real( sum_cerice_all(:,:,iz) ), &
                      & real( n_cerice_all(:,:,iz) ) )

      ! Cloud Number Concentration 1 [# cm-3]
      mean_numc1_all(:,:,iz) = &
         calc_ratio_2d( real( sum_numc1_all(:,:,iz) ), &
                      & real( n_numc1_all(:,:,iz) ) )

      ! Cloud Number Concentration 2 [# cm-3]
      mean_numc2_all(:,:,iz) = &
         calc_ratio_2d( real( sum_numc2_all(:,:,iz) ), &
                      & real( n_numc2_all(:,:,iz) ) )

      ! Cloud Number Concentration Liq [# cm-3]
      mean_numliq_all(:,:,iz) = &
         calc_ratio_2d( real( sum_numliq_all(:,:,iz) ), &
                      & real( n_numliq_all(:,:,iz) ) )

      ! Cloud Number Concentration Ice [# cm-3]
      mean_numice_all(:,:,iz) = &
         calc_ratio_2d( real( sum_numice_all(:,:,iz) ), &
                      & real( n_numice_all(:,:,iz) ) )

      enddo

   endif  ! ofinal

   !--- CFAD
   ! Monthly CFAD
   if ( ooutmon ) then
      do ilev = 1, dimlev
         do idbz = 1, dimdbz
            pdf_cfad( idbz,ilev,:,imon ) = &
               calc_ratio_1d( 100.0_sp * real( n_cfad_bin( idbz,ilev,:,imon ) ), &
                              real( n_cfad_obs( :,imon ) ) )
         enddo
      enddo
   endif

   ! Annual CFAD
   if ( ooutyr ) then
      do ilev = 1, dimlev
         do idbz = 1, dimdbz
            pdf_cfad( idbz,ilev,:,0 ) = &
               calc_ratio_1d( 100.0_sp * real( sum( n_cfad_bin( idbz,ilev,:,1:12 ), DIM=2 ) ), &
                              real( sum( n_cfad_obs( :,1:12 ), DIM=2 ) ) )
         enddo
      enddo

      ! update for accumulating whole analysis period
      n_cfad_obs_all( : ) = n_cfad_obs_all( : ) + sum( n_cfad_obs( :,1:12 ), DIM=2 )
      n_cfad_bin_all( :,:,: ) = n_cfad_bin_all( :,:,: ) + &
                                     sum( n_cfad_bin( :,:,:,1:12 ), DIM=4 )
   endif

   ! All-Term CFAD
   if ( ofinal ) then
      do ilev = 1, dimlev
         do idbz = 1, dimdbz
            pdf_cfad_all( idbz,ilev,: ) = &
               calc_ratio_1d( 100.0_sp * real( n_cfad_bin_all( idbz,ilev,: ) ), &
                              real( n_cfad_obs_all( : ) ) )
         enddo
      enddo
   endif

   !--- ISCCP-like Histogram
   ! Monthly ISCCP
   if ( ooutmon ) then
      smpl_isccp(:) = real( &
         sum( sum( n_isccp(1:dimcot,1:dimctp,:,imon), DIM=1 ), DIM=1 ), kind=sp )
      do ictp = 1, dimctp
         do icot = 1, dimcot
            isccp_fq( icot,ictp,:,imon ) = &
               calc_ratio_1d( 100.0_sp * real(n_isccp(icot,ictp,:,imon), kind=sp), smpl_isccp(:) )
            where ( n_isccp( icot,ictp,:,imon ) .eq. 0 )
               isccp_dv( icot,ictp,:,imon ) = undef
            else where
               isccp_dv( icot,ictp,:,imon ) = &
                  calc_ratio_1d( real( sum_isccpdv(icot,ictp,:,imon) ), &
                                 real( n_isccp(icot,ictp,:,imon), kind=sp ) )
            end where
         enddo
      enddo

      smpl_isccp_c4(:) = real( &
         sum( sum( n_isccp_c4(1:dimcot,1:dimctp,:,imon), DIM=1 ), DIM=1 ), kind=sp )
      do ictp = 1, dimctp
         do icot = 1, dimcot
            isccp_fq_c4( icot,ictp,:,imon ) = &
               calc_ratio_1d( 100.0_sp * real(n_isccp_c4(icot,ictp,:,imon), kind=sp), smpl_isccp_c4(:) )
            where ( n_isccpdv_c4( icot,ictp,:,imon ) .eq. 0 )
               isccp_dv_c4( icot,ictp,:,imon ) = undef
            else where
               isccp_dv_c4( icot,ictp,:,imon ) = &
                  calc_ratio_1d( real( sum_isccpdv_c4(icot,ictp,:,imon) ), &
                                 real( n_isccpdv_c4(icot,ictp,:,imon), kind=sp ) )
            end where
         enddo
      enddo
   endif

   ! Annual ISCCP
   if ( ooutyr ) then
      smpl_isccp(:) = real( &
         sum( &                                          ! month
              sum( &                                     ! ctp
                   sum( n_isccp(1:dimcot,1:dimctp,:,1:12), DIM=1 ), &  ! summing cot => (dimctp,nreg,nmonth)
              dim=1 ), &                                 ! summing ctp => (nreg,nmonth)
         dim=2 ), &                                      ! summing month => (nreg)
      kind=sp )
      do ictp = 1, dimctp
         do icot = 1, dimcot
            isccp_fq( icot,ictp,:,0 ) = &
               calc_ratio_1d( 100.0_sp * real( sum( n_isccp(icot,ictp,:,1:12), DIM=2 ), kind=sp ), &
                              smpl_isccp(:) )
            where ( sum( n_isccp( icot,ictp,:,1:12 ), DIM=2 ) .eq. 0 )
               isccp_dv( icot,ictp,:,0 ) = undef
            else where
               isccp_dv( icot,ictp,:,0 ) = &
                  calc_ratio_1d( real( sum( sum_isccpdv(icot,ictp,:,1:12), DIM=2 ) ), &
                                 real( sum( n_isccp(icot,ictp,:,1:12), DIM=2 ), kind=sp ) )
            end where
         enddo
      enddo

      smpl_isccp_c4(:) = real( &
         sum( &
              sum( &
                   sum( n_isccp_c4(1:dimcot,1:dimctp,:,1:12), DIM=1 ), &
              dim=1 ), &
         dim=2 ), &
      kind=sp )
      do ictp = 1, dimctp
         do icot = 1, dimcot
            isccp_fq_c4( icot,ictp,:,0 ) = &
               calc_ratio_1d( 100.0_sp * real( sum( n_isccp_c4(icot,ictp,:,1:12), DIM=2 ), kind=sp ), &
                              smpl_isccp_c4(:) )
            where ( sum( n_isccpdv_c4( icot,ictp,:,1:12 ), DIM=2 ) .eq. 0 )
               isccp_dv_c4( icot,ictp,:,0 ) = undef
            else where
               isccp_dv_c4( icot,ictp,:,0 ) = &
                  calc_ratio_1d( real( sum( sum_isccpdv_c4(icot,ictp,:,1:12), DIM=2 ) ), &
                                 real( sum( n_isccpdv_c4(icot,ictp,:,1:12), DIM=2 ), kind=sp ) )
            end where
         enddo
      enddo

      ! update for accumulating whole analysis period
      n_isccp_all( :,:,: ) = n_isccp_all( :,:,: ) + &
                                     sum( n_isccp( :,:,:,1:12 ), DIM=4 )
      sum_isccpdv_all( :,:,: ) = sum_isccpdv_all( :,:,: ) + &
                                     sum( sum_isccpdv( :,:,:,1:12 ), DIM=4 )
      n_isccp_c4_all( :,:,: ) = n_isccp_c4_all( :,:,: ) + &
                                     sum( n_isccp_c4( :,:,:,1:12 ), DIM=4 )
      n_isccpdv_c4_all( :,:,: ) = n_isccpdv_c4_all( :,:,: ) + &
                                     sum( n_isccpdv_c4( :,:,:,1:12 ), DIM=4 )
      sum_isccpdv_c4_all( :,:,: ) = sum_isccpdv_c4_all( :,:,: ) + &
                                     sum( sum_isccpdv_c4( :,:,:,1:12 ), DIM=4 )
   endif

   ! All-Term ISCCP
   if ( ofinal ) then
      smpl_isccp(:) = real( &
         sum( &
              sum( n_isccp_all(1:dimcot,1:dimctp,:), DIM=1 ), &
         DIM=1 ), &
      kind=sp )
      do ictp = 1, dimctp
         do icot = 1, dimcot
            isccp_fq_all( icot,ictp,: ) = &
               calc_ratio_1d( 100.0_sp * real( n_isccp_all( icot,ictp,: ), kind=sp ), &
                              smpl_isccp(:) )
            where ( n_isccp_all( icot,ictp,: ) .eq. 0 )
               isccp_dv_all( icot,ictp,: ) = undef
            else where
               isccp_dv_all( icot,ictp,: ) = &
                  calc_ratio_1d( real( sum_isccpdv_all( icot,ictp,: ) ), &
                                 real( n_isccp_all( icot,ictp,: ), kind=sp ) )
            end where
         enddo
      enddo

      smpl_isccp_c4(:) = real( &
         sum( &
              sum( n_isccp_c4_all(1:dimcot,1:dimctp,:), DIM=1 ), &
         DIM=1 ), &
      kind=sp )
      do ictp = 1, dimctp
         do icot = 1, dimcot
            isccp_fq_c4_all( icot,ictp,: ) = &
               calc_ratio_1d( 100.0_sp * real( n_isccp_c4_all( icot,ictp,: ), kind=sp ), &
                              smpl_isccp_c4(:) )
            where ( n_isccpdv_c4_all( icot,ictp,: ) .eq. 0 )
               isccp_dv_c4_all( icot,ictp,: ) = undef
            else where
               isccp_dv_c4_all( icot,ictp,: ) = &
                  calc_ratio_1d( real( sum_isccpdv_c4_all( icot,ictp,: ) ), &
                                 real( n_isccpdv_c4_all( icot,ictp,: ), kind=sp ) )
            end where
         enddo
      enddo
   endif

end subroutine calc_pdf

function calc_ratio_1d(numerator, denominator) result(ratio)
   implicit none
   real(sp), dimension(:), intent(in) :: numerator, denominator
   real(sp), dimension(:), allocatable :: ratio
   integer :: n

   if ( size(numerator) /= size(denominator) ) then
      write(error_unit,*) '### calc_ratio_1d: size mismatch!'
      stop
   end if

   n = size(denominator)
   allocate( ratio(n) )

   where ( denominator .gt. 0.0_sp )
      ratio = numerator / denominator
   else where
      ratio = undef
   end where

end function calc_ratio_1d

function calc_ratio_2d(numerator, denominator) result(ratio)
   implicit none
   real(sp), dimension(:,:), intent(in) :: numerator, denominator
   real(sp), dimension(:,:), allocatable :: ratio
   integer :: size1, size2
   integer, dimension(2) :: array_shape

   array_shape = shape(denominator)
   size1 = array_shape(1)
   size2 = array_shape(2)
   allocate( ratio(size1, size2) )

   where ( denominator .gt. 0.0_sp )
      ratio = numerator / denominator
   else where
      ratio = undef
   end where

end function calc_ratio_2d

subroutine increment_c4_count(count, count_name)
   implicit none
   integer(kind=4), intent(inout) :: count
   character(*), intent(in) :: count_name

   if ( count .eq. huge(count) ) then
      write(error_unit,*) "### C4 int32 count overflow: ", trim(count_name)
      stop 1
   endif
   count = count + 1_4
end subroutine increment_c4_count

end module output
