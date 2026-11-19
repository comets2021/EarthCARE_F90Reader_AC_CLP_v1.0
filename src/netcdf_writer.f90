module netcdf_writer
    use, intrinsic :: iso_fortran_env, only: error_unit
    use parameters
    use netcdf

    implicit none
    private check

contains

subroutine output_netcdf_cfad(outfname, gdbz, glev, pdf_cfad, &
                              aggregation_type, period_start, period_end)
    implicit none
    character(len=*), intent(in) :: outfname
    character(len=*), intent(in) :: aggregation_type, period_start, period_end
    real(sp), dimension(:), intent(in) :: gdbz(idbz_max + 1)  ! grid value (CFAD dBZ)
    real(sp), dimension(:), intent(in) :: glev(ilev_max + 1)  ! grid value (CFAD Lev)
    real(sp), dimension(:,:,:), intent(in) :: pdf_cfad
    character(len=6),   dimension(:), allocatable :: reg_name_list
    character(len=128), dimension(:), allocatable :: long_reg_name_list

    integer :: ncid
    integer :: file_type, dimt, idx, ireg
    integer, dimension(2) :: dim_id_list
    integer, dimension(1) :: dbz_dim_id, lev_dim_id
    integer :: dbz_var_id, lev_var_id
    integer, dimension(nreg) :: var_id_list

    ! set dimension of time
    dimt = 1

    ! select NetCDF format
    file_type = NF90_NETCDF4
    !file_type = NF90_CLASSIC_MODEL

    allocate(reg_name_list(nreg))
    allocate(long_reg_name_list(nreg))

    select case(nreg)
    case(3)
        reg_name_list = (/ "tropic", "midlat", "global" /)
        long_reg_name_list = (/ "Tropics", "Midlatitudes", "Global" /)
    case(16)
        reg_name_list = (/ "trwp", "itcz", "spcz", "nepa", "cast", "peru", "natl", "nami", &
                           "aust", "japn", "eqct", "easi", "soce", "arct", "anta", "glob" /)
        long_reg_name_list = (/ "Tropical Warm Pool (5S-15N,70E-135E)", &
                            "ITCZ (5N-15N,140E-140W)", &
                            "SPCZ (15S-5S,150E-130W)", &
                            "North East Pacific (25N-50N,160W-135W)", &
                            "California Stratocumulus (15N-35N,130W-110W)", &
                            "Peruvian (10S-30S,90W-70W)", &
                            "North Atlantic (30N-60N,45W-10W)", &
                            "Namibian (30S-0S,25W-15E)", &
                            "Australian (15S-35S,80E-115E)", &
                            "Japan (25N-50N,125E-150E)", &
                            "Eqt. Cold Tongue (5S-5N,130W-85W)", &
                            "Eastern Asia (20N-40N,105E-120E)", &
                            "Southern Ocean (40S-60S)", &
                            "Arctic (66N-90N)", &
                            "Antarctic (66S-90S)", &
                            "Global"/)
    end select

    ! create output file
    call check( nf90_create( trim(outfname), file_type, ncid ))

    ! define dimensions
    call check( nf90_def_dim( ncid, "dbz", dimdbz, dbz_dim_id(1) ) )
    call check( nf90_def_dim( ncid, "lev", dimlev, lev_dim_id(1) ) )

    ! define variables
    call define_var( ncid, var_id=dbz_var_id, &
                     dim_id=dbz_dim_id, var_name="dbz", var_type=NF90_FLOAT, &
                     units="dBZe", long_name="Reflectivity", &
                     standard_name="reflectivity", axis="X")

    call define_var( ncid, var_id=lev_var_id, &
                     dim_id=lev_dim_id, var_name="lev", var_type=NF90_FLOAT, &
                     units="km", long_name="Altitude", &
                     standard_name="altitude", axis="Y")

    dim_id_list = (/dbz_dim_id, lev_dim_id/)

    do ireg = 1, nreg
       idx = ireg
       call define_var( ncid, var_id=var_id_list(idx), dim_id=dim_id_list, &
                        var_name="cfad_" // trim(reg_name_list(ireg)), &
                        var_type=NF90_FLOAT, &
                        units="none", &
                        long_name=trim(long_reg_name_list(ireg)))
    enddo

    ! set global attributes
    call check( nf90_put_att(ncid,nf90_global,"title", &
                "EarthCARE Level3 (CFAD)"))
    call put_global_metadata(ncid, aggregation_type, period_start, period_end)

    call check( nf90_enddef( ncid ) )

    ! write variables
    call check( nf90_put_var(ncid, dbz_var_id, gdbz(2::2)) )
    call check( nf90_put_var(ncid, lev_var_id, glev(2::2)) )

    do ireg = 1, nreg
       idx = ireg
       call check( nf90_put_var( ncid, var_id_list(idx), pdf_cfad(:,:,ireg) ) )
    enddo

    deallocate(reg_name_list)
    deallocate(long_reg_name_list)

    call check( nf90_close(ncid) )

end subroutine output_netcdf_cfad

!subroutine output_netcdf_isccp(outfname, gcot, gctp, isccp_fq, isccp_dv)
subroutine output_netcdf_isccp(outfname, isccp_fq, isccp_dv, isccp_fq_c4, isccp_dv_c4, &
                               aggregation_type, period_start, period_end)
    implicit none
    character(len=*), intent(in) :: outfname
    character(len=*), intent(in) :: aggregation_type, period_start, period_end

    ! ISCCP bin *edges* (nonlinear). size = dim+1
!    real(sp), intent(in) :: gcot(dimcot+1)   ! Cloud Optical Thickness edges
!    real(sp), intent(in) :: gctp(dimctp+1)   ! Cloud-Top Pressure edges [hPa]

    ! ISCCP joint-hist output (frequency and Doppler velocity)
    real(sp), intent(in) :: isccp_fq(dimcot, dimctp, nreg)
    real(sp), intent(in) :: isccp_dv(dimcot, dimctp, nreg)
    real(sp), intent(in) :: isccp_fq_c4(dimcot, dimctp, nreg)
    real(sp), intent(in) :: isccp_dv_c4(dimcot, dimctp, nreg)

    character(len=6),   dimension(:), allocatable :: reg_name_list
    character(len=128), dimension(:), allocatable :: long_reg_name_list

    integer :: ncid
    integer :: file_type, idx, ireg, i
    integer, dimension(2) :: dim_id_list
    integer, dimension(1) :: cot_dim_id, ctp_dim_id, bnds_dim_id
    integer :: cot_var_id, ctp_var_id
    integer :: cot_bnds_var_id, ctp_bnds_var_id
    integer, dimension(2) :: cot_bnds_dimids, ctp_bnds_dimids
    integer, dimension(nreg) :: var_id_list_fq, var_id_list_dv
    integer, dimension(nreg) :: var_id_list_fq_c4, var_id_list_dv_c4
    real(sp) :: cot_center(dimcot), ctp_center(dimctp)
    real(sp) :: cot_bnds(dimcot,2), ctp_bnds(dimctp,2)

    allocate(reg_name_list(nreg))
    allocate(long_reg_name_list(nreg))

    select case(nreg)
    case(3)
        reg_name_list = (/ "tropic", "midlat", "global" /)
        long_reg_name_list = (/ "Tropics", "Midlatitudes", "Global" /)
    case(16)
        reg_name_list = (/ "trwp", "itcz", "spcz", "nepa", "cast", "peru", "natl", "nami", &
                           "aust", "japn", "eqct", "easi", "soce", "arct", "anta", "glob" /)
        long_reg_name_list = (/ "Tropical Warm Pool (5S-15N,70E-135E)", &
                            "ITCZ (5N-15N,140E-140W)", &
                            "SPCZ (15S-5S,150E-130W)", &
                            "North East Pacific (25N-50N,160W-135W)", &
                            "California Stratocumulus (15N-35N,130W-110W)", &
                            "Peruvian (10S-30S,90W-70W)", &
                            "North Atlantic (30N-60N,45W-10W)", &
                            "Namibian (30S-0S,25W-15E)", &
                            "Australian (15S-35S,80E-115E)", &
                            "Japan (25N-50N,125E-150E)", &
                            "Eqt. Cold Tongue (5S-5N,130W-85W)", &
                            "Eastern Asia (20N-40N,105E-120E)", &
                            "Southern Ocean (40S-60S)", &
                            "Arctic (66N-90N)", &
                            "Antarctic (66S-90S)", &
                            "Global" /)
    end select

    ! Build bin centers + bounds
    do i = 1, dimcot
       cot_bnds(i,1) = gcot(i)
       cot_bnds(i,2) = gcot(i+1)

       ! Recommended: geometric mean for tau bins if both edges > 0
       if (gcot(i) .gt. 0._sp .and. gcot(i+1) .gt. 0._sp) then
          cot_center(i) = sqrt( gcot(i) * gcot(i+1) )
       else
          cot_center(i) = 0.5_sp * (gcot(i) + gcot(i+1))
       endif
    enddo

    do i = 1, dimctp
       ctp_bnds(i,1) = gctp(i)
       ctp_bnds(i,2) = gctp(i+1)

       ! arithmetic mean is fine for pressure bins
       ctp_center(i) = 0.5_sp * (gctp(i) + gctp(i+1))
    enddo

    ! Create output file
    file_type = NF90_NETCDF4
    call check( nf90_create( trim(outfname), file_type, ncid ) )

    ! Define dimensions
    call check( nf90_def_dim( ncid, "cot",  dimcot, cot_dim_id(1) ) )
    call check( nf90_def_dim( ncid, "ctp",  dimctp, ctp_dim_id(1) ) )
    call check( nf90_def_dim( ncid, "bnds", 2,      bnds_dim_id(1) ) )

    ! Define coordinate variables (centers)
    call define_var( ncid, var_id=cot_var_id, &
                     dim_id=cot_dim_id, var_name="cot", var_type=NF90_FLOAT, &
                     units="1", long_name="Cloud Optical Thickness (bin center)", &
                     standard_name="cloud_optical_thickness", axis="X" )

    call define_var( ncid, var_id=ctp_var_id, &
                     dim_id=ctp_dim_id, var_name="ctp", var_type=NF90_FLOAT, &
                     units="hPa", long_name="Cloud-Top Pressure (bin center)", &
                     standard_name="air_pressure", axis="Y" )
    call check( nf90_put_att(ncid, ctp_var_id, "positive", "down") )

    ! Define bounds variables (N x 2)
    cot_bnds_dimids = (/ cot_dim_id(1), bnds_dim_id(1) /)
    ctp_bnds_dimids = (/ ctp_dim_id(1), bnds_dim_id(1) /)

    call define_var( ncid, var_id=cot_bnds_var_id, &
                     dim_id=cot_bnds_dimids, var_name="cot_bnds", var_type=NF90_FLOAT, &
                     units="1", long_name="Cloud Optical Thickness bin bounds" )

    call define_var( ncid, var_id=ctp_bnds_var_id, &
                     dim_id=ctp_bnds_dimids, var_name="ctp_bnds", var_type=NF90_FLOAT, &
                     units="hPa", long_name="Cloud-Top Pressure bin bounds" )

    ! Link bounds to coordinates
    call check( nf90_put_att(ncid, cot_var_id, "bounds", "cot_bnds") )
    call check( nf90_put_att(ncid, ctp_var_id, "bounds", "ctp_bnds") )

    ! Define histogram variables
    dim_id_list = (/ cot_dim_id(1), ctp_dim_id(1) /)

    do ireg = 1, nreg
       idx = ireg
       call define_var( ncid, var_id=var_id_list_fq(idx), dim_id=dim_id_list, &
                        var_name="isccp_fq_" // trim(reg_name_list(ireg)), &
                        var_type=NF90_FLOAT, units="%", &
                        long_name=trim(long_reg_name_list(ireg)) // " frequency" )
       call define_var( ncid, var_id=var_id_list_dv(idx), dim_id=dim_id_list, &
                        var_name="isccp_dv_" // trim(reg_name_list(ireg)), &
                        var_type=NF90_FLOAT, units="m s-1", &
                        long_name=trim(long_reg_name_list(ireg)) // " mean Doppler velocity" )
    enddo

    do ireg = 1, nreg
       idx = ireg
       call define_var( ncid, var_id=var_id_list_fq_c4(idx), dim_id=dim_id_list, &
                        var_name="isccp_fq_c4_" // trim(reg_name_list(ireg)), &
                        var_type=NF90_FLOAT, units="%", &
                        long_name=trim(long_reg_name_list(ireg)) // &
                                  " COT-C4 CTP joint frequency" )
       call check( nf90_put_att(ncid, var_id_list_fq_c4(idx), "cloud_top_definition", &
                   "Highest valid C4 Cloud Mask (CPR or ATLID) cloud bin" ) )
       call check( nf90_put_att(ncid, var_id_list_fq_c4(idx), "population", &
                   "Profiles with valid COT and C4 cloud-top pressure in histogram range; Doppler not required" ) )
       call define_var( ncid, var_id=var_id_list_dv_c4(idx), dim_id=dim_id_list, &
                        var_name="isccp_dv_c4_" // trim(reg_name_list(ireg)), &
                        var_type=NF90_FLOAT, units="m s-1", &
                        long_name=trim(long_reg_name_list(ireg)) // &
                                  " mean Doppler velocity in COT-C4 CTP cells" )
       call check( nf90_put_att(ncid, var_id_list_dv_c4(idx), "cloud_top_definition", &
                   "Highest valid C4 Cloud Mask (CPR or ATLID) cloud bin" ) )
       call check( nf90_put_att(ncid, var_id_list_dv_c4(idx), "population", &
                   "Profiles in the frequency population with a valid profile-mean Doppler velocity" ) )
    enddo

    ! Global attributes
    call check( nf90_put_att(ncid, nf90_global, "title", &
                "EarthCARE Level3 (ISCCP-like Joint Histogram)" ) )
    call put_global_metadata(ncid, aggregation_type, period_start, period_end)

    call check( nf90_enddef( ncid ) )

    ! Write coordinate variables
    call check( nf90_put_var(ncid, cot_var_id, cot_center) )
    call check( nf90_put_var(ncid, ctp_var_id, ctp_center) )

    call check( nf90_put_var(ncid, cot_bnds_var_id, cot_bnds) )
    call check( nf90_put_var(ncid, ctp_bnds_var_id, ctp_bnds) )

    ! Write data variables
    do ireg = 1, nreg
       idx = ireg
       call check( nf90_put_var( ncid, var_id_list_fq(idx), isccp_fq(:,:,ireg) ) )
       call check( nf90_put_var( ncid, var_id_list_dv(idx), isccp_dv(:,:,ireg) ) )
       call check( nf90_put_var( ncid, var_id_list_fq_c4(idx), isccp_fq_c4(:,:,ireg) ) )
       call check( nf90_put_var( ncid, var_id_list_dv_c4(idx), isccp_dv_c4(:,:,ireg) ) )
    enddo

    deallocate(reg_name_list)
    deallocate(long_reg_name_list)

    call check( nf90_close(ncid) )

end subroutine output_netcdf_isccp

subroutine output_netcdf_map2d( &
    glon, glat, &
    n_obs2d,  pdf_clr,  pdf_cld,  pdf_drz, &
    pdf_prc,  mean_lwp, mean_iwp, &
    mean_cot, mean_cth, mean_cth_c4, mean_ctp_c4, &
    outfname, output_type, aggregation_type, period_start, period_end)

    implicit none
    real(sp), dimension(:), intent(in) :: glon
    real(sp), dimension(:), intent(in) :: glat
    real(sp), dimension(:,:), intent(in) :: n_obs2d
    real(sp), dimension(:,:), intent(in) :: pdf_clr,  pdf_cld,  pdf_drz,  pdf_prc
    real(sp), dimension(:,:), intent(in) :: mean_lwp, mean_iwp, mean_cot, mean_cth
    real(sp), dimension(:,:), intent(in) :: mean_cth_c4, mean_ctp_c4
    character(len=*), intent(in) :: outfname, output_type
    character(len=*), intent(in) :: aggregation_type, period_start, period_end

    integer, dimension(0:nmon) :: time

    integer :: ncid
    integer, dimension(1) :: lon_dim_id, lat_dim_id
    integer, dimension(nvar2d+2) :: var_id_list
    character(len=20), dimension(nvar2d+2) :: var_name_list, units_list
    character(len=128), dimension(nvar2d+2) :: long_name_list
    integer :: lon_var_id, lat_var_id
    integer, dimension(2) :: dim_id_list
    integer :: file_type, mon_ini, mon_fin, dimt, ivar

    time = (/0:nmon/)
    select case(output_type)
    case("12_months")
        mon_ini = 1
        mon_fin = 12
        dimt = 12
    case("annual")
        mon_ini = 1
        mon_fin = 1
        dimt = 1
    case("monthly")
        mon_ini = 1
        mon_fin = 1
        dimt = 1
    end select

    ! select NetCDF format
    file_type = NF90_NETCDF4
    !file_type = NF90_CLASSIC_MODEL

    ! create output file
    call check( nf90_create( trim(outfname), file_type, ncid ))

    ! define dimensions
    call check( nf90_def_dim( ncid, "lon", dimx, lon_dim_id(1) ) )
    call check( nf90_def_dim( ncid, "lat", dimy, lat_dim_id(1) ) )

    ! define variables
    call define_var( ncid, var_id=lon_var_id, &
                     dim_id=lon_dim_id, var_name="lon", var_type=NF90_FLOAT, &
                     units="degrees_east", long_name="Longitude", &
                     standard_name="longitude", axis="X")

    call define_var( ncid, var_id=lat_var_id, &
                     dim_id=lat_dim_id, var_name="lat", var_type=NF90_FLOAT, &
                     units="degrees_north", long_name="Latitude", &
                     standard_name="latitude", axis="Y")

    dim_id_list = (/lon_dim_id, lat_dim_id/)

    var_name_list = (/ "n_obs2d", &
                       "pdf_clr",  "pdf_cld",  "pdf_drz",  "pdf_prc",  &
                       "lwp",      "iwp",      "cot",      "cth",      &
                       "cth_c4",   "ctp_c4"                                /)

    units_list = (/ "#", &
                    "%",     "%",     "%",     "%", &
                    "g m-2", "g m-2", "none",  "m", &
                    "m",     "hPa"                       /)

    long_name_list = (/ "# of Observation", &
                        "PDF of Clear-sky (no C4 cloudy bin)", &
                        "PDF of Non-Precipitating Clouds (-35 <= in-cloud Zmax < -15 dBZ)", &
                        "PDF of Drizzling Clouds (-15 <= in-cloud Zmax < 0 dBZ)", &
                        "PDF of Precipitating Clouds (0 <= in-cloud Zmax)", &
                        "LWP (cloud_water_path)", &
                        "IWP (cloud_ice_water_path)", &
                        "COT (optical_thickness)", &
                        "CPR echo-top height (reflectivity >= -35 dBZ)", &
                        "C4 Cloud-Top Height (CPR or ATLID)", &
                        "C4 Cloud-Top Pressure (CPR or ATLID)" /)

    do ivar = 1, nvar2d+2
       call define_var( ncid, var_id=var_id_list(ivar), &
                        dim_id=dim_id_list, var_name=var_name_list(ivar), var_type=NF90_FLOAT, &
                        units=units_list(ivar), long_name=long_name_list(ivar))
    enddo

    ! set global attributes
    call check( nf90_put_att(ncid,nf90_global,"title", &
                "EarthCARE Level3 (Global 2D Map)"))
    call put_global_metadata(ncid, aggregation_type, period_start, period_end)

    call check( nf90_enddef( ncid ) )

    ! write variables
    call check( nf90_put_var(ncid, lon_var_id, glon(2::2)) )
    call check( nf90_put_var(ncid, lat_var_id, glat(2::2)) )

    call check( nf90_put_var(ncid, var_id_list( 1), n_obs2d     ) )
    call check( nf90_put_var(ncid, var_id_list( 2), pdf_clr     ) )
    call check( nf90_put_var(ncid, var_id_list( 3), pdf_cld     ) )
    call check( nf90_put_var(ncid, var_id_list( 4), pdf_drz     ) )
    call check( nf90_put_var(ncid, var_id_list( 5), pdf_prc     ) )
    call check( nf90_put_var(ncid, var_id_list( 6), mean_lwp    ) )
    call check( nf90_put_var(ncid, var_id_list( 7), mean_iwp    ) )
    call check( nf90_put_var(ncid, var_id_list( 8), mean_cot    ) )
    call check( nf90_put_var(ncid, var_id_list( 9), mean_cth    ) )
    call check( nf90_put_var(ncid, var_id_list(10), mean_cth_c4 ) )
    call check( nf90_put_var(ncid, var_id_list(11), mean_ctp_c4 ) )

    ! close file
    call check( nf90_close(ncid) )

end subroutine output_netcdf_map2d

subroutine output_netcdf_map3d( &
    glon, glat, gplev, &
    n_obs3d, temp, shum, rhum, &
    cldavel, clddvel, cldtvel1, cldtvel2, cldtvelt, &
    cwc1, cwc2, ciwc, cer1, cer2, cerliq, cerice, &
    numc1, numc2, numliq, numice, &
    outfname, output_type, aggregation_type, period_start, period_end)

    implicit none
    real(sp), dimension(:), intent(in) :: glon
    real(sp), dimension(:), intent(in) :: glat
    real(sp), dimension(:), intent(in) :: gplev
    real(sp), dimension(:,:,:), intent(in) :: n_obs3d
    real(sp), dimension(:,:,:), intent(in) :: temp, shum, rhum
    real(sp), dimension(:,:,:), intent(in) :: cldavel, clddvel, cldtvel1, cldtvel2, cldtvelt
    real(sp), dimension(:,:,:), intent(in) :: cwc1, cwc2, ciwc, cer1, cer2, cerliq, cerice
    real(sp), dimension(:,:,:), intent(in) :: numc1, numc2, numliq, numice
    character(len=*), intent(in) :: outfname, output_type
    character(len=*), intent(in) :: aggregation_type, period_start, period_end

    integer, dimension(0:nmon) :: time

    integer :: ncid
    integer, dimension(1) :: lon_dim_id, lat_dim_id, plev_dim_id
    integer, dimension(nvar3d) :: var_id_list
    character(len=20), dimension(nvar3d) :: var_name_list, units_list
    character(len=128), dimension(nvar3d) :: long_name_list
    integer :: lon_var_id, lat_var_id, plev_var_id
    integer, dimension(3) :: dim_id_list
    integer :: file_type, mon_ini, mon_fin, dimt, ivar

    time = (/0:nmon/)
    select case(output_type)
    case("12_months")
        mon_ini = 1
        mon_fin = 12
        dimt = 12
    case("annual")
        mon_ini = 1
        mon_fin = 1
        dimt = 1
    case("monthly")
        mon_ini = 1
        mon_fin = 1
        dimt = 1
    end select

    ! select NetCDF format
    file_type = NF90_NETCDF4
    !file_type = NF90_CLASSIC_MODEL

    ! create output file
    call check( nf90_create( trim(outfname), file_type, ncid ))

    ! define dimensions
    call check( nf90_def_dim( ncid, "lon",  dimx, lon_dim_id(1) ) )
    call check( nf90_def_dim( ncid, "lat",  dimy, lat_dim_id(1) ) )
    call check( nf90_def_dim( ncid, "plev", dimz, plev_dim_id(1) ) )

    ! define variables
    call define_var( ncid, var_id=lon_var_id, &
                     dim_id=lon_dim_id, var_name="lon", var_type=NF90_FLOAT, &
                     units="degrees_east", long_name="Longitude", &
                     standard_name="longitude", axis="X")

    call define_var( ncid, var_id=lat_var_id, &
                     dim_id=lat_dim_id, var_name="lat", var_type=NF90_FLOAT, &
                     units="degrees_north", long_name="Latitude", &
                     standard_name="latitude", axis="Y")

    call define_var( ncid, var_id=plev_var_id, &
                     dim_id=plev_dim_id, var_name="plev", var_type=NF90_FLOAT, &
                     units="hPa", long_name="Pressure", &
                     standard_name="air_pressure", axis="Z", &
                     positive="down")

    dim_id_list = (/ lon_dim_id(1), lat_dim_id(1), plev_dim_id(1) /)

    var_name_list = (/ "n_obs3d", &
                       "temp", "shum", "rhum", &
                       "cldavel", "clddvel", "cldtvel1", "cldtvel2", "cldtvelt", &
                       "cwc1", "cwc2", "ciwc", "cer1", "cer2", "cerliq", "cerice", &
                       "numc1", "numc2", "numliq", "numice" /)

    units_list = (/ "#", &
                    "K", "kg kg -1", "%",   &
                    "m s-1", "m s-1", "m s-1",  "m s-1", "m s-1",  &
                    "g m-3", "g m-3", "g m-3",  "micron", "micron", "micron", "micron",  &
                    "# cm-3", "# cm-3", "# cm-3", "# cm-3" /)

    long_name_list = (/ "# of Observation", &
                        "Temperature", &
                        "Specific Humidity", &
                        "Relative Humidity", &
                        "Air Vertical Velocity", &
                        "Cloud Doppler Velocity", &
                        "Cloud Terminal Velocity 1", &
                        "Cloud Terminal Velocity 2", &
                        "Total Cloud Terminal Velocity", &
                        "Cloud Water Content 1", &
                        "Cloud Water Content 2", &
                        "Ice Water Content", &
                        "Cloud Effective Radius 1", &
                        "Cloud Effective Radius 2", &
                        "Cloud Effective Radius Liq", &
                        "Cloud Effective Radius Ice", &
                        "Cloud Number Concentration 1", &
                        "Cloud Number Concentration 2", &
                        "Cloud Number Concentration Liq", &
                        "Cloud Number Concentration Ice" /)

    do ivar = 1, nvar3d
       call define_var( ncid, var_id=var_id_list(ivar), &
                        dim_id=dim_id_list, var_name=var_name_list(ivar), var_type=NF90_FLOAT, &
                        units=units_list(ivar), long_name=long_name_list(ivar))
    enddo

    ! set global attributes
    call check( nf90_put_att(ncid,nf90_global,"title", &
                "EarthCARE Level3 (Global 3D Map)"))
    call put_global_metadata(ncid, aggregation_type, period_start, period_end)

    call check( nf90_enddef( ncid ) )

    ! write variables
    call check( nf90_put_var(ncid, lon_var_id,  glon(2::2)) )
    call check( nf90_put_var(ncid, lat_var_id,  glat(2::2)) )
    !call check( nf90_put_var(ncid, plev_var_id, gplev(2::2)) )
    call check( nf90_put_var(ncid, plev_var_id, gplev(2:dimz+1)) )

    call check( nf90_put_var(ncid, var_id_list( 1), n_obs3d     ) )
    call check( nf90_put_var(ncid, var_id_list( 2), temp        ) )
    call check( nf90_put_var(ncid, var_id_list( 3), shum        ) )
    call check( nf90_put_var(ncid, var_id_list( 4), rhum        ) )
    call check( nf90_put_var(ncid, var_id_list( 5), cldavel     ) )
    call check( nf90_put_var(ncid, var_id_list( 6), clddvel     ) )
    call check( nf90_put_var(ncid, var_id_list( 7), cldtvel1    ) )
    call check( nf90_put_var(ncid, var_id_list( 8), cldtvel2    ) )
    call check( nf90_put_var(ncid, var_id_list( 9), cldtvelt    ) )
    call check( nf90_put_var(ncid, var_id_list(10), cwc1        ) )
    call check( nf90_put_var(ncid, var_id_list(11), cwc2        ) )
    call check( nf90_put_var(ncid, var_id_list(12), ciwc        ) )
    call check( nf90_put_var(ncid, var_id_list(13), cer1        ) )
    call check( nf90_put_var(ncid, var_id_list(14), cer2        ) )
    call check( nf90_put_var(ncid, var_id_list(15), cerliq      ) )
    call check( nf90_put_var(ncid, var_id_list(16), cerice      ) )
    call check( nf90_put_var(ncid, var_id_list(17), numc1       ) )
    call check( nf90_put_var(ncid, var_id_list(18), numc2       ) )
    call check( nf90_put_var(ncid, var_id_list(19), numliq      ) )
    call check( nf90_put_var(ncid, var_id_list(20), numice      ) )

    ! close file
    call check( nf90_close(ncid) )

end subroutine output_netcdf_map3d

subroutine output_netcdf_cld3d( &
    glon, glat, gplev, &
    temp,              &
    n_obs3d,           &
    n_cld3d_0, n_cld3d_1, n_cld3d_2, n_cld3d_3, n_cld3d_4, &
    n_cld3d_5, n_cld3d_6, n_cld3d_7, n_cld3d_8, n_cld3d_9, &
    n_cld3d_10, n_cld3d_11, n_cld3d_12, n_cld3d_13, n_cld3d_14, &
    outfname, output_type, aggregation_type, period_start, period_end)

    implicit none
    real(sp), dimension(:), intent(in) :: glon
    real(sp), dimension(:), intent(in) :: glat
    real(sp), dimension(:), intent(in) :: gplev
    real(sp), dimension(:,:,:), intent(in) :: temp
    real(sp), dimension(:,:,:), intent(in) :: n_obs3d
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_0
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_1
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_2
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_3
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_4
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_5
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_6
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_7
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_8
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_9
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_10
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_11
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_12
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_13
    real(sp), dimension(:,:,:), intent(in) :: n_cld3d_14
    character(len=*), intent(in) :: outfname, output_type
    character(len=*), intent(in) :: aggregation_type, period_start, period_end

    integer, dimension(0:nmon) :: time

    integer :: ncid
    integer, dimension(1) :: lon_dim_id, lat_dim_id, plev_dim_id
    integer, dimension(nvar3dcld) :: var_id_list
    character(len=20), dimension(nvar3dcld) :: var_name_list, units_list
    character(len=128), dimension(nvar3dcld) :: long_name_list
    integer :: lon_var_id, lat_var_id, plev_var_id
    integer, dimension(3) :: dim_id_list
    integer :: file_type, mon_ini, mon_fin, dimt, ivar

    time = (/0:nmon/)
    select case(output_type)
    case("12_months")
        mon_ini = 1
        mon_fin = 12
        dimt = 12
    case("annual")
        mon_ini = 1
        mon_fin = 1
        dimt = 1
    case("monthly")
        mon_ini = 1
        mon_fin = 1
        dimt = 1
    end select

    ! select NetCDF format
    file_type = NF90_NETCDF4
    !file_type = NF90_CLASSIC_MODEL

    ! create output file
    call check( nf90_create( trim(outfname), file_type, ncid ))

    ! define dimensions
    call check( nf90_def_dim( ncid, "lon",  dimx, lon_dim_id(1) ) )
    call check( nf90_def_dim( ncid, "lat",  dimy, lat_dim_id(1) ) )
    call check( nf90_def_dim( ncid, "plev", dimz, plev_dim_id(1) ) )

    ! define variables
    call define_var( ncid, var_id=lon_var_id, &
                     dim_id=lon_dim_id, var_name="lon", var_type=NF90_FLOAT, &
                     units="degrees_east", long_name="Longitude", &
                     standard_name="longitude", axis="X")

    call define_var( ncid, var_id=lat_var_id, &
                     dim_id=lat_dim_id, var_name="lat", var_type=NF90_FLOAT, &
                     units="degrees_north", long_name="Latitude", &
                     standard_name="latitude", axis="Y")

    call define_var( ncid, var_id=plev_var_id, &
                     dim_id=plev_dim_id, var_name="plev", var_type=NF90_FLOAT, &
                     units="hPa", long_name="Pressure", &
                     standard_name="air_pressure", axis="Z", &
                     positive="down")

    dim_id_list = (/ lon_dim_id(1), lat_dim_id(1), plev_dim_id(1) /)

    var_name_list = (/ "temp",       "n_obs3d",                  &
                       "n_cld3d_0",  "n_cld3d_1",  "n_cld3d_2",  &
                       "n_cld3d_3",  "n_cld3d_4",  "n_cld3d_5",  &
                       "n_cld3d_6",  "n_cld3d_7",  "n_cld3d_8",  &
                       "n_cld3d_9",  "n_cld3d_10", "n_cld3d_11", &
                       "n_cld3d_12", "n_cld3d_13", "n_cld3d_14"  /)

    units_list = (/ "K", "#",      &
                    "#", "#", "#", &
                    "#", "#", "#", &
                    "#", "#", "#", &
                    "#", "#", "#", &
                    "#", "#", "#"  /)

    long_name_list = (/ "Temperature",              &
                        "# of Observation",         &
                        "# of cloudy-sky",          &  ! NOTE: not "clear-sky" here
                        "# of warm water",          &
                        "# of supercooled water",   &
                        "# of 3D ice",              &
                        "# of 2D plate ice",        &
                        "# of 3D+2D ice mixture",   &
                        "# of liquid drizzle",      &
                        "# of mixed-phase drizzle", &
                        "# of rain",                &
                        "# of snow",                &
                        "# of water and drizzle",   &
                        "# of water and rain",      &
                        "# of mixed phase",         &
                        "# of unknown",             &
                        "# of melting layer" /)

    do ivar = 1, nvar3dcld
       call define_var( ncid, var_id=var_id_list(ivar), &
                        dim_id=dim_id_list, var_name=var_name_list(ivar), var_type=NF90_FLOAT, &
                        units=units_list(ivar), long_name=long_name_list(ivar))
    enddo

    ! set global attributes
    call check( nf90_put_att(ncid,nf90_global,"title", &
                "EarthCARE Level3 (Global 3D Cloud Fraction)"))
    call put_global_metadata(ncid, aggregation_type, period_start, period_end)

    call check( nf90_enddef( ncid ) )

    ! write variables
    call check( nf90_put_var(ncid, lon_var_id,  glon(2::2)) )
    call check( nf90_put_var(ncid, lat_var_id,  glat(2::2)) )
    !call check( nf90_put_var(ncid, plev_var_id, gplev(2::2)) )
    call check( nf90_put_var(ncid, plev_var_id, gplev(2:dimz+1)) )

    call check( nf90_put_var(ncid, var_id_list( 1), temp        ) )
    call check( nf90_put_var(ncid, var_id_list( 2), n_obs3d     ) )
    call check( nf90_put_var(ncid, var_id_list( 3), n_cld3d_0   ) )
    call check( nf90_put_var(ncid, var_id_list( 4), n_cld3d_1   ) )
    call check( nf90_put_var(ncid, var_id_list( 5), n_cld3d_2   ) )
    call check( nf90_put_var(ncid, var_id_list( 6), n_cld3d_3   ) )
    call check( nf90_put_var(ncid, var_id_list( 7), n_cld3d_4   ) )
    call check( nf90_put_var(ncid, var_id_list( 8), n_cld3d_5   ) )
    call check( nf90_put_var(ncid, var_id_list( 9), n_cld3d_6   ) )
    call check( nf90_put_var(ncid, var_id_list(10), n_cld3d_7   ) )
    call check( nf90_put_var(ncid, var_id_list(11), n_cld3d_8   ) )
    call check( nf90_put_var(ncid, var_id_list(12), n_cld3d_9   ) )
    call check( nf90_put_var(ncid, var_id_list(13), n_cld3d_10  ) )
    call check( nf90_put_var(ncid, var_id_list(14), n_cld3d_11  ) )
    call check( nf90_put_var(ncid, var_id_list(15), n_cld3d_12  ) )
    call check( nf90_put_var(ncid, var_id_list(16), n_cld3d_13  ) )
    call check( nf90_put_var(ncid, var_id_list(17), n_cld3d_14  ) )

    ! close file
    call check( nf90_close(ncid) )

end subroutine output_netcdf_cld3d

subroutine output_netcdf_c4_3d( &
    glon, glat, gplev, &
    n_valid, n_cloud, n_mask_invalid, n_type, n_category, &
    n_type_qc, n_category_invalid, outfname, &
    aggregation_type, period_start, period_end )
    implicit none
    real(sp), intent(in) :: glon(:), glat(:), gplev(:)
    integer(kind=4), intent(in) :: n_valid(dimx,dimy,dimz)
    integer(kind=4), intent(in) :: n_cloud(dimx,dimy,dimz)
    integer(kind=4), intent(in) :: n_mask_invalid(dimx,dimy,dimz)
    integer(kind=4), intent(in) :: n_type(dimx,dimy,dimz,ntype_c4)
    integer(kind=4), intent(in) :: n_category(dimx,dimy,dimz,ncategory_c4)
    integer(kind=4), intent(in) :: n_type_qc(dimx,dimy,dimz,ntype_qc_c4)
    integer(kind=4), intent(in) :: n_category_invalid(dimx,dimy,dimz)
    character(len=*), intent(in) :: outfname
    character(len=*), intent(in) :: aggregation_type, period_start, period_end

    integer :: ncid, k
    integer :: lon_dim_id, lat_dim_id, plev_dim_id
    integer :: type_dim_id, category_dim_id, type_qc_dim_id
    integer :: lon_var_id, lat_var_id, plev_var_id
    integer :: type_coord_id, category_coord_id, type_qc_coord_id
    integer :: valid_id, cloud_id, clear_id, mask_invalid_id, cloud_fraction_id
    integer :: type_id, type_assigned_id, type_unassigned_id
    integer :: type_occurrence_id, type_relative_id
    integer :: category_id, category_assigned_id, category_unassigned_id
    integer :: category_occurrence_id, category_relative_id
    integer :: type_qc_id, category_invalid_id
    integer, dimension(3) :: dim3
    integer, dimension(4) :: type_dim4, category_dim4, type_qc_dim4
    integer(kind=4) :: type_codes(ntype_c4), category_codes(ncategory_c4)
    integer(kind=4) :: type_qc_codes(ntype_qc_c4)
    integer(kind=4), allocatable :: derived(:,:,:)
    integer(kind=4), allocatable :: type_assigned(:,:,:)
    integer(kind=4), allocatable :: category_assigned(:,:,:)
    real(sp), allocatable :: frequency_slab(:,:,:,:)

    type_codes = (/ (k, k=1,ntype_c4) /)
    category_codes = (/ (k, k=1,ncategory_c4) /)
    type_qc_codes = (/ 15_4, 16_4, 17_4, -999_4 /)

    allocate( derived(dimx,dimy,dimz) )
    allocate( type_assigned(dimx,dimy,dimz) )
    allocate( category_assigned(dimx,dimy,dimz) )
    allocate( frequency_slab(dimx,dimy,dimz,1) )

    type_assigned = sum(n_type, DIM=4)
    category_assigned = sum(n_category, DIM=4)
    if ( any(n_cloud .gt. n_valid) .or. &
       & any(type_assigned .gt. n_cloud) .or. &
       & any(category_assigned .gt. n_cloud) ) then
       write(error_unit,*) "### Invalid C4 count relationship before NetCDF output"
       stop 1
    endif

    call check( nf90_create(trim(outfname), NF90_NETCDF4, ncid) )
    call check( nf90_def_dim(ncid, "lon", dimx, lon_dim_id) )
    call check( nf90_def_dim(ncid, "lat", dimy, lat_dim_id) )
    call check( nf90_def_dim(ncid, "plev", dimz, plev_dim_id) )
    call check( nf90_def_dim(ncid, "particle_type", ntype_c4, type_dim_id) )
    call check( nf90_def_dim(ncid, "particle_category", ncategory_c4, category_dim_id) )
    call check( nf90_def_dim(ncid, "type_qc_class", ntype_qc_c4, type_qc_dim_id) )

    call define_var(ncid, lon_var_id, (/lon_dim_id/), "lon", NF90_FLOAT, &
                    "degrees_east", "Longitude", "longitude", "X")
    call define_var(ncid, lat_var_id, (/lat_dim_id/), "lat", NF90_FLOAT, &
                    "degrees_north", "Latitude", "latitude", "Y")
    call define_var(ncid, plev_var_id, (/plev_dim_id/), "plev", NF90_FLOAT, &
                    "hPa", "Pressure", "air_pressure", "Z", "down")

    call check( nf90_def_var(ncid, "particle_type", NF90_INT, &
                             (/type_dim_id/), type_coord_id) )
    call check( nf90_put_att(ncid, type_coord_id, "units", "1") )
    call check( nf90_put_att(ncid, type_coord_id, "long_name", &
               "JAXA AC_CLP particle type code") )
    call check( nf90_put_att(ncid, type_coord_id, "flag_values", type_codes) )
    call check( nf90_put_att(ncid, type_coord_id, "flag_meanings", &
               "warm_water supercooled_water 3D_ice 2D_plate_ice " // &
               "3D_and_2D_ice_mixture liquid_drizzle mixed_phase_drizzle rain snow " // &
               "water_and_drizzle water_and_rain mixed_phase unknown melting_layer") )

    call check( nf90_def_var(ncid, "particle_category", NF90_INT, &
                             (/category_dim_id/), category_coord_id) )
    call check( nf90_put_att(ncid, category_coord_id, "units", "1") )
    call check( nf90_put_att(ncid, category_coord_id, "long_name", &
               "JAXA AC_CLP particle category code") )
    call check( nf90_put_att(ncid, category_coord_id, "flag_values", category_codes) )
    call check( nf90_put_att(ncid, category_coord_id, "flag_meanings", &
               "2D_plates 2D_columns bullet_rosette_or_3D_types droxtals " // &
               "2D_and_3D_Voronois supercooled_water warm_water liquid_drizzle rain " // &
               "water_and_liquid_drizzle water_and_rain unknown") )

    call check( nf90_def_var(ncid, "type_qc_class", NF90_INT, &
                             (/type_qc_dim_id/), type_qc_coord_id) )
    call check( nf90_put_att(ncid, type_qc_coord_id, "units", "1") )
    call check( nf90_put_att(ncid, type_qc_coord_id, "long_name", &
               "Particle type Core QC class") )
    call check( nf90_put_att(ncid, type_qc_coord_id, "flag_values", type_qc_codes) )
    call check( nf90_put_att(ncid, type_qc_coord_id, "flag_meanings", &
               "non_cloud_echo_1 fully_attenuated non_cloud_echo_2 other_unsupported") )

    dim3 = (/lon_dim_id, lat_dim_id, plev_dim_id/)
    type_dim4 = (/lon_dim_id, lat_dim_id, plev_dim_id, type_dim_id/)
    category_dim4 = (/lon_dim_id, lat_dim_id, plev_dim_id, category_dim_id/)
    type_qc_dim4 = (/lon_dim_id, lat_dim_id, plev_dim_id, type_qc_dim_id/)

    call define_c4_data_var(ncid, valid_id, dim3, "n_c4_valid", NF90_INT, "1", &
         "Number of native bins with C4 mask 0 or 1", "not_applicable", &
         "valid_geolocation_finite_height_and_0_lt_pressure_lt_1100_hPa", &
         "C4 mask equals 0 or 1")
    call define_c4_data_var(ncid, cloud_id, dim3, "n_c4_cloud", NF90_INT, "1", &
         "Number of native bins with C4 cloud mask", "not_applicable", &
         "valid_geolocation_finite_height_and_0_lt_pressure_lt_1100_hPa", &
         "C4 mask equals 1")
    call define_c4_data_var(ncid, clear_id, dim3, "n_c4_clear", NF90_INT, "1", &
         "Number of native bins with C4 clear mask", "not_applicable", &
         "valid_geolocation_finite_height_and_0_lt_pressure_lt_1100_hPa", &
         "C4 mask equals 0; derived as n_c4_valid minus n_c4_cloud")
    call define_c4_data_var(ncid, mask_invalid_id, dim3, "n_c4_mask_invalid", &
         NF90_INT, "1", "Number of native bins with invalid C4 mask", &
         "not_applicable", &
         "valid_geolocation_finite_height_and_0_lt_pressure_lt_1100_hPa", &
         "C4 mask is neither 0 nor 1")
    call define_c4_data_var(ncid, cloud_fraction_id, dim3, "c4_cloud_fraction", &
         NF90_FLOAT, "%", "Native-bin-based C4 cloud fraction", "n_c4_valid", &
         "valid_geolocation_finite_height_pressure_and_C4_mask_0_or_1", &
         "100 times n_c4_cloud divided by n_c4_valid; no CPR dBZ threshold")

    call define_c4_data_var(ncid, type_id, type_dim4, "n_c4_type", NF90_INT, "1", &
         "Number of C4 cloud native bins by particle type", "not_applicable", &
         "valid_geolocation_height_pressure_C4_mask_1_and_Type_1_to_14", &
         "C4 mask equals 1 and particle type equals coordinate code")
    call define_c4_data_var(ncid, type_assigned_id, dim3, "n_c4_type_assigned", &
         NF90_INT, "1", "Number of C4 cloud bins assigned Type 1 to 14", &
         "not_applicable", "valid_geolocation_height_pressure_and_C4_mask_1", &
         "sum of n_c4_type over particle_type")
    call define_c4_data_var(ncid, type_unassigned_id, dim3, "n_c4_type_unassigned", &
         NF90_INT, "1", "Number of C4 cloud bins without Type 1 to 14", &
         "not_applicable", "valid_geolocation_height_pressure_and_C4_mask_1", &
         "n_c4_cloud minus n_c4_type_assigned")
    call define_c4_data_var(ncid, type_occurrence_id, type_dim4, &
         "c4_type_occurrence_frequency", NF90_FLOAT, "%", &
         "C4 particle type occurrence frequency", "n_c4_valid", &
         "valid_geolocation_height_pressure_and_C4_mask_0_or_1", &
         "100 times n_c4_type divided by n_c4_valid")
    call define_c4_data_var(ncid, type_relative_id, type_dim4, &
         "c4_type_relative_frequency", NF90_FLOAT, "%", &
         "C4 particle type relative frequency among assigned types", &
         "n_c4_type_assigned", "C4_cloud_bins_with_Type_1_to_14", &
         "100 times n_c4_type divided by n_c4_type_assigned")

    call define_c4_data_var(ncid, category_id, category_dim4, "n_c4_category", &
         NF90_INT, "1", "Number of C4 cloud native bins by particle category", &
         "not_applicable", &
         "valid_geolocation_height_pressure_C4_mask_1_and_Category_1_to_12", &
         "C4 mask equals 1 and particle category equals coordinate code")
    call define_c4_data_var(ncid, category_assigned_id, dim3, &
         "n_c4_category_assigned", NF90_INT, "1", &
         "Number of C4 cloud bins assigned Category 1 to 12", "not_applicable", &
         "valid_geolocation_height_pressure_and_C4_mask_1", &
         "sum of n_c4_category over particle_category")
    call define_c4_data_var(ncid, category_unassigned_id, dim3, &
         "n_c4_category_unassigned", NF90_INT, "1", &
         "Number of C4 cloud bins without Category 1 to 12", "not_applicable", &
         "valid_geolocation_height_pressure_and_C4_mask_1", &
         "n_c4_cloud minus n_c4_category_assigned")
    call define_c4_data_var(ncid, category_occurrence_id, category_dim4, &
         "c4_category_occurrence_frequency", NF90_FLOAT, "%", &
         "C4 particle category occurrence frequency", "n_c4_valid", &
         "valid_geolocation_height_pressure_and_C4_mask_0_or_1", &
         "100 times n_c4_category divided by n_c4_valid")
    call define_c4_data_var(ncid, category_relative_id, category_dim4, &
         "c4_category_relative_frequency", NF90_FLOAT, "%", &
         "C4 particle category relative frequency among assigned categories", &
         "n_c4_category_assigned", "C4_cloud_bins_with_Category_1_to_12", &
         "100 times n_c4_category divided by n_c4_category_assigned")

    call define_c4_data_var(ncid, type_qc_id, type_qc_dim4, "n_c4_type_qc", &
         NF90_INT, "1", "Particle Type Core QC native-bin count", &
         "not_applicable", &
         "all_valid_geolocation_finite_height_and_0_lt_pressure_lt_1100_hPa_bins", &
         "Type 15, 16, 17, or other unsupported raw Type; classes are exclusive")
    call define_c4_data_var(ncid, category_invalid_id, dim3, &
         "n_c4_category_invalid", NF90_INT, "1", &
         "Particle Category invalid or unassigned native-bin count", "not_applicable", &
         "all_valid_geolocation_finite_height_and_0_lt_pressure_lt_1100_hPa_bins", &
         "Category outside 0 to 12 after input conversion; raw missing classes are merged")

    call check( nf90_put_att(ncid, nf90_global, "title", &
               "EarthCARE Level3 C4 3D Classification Diagnostics") )
    call check( nf90_put_att(ncid, nf90_global, "reader_version", "1.0") )
    call check( nf90_put_att(ncid, nf90_global, "input_product", &
               "JAXA EarthCARE 2B AC_CLP") )
    call check( nf90_put_att(ncid, nf90_global, "input_product_version", "vCb") )
    call check( nf90_put_att(ncid, nf90_global, "aggregation_type", &
               trim(aggregation_type)) )
    call check( nf90_put_att(ncid, nf90_global, "processing_period_start_month", &
               trim(period_start)) )
    call check( nf90_put_att(ncid, nf90_global, "processing_period_end_month", &
               trim(period_end)) )
    call check( nf90_put_att(ncid, nf90_global, "cloud_mask_definition", &
               "cloud_mask_cpr_atlid_10km: 0=clear, 1=cloud, other=invalid") )
    call check( nf90_put_att(ncid, nf90_global, "sampling_method", &
               "native vertical bins on the 40 pressure levels") )
    call check( nf90_put_att(ncid, nf90_global, "reflectivity_selection", &
               "No CPR dBZ threshold is applied to C4 diagnostics") )
    call check( nf90_put_att(ncid, nf90_global, "geolocation_selection", &
               "finite in-range latitude/longitude; zero coordinates are not excluded") )
    call check( nf90_enddef(ncid) )

    call check( nf90_put_var(ncid, lon_var_id, glon(2::2)) )
    call check( nf90_put_var(ncid, lat_var_id, glat(2::2)) )
    call check( nf90_put_var(ncid, plev_var_id, gplev(2:dimz+1)) )
    call check( nf90_put_var(ncid, type_coord_id, type_codes) )
    call check( nf90_put_var(ncid, category_coord_id, category_codes) )
    call check( nf90_put_var(ncid, type_qc_coord_id, type_qc_codes) )

    call check( nf90_put_var(ncid, valid_id, n_valid) )
    call check( nf90_put_var(ncid, cloud_id, n_cloud) )
    derived = n_valid - n_cloud
    call check( nf90_put_var(ncid, clear_id, derived) )
    call check( nf90_put_var(ncid, mask_invalid_id, n_mask_invalid) )
    frequency_slab(:,:,:,1) = undef
    where (n_valid .gt. 0)
       frequency_slab(:,:,:,1) = 100.0_sp * real(n_cloud,sp) / real(n_valid,sp)
    end where
    call check( nf90_put_var(ncid, cloud_fraction_id, frequency_slab(:,:,:,1)) )

    call check( nf90_put_var(ncid, type_id, n_type) )
    call check( nf90_put_var(ncid, type_assigned_id, type_assigned) )
    derived = n_cloud - type_assigned
    call check( nf90_put_var(ncid, type_unassigned_id, derived) )
    do k = 1, ntype_c4
       frequency_slab(:,:,:,1) = undef
       where (n_valid .gt. 0)
          frequency_slab(:,:,:,1) = 100.0_sp * real(n_type(:,:,:,k),sp) / &
                                    real(n_valid,sp)
       end where
       call check( nf90_put_var(ncid, type_occurrence_id, frequency_slab, &
                   start=(/1,1,1,k/), count=(/dimx,dimy,dimz,1/)) )
       frequency_slab(:,:,:,1) = undef
       where (type_assigned .gt. 0)
          frequency_slab(:,:,:,1) = 100.0_sp * real(n_type(:,:,:,k),sp) / &
                                    real(type_assigned,sp)
       end where
       call check( nf90_put_var(ncid, type_relative_id, frequency_slab, &
                   start=(/1,1,1,k/), count=(/dimx,dimy,dimz,1/)) )
    enddo

    call check( nf90_put_var(ncid, category_id, n_category) )
    call check( nf90_put_var(ncid, category_assigned_id, category_assigned) )
    derived = n_cloud - category_assigned
    call check( nf90_put_var(ncid, category_unassigned_id, derived) )
    do k = 1, ncategory_c4
       frequency_slab(:,:,:,1) = undef
       where (n_valid .gt. 0)
          frequency_slab(:,:,:,1) = 100.0_sp * real(n_category(:,:,:,k),sp) / &
                                    real(n_valid,sp)
       end where
       call check( nf90_put_var(ncid, category_occurrence_id, frequency_slab, &
                   start=(/1,1,1,k/), count=(/dimx,dimy,dimz,1/)) )
       frequency_slab(:,:,:,1) = undef
       where (category_assigned .gt. 0)
          frequency_slab(:,:,:,1) = 100.0_sp * real(n_category(:,:,:,k),sp) / &
                                    real(category_assigned,sp)
       end where
       call check( nf90_put_var(ncid, category_relative_id, frequency_slab, &
                   start=(/1,1,1,k/), count=(/dimx,dimy,dimz,1/)) )
    enddo

    call check( nf90_put_var(ncid, type_qc_id, n_type_qc) )
    call check( nf90_put_var(ncid, category_invalid_id, n_category_invalid) )
    call check( nf90_close(ncid) )

    deallocate(derived, type_assigned, category_assigned, frequency_slab)
end subroutine output_netcdf_c4_3d

subroutine define_c4_data_var(ncid, var_id, dim_id, var_name, var_type, &
                              units, long_name, denominator, population, criteria)
    implicit none
    integer, intent(in) :: ncid, var_type
    integer, intent(out) :: var_id
    integer, intent(in) :: dim_id(:)
    character(len=*), intent(in) :: var_name, units, long_name
    character(len=*), intent(in) :: denominator, population, criteria
    integer, allocatable :: chunks(:)

    call check( nf90_def_var(ncid, var_name, var_type, dim_id, var_id) )
    call check( nf90_put_att(ncid, var_id, "units", units) )
    call check( nf90_put_att(ncid, var_id, "long_name", long_name) )
    call check( nf90_put_att(ncid, var_id, "denominator", denominator) )
    call check( nf90_put_att(ncid, var_id, "sampling_population", population) )
    call check( nf90_put_att(ncid, var_id, "cloud_mask_definition", &
               "0=clear, 1=cloud, other=invalid") )
    call check( nf90_put_att(ncid, var_id, "relevant_selection_criteria", criteria) )
    if (var_type .eq. NF90_INT) then
       call check( nf90_def_var_fill(ncid, var_id, 0, undef_i4) )
    else
       call check( nf90_def_var_fill(ncid, var_id, 0, undef) )
    endif

    allocate(chunks(size(dim_id)))
    if (size(dim_id) .eq. 3) then
       chunks = (/36,18,10/)
    elseif (size(dim_id) .eq. 4) then
       chunks = (/36,18,10,1/)
    else
       write(error_unit,*) "### Unsupported C4 NetCDF rank: ", size(dim_id)
       stop 1
    endif
    call check( nf90_def_var_chunking(ncid, var_id, NF90_CHUNKED, chunks) )
    call check( nf90_def_var_deflate(ncid, var_id, 1, 1, 1) )
    deallocate(chunks)
end subroutine define_c4_data_var

subroutine put_global_metadata(ncid, aggregation_type, period_start, period_end)
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: aggregation_type, period_start, period_end

    call check( nf90_put_att(ncid, nf90_global, "reader_version", "1.0") )
    call check( nf90_put_att(ncid, nf90_global, "input_product", "JAXA EarthCARE 2B AC_CLP") )
    call check( nf90_put_att(ncid, nf90_global, "input_product_version", "vCb") )
    call check( nf90_put_att(ncid, nf90_global, "aggregation_type", trim(aggregation_type)) )
    call check( nf90_put_att(ncid, nf90_global, "processing_period_start_month", trim(period_start)) )
    call check( nf90_put_att(ncid, nf90_global, "processing_period_end_month", trim(period_end)) )
end subroutine put_global_metadata

! define netcdf variables and put attributions
subroutine define_var(&
    ncid, var_id, &
    dim_id, var_name, var_type, &
    units, long_name, &
    standard_name, axis, positive)

    implicit none
    integer, intent(out) :: var_id
    integer, intent(in) :: ncid
    integer, dimension(:), intent(in) :: dim_id
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: var_name, units, long_name
    character(len=*), intent(in), optional :: standard_name, axis, positive

    call check( nf90_def_var( ncid, var_name, var_type, dim_id, var_id) )

    ! put attribute
    call check( nf90_put_att( ncid, var_id, "units", units ))
    call check( nf90_put_att( ncid, var_id, "long_name", long_name ))

    if ( present( standard_name ) ) then
        call check( nf90_put_att( ncid, var_id, "standard_name", standard_name ))
    endif

    if ( present( axis ) ) then
        call check( nf90_put_att( ncid, var_id, "axis", axis ))
    endif
    if ( present( positive ) ) then
        call check( nf90_put_att( ncid, var_id, "positive", positive ))
    endif

    ! specify undefined values
    call check( nf90_def_var_fill( ncid, var_id, 0, undef ))

end subroutine define_var

! check operations
subroutine check(status)
    implicit none
    integer, intent(in) :: status

    if(status /= nf90_noerr) then
       write(error_unit,*) trim(nf90_strerror(status))
       stop "NetCDF Error"
    end if
end subroutine check

end module netcdf_writer
