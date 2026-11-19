module summary_processor
    use parameters
    use numbers, only: n_c4_valid, n_c4_cloud, n_c4_mask_invalid, &
                       n_c4_type, n_c4_category, n_c4_type_qc, &
                       n_c4_category_invalid, n_c4_valid_all, n_c4_cloud_all, &
                       n_c4_mask_invalid_all, n_c4_type_all, n_c4_category_all, &
                       n_c4_type_qc_all, n_c4_category_invalid_all
    use netcdf_writer
    use pdfs

contains

subroutine output_summary( &
    ooutyr, ooutmon, ofinal, &
    ym, year, imon, &
    annual_start_month, processing_start_month, &
    glon, glat, gplev, gdbz, glev, &
    n_obs2d, n_obs2d_all, &
    n_obs3d, n_obs3d_all, &
    n_cld3d, n_cld3d_all  )

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
    real(sp), dimension(:, :, 0:), intent(in) :: n_obs2d
    real(sp), dimension(:, :), intent(in) :: n_obs2d_all
    real(sp), dimension(:, :, :, 0:), intent(in) :: n_obs3d
    real(sp), dimension(:, :, :), intent(in) :: n_obs3d_all
    real(sp), dimension(:, :, :, 0:, 0:), intent(in) :: n_cld3d
    real(sp), dimension(:, :, :, 0:), intent(in) :: n_cld3d_all
    character(len = 128) :: outfname
    character(len = 7) :: current_month

    current_month = ym(1:4) // "-" // ym(5:6)

    !! Monthly Output
    if ( ooutmon ) then

        outfname = trim(outpath_lev3) // "/map2d/ecare_map2d_" // trim(ym) // ".nc"
        call output_netcdf_map2d( &
            glon,               glat,               n_obs2d(:,:,imon),    & ! 01-03
            pdf_clr(:,:,imon),  pdf_cld(:,:,imon),  pdf_drz(:,:,imon),    & ! 04-06
            pdf_prc(:,:,imon),  mean_lwp(:,:,imon), mean_iwp(:,:,imon),   & ! 07-09
            mean_cot(:,:,imon), mean_cth(:,:,imon),                       & ! 10-11
            mean_cth_c4(:,:,imon), mean_ctp_c4(:,:,imon),                 & ! 12-13
            outfname, output_type='monthly', &
            aggregation_type='monthly', period_start=current_month, period_end=current_month )

        outfname = trim(outpath_lev3) // "/map3d/ecare_map3d_" // trim(ym) // ".nc"
        call output_netcdf_map3d( &
             glon,                      glat,                      gplev,                     & ! 01-03
             n_obs3d(:,:,:,imon),                                                             & ! 04
             mean_temp(:,:,:,imon),     mean_shum(:,:,:,imon),     mean_rhum(:,:,:,imon),     & ! 05-07
             mean_cldavel(:,:,:,imon),  mean_clddvel(:,:,:,imon),                             & ! 08-09
             mean_cldtvel1(:,:,:,imon), mean_cldtvel2(:,:,:,imon), mean_cldtvelt(:,:,:,imon), & ! 10-12
             mean_cwc1(:,:,:,imon),     mean_cwc2(:,:,:,imon),     mean_ciwc(:,:,:,imon),     & ! 13-15
             mean_cer1(:,:,:,imon),     mean_cer2(:,:,:,imon),                                & ! 16-17
             mean_cerliq(:,:,:,imon),   mean_cerice(:,:,:,imon),                              & ! 18-19
             mean_numc1(:,:,:,imon),    mean_numc2(:,:,:,imon),                               & ! 20-21
             mean_numliq(:,:,:,imon),   mean_numice(:,:,:,imon),                              & ! 22-23
             outfname, output_type='monthly', &
             aggregation_type='monthly', period_start=current_month, period_end=current_month )

        outfname = trim(outpath_lev3) // "/map3d/ecare_cld3d_" // trim(ym) // ".nc"
        call output_netcdf_cld3d( &
             glon,                      glat,                      gplev,                     & ! 01-03
             mean_temp(:,:,:,imon),                                                           & ! 04
             n_obs3d(:,:,:,imon),                                                             & ! 05
             n_cld3d(:,:,:, 0,imon),    n_cld3d(:,:,:, 1,imon),    n_cld3d(:,:,:, 2,imon),    & ! 06-08
             n_cld3d(:,:,:, 3,imon),    n_cld3d(:,:,:, 4,imon),    n_cld3d(:,:,:, 5,imon),    & ! 09-11
             n_cld3d(:,:,:, 6,imon),    n_cld3d(:,:,:, 7,imon),    n_cld3d(:,:,:, 8,imon),    & ! 12-14
             n_cld3d(:,:,:, 9,imon),    n_cld3d(:,:,:,10,imon),    n_cld3d(:,:,:,11,imon),    & ! 15-17
             n_cld3d(:,:,:,12,imon),    n_cld3d(:,:,:,13,imon),    n_cld3d(:,:,:,14,imon),    & ! 18-20
             outfname, output_type='monthly', &
             aggregation_type='monthly', period_start=current_month, period_end=current_month )

        outfname = trim(outpath_lev3) // "/cfad/ecare_cfad_" // trim(ym) // ".nc"
        call output_netcdf_cfad( outfname, gdbz, glev, &
                                 pdf_cfad( 1:dimdbz, 1:dimlev, 1:nreg, imon ), &
                                 'monthly', current_month, current_month )

        outfname = trim(outpath_lev3) // "/isccp/ecare_isccp_" // trim(ym) // ".nc"
        !call output_netcdf_isccp( outfname, gcot, gctp, &
        call output_netcdf_isccp( outfname, &
                                  isccp_fq( 1:dimcot, 1:dimctp, 1:nreg, imon ), &
                                  isccp_dv( 1:dimcot, 1:dimctp, 1:nreg, imon ), &
                                  isccp_fq_c4( 1:dimcot, 1:dimctp, 1:nreg, imon ), &
                                  isccp_dv_c4( 1:dimcot, 1:dimctp, 1:nreg, imon ), &
                                  'monthly', current_month, current_month  )

        outfname = trim(outpath_lev3) // "/map3d/ecare_c4_3d_" // trim(ym) // ".nc"
        call output_netcdf_c4_3d( glon, glat, gplev, &
             n_c4_valid(:,:,:,imon), n_c4_cloud(:,:,:,imon), &
             n_c4_mask_invalid(:,:,:,imon), n_c4_type(:,:,:,:,imon), &
             n_c4_category(:,:,:,:,imon), n_c4_type_qc(:,:,:,:,imon), &
             n_c4_category_invalid(:,:,:,imon), outfname, &
             'monthly', current_month, current_month )

    endif

    !! Annual Output
    if ( ooutyr ) then

        outfname = trim(outpath_lev3) // "/map2d/ecare_map2d_" // trim(year) // ".nc"
        call output_netcdf_map2d( &
            glon,               glat,               n_obs2d(:,:,0),    & ! 01-03
            pdf_clr(:,:,0),     pdf_cld(:,:,0),     pdf_drz(:,:,0),    & ! 04-06
            pdf_prc(:,:,0),     mean_lwp(:,:,0),    mean_iwp(:,:,0),   & ! 07-09
            mean_cot(:,:,0),    mean_cth(:,:,0),                       & ! 10-11
            mean_cth_c4(:,:,0), mean_ctp_c4(:,:,0),                    & ! 12-13
            outfname,           output_type='annual', &
            aggregation_type='annual', period_start=annual_start_month, period_end=current_month )

        outfname = trim(outpath_lev3) // "/map3d/ecare_map3d_" // trim(year) // ".nc"
        call output_netcdf_map3d( &
             glon,                   glat,                   gplev,                   & ! 01-03
             n_obs3d(:,:,:,0),                                                        & ! 04
             mean_temp(:,:,:,0),     mean_shum(:,:,:,0),     mean_rhum(:,:,:,0),      & ! 05-07
             mean_cldavel(:,:,:,0),  mean_clddvel(:,:,:,0),                           & ! 08-09
             mean_cldtvel1(:,:,:,0), mean_cldtvel2(:,:,:,0), mean_cldtvelt(:,:,:,0),  & ! 10-12
             mean_cwc1(:,:,:,0),     mean_cwc2(:,:,:,0),     mean_ciwc(:,:,:,0),      & ! 13-15
             mean_cer1(:,:,:,0),     mean_cer2(:,:,:,0),                              & ! 16-17
             mean_cerliq(:,:,:,0),   mean_cerice(:,:,:,0),                            & ! 18-19
             mean_numc1(:,:,:,0),    mean_numc2(:,:,:,0),                             & ! 20-21
             mean_numliq(:,:,:,0),   mean_numice(:,:,:,0),                            & ! 22-23
             outfname,               output_type='annual', &
             aggregation_type='annual', period_start=annual_start_month, period_end=current_month )

        outfname = trim(outpath_lev3) // "/map3d/ecare_cld3d_" // trim(year) // ".nc"
        call output_netcdf_cld3d( &
             glon,                   glat,                   gplev,                   & ! 01-03
             mean_temp(:,:,:,0),                                                      & ! 04
             n_obs3d(:,:,:,0),                                                        & ! 05
             n_cld3d(:,:,:, 0,0),    n_cld3d(:,:,:, 1,0),    n_cld3d(:,:,:, 2,0),     & ! 06-08
             n_cld3d(:,:,:, 3,0),    n_cld3d(:,:,:, 4,0),    n_cld3d(:,:,:, 5,0),     & ! 09-11
             n_cld3d(:,:,:, 6,0),    n_cld3d(:,:,:, 7,0),    n_cld3d(:,:,:, 8,0),     & ! 12-14
             n_cld3d(:,:,:, 9,0),    n_cld3d(:,:,:,10,0),    n_cld3d(:,:,:,11,0),     & ! 15-17
             n_cld3d(:,:,:,12,0),    n_cld3d(:,:,:,13,0),    n_cld3d(:,:,:,14,0),     & ! 18-20
             outfname,               output_type='annual', &
             aggregation_type='annual', period_start=annual_start_month, period_end=current_month )

        outfname = trim(outpath_lev3) // "/cfad/ecare_cfad_" // trim(year) // ".nc"
        call output_netcdf_cfad( outfname, gdbz, glev, &
                                 pdf_cfad( 1:dimdbz, 1:dimlev, 1:nreg, 0 ), &
                                 'annual', annual_start_month, current_month )

        outfname = trim(outpath_lev3) // "/isccp/ecare_isccp_" // trim(year) // ".nc"
        call output_netcdf_isccp( outfname, &
                                  isccp_fq( 1:dimcot, 1:dimctp, 1:nreg, 0 ), &
                                  isccp_dv( 1:dimcot, 1:dimctp, 1:nreg, 0 ), &
                                  isccp_fq_c4( 1:dimcot, 1:dimctp, 1:nreg, 0 ), &
                                  isccp_dv_c4( 1:dimcot, 1:dimctp, 1:nreg, 0 ), &
                                  'annual', annual_start_month, current_month  )

        outfname = trim(outpath_lev3) // "/map3d/ecare_c4_3d_" // trim(year) // ".nc"
        call output_netcdf_c4_3d( glon, glat, gplev, &
             n_c4_valid(:,:,:,0), n_c4_cloud(:,:,:,0), &
             n_c4_mask_invalid(:,:,:,0), n_c4_type(:,:,:,:,0), &
             n_c4_category(:,:,:,:,0), n_c4_type_qc(:,:,:,:,0), &
             n_c4_category_invalid(:,:,:,0), outfname, &
             'annual', annual_start_month, current_month )
    endif

    !! All-Term Output
    if ( ofinal ) then

        outfname = trim(outpath_lev3) // "/map2d/ecare_map2d_all.nc"
        call output_netcdf_map2d( &
            glon,               glat,               n_obs2d_all(:,:),    & ! 01-03
            pdf_clr_all(:,:),   pdf_cld_all(:,:),   pdf_drz_all(:,:),    & ! 04-06
            pdf_prc_all(:,:),   mean_lwp_all(:,:),  mean_iwp_all(:,:),   & ! 07-09
            mean_cot_all(:,:),  mean_cth_all(:,:),                       & ! 10-11
            mean_cth_c4_all(:,:), mean_ctp_c4_all(:,:),                  & ! 12-13
            outfname,           output_type='all', &
            aggregation_type='all-period', period_start=processing_start_month, period_end=current_month )

        outfname = trim(outpath_lev3) // "/map3d/ecare_map3d_all.nc"
        call output_netcdf_map3d( &
             glon,                     glat,                     gplev,                    & ! 01-03
             n_obs3d_all(:,:,:),                                                           & ! 04
             mean_temp_all(:,:,:),     mean_shum_all(:,:,:),     mean_rhum_all(:,:,:),     & ! 05-07
             mean_cldavel_all(:,:,:),  mean_clddvel_all(:,:,:),                            & ! 08-09
             mean_cldtvel1_all(:,:,:), mean_cldtvel2_all(:,:,:), mean_cldtvelt_all(:,:,:), & ! 10-12
             mean_cwc1_all(:,:,:),     mean_cwc2_all(:,:,:),     mean_ciwc_all(:,:,:),     & ! 13-15
             mean_cer1_all(:,:,:),     mean_cer2_all(:,:,:),                               & ! 16-17
             mean_cerliq_all(:,:,:),   mean_cerice_all(:,:,:),                             & ! 18-19
             mean_numc1_all(:,:,:),    mean_numc2_all(:,:,:),                              & ! 20-21
             mean_numliq_all(:,:,:),   mean_numice_all(:,:,:),                             & ! 22-23
             outfname,                 output_type='all', &
             aggregation_type='all-period', period_start=processing_start_month, period_end=current_month )

        outfname = trim(outpath_lev3) // "/map3d/ecare_cld3d_all.nc"
        call output_netcdf_cld3d( &
             glon,                     glat,                     gplev,                    & ! 01-03
             mean_temp_all(:,:,:),                                                         & ! 04
             n_obs3d_all(:,:,:),                                                           & ! 05
             n_cld3d_all(:,:,:, 0),    n_cld3d_all(:,:,:, 1),    n_cld3d_all(:,:,:, 2),    & ! 06-08
             n_cld3d_all(:,:,:, 3),    n_cld3d_all(:,:,:, 4),    n_cld3d_all(:,:,:, 5),    & ! 09-11
             n_cld3d_all(:,:,:, 6),    n_cld3d_all(:,:,:, 7),    n_cld3d_all(:,:,:, 8),    & ! 12-14
             n_cld3d_all(:,:,:, 9),    n_cld3d_all(:,:,:,10),    n_cld3d_all(:,:,:,11),    & ! 15-17
             n_cld3d_all(:,:,:,12),    n_cld3d_all(:,:,:,13),    n_cld3d_all(:,:,:,14),    & ! 18-20
             outfname,                 output_type='all', &
             aggregation_type='all-period', period_start=processing_start_month, period_end=current_month )

        outfname = trim(outpath_lev3) // "/cfad/ecare_cfad_all.nc"
        call output_netcdf_cfad( outfname, gdbz, glev, &
                                 pdf_cfad_all( 1:dimdbz, 1:dimlev, 1:nreg ), &
                                 'all-period', processing_start_month, current_month )

        outfname = trim(outpath_lev3) // "/isccp/ecare_isccp_all.nc"
        call output_netcdf_isccp( outfname, &
                                  isccp_fq_all( 1:dimcot, 1:dimctp, 1:nreg ), &
                                  isccp_dv_all( 1:dimcot, 1:dimctp, 1:nreg ), &
                                  isccp_fq_c4_all( 1:dimcot, 1:dimctp, 1:nreg ), &
                                  isccp_dv_c4_all( 1:dimcot, 1:dimctp, 1:nreg ), &
                                  'all-period', processing_start_month, current_month  )

        outfname = trim(outpath_lev3) // "/map3d/ecare_c4_3d_all.nc"
        call output_netcdf_c4_3d( glon, glat, gplev, &
             n_c4_valid_all, n_c4_cloud_all, n_c4_mask_invalid_all, &
             n_c4_type_all, n_c4_category_all, n_c4_type_qc_all, &
             n_c4_category_invalid_all, outfname, &
             'all-period', processing_start_month, current_month )
    endif

end subroutine output_summary

end module summary_processor
