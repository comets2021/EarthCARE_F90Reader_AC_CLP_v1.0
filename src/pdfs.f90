module pdfs
    use parameters

    implicit none

    ! Monthly 2D
    real(sp), dimension( dimx,dimy,0:12 ), save  :: mean_lwp            ! LWP [g m-2]
    real(sp), dimension( dimx,dimy,0:12 ), save  :: mean_iwp            ! IWP [g m-2]
    real(sp), dimension( dimx,dimy,0:12 ), save  :: mean_cot            ! COT 
    real(sp), dimension( dimx,dimy,0:12 ), save  :: mean_cth            ! CTH [m]
    real(sp), dimension( dimx,dimy,0:12 ), save  :: mean_cth_c4         ! C4 CTH [m]
    real(sp), dimension( dimx,dimy,0:12 ), save  :: mean_ctp_c4         ! C4 CTP [hPa]
    real(sp), dimension( dimx,dimy,0:12 ), save  :: pdf_clr             ! PDF of Clear-sky
    real(sp), dimension( dimx,dimy,0:12 ), save  :: pdf_cld             ! PDF of -35 < dBZ < -15
    real(sp), dimension( dimx,dimy,0:12 ), save  :: pdf_drz             ! PDF of -15 < dBZ < 0
    real(sp), dimension( dimx,dimy,0:12 ), save  :: pdf_prc             ! PDF of   0 < dBZ
    ! Monthly 3D
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_temp      ! Temperature [K]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_shum      ! Specific Humidity [kg kg-1]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_rhum      ! Relative Humidity [%]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_cldavel   ! Air Vertical Velocity [m s-1]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_clddvel   ! Cloud Doppler Velocity [m s-1]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_cldtvel1  ! Cloud Terminal Velocity 1 [m s-1]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_cldtvel2  ! Cloud Terminal Velocity 2 [m s-1]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_cldtvelt  ! Total Cloud Terminal Velocity [m s-1]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_cwc1      ! Cloud Water Content 1 [g m-3]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_cwc2      ! Cloud Water Content 2 [g m-3]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_ciwc      ! Ice Water Content [g m-3]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_cer1      ! Cloud Effective Radius 1 [micron]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_cer2      ! Cloud Effective Radius 2 [micron]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_cerliq    ! Cloud Effective Radius Liq [micron]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_cerice    ! Cloud Effective Radius Ice [micron]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_numc1     ! Cloud Number Concentration 1 [# cm-3]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_numc2     ! Cloud Number Concentration 2 [# cm-3]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_numliq    ! Cloud Number Concentration Liq [# cm-3]
    real(sp), dimension( dimx,dimy,dimz,0:12 ), save  :: mean_numice    ! Cloud Number Concentration Ice [# cm-3]
    real(sp), dimension( dimdbz,dimlev,nreg,0:12 ), save :: pdf_cfad    ! CFAD
    real(sp), dimension( dimcot,dimctp,nreg,0:12 ), save :: isccp_fq    ! ISCCP Frequency
    real(sp), dimension( dimcot,dimctp,nreg,0:12 ), save :: isccp_dv    ! ISCCP Doppler Velo.
    real(sp), dimension( dimcot,dimctp,nreg,0:12 ), save :: isccp_fq_c4 ! C4 ISCCP Frequency
    real(sp), dimension( dimcot,dimctp,nreg,0:12 ), save :: isccp_dv_c4 ! C4 ISCCP Doppler Velo.

    ! Annual 2D
    real(sp), dimension( dimx,dimy ), save  :: mean_lwp_all             ! LWP [g m-2]
    real(sp), dimension( dimx,dimy ), save  :: mean_iwp_all             ! IWP [g m-2]
    real(sp), dimension( dimx,dimy ), save  :: mean_cot_all             ! COT
    real(sp), dimension( dimx,dimy ), save  :: mean_cth_all             ! CTH [m]
    real(sp), dimension( dimx,dimy ), save  :: mean_cth_c4_all          ! C4 CTH [m]
    real(sp), dimension( dimx,dimy ), save  :: mean_ctp_c4_all          ! C4 CTP [hPa]
    real(sp), dimension( dimx,dimy ), save  :: pdf_clr_all              ! PDF of Clear-sky
    real(sp), dimension( dimx,dimy ), save  :: pdf_cld_all              ! PDF of -35 < dBZ < -15
    real(sp), dimension( dimx,dimy ), save  :: pdf_drz_all              ! PDF of -15 < dBZ < 0
    real(sp), dimension( dimx,dimy ), save  :: pdf_prc_all              ! PDF of   0 < dBZ
    ! Annual 3D
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_temp_all       ! Temperature [K]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_shum_all       ! Specific Humidity [kg kg-1]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_rhum_all       ! Relative Humidity [%]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_cldavel_all    ! Air Vertical Velocity [m s-1]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_clddvel_all    ! Cloud Doppler Velocity [m s-1]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_cldtvel1_all   ! Cloud Terminal Velocity 1 [m s-1]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_cldtvel2_all   ! Cloud Terminal Velocity 2 [m s-1]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_cldtvelt_all   ! Total Cloud Terminal Velocity [m s-1]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_cwc1_all       ! Cloud Water Content 1 [g m-3]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_cwc2_all       ! Cloud Water Content 2 [g m-3]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_ciwc_all       ! Ice Water Content [g m-3]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_cer1_all       ! Cloud Effective Radius 1 [micron]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_cer2_all       ! Cloud Effective Radius 2 [micron]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_cerliq_all     ! Cloud Effective Radius Liq [micron]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_cerice_all     ! Cloud Effective Radius Ice [micron]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_numc1_all      ! Cloud Number Concentration 1 [# cm-3]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_numc2_all      ! Cloud Number Concentration 2 [# cm-3]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_numliq_all     ! Cloud Number Concentration Liq [# cm-3]
    real(sp), dimension( dimx,dimy,dimz ), save  :: mean_numice_all     ! Cloud Number Concentration Ice [# cm-3]
    real(sp), dimension( dimdbz,dimlev,nreg ), save :: pdf_cfad_all     ! CFAD
    real(sp), dimension( dimcot,dimctp,nreg ), save :: isccp_fq_all     ! ISCCP Frequency
    real(sp), dimension( dimcot,dimctp,nreg ), save :: isccp_dv_all     ! ISCCP Doppler Velo.
    real(sp), dimension( dimcot,dimctp,nreg ), save :: isccp_fq_c4_all  ! C4 ISCCP Frequency
    real(sp), dimension( dimcot,dimctp,nreg ), save :: isccp_dv_c4_all  ! C4 ISCCP Doppler Velo.

end module pdfs
