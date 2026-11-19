module numbers
    use parameters

    implicit none

    ! monthly
    real(sp),   dimension( dimx,dimy,0:12 ),       save :: n_obs2d          ! # of profiles
    real(sp),   dimension( dimx,dimy,dimz,0:12 ),  save :: n_obs3d          ! # of observations
    real(sp),   dimension( dimx,dimy,dimz,0:ntype,0:12 ), save :: n_cld3d   ! # of clouds
    integer(kind=4), dimension( dimx,dimy,dimz,0:12 ), save :: n_c4_valid
    integer(kind=4), dimension( dimx,dimy,dimz,0:12 ), save :: n_c4_cloud
    integer(kind=4), dimension( dimx,dimy,dimz,0:12 ), save :: n_c4_mask_invalid
    integer(kind=4), dimension( dimx,dimy,dimz,ntype_c4,0:12 ), save :: n_c4_type
    integer(kind=4), dimension( dimx,dimy,dimz,ncategory_c4,0:12 ), save :: n_c4_category
    integer(kind=4), dimension( dimx,dimy,dimz,ntype_qc_c4,0:12 ), save :: n_c4_type_qc
    integer(kind=4), dimension( dimx,dimy,dimz,0:12 ), save :: n_c4_category_invalid
    integer,    dimension( dimx,dimy,0:12 ),       save :: n_clr            ! C4 mask: no cloudy bin
    integer,    dimension( dimx,dimy,0:12 ),       save :: n_cld            ! -35 <= in-cloud Zmax < -15
    integer,    dimension( dimx,dimy,0:12 ),       save :: n_drz            ! -15 <= in-cloud Zmax < 0
    integer,    dimension( dimx,dimy,0:12 ),       save :: n_prc            ! 0 <= in-cloud Zmax
    integer,    dimension( dimx,dimy,0:12 ),       save :: n_lwp            ! LWP
    real(dp),   dimension( dimx,dimy,0:12 ),       save :: sum_lwp          ! LWP
    integer,    dimension( dimx,dimy,0:12 ),       save :: n_iwp            ! IWP
    real(dp),   dimension( dimx,dimy,0:12 ),       save :: sum_iwp          ! IWP
    integer,    dimension( dimx,dimy,0:12 ),       save :: n_cot            ! COT
    real(dp),   dimension( dimx,dimy,0:12 ),       save :: sum_cot          ! COT
    integer,    dimension( dimx,dimy,0:12 ),       save :: n_cth            ! CTH
    real(dp),   dimension( dimx,dimy,0:12 ),       save :: sum_cth          ! CTH
    integer,    dimension( dimx,dimy,0:12 ),       save :: n_cth_c4         ! C4 CTH
    real(dp),   dimension( dimx,dimy,0:12 ),       save :: sum_cth_c4       ! C4 CTH
    integer,    dimension( dimx,dimy,0:12 ),       save :: n_ctp_c4         ! C4 CTP
    real(dp),   dimension( dimx,dimy,0:12 ),       save :: sum_ctp_c4       ! C4 CTP
    integer,    dimension( dimx,dimy,nctgr,0:12 ), save :: n_pctgr          ! Particle Category
    integer,    dimension( dimx,dimy,ntype,0:12 ), save :: n_ptype          ! Particle Type
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_temp         ! Temperature
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_shum         ! Spc. Humid.
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_rhum         ! Rel. Humid.
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_cldavel        ! Cloud Air Velo
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_cldavel      ! Cloud Air Velo
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_clddvel        ! Doppler Velo
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_clddvel      ! Doppler Velo
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_cldtvel1       ! Term. Velo 1
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_cldtvel1     ! Term. Velo 1
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_cldtvel2       ! Term. Velo 2
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_cldtvel2     ! Term. Velo 2
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_cldtvelt       ! Sedim. Velo
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_cldtvelt     ! Sedim. Velo
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_cwc1           ! CWC1
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_cwc1         ! CWC1
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_cwc2           ! CWC2
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_cwc2         ! CWC2
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_ciwc           ! CIWC
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_ciwc         ! CIWC
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_cer1           ! CER1
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_cer1         ! CER1
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_cer2           ! CER2
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_cer2         ! CER2
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_cerliq         ! CERLIQ
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_cerliq       ! CERLIQ
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_cerice         ! CERICE
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_cerice       ! CERICE
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_numc1          ! NUMC1
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_numc1        ! NUMC1
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_numc2          ! NUMC2
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_numc2        ! NUMC2
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_numliq         ! NUMLIQ
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_numliq       ! NUMLIQ
    integer,    dimension( dimx,dimy,dimz,0:12 ),  save :: n_numice         ! NUMICE
    real(dp),   dimension( dimx,dimy,dimz,0:12 ),  save :: sum_numice       ! NUMICE
    integer(8), dimension( nreg,0:12 ),               save :: n_cfad_obs    ! CFAD All Samples
    integer(8), dimension( dimdbz,dimlev,nreg,0:12 ), save :: n_cfad_bin    ! Binned by dBZ/Height
    integer(8), dimension( dimcot,dimctp,nreg,0:12 ), save :: n_isccp       ! ISCCP Hist.
    real(dp),   dimension( dimcot,dimctp,nreg,0:12 ), save :: sum_isccpdv   ! ISCCP Hist.
    integer(8), dimension( dimcot,dimctp,nreg,0:12 ), save :: n_isccp_c4    ! C4 ISCCP Frequency
    integer(8), dimension( dimcot,dimctp,nreg,0:12 ), save :: n_isccpdv_c4  ! C4 ISCCP Doppler
    real(dp),   dimension( dimcot,dimctp,nreg,0:12 ), save :: sum_isccpdv_c4 ! C4 ISCCP Doppler

    ! annual
    real(sp),   dimension( dimx,dimy ),            save :: n_obs2d_all      ! # of profiles
    real(sp),   dimension( dimx,dimy,dimz ),       save :: n_obs3d_all      ! # of observations
    real(sp),   dimension( dimx,dimy,dimz,0:ntype ), save :: n_cld3d_all    ! # of clouds
    integer(kind=4), dimension( dimx,dimy,dimz ), save :: n_c4_valid_all
    integer(kind=4), dimension( dimx,dimy,dimz ), save :: n_c4_cloud_all
    integer(kind=4), dimension( dimx,dimy,dimz ), save :: n_c4_mask_invalid_all
    integer(kind=4), dimension( dimx,dimy,dimz,ntype_c4 ), save :: n_c4_type_all
    integer(kind=4), dimension( dimx,dimy,dimz,ncategory_c4 ), save :: n_c4_category_all
    integer(kind=4), dimension( dimx,dimy,dimz,ntype_qc_c4 ), save :: n_c4_type_qc_all
    integer(kind=4), dimension( dimx,dimy,dimz ), save :: n_c4_category_invalid_all
    real(dp),   dimension( dimx,dimy ),            save :: n_clr_all        ! C4 mask: no cloudy bin
    real(dp),   dimension( dimx,dimy ),            save :: n_cld_all        ! -35 <= in-cloud Zmax < -15
    real(dp),   dimension( dimx,dimy ),            save :: n_drz_all        ! -15 <= in-cloud Zmax < 0
    real(dp),   dimension( dimx,dimy ),            save :: n_prc_all        ! 0 <= in-cloud Zmax
    integer,    dimension( dimx,dimy ),            save :: n_lwp_all        ! LWP
    real(dp),   dimension( dimx,dimy ),            save :: sum_lwp_all      ! LWP
    integer,    dimension( dimx,dimy ),            save :: n_iwp_all        ! IWP
    real(dp),   dimension( dimx,dimy ),            save :: sum_iwp_all      ! IWP
    integer,    dimension( dimx,dimy ),            save :: n_cot_all        ! COT
    real(dp),   dimension( dimx,dimy ),            save :: sum_cot_all      ! COT
    integer,    dimension( dimx,dimy ),            save :: n_cth_all        ! CTH
    real(dp),   dimension( dimx,dimy ),            save :: sum_cth_all      ! CTH
    integer,    dimension( dimx,dimy ),            save :: n_cth_c4_all     ! C4 CTH
    real(dp),   dimension( dimx,dimy ),            save :: sum_cth_c4_all   ! C4 CTH
    integer,    dimension( dimx,dimy ),            save :: n_ctp_c4_all     ! C4 CTP
    real(dp),   dimension( dimx,dimy ),            save :: sum_ctp_c4_all   ! C4 CTP
    integer,    dimension( dimx,dimy,nctgr ),      save :: n_pctgr_all      ! Particle Category
    integer,    dimension( dimx,dimy,ntype ),      save :: n_ptype_all      ! Particle Type
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_temp_all     ! Temperature
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_shum_all     ! Spc. Humid.
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_rhum_all     ! Rel. Humid.
    integer,    dimension( dimx,dimy,dimz ),       save :: n_cldavel_all    ! Cloud Air Velo
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_cldavel_all  ! Cloud Air Velo
    integer,    dimension( dimx,dimy,dimz ),       save :: n_clddvel_all    ! Doppler Velo
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_clddvel_all  ! Doppler Velo
    integer,    dimension( dimx,dimy,dimz ),       save :: n_cldtvel1_all   ! Term. Velo 1
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_cldtvel1_all ! Term. Velo 1
    integer,    dimension( dimx,dimy,dimz ),       save :: n_cldtvel2_all   ! Term. Velo 2
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_cldtvel2_all ! Term. Velo 2
    integer,    dimension( dimx,dimy,dimz ),       save :: n_cldtvelt_all   ! Sedim. Velo
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_cldtvelt_all ! Sedim. Velo
    integer,    dimension( dimx,dimy,dimz ),       save :: n_cwc1_all       ! CWC1
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_cwc1_all     ! CWC1
    integer,    dimension( dimx,dimy,dimz ),       save :: n_cwc2_all       ! CWC2
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_cwc2_all     ! CWC2
    integer,    dimension( dimx,dimy,dimz ),       save :: n_ciwc_all       ! CIWC
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_ciwc_all     ! CIWC
    integer,    dimension( dimx,dimy,dimz ),       save :: n_cer1_all       ! CER1
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_cer1_all     ! CER1
    integer,    dimension( dimx,dimy,dimz ),       save :: n_cer2_all       ! CER2
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_cer2_all     ! CER2
    integer,    dimension( dimx,dimy,dimz ),       save :: n_cerliq_all     ! CERLIQ
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_cerliq_all   ! CERLIQ
    integer,    dimension( dimx,dimy,dimz ),       save :: n_cerice_all     ! CERICE
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_cerice_all   ! CERICE
    integer,    dimension( dimx,dimy,dimz ),       save :: n_numc1_all      ! NUMC1
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_numc1_all    ! NUMC1
    integer,    dimension( dimx,dimy,dimz ),       save :: n_numc2_all      ! NUMC2
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_numc2_all    ! NUMC2
    integer,    dimension( dimx,dimy,dimz ),       save :: n_numliq_all     ! NUMLIQ
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_numliq_all   ! NUMLIQ
    integer,    dimension( dimx,dimy,dimz ),       save :: n_numice_all     ! NUMICE
    real(dp),   dimension( dimx,dimy,dimz ),       save :: sum_numice_all   ! NUMICE
    integer(8), dimension( nreg ),                 save :: n_cfad_obs_all   ! CFAD All Samples
    integer(8), dimension( dimdbz,dimlev,nreg ),   save :: n_cfad_bin_all   ! Binning by dBZ/Height
    integer(8), dimension( dimcot,dimctp,nreg ),   save :: n_isccp_all      ! ISCCP Hist.
    real(dp),   dimension( dimcot,dimctp,nreg ),   save :: sum_isccpdv_all  ! ISCCP Hist.
    integer(8), dimension( dimcot,dimctp,nreg ),   save :: n_isccp_c4_all   ! C4 ISCCP Frequency
    integer(8), dimension( dimcot,dimctp,nreg ),   save :: n_isccpdv_c4_all ! C4 ISCCP Doppler
    real(dp),   dimension( dimcot,dimctp,nreg ),   save :: sum_isccpdv_c4_all ! C4 ISCCP Doppler
 
end module numbers
