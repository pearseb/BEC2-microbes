!zhc: This file includes parameters which were originally defined in ecosys_parms.F90 from WHOI

        real c1, c0, c2,c1000,p5,spd,dps,t0_kelvin
         parameter ( c1=1., c0=0.0,c2=2., 
     &  c1000=1000.,p5=0.5,
     &  spd = 86400.0,  dps = c1 / spd ,   
     &  t0_kelvin= 273.16)
       
        real  parm_Red_D_C_P,parm_Red_P_C_P,parm_Red_D_C_N,
     &  parm_Red_P_C_N,parm_Red_D_C_O2,parm_Red_P_C_O2,parm_Red_Fe_C 
        parameter ( parm_Red_D_C_P  = 117.0,
     &  parm_Red_P_C_P  = 117.0,     
     &  parm_Red_D_C_N  = 117.0 / 16.0,  
     &  parm_Red_P_C_N  = 117.0 / 16.0,   
     &  parm_Red_D_C_O2 = 117.0 / 170.0, 
     &  parm_Red_P_C_O2 = 117.0 / 170.0,
     &  parm_Red_Fe_C   = 3.0e-6)

       REAL ::
     &  parm_Fe_bioavail,       ! fraction of Fe flux that is bioavailable
     &  parm_o2_min,            ! lower limit of O2 for prod  consumption (mmol/m^3)
     &  parm_kappa_nitrif,      ! nitrification inverse time constant (1/sec)
     &  parm_nitrif_par_lim,    ! PAR limit for nitrif. (W/m^2)
     &  parm_POC_flux_ref,      ! reference POC flux (mmol C/m^2/sec)
     &  parm_z_umax_0,          ! max. zoopl growth rate on sphyto at tref (1/sec)
     &  parm_diat_umax_0,       ! max. zoopl growth rate on diatoms at tref (1/sec)
     &  parm_z_mort_0,          ! zoopl linear mort rate (1/sec)
     &  parm_z_mort2_0,         ! zoopl quad mort rate, higher trophic level grazing (1/sec/((mmol C/m3))
     &  parm_sd_remin_0,        ! small detrital remineralization rate (1/sec)
     &  parm_sp_kNO3,           ! small phyto nitrate uptake half saturation coeff. (mmol N/m3)
     &  parm_diat_kNO3,         ! diatom nitrate uptake half saturation coeff. (mmol N/m3)
     &  parm_sp_kNH4,           ! small phyto ammonium uptake half saturation coeff. (mmol N/m3)
     &  parm_diat_kNH4,         ! diatom ammonium uptake half saturation coeff. (mmol N/m3)
     &  parm_sp_kFe,            ! small phyto iron uptake half saturation coefficient (mmol Fe/m3)
     &  parm_diat_kFe,          ! diatom iron uptake half saturation coefficient (mmol Fe/m3)
     &  parm_diat_kSiO3,        ! diatom si uptake half saturation coefficient (mmol SiO3/m3)
     &  parm_sp_kPO4,           ! small phyto PO4 uptake (mmol P/m^3)
     &  parm_diat_kPO4,         ! diatom PO4 uptate (mmol P/m^3)
     &  parm_z_grz,             ! grazing coefficient for small phyto (mmol C/m^3)
     &  parm_alphaChl,          ! Chl. specific initial slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
     &  parm_labile_ratio,      ! portion of loss to DOC that routed directly to DIC (non-dimensional)
     &  parm_alphaDiaz,         ! chl. spec. init. slope of P_I curve for diazotrophs
     &  parm_diaz_umax_0,       !  max. zoopl growth rate on diazotrophs at tre
     &  gQsi_0, gQsi_coef, gQsi_max    ! Elemental Ratios for Growth

! JDS: parameters that were defined within ecosys_bec.F are now included here:

       REAL ::
     &  PCref,           ! max phyto C-specific growth rate at tref (GD98) (1/sec)
     &  sp_mort,         ! small phyto non-grazing death rate (1/sec)
     &  sp_mort2,        ! small phyto quad mort rate, agg (1/sec/((mmol C/m3))
     &  diat_mort,       ! diatom non-grazing death rate (1/sec)
     &  diat_mort2,      ! diatom quad mort rate, agg/sinking (1/sec/((mmol C/m3))
     &  z_ingest,        ! zoo ingestion coefficient (non-dim)
     &  Q,               ! N/C ratio (mmol N/mmol C) of all phyto and zoopl
     &  Qp,              ! P/C ratio (mmol P/mmol C) of small phyto, diatoms, zoo
     &  thres_z1,        ! threshold = C_loss_thres for z shallower than this (m) ! HF: used to be in cm
     &  thres_z2,        ! threshold = 0 for z deeper than this (m) ! HF: used to be in cm
     &  PCrefDiaz,       ! max Diaz C-specific growth rate at tref (GD98) (1/sec)
     &  Qp_diaz,          ! diazotroph P/C ratio
     &  diaz_mort,        ! diazotroph non-grazing death rate (1/sec)
     &  diaz_kPO4,        ! diazotroph half-saturation const. for P uptake
     &  diaz_kFe,         ! diazotroph half-saturation const. for P uptake
     &  Qfe_zoo,           ! zooplankton fe/C ratio
     &  parm_diss,
     &  parm_gamma,
     &  Q10_POC,Q10_growth,
     &  parm_decay_hard,
     &  parm_nitrif_o2

       common/eco_para/parm_Fe_bioavail,  
     &   parm_o2_min,  
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref, 
     &   parm_z_umax_0, 
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0, 
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe, 
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4, 
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl, 
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max,
! JDS added Parameters
     &   PCref,sp_mort,sp_mort2,diat_mort,diat_mort2,z_ingest,
     &   Q,Qp,thres_z1,thres_z2,PCrefDiaz,Qp_diaz,diaz_mort,
     &   diaz_kPO4,diaz_kFe,Qfe_zoo,parm_diss,parm_gamma,
     &   Q10_POC,Q10_growth,parm_decay_hard,parm_nitrif_o2

       namelist/biology_param_bec_nml/parm_Fe_bioavail, 
     &   parm_o2_min, 
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref,
     &   parm_z_umax_0,
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0,
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe,
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4,
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl,
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max,
! JDS added Parameters
     &   PCref,sp_mort,sp_mort2,diat_mort,diat_mort2,z_ingest,
     &   Q,Qp,thres_z1,thres_z2,PCrefDiaz,Qp_diaz,diaz_mort,
     &   diaz_kPO4,diaz_kFe,Qfe_zoo,parm_diss,parm_gamma,
     &   Q10_POC,Q10_growth,parm_decay_hard,parm_nitrif_o2
