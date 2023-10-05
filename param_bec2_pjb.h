! Parameters related to BEC model
      real c1, c0, c2,c10,c1000,p5,spd,dps,t0_kelvin,yps,mpercm
       parameter ( c1=1., c0=0.0,c2=2., 
     & c10=10.,c1000=1000.,p5=0.5,
     & spd = 86400.0,  
     & dps = c1 / spd,         ! number of days per second
     & yps = c1 / (365.0*spd), ! number of years in a second
     & mpercm = .01,        ! meters per cm
     & t0_kelvin= 273.16)

! Parameters related to autotrophs: how many, how they are ordered
       integer, parameter :: autotroph_cnt=3
#ifdef BEC_COCCO
     &          +1
#endif
       integer, parameter ::
     &   sp_ind   = 1,  ! small phytoplankton
     &   diat_ind = 2,  ! diatoms
     &   diaz_ind = 3   ! diazotrophs
#ifdef BEC_COCCO
     &  ,cocco_ind = 4  ! Coccolithophores
#endif

!
! The following arrays contain one parameter for all of the 3 or 4 autotrophs, in the
! following order:
!  1 --> small phytoplankton
!  2 --> diatoms
!  3 --> diazotrophs
!  4 --> coccolithophores
!
       character*24 sname(autotroph_cnt)     ! short name of each autotroph
       character*80 lname(autotroph_cnt)     ! long name of each autotroph
       integer
     &   Chl_ind(autotroph_cnt),             ! tracer indices for Chl, C, Fe content
     &   C_ind(autotroph_cnt), 
     &   Fe_ind(autotroph_cnt),
     &   Si_ind(autotroph_cnt),
     &   CaCO3_ind(autotroph_cnt)            ! tracer indices for Si, CaCO3 content
       real
     &   kFe(autotroph_cnt),                 ! nutrient uptake half-sat constants
     &   kPO4(autotroph_cnt),
     &   kDOP(autotroph_cnt),
     &   kNO3(autotroph_cnt),
     &   kNO2(autotroph_cnt),
     &   kNH4(autotroph_cnt),
     &   kSiO3(autotroph_cnt),
     &   Qp(autotroph_cnt),                  ! P/C ratio
     &   gQfe_0(autotroph_cnt),
     &   gQfe_min(autotroph_cnt),            ! initial and minimum fe/C ratios
     &   alphaPI(autotroph_cnt),             ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
     &   PCref(autotroph_cnt),               ! max C-spec. grth rate at tref (1/sec)
     &   thetaN_max(autotroph_cnt),          ! max thetaN (Chl/N) (mg Chl/mmol N)
     &   loss_thres(autotroph_cnt),
     &   loss_thres2(autotroph_cnt),         ! conc. where losses go to zero
     &   temp_thres(autotroph_cnt),          ! Temp. where concentration threshold and photosynth. rate drops
     &   mort(autotroph_cnt),                ! linear mortality rates (1/sec), (1/sec/((mmol C/m3))
     &   mort2(autotroph_cnt),               ! quadratic mortality rates (1/sec), (1/sec/((mmol C/m3))
     &   agg_rate_max(autotroph_cnt),        ! max agg. rate (1/d)
     &   agg_rate_min(autotroph_cnt),        ! min agg. rate (1/d)
     &   z_umax_0(autotroph_cnt),            ! max zoo growth rate at tref (1/sec)
     &   z_grz(autotroph_cnt),               ! grazing coef. (mmol C/m^3)
     &   graze_zoo(autotroph_cnt),           ! routing of grazed term, remainder goes to dic
     &   graze_poc(autotroph_cnt),
     &   graze_doc(autotroph_cnt),    
     &   loss_poc(autotroph_cnt),            ! routing of loss term
     &   f_zoo_detr(autotroph_cnt)           ! fraction of zoo losses to detrital
       integer grazee_ind(autotroph_cnt)     ! which grazee category does autotroph belong to
       logical Nfixer(autotroph_cnt)         ! flag set to true for autotrophs that fix N2
       logical imp_calcifier(autotroph_cnt)  ! flag set to true if autotroph implicitly handles calcification
       logical exp_calcifier(autotroph_cnt)  ! flag set to true if autotroph explicitly handles calcification
       
       common /ecosys_bec1/ Chl_ind, C_ind, Fe_ind, Si_ind, CaCO3_ind, grazee_ind
     &        , Nfixer, imp_calcifier, exp_calcifier, sname, lname
       common /ecosys_bec_reals/ kFe, kPO4, kDOP, kNO3, kNO2, kNH4, kSiO3,
     &        Qp, gQfe_0, gQfe_min, alphaPI, PCref, thetaN_max, loss_thres, loss_thres2,
     &        temp_thres, mort, mort2, agg_rate_max, agg_rate_min, z_umax_0, z_grz,
     &        graze_zoo, graze_poc, graze_doc, loss_poc, f_zoo_detr

  !-----------------------------------------------------------------------------
  !   Redfield Ratios, dissolved & particulate
  !-----------------------------------------------------------------------------

       real parm_Red_D_C_P, parm_Red_D_N_P, parm_Red_D_O2_P, parm_Remin_D_O2_P,
     &   parm_Red_P_C_P, parm_Red_D_C_N, parm_Red_P_C_N, parm_Red_D_C_O2,
     &   parm_Remin_D_C_O2, parm_Red_P_C_O2, parm_Red_Fe_C, parm_Red_D_C_O2_diaz
# ifdef Ncycle_SY
     &   ,parm_Red_D_C_O2_NO2V
# endif
       parameter(
     &   parm_Red_D_C_P  = 117.0,                  ! carbon:phosphorus
     &   parm_Red_D_N_P  =  16.0,                  ! nitrogen:phosphorus
     &   parm_Red_D_O2_P = 170.0,                  ! oxygen:phosphorus
     &   parm_Remin_D_O2_P = 138.0,                ! oxygen:phosphorus
     &   parm_Red_P_C_P  = parm_Red_D_C_P,                  ! carbon:phosphorus
     &   parm_Red_D_C_N  = parm_Red_D_C_P/parm_Red_D_N_P,   ! carbon:nitrogen
     &   parm_Red_P_C_N  = parm_Red_D_C_N,                  ! carbon:nitrogen
     &   parm_Red_D_C_O2 = parm_Red_D_C_P/parm_Red_D_O2_P,  ! carbon:oxygen for HNO3 uptake (assuming OM is C117H297O85N16P)
     &   parm_Remin_D_C_O2 = parm_Red_D_C_P/parm_Remin_D_O2_P,  ! carbon:oxygen for NH3 uptake (assuming OM is C117H297O85N16P)
     &   parm_Red_P_C_O2 = parm_Red_D_C_O2,                 ! carbon:oxygen
     &   parm_Red_Fe_C   = 3.0e-6,                 ! iron:carbon
     &   parm_Red_D_C_O2_diaz = parm_Red_D_C_P/150.0   ! carbon:oxygen for diazotrophs
# ifdef Ncycle_SY
     &   ,parm_Red_D_C_O2_NO2V = parm_Red_D_C_P/162.0   ! carbon:oxygen for HNO2 uptake (assuming OM is C117H297O85N16P)
# endif 
     & )

  !----------------------------------------------------------------------------
  !   ecosystem parameters accessible via namelist input
  !----------------------------------------------------------------------------

       real
     &   parm_Fe_bioavail,       ! fraction of Fe flux that is bioavailable
     &   parm_o2_min,            ! min O2 needed for prod & consump. (nmol/cm^3)
     &   parm_o2_min_delta,      ! width of min O2 range (nmol/cm^3)
     &   parm_lowo2_remin_factor, ! remineralization slowdown factor at low o2
     &   parm_kappa_nitrif,      ! nitrification inverse time constant (1/sec)
# ifdef TDEP_REMIN
     &   parm_ktfunc_soft, ! parameter for the temperature dependance of remin on temp (Laufkoeuter 2017)
# endif
# ifdef Ncycle_SY
     &   parm_kao,        ! max ammonium oxidation rate (1/s)
     &   parm_kno,        ! max nitrite oxidation rate (1/s)
     &   parm_ko2_ao,     ! Michaelis Menton O2 constant for ammonium oxidation (mmol m-3)
     &   parm_knh4_ao,    ! Michaelis Menton NH4 constant for ammonium oxidation (mmol m-3)
     &   parm_ko2_no,     ! Michaelis Menton O2 constant for nitrite oxidation (mmol m-3)
     &   parm_kno2_no,    ! Michaelis Menton NO2 constant for nitrite oxidation (mmol m-3)
     &   parm_kno3_den1,       ! no3 half saturation constant for denitrification 1 (no3-> no2, mmol/m^3)
     &   parm_kno2_den2,       ! no2 half saturation constant for denitrification 2 (no2-> n2o, mmol/m^3)
     &   parm_kn2o_den3,       ! n2o half saturation constant for denitrification 3 (n2o-> n2, mmol/m^3)
     &   parm_ko2_oxic,       ! half saturation constant for oxygen consumption during oxic remin (mmol/m^3)	
     &   parm_ko2_den1,       ! exponential decay constant for denitrification 1 (NO3-> NO2, mmol/m^3)
     &   parm_ko2_den2,       ! exponential decay constant for denitrification 2 (NO2-> N2O, mmol/m^3)
     &   parm_ko2_den3,       ! exponential decay constant for denitrification 3 (N2O-> N2, mmol/m^3)
     &   parm_koxic,       ! maximum oxic remin specific rate (mmol C/m^3/s)
     &   parm_kden1,       ! maximum denitrification 1 specific rate (mmol C/m^3/s)  
     &   parm_kden2,       ! maximum denitrification 2 specific rate (mmol C/m^3/s)
     &   parm_kden3,       ! maximum denitrification 3 specific rate (mmol C/m^3/s)
     &   parm_kax,         ! maximum anaerobic ammonium oxidation specific rate (mmol N/m^3/s) 
     &   parm_knh4_ax,     ! NH4 half saturation constant for anammox (mmol/m^3)
     &   parm_kno2_ax,    ! NO2 half saturation constant for anammox (mmol/m^3)
     &   parm_ko2_ax,      ! exponential decay constant for anammox (mmol/m^3)
     &   r_no2tonh4_ax, ! ratio of N consumed from NO2 vs NH4 during anammox (unitless)
     &   parm_n2o_ji_a, ! n2o yield constant (Ji et al.  2015)
     &   parm_n2o_ji_b, ! n2o yield constant (Ji et al.  2015)
     &   parm_n2o_gor_a, ! n2o yield constant (Goreau et al. 1980)
     &   parm_n2o_gor_b, ! n2o yield constant (Goreau et al. 1980)
# endif
# ifdef N2O_NEV
     &   parm_n2o_nev_a1, ! n2o production constant (Nevison et al. 2oo3)
     &   parm_n2o_nev_a2, ! n2o production constant (Nevison et al. 2oo3)
     &   O2_crit_nev, ! [O2] at which we switch to N2O consumption (Nevison et al. 2003)
     &   N2O_cons_tau_nev, ! n2o consumption timescale [1/s] (Cornejo and Farias 2007)
     &   z_scale_nev, ! Depth scale for N2O production (Nevison et al. 2oo3)
# endif
# ifdef EXPLICIT_MICROBES
     &   parm_muaoa,      ! max growth rate of AOA (per day)
     &   parm_munob,      ! max growth rate of NOB (per day)
     &   parm_muaox,      ! max growth rate of AOX (per day)
     &   parm_muaer,      ! max growth rate of AER (per day)
     &   parm_munar,      ! max growth rate of NAR (per day)
     &   parm_munai,      ! max growth rate of NAI (per day)
     &   parm_munir,      ! max growth rate of NIR (per day)
     &   parm_munio,      ! max growth rate of NIO (per day)
     &   parm_munos,      ! max growth rate of NOS (per day)
     &   parm_munao,      ! max growth rate of NOS (per day)
     &   parm_Vmax_doc,   ! max uptake rate of DOC into cell (per day)
     &   parm_Vmax_no3,   ! max uptake rate of NO3 into cell (per day)
     &   parm_Vmax_no2,   ! max uptake rate of NO2 into cell (per day)
     &   parm_yaoa_nh4,   ! biomass yield per NH4 for AOA (mol bioC per mol NH4)
     &   parm_yaoa_oxy,   ! biomass yield per O2 for AOA  (mol bioC per mol O2)
     &   parm_ynob_no2,   ! biomass yield per NO2 for NOB (mol bioC per mol NO2)
     &   parm_ynob_oxy,   ! biomass yield per O2 for NOB  (mol bioC per mol O2)
     &   parm_yaox_nh4,   ! biomass yield per NH4 for AOX (mol bioC per mol NH4)
     &   parm_yaox_no2,   ! biomass yield per NO2 for AOX (mol bioC per mol NO2)
     &   parm_yaer_doc,   ! biomass yield per DOC for AER (mol bioC per mol DOC)
     &   parm_yaer_oxy,   ! biomass yield per O2 for AER  (mol bioC per mol O2)
     &   parm_ynar_aer,   ! biomass yield per DOC for NAR (mol bioC per mol DOC)
     &   parm_ynar_oxy,   ! biomass yield per O2 for NAR  (mol bioC per mol O2)
     &   parm_ynai_aer,   ! biomass yield per DOC for NAI (mol bioC per mol DOC)
     &   parm_ynai_oxy,   ! biomass yield per O2 for NAI  (mol bioC per mol O2)
     &   parm_ynir_aer,   ! biomass yield per DOC for NIR (mol bioC per mol DOC)
     &   parm_ynir_oxy,   ! biomass yield per O2 for NIR  (mol bioC per mol O2)
     &   parm_ynio_aer,   ! biomass yield per DOC for NIO (mol bioC per mol DOC)
     &   parm_ynio_oxy,   ! biomass yield per O2 for NIO  (mol bioC per mol O2)
     &   parm_ynos_aer,   ! biomass yield per DOC for NOS (mol bioC per mol DOC)
     &   parm_ynos_oxy,   ! biomass yield per O2 for NOS  (mol bioC per mol O2)
     &   parm_ynao_aer,   ! biomass yield per DOC for NAO (mol bioC per mol DOC)
     &   parm_ynao_oxy,   ! biomass yield per O2 for NAO  (mol bioC per mol O2)
     &   parm_ynar_ana,   ! biomass yield per DOC for NAR (mol bioC per mol DOC)
     &   parm_ynar_no3,   ! biomass yield per NO3 for NAR (mol bioC per mol NO3)
     &   parm_ynir_ana,   ! biomass yield per DOC for NIR (mol bioC per mol DOC)
     &   parm_ynir_no2,   ! biomass yield per NO2 for NIR (mol bioC per mol NO2)
     &   parm_ynos_ana,   ! biomass yield per DOC for NOS (mol bioC per mol DOC)
     &   parm_ynos_n2o,   ! biomass yield per N2O for NOS (mol bioC per mol N2O)
     &   parm_ynai_ana,   ! biomass yield per DOC for NAI (mol bioC per mol DOC)
     &   parm_ynai_no3,   ! biomass yield per NO3 for NAI (mol bioC per mol NO3)
     &   parm_ynio_ana,   ! biomass yield per DOC for NIO (mol bioC per mol DOC)
     &   parm_ynio_no2,   ! biomass yield per NO2 for NIO (mol bioC per mol NO2)
     &   parm_ynao_ana,   ! biomass yield per DOC for NAO (mol bioC per mol DOC)
     &   parm_ynao_no3,   ! biomass yield per NO3 for NAO (mol bioC per mol NO3)
     &   parm_paox_no3,   ! production of NO3 per biomass (mol NO3 per mol bioC)
     &   parm_aoa_po2,    ! diffusive uptake limit of O2 by AOA (m3 / mmol C / day)
     &   parm_nob_po2,    ! diffusive uptake limit of O2 by NOB (m3 / mmol C / day)
     &   parm_aer_po2,    ! diffusive uptake limit of O2 by AER (m3 / mmol C / day)
     &   parm_nar_po2,    ! diffusive uptake limit of O2 by NAR (m3 / mmol C / day)
     &   parm_nai_po2,    ! diffusive uptake limit of O2 by NAI (m3 / mmol C / day)
     &   parm_nir_po2,    ! diffusive uptake limit of O2 by NIR (m3 / mmol C / day)
     &   parm_nio_po2,    ! diffusive uptake limit of O2 by NIO (m3 / mmol C / day)
     &   parm_nos_po2,    ! diffusive uptake limit of O2 by NOS (m3 / mmol C / day)
     &   parm_nao_po2,    ! diffusive uptake limit of O2 by NAO (m3 / mmol C / day)
     &   parm_nos_pn2o,   ! diffusive uptake limit of N2O by NOS (m3 / mmol C / day)
     &   parm_aoa_CN,     ! C:N ratio of AOA biomass
     &   parm_nob_CN,     ! C:N ratio of NOB biomass
     &   parm_aox_CN,     ! C:N ratio of AOX biomass
     &   parm_aer_CN,     ! C:N ratio of AER biomass
     &   parm_nar_CN,     ! C:N ratio of NAR biomass
     &   parm_nai_CN,     ! C:N ratio of NAI biomass
     &   parm_nir_CN,     ! C:N ratio of NIR biomass
     &   parm_nio_CN,     ! C:N ratio of NIO biomass
     &   parm_nos_CN,     ! C:N ratio of NOS biomass
     &   parm_nao_CN,     ! C:N ratio of NAO biomass
     &   parm_aoa_CP,     ! C:P ratio of AOA biomass
     &   parm_nob_CP,     ! C:P ratio of NOB biomass
     &   parm_aox_CP,     ! C:P ratio of AOX biomass
     &   parm_aer_CP,     ! C:P ratio of AER biomass
     &   parm_nar_CP,     ! C:P ratio of NAR biomass
     &   parm_nai_CP,     ! C:P ratio of NAI biomass
     &   parm_nir_CP,     ! C:P ratio of NIR biomass
     &   parm_nio_CP,     ! C:P ratio of NIO biomass
     &   parm_nos_CP,     ! C:P ratio of NOS biomass
     &   parm_nao_CP,     ! C:P ratio of NAO biomass
     &   parm_aoa_CFe,    ! C:Fe ratio of AOA biomass
     &   parm_nob_CFe,    ! C:Fe ratio of NOB biomass
     &   parm_aox_CFe,    ! C:Fe ratio of AOX biomass
     &   parm_aer_CFe,    ! C:Fe ratio of AER biomass
     &   parm_nar_CFe,    ! C:Fe ratio of NAR biomass
     &   parm_nai_CFe,    ! C:Fe ratio of NAI biomass
     &   parm_nir_CFe,    ! C:Fe ratio of NIR biomass
     &   parm_nio_CFe,    ! C:Fe ratio of NIO biomass
     &   parm_nos_CFe,    ! C:Fe ratio of NOS biomass
     &   parm_nao_CFe,    ! C:Fe ratio of NAO biomass
     &   parm_kaoa_oxy,   ! half-saturation constant for O2 uptake by AOA
     &   parm_kaoa_nh4,   ! half-saturation constant for NH4 uptake by AOA
     &   parm_kaoa_po4,   ! half-saturation constant for PO4 uptake by AOA
     &   parm_kaoa_fer,   ! half-saturation constant for Iron uptake by AOA
     &   parm_kaoa_par,   ! half-saturation constant for light lim of AOA
     &   parm_knob_oxy,   ! half-saturation constant for O2 uptake by NOB
     &   parm_knob_no2,   ! half-saturation constant for NO2 uptake by NOB
     &   parm_knob_po4,   ! half-saturation constant for PO4 uptake by NOB
     &   parm_knob_fer,   ! half-saturation constant for Iron uptake by NOB
     &   parm_kaox_nh4,   ! half-saturation constant for NH4 uptake by AOX
     &   parm_kaox_no2,   ! half-saturation constant for NO2 uptake by AOX
     &   parm_kaox_po4,   ! half-saturation constant for PO4 uptake by AOX
     &   parm_kaox_fer,   ! half-saturation constant for Iron uptake by AOX
     &   parm_kaer_oxy,   ! half-saturation constant for O2 uptake by AER
     &   parm_kaer_doc,   ! half-saturation constant for DOC uptake by AER
     &   parm_kaer_docr,  ! half-saturation constant for DOCr uptake by AER
     &   parm_kaer_po4,   ! half-saturation constant for P uptake by AER
     &   parm_kaer_fer,   ! half-saturation constant for Iron uptake by AER
     &   parm_knar_oxy,   ! half-saturation constant for O2 uptake by NAR
     &   parm_knar_doc,   ! half-saturation constant for DOC uptake by NAR
     &   parm_knar_docr,  ! half-saturation constant for DOCr uptake by NAR
     &   parm_knar_no3,   ! half-saturation constant for NO3 uptake by NAR
     &   parm_knar_po4,   ! half-saturation constant for P uptake by NAR
     &   parm_knar_fer,   ! half-saturation constant for Iron uptake by NAR
     &   parm_knai_oxy,   ! half-saturation constant for O2 uptake by NAI
     &   parm_knai_doc,   ! half-saturation constant for DOC uptake by NAI
     &   parm_knai_docr,  ! half-saturation constant for DOCr uptake by NAI
     &   parm_knai_no3,   ! half-saturation constant for NO3 uptake by NAI
     &   parm_knai_po4,   ! half-saturation constant for P uptake by NAI
     &   parm_knai_fer,   ! half-saturation constant for Iron uptake by NAI
     &   parm_knir_oxy,   ! half-saturation constant for O2 uptake by NIR
     &   parm_knir_doc,   ! half-saturation constant for DOC uptake by NIR
     &   parm_knir_docr,  ! half-saturation constant for DOCr uptake by NIR
     &   parm_knir_no2,   ! half-saturation constant for NO2 uptake by NIR
     &   parm_knir_po4,   ! half-saturation constant for PO4 uptake by NIR
     &   parm_knir_fer,   ! half-saturation constant for Iron uptake by NIR
     &   parm_knio_oxy,   ! half-saturation constant for O2 uptake by NIO
     &   parm_knio_doc,   ! half-saturation constant for DOC uptake by NIO
     &   parm_knio_docr,  ! half-saturation constant for DOCr uptake by NIO
     &   parm_knio_no2,   ! half-saturation constant for NO2 uptake by NIO
     &   parm_knio_po4,   ! half-saturation constant for PO4 uptake by NIO
     &   parm_knio_fer,   ! half-saturation constant for Iron uptake by NIO
     &   parm_knos_oxy,   ! half-saturation constant for O2 uptake by NOS
     &   parm_knos_doc,   ! half-saturation constant for DOC uptake by NOS
     &   parm_knos_docr,  ! half-saturation constant for DOCr uptake by NOS
     &   parm_knos_n2o,   ! half-saturation constant for N2O uptake by NOS
     &   parm_knos_po4,   ! half-saturation constant for PO4 uptake by NOS
     &   parm_knos_fer,   ! half-saturation constant for Iron uptake by NOS
     &   parm_knao_oxy,   ! half-saturation constant for O2 uptake by NAO
     &   parm_knao_doc,   ! half-saturation constant for DOC uptake by NAO
     &   parm_knao_docr,  ! half-saturation constant for DOCr uptake by NAO
     &   parm_knao_no3,   ! half-saturation constant for NO3 uptake by NAO
     &   parm_knao_po4,   ! half-saturation constant for P uptake by NAO
     &   parm_knao_fer,   ! half-saturation constant for Iron uptake by NAO
     &   parm_che_bmin,   ! minimum C biomass beneath which losses do not occur (mmol C / m3)
     &   parm_het_bmin,   ! minimum C biomass beneath which losses do not occur (mmol C / m3)
     &   parm_szoo_mumax, ! maximum growth rate of small zooplankton
     &   parm_szoo_lmort, ! linear mortality rate for small zooplankton
     &   parm_szoo_qmort, ! quadratic mortality rate for small zooplankton
     &   parm_szoo_labil, ! fraction of small zooplankton loss routed to DIC
     &   parm_szoo_fdetr, ! fraction of small zooplankton loss routed to POC
     &   parm_szoo_losst, ! concentration of small zooplankton beneath which losses don't occur
     &   parm_szoo_bmin,  ! minimum C biomass beneath which losses do not occur (mmol C / m3)
     &   parm_kzoo_aoa,   ! half-saturation constant for zoo grazing on AOA
     &   parm_kzoo_nob,   ! half-saturation constant for zoo grazing on NOB
     &   parm_kzoo_aox,   ! half-saturation constant for zoo grazing on AOX
     &   parm_kszoo_aer,  ! half-saturation constant for szoo grazing on AER
     &   parm_kszoo_nar,  ! half-saturation constant for szoo grazing on NAR
     &   parm_kszoo_nai,  ! half-saturation constant for szoo grazing on NAI
     &   parm_kszoo_nir,  ! half-saturation constant for szoo grazing on NIR
     &   parm_kszoo_nio,  ! half-saturation constant for szoo grazing on NIO
     &   parm_kszoo_nos,  ! half-saturation constant for szoo grazing on NOS
     &   parm_kszoo_nao,  ! half-saturation constant for szoo grazing on NAO
     &   parm_aer_qmort,  ! quadratic mortality rate for AER
     &   parm_aoa_qmort,  ! quadratic mortality rate for AOA
     &   parm_nob_qmort,  ! quadratic mortality rate for NOB
     &   parm_aox_qmort,  ! quadratic mortality rate for AOX
     &   parm_nar_qmort,  ! quadratic mortality rate for NAR
     &   parm_nai_qmort,  ! quadratic mortality rate for NAI
     &   parm_nir_qmort,  ! quadratic mortality rate for NIR
     &   parm_nio_qmort,  ! quadratic mortality rate for NIO
     &   parm_nos_qmort,  ! quadratic mortality rate for NOS
     &   parm_nao_qmort,  ! quadratic mortality rate for NAO
     &   parm_aer_lmort,  ! linear mortality rate for AER
     &   parm_aoa_lmort,  ! linear mortality rate for AOA
     &   parm_nob_lmort,  ! linear mortality rate for NOB
     &   parm_aox_lmort,  ! linear mortality rate for AOX
     &   parm_nar_lmort,  ! linear mortality rate for NAR
     &   parm_nai_lmort,  ! linear mortality rate for NAI
     &   parm_nir_lmort,  ! linear mortality rate for NIR
     &   parm_nio_lmort,  ! linear mortality rate for NIO
     &   parm_nos_lmort,  ! linear mortality rate for NOS
     &   parm_nao_lmort,  ! linear mortality rate for NAO
     &   parm_che_mortdoc,  ! fraction of mortality of chemoautotrophs routed to DOC
     &   parm_che_mortpoc,  ! fraction of mortality of chemoautotrophs routed to POC
     &   parm_het_mortdoc,  ! fraction of mortality of heterotrophs routed to DOC
     &   parm_het_mortpoc,  ! fraction of mortality of heterotrophs routed to POC
# endif
     &   parm_nitrif_par_lim,    ! PAR limit for nitrif. (W/m^2)
     &   parm_z_mort_0,          ! zoo linear mort rate (1/sec)
     &   parm_z_mort2_0,         ! zoo quad mort rate (1/sec/((mmol C/m3))
     &   parm_labile_ratio,      ! fraction of loss to DOC that routed directly to DIC (non-dimensional)
     &   parm_POMbury,           ! scale factor for burial of POC, PON, and POP
     &   parm_BSIbury,           ! scale factor burial of bSi
     &   parm_Fe_scavenge_rate0, ! base scavenging rate
     &   parm_f_prod_sp_CaCO3,   ! fraction of sp prod. as CaCO3 prod.
     &   parm_POC_diss,          ! base POC diss len scale
     &   parm_SiO2_diss,         ! base SiO2 diss len scale
     &   parm_CaCO3_diss         ! base CaCO3 diss len scale

       real 
     &   parm_scalelen_z(4),     ! depths of prescribed scalelen values
     &   parm_scalelen_vals(4)   ! prescribed scalelen values

       common /ecosys_bec2/ parm_Fe_bioavail, parm_o2_min, parm_o2_min_delta, parm_lowo2_remin_factor, parm_kappa_nitrif,
     &   parm_nitrif_par_lim, parm_z_mort_0, parm_z_mort2_0, parm_labile_ratio, parm_POMbury,
     &   parm_BSIbury, parm_Fe_scavenge_rate0, parm_f_prod_sp_CaCO3, parm_POC_diss,
     &   parm_SiO2_diss, parm_CaCO3_diss,
     &   parm_scalelen_z, parm_scalelen_vals
# ifdef TDEP_REMIN
     &   , parm_ktfunc_soft
# endif
# ifdef Ncycle_SY
     &   , parm_kao, parm_kno, parm_ko2_ao, parm_knh4_ao, parm_ko2_no, parm_kno2_no, parm_kno3_den1,
     &   parm_kno2_den2, parm_kn2o_den3, parm_ko2_oxic, parm_ko2_den1, parm_ko2_den2, parm_ko2_den3,
     &   parm_koxic, parm_kden1, parm_kden2, parm_kden3, parm_kax, parm_knh4_ax, 
     &   parm_kno2_ax, parm_ko2_ax, r_no2tonh4_ax, parm_n2o_ji_a, parm_n2o_ji_b, parm_n2o_gor_a,
     &   parm_n2o_gor_b
# endif
# ifdef N2O_NEV
     &   ,parm_n2o_nev_a1, parm_n2o_nev_a2, O2_crit_nev, N2O_cons_tau_nev, z_scale_nev
# endif
# ifdef EXPLICIT_MICROBES
     &   , parm_muaer, parm_muaoa, parm_munob, parm_muaox,
     &   parm_munar, parm_munai, parm_munir, parm_munio, parm_munos, parm_munao,
     &   parm_Vmax_doc, parm_Vmax_no3, parm_Vmax_no2,
     &   parm_yaer_doc, parm_yaer_oxy, 
     &   parm_yaoa_nh4, parm_yaoa_oxy, parm_ynob_no2, parm_ynob_oxy, parm_yaox_nh4, parm_yaox_no2,
     &   parm_ynar_aer, parm_ynar_oxy, 
     &   parm_ynai_aer, parm_ynai_oxy, 
     &   parm_ynir_aer, parm_ynir_oxy, 
     &   parm_ynio_aer, parm_ynio_oxy, 
     &   parm_ynos_aer, parm_ynos_oxy, 
     &   parm_ynao_aer, parm_ynao_oxy, 
     &   parm_ynar_ana, parm_ynar_no3,
     &   parm_ynai_ana, parm_ynai_no3,
     &   parm_ynir_ana, parm_ynir_no2,
     &   parm_ynio_ana, parm_ynio_no2,
     &   parm_ynos_ana, parm_ynos_n2o,
     &   parm_ynao_ana, parm_ynao_no3,
     &   parm_paox_no3, parm_aer_po2, parm_aoa_po2, parm_nob_po2,
     &   parm_nar_po2, parm_nai_po2, parm_nir_po2, parm_nio_po2, parm_nos_po2, parm_nao_po2, parm_nos_pn2o,
     &   parm_aer_CN, parm_aoa_CN, parm_nob_CN, parm_aox_CN,
     &   parm_nar_CN, parm_nai_CN, parm_nir_CN, parm_nio_CN, parm_nos_CN, parm_nao_CN,
     &   parm_aer_CP, parm_aoa_CP, parm_nob_CP, parm_aox_CP,
     &   parm_nar_CP, parm_nai_CP, parm_nir_CP, parm_nio_CP, parm_nos_CP, parm_nao_CP,
     &   parm_aer_CFe, parm_aoa_CFe, parm_nob_CFe, parm_aox_CFe,
     &   parm_nar_CFe, parm_nai_CFe, parm_nir_CFe, parm_nio_CFe, parm_nos_CFe, parm_nao_CFe,
     &   parm_kaer_oxy, parm_kaer_doc, parm_kaer_docr, parm_kaer_po4, parm_kaer_fer,
     &   parm_kaoa_oxy, parm_kaoa_nh4, parm_kaoa_po4, parm_kaoa_fer, parm_kaoa_par,
     &   parm_knob_oxy, parm_knob_no2, parm_knob_po4, parm_knob_fer,
     &   parm_kaox_nh4, parm_kaox_no2, parm_kaox_po4, parm_kaox_fer,
     &   parm_knar_oxy, parm_knar_doc, parm_knar_docr, parm_knar_no3,
     &   parm_knar_po4, parm_knar_fer,
     &   parm_knai_oxy, parm_knai_doc, parm_knai_docr, parm_knai_no3,
     &   parm_knai_po4, parm_knai_fer,
     &   parm_knir_oxy, parm_knir_doc, parm_knir_docr, parm_knir_no2,
     &   parm_knir_po4, parm_knir_fer,
     &   parm_knio_oxy, parm_knio_doc, parm_knio_docr, parm_knio_no2,
     &   parm_knio_po4, parm_knio_fer,
     &   parm_knos_oxy, parm_knos_doc, parm_knos_docr, parm_knos_n2o,
     &   parm_knos_po4, parm_knos_fer,
     &   parm_knao_oxy, parm_knao_doc, parm_knao_docr, parm_knao_no3,
     &   parm_knao_po4, parm_knao_fer,
     &   parm_che_bmin, parm_het_bmin,
     &   parm_szoo_mumax, parm_szoo_lmort, parm_szoo_qmort,
     &   parm_szoo_labil, parm_szoo_fdetr, parm_szoo_losst, parm_szoo_bmin,
     &   parm_kzoo_aoa, parm_kzoo_nob, parm_kzoo_aox, parm_kszoo_aer,
     &   parm_kszoo_nar, parm_kszoo_nai, parm_kszoo_nir, parm_kszoo_nio, parm_kszoo_nos, parm_kszoo_nao,
     &   parm_aoa_qmort, parm_nob_qmort, parm_aox_qmort, parm_aer_qmort, 
     &   parm_nar_qmort, parm_nai_qmort, parm_nir_qmort, parm_nio_qmort, parm_nos_qmort, parm_nao_qmort,
     &   parm_aoa_lmort, parm_nob_lmort, parm_aox_lmort, parm_aer_lmort,
     &   parm_nar_lmort, parm_nai_lmort, parm_nir_lmort, parm_nio_lmort, parm_nos_lmort, parm_nao_lmort,
     &   parm_che_mortdoc, parm_che_mortpoc,
     &   parm_het_mortdoc, parm_het_mortpoc
# endif

  !---------------------------------------------------------------------
  !     Misc. Rate constants
  !---------------------------------------------------------------------
  !
  ! DL: dust_fescav_scale changed from 1e9 to 1e10 since dust fluxes in ROMS
  ! are in kg dust/(m^2 s) = 0.1 g dust/(cm^2 s)
  !---------------------------------------------------------------------

       real fe_scavenge_thres1, dust_fescav_scale, fe_max_scale2
       parameter(
     &   fe_scavenge_thres1 = 0.6e-3,   ! upper thres. for Fe scavenging (mmol/m^3)
     &   dust_fescav_scale  = 1.0e10,   ! dust scavenging scale factor (was 1e9 in CESM)
     &   fe_max_scale2      = 1200.0    ! unitless scaling coeff.
     & )

  !---------------------------------------------------------------------
  !     Compute iron remineralization and flux out. In CESM units
  !     dust remin gDust = 0.035 gFe      mol Fe     1e9 nmolFe
  !                        --------- *  ---------- * ----------
  !                         gDust       55.847 gFe     molFe
  !
  !     dust_to_Fe          conversion - dust to iron (CESM: nmol Fe/g Dust) 
  !---------------------------------------------------------------------
  !
  ! DL: in ROMS we have to convert kg dust -> mmol Fe, so the above calculation
  ! becomes:
  !
  !                    0.035 kg Fe         mol Fe         1e3 mmolFe              mmol Fe
  ! dust remin gDust = ---------    *  ---------------- * ----------  =  626.712  -------
  !                     kg dust        55.847e-3 kg Fe      molFe                 kg dust

       real dust_to_Fe
       parameter(dust_to_Fe=0.035/55.847*1.0e6)  ! mmol Fe/kg dust

  !----------------------------------------------------------------------------
  !     Partitioning of phytoplankton growth, grazing and losses
  !
  !     All f_* variables are fractions and are non-dimensional
  !----------------------------------------------------------------------------

       real caco3_poc_min, spc_poc_fac, f_graze_sp_poc_lim,
     &   f_photosp_CaCO3, f_graze_CaCO3_remin, f_graze_si_remin
       parameter(
     &   caco3_poc_min    = 0.4,  ! minimum proportionality between 
     &                            !   QCaCO3 and grazing losses to POC 
     &                            !   (mmol C/mmol CaCO3)
     &   spc_poc_fac      = 0.14, ! small phyto grazing factor (1/mmolC)
     &   f_graze_sp_poc_lim = 0.36, 
     &   f_photosp_CaCO3  = 0.4,  ! proportionality between small phyto 
     &                            ! production and CaCO3 production
     &   f_graze_CaCO3_remin = 0.33, ! fraction of spCaCO3 grazing 
     &                               !          which is remin
     &   f_graze_si_remin    = 0.35  ! fraction of diatom Si grazing 
     &                               !          which is remin
     & )

  !----------------------------------------------------------------------------
  !     fixed ratios
  !----------------------------------------------------------------------------

       real r_Nfix_photo
       parameter(r_Nfix_photo=1.25)         ! N fix relative to C fix (non-dim)

  !-----------------------------------------------------------------------
  !     SET FIXED RATIOS for N/C, P/C, SiO3/C, Fe/C
  !     assumes C/N/P of 117/16/1 based on Anderson and Sarmiento, 1994
  !     for diazotrophs a N/P of 45 is assumed based on Letelier & Karl, 1998
  !-----------------------------------------------------------------------

       real Q, Qp_zoo_pom, Qfe_zoo, gQsi_0, gQsi_max, gQsi_min, QCaCO3_max,
     &   denitrif_C_N, denitrif_NO3_C, denitrif_NO2_C, denitrif_N2O_C
       parameter(
     &   Q             = 0.137,   !N/C ratio (mmol/mmol) of phyto & zoo
     &   Qp_zoo_pom    = 0.00855, !P/C ratio (mmol/mmol) zoo & pom
     &   Qfe_zoo       = 3.0e-6,  !zooplankton fe/C ratio
     &   gQsi_0        = 0.137,   !initial Si/C ratio
     &   gQsi_max      = 0.8,   !max Si/C ratio
     &   gQsi_min      = 0.0429,  !min Si/C ratio
     &   QCaCO3_max    = 0.4,     !max QCaCO3
     &   ! carbon:nitrogen ratio for denitrification
     &   denitrif_C_N  = parm_Red_D_C_P/136.0,
     &   denitrif_NO3_C  = 472.0 / 2.0 / 106.0, ! need to comment on that and check 
     &   denitrif_NO2_C  = 472.0 / 2.0 / 106.0,
     &   denitrif_N2O_C  = 472.0 / 2.0 / 106.0
     & )

  !----------------------------------------------------------------------------
  !     loss term threshold parameters, chl:c ratios
  !----------------------------------------------------------------------------

       real thres_z1, thres_z2, loss_thres_zoo, CaCO3_temp_thres1,
     &   CaCO3_temp_thres2, CaCO3_sp_thres
       parameter(
     &   thres_z1          = 100.0,  ! threshold = C_loss_thres for z shallower than this (m)
     &   thres_z2          = 150.0,  ! threshold = 0 for z deeper than this (m)
     &   loss_thres_zoo    = 0.06,    ! zoo conc. where losses go to zero
     &   CaCO3_temp_thres1 = 6.0,      ! upper temp threshold for CaCO3 prod
     &   CaCO3_temp_thres2 = -2.0,     ! lower temp threshold
     &   CaCO3_sp_thres    = 2.5       ! bloom condition thres (mmolC/m3)
     & )

  !---------------------------------------------------------------------
  !     fraction of incoming shortwave assumed to be PAR
  !---------------------------------------------------------------------

       real f_qsw_par
       parameter(f_qsw_par = 0.45)   ! PAR fraction

  !---------------------------------------------------------------------
  !     Temperature parameters
  !---------------------------------------------------------------------

       real Tref, Q_10
       parameter(
     &   Tref = 30.0,   ! reference temperature (C)
     &   Q_10 = 1.7     ! factor for temperature dependence (non-dim)
     & )

  !---------------------------------------------------------------------
  !  DOM parameters for refractory components and DOP uptake
  !---------------------------------------------------------------------

       real DOC_reminR, DON_reminR, DOFe_reminR, DOP_reminR, DONr_reminR, DOPr_reminR,
     &      DONrefract, DOPrefract
# ifdef EXPLICIT_MICROBES
     &      , DOCrefract
# endif
       parameter(
     &   DOC_reminR  = (c1/(365.0*15.0)) * dps,         ! rate for semi-labile DOC 1/15years
     &   DON_reminR  = (c1/(365.0*15.0)) * dps,         ! rate for semi-labile DON 1/15years
     &   DOFe_reminR = (c1/(365.0*9.0)) * dps,         ! rate for semi-labile DOFe 1/9years
     &   DOP_reminR  = (c1/(365.0*60.0)) * dps,         ! rate for semi-labile DOP 1/60years  
     &   DONr_reminR = (c1/(365.0*9500.0)) * dps,   ! timescale for refrac DON 1/9500yrs
     &   DOPr_reminR = (c1/(365.0*16000.0)) * dps,   ! timescale for refrac DOP 1/16000yrs
# ifdef EXPLICIT_MICROBES
     &   DOCrefract = 0.3,                      ! fraction of DOC to refractory pool (Guo et al 2023 GRL)
     &   DONrefract = 0.3,                      ! fraction of DON to refractory pool
     &   DOPrefract = 0.3                      ! fraction of DOP to refractory pool
# else
     &   DONrefract = 0.0115,                      ! fraction of DON to refractory pool
     &   DOPrefract = 0.003                      ! fraction of DOP to refractory pool
# endif
     & )

  !---------------------------------------------------------------------
  !  Threshold for PAR used in computation of pChl:
  !  Introduced by CN in April 2015 (BEC blew up in SO setup otherwise)
  !---------------------------------------------------------------------

       real PAR_thres_pChl
       parameter(
     &   PAR_thres_pChl = 1e-10
     & )

   !---------------------------------------------------------------------
   !  Thresholds below which denitrification will be reduced. If NO3 conc
   !  is below the values set here, then the denitrification rate will be
   !  multiplied by the inverse of the parameter value and by the NO3 conc.
   !  This helps reduce negative concentrations somewhat. Introduced by
   !  Cara Nissen in April 2015.
   !---------------------------------------------------------------------

       real parm_denitrif_NO3_limit, parm_sed_denitrif_NO3_limit
       parameter(
     &   parm_denitrif_NO3_limit = 5.0,   ! threshold for reducing water column denitrification
     &   parm_sed_denitrif_NO3_limit = 5.0 ! threshold for reducing sediment denitrification
     & )

#ifdef USE_EXPLICIT_VSINK
   !---------------------------------------------------------------------
   !  Vertical sinking speeds used for explicit sinking (in m/day)
   !  Remineralization rates (in day^{-1})
   !
   !  MF:
   !  hard/soft components are both mineral associated (except for POC) 
   !  and thus have the same sinking speed, but different remineralization rates.
   !  Set values to reflect what we had in implicit sinking. Iron desorption
   !  following Moore and Braucher 2008.
   ! 
   !  With sinking speeds (w) and remin rates (k) I tried to represent more or
   !  less what the implicit sinking parameter for dissolution length (D=w/k)
   !  is set to: D(SiO2)=360, D(CaCO3)=400, D(POC)=82
   !  NOTE: This only applies to soft components! Hard components are more resistant
   !        to dissolution and should have approx. D=40000m (Moore 2013, Lima 2014)
   !
   !---------------------------------------------------------------------
       real wDustHard, wPOCHard, wPCaCO3Hard, wPSiO2Hard, wPIronHard
       real wDustSoft, wPOCSoft, wPCaCO3Soft, wPSiO2Soft, wPIronSoft
       parameter( wDustHard=50.0,  wPOCHard=10., wPCaCO3Hard=50., wPSiO2Hard=50., wPIronHard=50. )
       parameter( wDustSoft=50.0,  wPOCSoft=10., wPCaCO3Soft=50., wPSiO2Soft=50., wPIronSoft=50. )
       real param_dusthard_remin, param_pochard_remin, param_pcaco3hard_remin,
     &   param_psio2hard_remin, param_pironhard_remin
       real param_dustsoft_remin, param_pocsoft_remin, param_pcaco3soft_remin,
     &   param_psio2soft_remin, param_pironsoft_remin
       parameter( param_dusthard_remin=8.0e-5, param_pochard_remin=0.0005, param_pcaco3hard_remin=0.00125,
     &   param_psio2hard_remin=0.00125, param_pironhard_remin=0.00125 )
       parameter( param_dustsoft_remin=8.0e-5, param_pocsoft_remin=0.125, param_pcaco3soft_remin=0.125,
     &   param_psio2soft_remin=0.125, param_pironsoft_remin=0.125 )
#endif
