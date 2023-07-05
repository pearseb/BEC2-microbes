! SY : March 4, decrease Zumax on diatoms
#include "cppdefs.h"
#ifdef BIOLOGY_BEC2
      subroutine ecosys2_init()
      implicit none
 
#include "param.h" 
#include "param_bec2.h"
#include "ecosys_bec2.h"
#include "grid.h"

      integer auto_ind
      ! Variables used for namelist parameter input:
      integer status,  lvar,itrc,lenstr
#ifdef BIOPAR_NAMELIST
      character*100 filename_biopar_bec_nml, filename_biopar_out_bec_nml

      ! Define namelist for ecosystem parameters:
      namelist /ecosys_parms_nml/ 
     &    parm_Fe_bioavail,
     &    parm_o2_min,
     &    parm_o2_min_delta,
     &    parm_lowo2_remin_factor,
     &    parm_kappa_nitrif,
# ifdef TDEP_REMIN
     &    parm_ktfunc_soft,
# endif
# ifdef Ncycle_SY
     &    parm_kao,
     &    parm_kno,
     &    parm_ko2_ao,          
     &    parm_knh4_ao,        
     &    parm_ko2_no,          
     &    parm_kno2_no,         
     &    parm_kno3_den1,      
     &    parm_kno2_den2,       
     &    parm_kn2o_den3,       
     &    parm_ko2_den1,       
     &    parm_ko2_oxic,       
     &    parm_ko2_den2,       
     &    parm_ko2_den3,    
     &    parm_koxic,       
     &    parm_kden1,       
     &    parm_kden2,       
     &    parm_kden3,       
     &    parm_kax,       
     &    parm_knh4_ax,       
     &    parm_kno2_ax,
     &    parm_ko2_ax,       
     &    r_no2tonh4_ax
     &    parm_n2o_ji_a,
     &    parm_n2o_ji_b,
     &    parm_n2o_gor_a,
     &    parm_n2o_gor_b,
# endif
# ifdef N2O_NEV
     &    parm_n2o_nev_a1,
     &    parm_n2o_nev_a2,
     &    N2O_cons_tau_nev,
     &    O2_crit_nev,
     &    z_scale_nev,
# endif
# ifdef EXPLICIT_MICROBES
     &    parm_muaoa,
     &    parm_munob,
     &    parm_muaox,
     &    parm_munar,
     &    parm_munai,
     &    parm_munir,
     &    parm_munio,
     &    parm_munos,
     &    parm_yaoa_nh4,
     &    parm_yaoa_oxy,
     &    parm_ynob_no2,
     &    parm_ynob_oxy,
     &    parm_yaox_nh4,
     &    parm_yaox_no2,
     &    parm_paox_no3,
     &    parm_ynar_aer,
     &    parm_ynar_oxy,
     &    parm_ynar_ana,
     &    parm_ynar_no3,
     &    parm_ynai_aer,
     &    parm_ynai_oxy,
     &    parm_ynai_ana,
     &    parm_ynai_no3,
     &    parm_ynir_aer,
     &    parm_ynir_oxy,
     &    parm_ynir_ana,
     &    parm_ynir_no2,
     &    parm_ynio_aer,
     &    parm_ynio_oxy,
     &    parm_ynio_ana,
     &    parm_ynio_no2,
     &    parm_ynos_aer,
     &    parm_ynos_oxy,
     &    parm_ynos_ana,
     &    parm_ynos_n2o,
     &    parm_aoa_po2,
     &    parm_nob_po2,
     &    parm_nar_po2,
     &    parm_nai_po2,
     &    parm_nir_po2,
     &    parm_nio_po2,
     &    parm_nos_po2,
     &    parm_nos_pn2o,
     &    parm_aoa_CN,
     &    parm_nob_CN,
     &    parm_aox_CN,
     &    parm_nar_CN,
     &    parm_nai_CN,
     &    parm_nir_CN,
     &    parm_nio_CN,
     &    parm_nos_CN,
     &    parm_aoa_CP,
     &    parm_nob_CP,
     &    parm_aox_CP,
     &    parm_nar_CP,
     &    parm_nai_CP,
     &    parm_nir_CP,
     &    parm_nio_CP,
     &    parm_nos_CP,
     &    parm_aoa_CFe,
     &    parm_nob_CFe,
     &    parm_aox_CFe,
     &    parm_nar_CFe,
     &    parm_nai_CFe,
     &    parm_nir_CFe,
     &    parm_nio_CFe,
     &    parm_nos_CFe,
     &    parm_kaoa_oxy,
     &    parm_kaoa_nh4,
     &    parm_kaoa_po4,
     &    parm_kaoa_fer,
     &    parm_kaoa_par,
     &    parm_knob_oxy,
     &    parm_knob_no2,
     &    parm_knob_po4,
     &    parm_knob_fer,
     &    parm_kaox_nh4,
     &    parm_kaox_no2,
     &    parm_kaox_po4,
     &    parm_kaox_fer,
     &    parm_knar_oxy,
     &    parm_knar_doc,
     &    parm_knar_docr,
     &    parm_knar_no3,
     &    parm_knar_po4,
     &    parm_knar_fer,
     &    parm_knai_oxy,
     &    parm_knai_doc,
     &    parm_knai_docr,
     &    parm_knai_no3,
     &    parm_knai_po4,
     &    parm_knai_fer,
     &    parm_knir_oxy,
     &    parm_knir_doc,
     &    parm_knir_docr,
     &    parm_knir_no2,
     &    parm_knir_po4,
     &    parm_knir_fer,
     &    parm_knio_oxy,
     &    parm_knio_doc,
     &    parm_knio_docr,
     &    parm_knio_no2,
     &    parm_knio_po4,
     &    parm_knio_fer,
     &    parm_knos_oxy,
     &    parm_knos_doc,
     &    parm_knos_docr,
     &    parm_knos_n2o,
     &    parm_knos_po4,
     &    parm_knos_fer,
     &    parm_che_bmin,
     &    parm_het_bmin,
     &    parm_szoo_mumax,
     &    parm_szoo_lmort,
     &    parm_szoo_qmort,
     &    parm_szoo_labil,
     &    parm_szoo_fdetr,
     &    parm_szoo_losst,
     &    parm_szoo_bmin,
     &    parm_kzoo_aoa,
     &    parm_kzoo_nob,
     &    parm_kzoo_aox,
     &    parm_kszoo_nar,
     &    parm_kszoo_nai,
     &    parm_kszoo_nir,
     &    parm_kszoo_nio,
     &    parm_kszoo_nos,
     &    parm_aoa_qmort,
     &    parm_nob_qmort,
     &    parm_aox_qmort,
     &    parm_nar_qmort,
     &    parm_nai_qmort,
     &    parm_nir_qmort,
     &    parm_nio_qmort,
     &    parm_nos_qmort,
     &    parm_aoa_lmort,
     &    parm_nob_lmort,
     &    parm_aox_lmort,
     &    parm_nar_lmort,
     &    parm_nai_lmort,
     &    parm_nir_lmort,
     &    parm_nio_lmort,
     &    parm_nos_lmort,
     &    parm_che_mortdoc,
     &    parm_che_mortpoc,
     &    parm_het_mortdoc,
     &    parm_het_mortpoc,
# endif
     &    parm_nitrif_par_lim,
     &    parm_z_mort_0,
     &    parm_z_mort2_0,
     &    parm_labile_ratio,
     &    parm_POMbury,
     &    parm_BSIbury,
     &    parm_Fe_scavenge_rate0,
     &    parm_f_prod_sp_CaCO3,
     &    parm_POC_diss,
     &    parm_SiO2_diss,
     &    parm_CaCO3_diss,
     &    parm_scalelen_z,
     &    parm_scalelen_vals,
     &    kFe, kPO4, kDOP, kNO3, kNH4, kSiO3, Qp, gQfe_0, gQfe_min, 
     &    alphaPI, PCref,
     &    thetaN_max, loss_thres, loss_thres2, temp_thres, mort, mort2, 
     &    agg_rate_max,
     &    agg_rate_min, z_umax_0, z_grz, graze_zoo, graze_poc, 
     &    graze_doc, loss_poc, f_zoo_detr


      filename_biopar_bec_nml = 'biopar_bec.nml'
      filename_biopar_out_bec_nml = 'biopar_out_bec.nml'
#endif

      !---------------------------------------------------------------------------
      !   default namelist settings
      !---------------------------------------------------------------------------

!      parm_Fe_bioavail       = 0.02 ! SY: switched to soluble iron forcing
      parm_Fe_bioavail       = 1.0 ! SY: switched to soluble iron forcing
      !parm_o2_min            = 4.0  ! mmol/m^3 = nmol/cm^3
      parm_o2_min            = 1.0  ! mmol/m^3 = nmol/cm^3 (MF 11.11.2016: Taking down DENITRIF following discussion with SY (Babbin et al. 2014))
      parm_lowo2_remin_factor = 3.3
      parm_o2_min_delta      = 2.0  ! mmol/m^3 = nmol/cm^3
      parm_kappa_nitrif      = 0.06 * dps ! (= 1/( days))
      parm_nitrif_par_lim    = 1.0
      parm_z_mort_0          = 0.1 * dps
      parm_z_mort2_0         = 0.4 * dps
      parm_labile_ratio      = 0.85
      parm_POMbury           = 1.0 ! x1 default
# ifdef TDEP_REMIN
      parm_ktfunc_soft       = 0.055
# endif
# ifdef Ncycle_SY
      parm_kao = 0.0500 * dps 
      parm_kno = 0.0500 * dps
      parm_ko2_ao = 0.333    
      parm_knh4_ao = 0.305   
      parm_ko2_no = 0.778    
      parm_kno2_no = 0.509   
      parm_kno3_den1 = 1     
      parm_kno2_den2 = 0.01
      parm_kn2o_den3 = 0.159
      parm_ko2_oxic = 1
      parm_ko2_den1 = 6
      parm_ko2_den2 = 2.3  ! from 2.969 in tune5
      parm_ko2_den3 = 0.506
      parm_koxic = 0.08 * dps ! CONSTANT
      parm_kden1 = 0.0160 * dps ! from 0.205 in tune5
      parm_kden2 = 0.008 * dps
      parm_kden3 = 0.0496 * dps
      parm_kax = 0.441 * dps
      parm_knh4_ax = 1
      parm_kno2_ax = 1
      parm_ko2_ax = 6
      r_no2tonh4_ax = 1.00 ! Strous et al. 1998 (Some NO2 is oxidized to NO3 via CO2)
      parm_n2o_ji_a = 0.3
      parm_n2o_ji_b = 0.1
      parm_n2o_gor_a = 0.2
      parm_n2o_gor_b = -0.0004
# endif
# ifdef N2O_NEV
      parm_n2o_nev_a1 = 0.2
      parm_n2o_nev_a2 = 0.08
      O2_crit_nev = 3.0
      N2O_cons_tau_nev = 0.5 * dps
      z_scale_nev = 3000.0
#endif
# ifdef EXPLICIT_MICROBES
      parm_muaoa = 1.0 * dps
      parm_munob = 2.0 * dps
      parm_muaox = 0.5 * dps
      parm_munar = 1.0 * dps
      parm_munai = 1.0 * dps
      parm_munir = 1.0 * dps
      parm_munio = 1.0 * dps
      parm_munos = 1.0 * dps
      parm_yaoa_nh4 = 11.0
      parm_yaoa_oxy = 15.5
      parm_ynob_no2 = 27.8
      parm_ynob_oxy = 12.9
      parm_yaox_nh4 = 14.0
      parm_yaox_no2 = 16.2
      parm_paox_no3 = 2.0
      parm_ynar_aer = 5.9  ! 5.9 (default)
      parm_ynar_oxy = 5.4  ! 5.4 (default)
      parm_ynar_ana = 8.6  ! 6.5 (default)  ; 8.6 when y(NAR) < y(NIR) < y(NOS)
      parm_ynar_no3 = 16.8  ! 12.3 (default); 16.8 when y(NAR) < y(NIR) < y(NOS)
      parm_ynir_aer = 5.9  ! 5.9 (default)
      parm_ynir_oxy = 5.4  ! 5.4 (default)
      parm_ynir_ana = 6.5  ! 6.5 (default)  ; 6.5 when y(NAR) < y(NIR) < y(NOS)
      parm_ynir_no2 = 12.1  ! 12.3 (default); 12.1 when y(NAR) < y(NIR) < y(NOS)
      parm_ynos_aer = 5.9  ! 5.9 (default)
      parm_ynos_oxy = 5.4  ! 5.4 (default)
      parm_ynos_ana = 4.3  ! 6.5 (default)  ; 4.3 when y(NAR) < y(NIR) < y(NOS)
      parm_ynos_n2o = 14.8  ! 24.6 (default); 14.8 when y(NAR) < y(NIR) < y(NOS)
      parm_ynai_aer = 5.9  ! 5.9 (default)
      parm_ynai_oxy = 5.4  ! 5.4 (default)
      parm_ynai_ana = 9.3  ! 6.5 (default)  ; 9.3 when y(NAR) < y(NIR) < y(NOS)
      parm_ynai_no3 = 9.1  ! 12.3 (default); 9.1 when y(NAR) < y(NIR) < y(NOS)
      parm_ynio_aer = 5.9  ! 5.9 (default)
      parm_ynio_oxy = 5.4  ! 5.4 (default)
      parm_ynio_ana = 8.2  ! 6.5 (default)  ; 8.2 when y(NAR) < y(NIR) < y(NOS)
      parm_ynio_no2 = 10.6  ! 12.3 (default); 10.6 when y(NAR) < y(NIR) < y(NOS)
      parm_aoa_po2 = 293.0 * dps
      parm_nob_po2 = 71.0 * dps
      parm_nar_po2 = 450.0 * dps  ! 450.0 (default)
      parm_nai_po2 = 450.0 * dps  ! 450.0 (default)
      parm_nir_po2 = 450.0 * dps  ! 450.0 (default)
      parm_nio_po2 = 450.0 * dps  ! 450.0 (default)
      parm_nos_po2 = 450.0 * dps  ! 450.0 (default)
      parm_nos_pn2o = 452.0 * dps
      parm_aoa_CN = 4.0    ! Important ! C:N(chemos) <= C:N(zoo) = 7.3
      parm_nob_CN = 4.0
      parm_aox_CN = 5.0
      parm_nar_CN = 4.5      ! Important ! C:N(chemos) <= C:N(zoo) 
      parm_nai_CN = 4.5      ! Important ! C:N(chemos) <= C:N(zoo) 
      parm_nir_CN = 4.5
      parm_nio_CN = 4.5
      parm_nos_CN = 4.5
      parm_aoa_CP = 55.0   ! Shafiee et al. 2019, Meador et al. 2020
      parm_nob_CP = 55.0   ! Important ! C:P(chemos) <= C:P(zoo) = 117 
      parm_aox_CP = 55.0
      parm_nar_CP = 35.0   ! Important ! C:P(chemos) <= C:P(zoo)
      parm_nai_CP = 35.0   ! Important ! C:P(chemos) <= C:P(zoo)
      parm_nir_CP = 35.0   ! White et al. 2019; Zimmerman et al. 2014
      parm_nio_CP = 35.0   ! White et al. 2019; Zimmerman et al. 2014
      parm_nos_CP = 35.0
      parm_aoa_CFe = 1.0/20.0e-6  ! Shafiee et al. 2019, Meador et al. 2020
      parm_nob_CFe = 1.0/20.0e-6  ! Important ! C:Fe(chemos) <= C:Fe(zoo) = 1/3.0e-6 
      parm_aox_CFe = 1.0/20.0e-6
      parm_nar_CFe = 1.0/20.0e-6  ! Important ! C:Fe(chemos) <= C:Fe(zoo)
      parm_nai_CFe = 1.0/20.0e-6  ! Important ! C:Fe(chemos) <= C:Fe(zoo)
      parm_nir_CFe = 1.0/20.0e-6  ! Fourquez et al. 2020 says 1/3.6e-6
      parm_nio_CFe = 1.0/20.0e-6  ! Fourquez et al. 2020 says 1/3.6e-6
      parm_nos_CFe = 1.0/20.0e-6
      parm_kaoa_oxy = 0.3  ! 333 nM at oxic-anoxic interface (Bristow et al. 2016)
      parm_kaoa_nh4 = 0.1
      parm_kaoa_po4 = 0.0
      parm_kaoa_fer = 0.36
      parm_kaoa_par = 0.25 ! 50% reduction in NH4ox at 25% surface PAR (Travis et al. (2023) Biogeosciences) 
      parm_knob_oxy = 0.7  ! 778 nM at oxic-anoxic interface (Bristow et al. 2016)
      parm_knob_no2 = 0.1
      parm_knob_po4 = 0.0
      parm_knob_fer = 0.0
      parm_kaox_nh4 = 0.45
      parm_kaox_no2 = 0.45
      parm_kaox_po4 = 0.0
      parm_kaox_fer = 0.0
      parm_knar_oxy = 0.05  ! 10 to 200 nM at extremely low oxygen concentrations (Tiano et al., 2014)  
      parm_knar_doc = 10.0
      parm_knar_docr = 1000.0
      parm_knar_no3 = 15.0
      parm_knar_po4 = 0.0
      parm_knar_fer = 0.0
      parm_knai_oxy = 0.05  ! 10 to 200 nM at extremely low oxygen concentrations (Tiano et al., 2014)  
      parm_knai_doc = 10.0
      parm_knai_docr = 1000.0
      parm_knai_no3 = 15.0
      parm_knai_po4 = 0.0
      parm_knai_fer = 0.0
      parm_knir_oxy = 0.05  ! 10 to 200 nM at extremely low oxygen concentrations (Tiano et al., 2014)  
      parm_knir_doc = 10.0
      parm_knir_docr = 1000.0
      parm_knir_no2 = 10.0
      parm_knir_po4 = 0.0
      parm_knir_fer = 0.0
      parm_knio_oxy = 0.05  ! 10 to 200 nM at extremely low oxygen concentrations (Tiano et al., 2014)  
      parm_knio_doc = 10.0
      parm_knio_docr = 1000.0
      parm_knio_no2 = 10.0
      parm_knio_po4 = 0.0
      parm_knio_fer = 0.0
      parm_knos_oxy = 0.05  ! 10 to 200 nM at extremely low oxygen concentrations (Tiano et al., 2014)  
      parm_knos_doc = 10.0
      parm_knos_docr = 1000.0
      parm_knos_n2o = 0.3  ! based on (Sun et al., 2021) in ODZ k = 0.3 µM N2O
      parm_knos_po4 = 0.0
      parm_knos_fer = 0.0
      parm_che_bmin = 1e-4   ! Equivalent to 30-100 cells per mL for cell C content of 40-10 fg C
      parm_het_bmin = 1e-4   ! Equivalent to 80-240 cells per mL for cell C content of 15-5 fg C
      parm_szoo_mumax = 3.3 * dps
      parm_szoo_lmort = 0.05 * dps
      parm_szoo_qmort = 0.5 * dps
      parm_szoo_labil = 0.85
      parm_szoo_fdetr = 0.25
      parm_szoo_losst = 0.01
      parm_szoo_bmin = 1e-3
      parm_kzoo_aoa = 2.4
      parm_kzoo_nob = 2.4
      parm_kzoo_aox = 2.4
      parm_kszoo_nar = 6.2
      parm_kszoo_nai = 6.2
      parm_kszoo_nir = 6.2
      parm_kszoo_nio = 6.2
      parm_kszoo_nos = 6.2
      parm_aoa_qmort = 0.1 * dps
      parm_nob_qmort = 0.1 * dps
      parm_aox_qmort = 0.1 * dps
      parm_nar_qmort = 0.1 * dps
      parm_nai_qmort = 0.1 * dps
      parm_nir_qmort = 0.1 * dps
      parm_nio_qmort = 0.1 * dps
      parm_nos_qmort = 0.1 * dps
      parm_aoa_lmort = 0.025 * parm_muaoa  ! average of % heterotrophic cells infected by viruses
      parm_nob_lmort = 0.025 * parm_munob  ! (Proctor & Fuhrman 1990 Nature)
      parm_aox_lmort = 0.025 * parm_muaox
      parm_nar_lmort = 0.025 * parm_munar
      parm_nai_lmort = 0.025 * parm_munai
      parm_nir_lmort = 0.025 * parm_munir
      parm_nio_lmort = 0.025 * parm_munio
      parm_nos_lmort = 0.025 * parm_munos
      parm_che_mortdoc = 0.5
      parm_che_mortpoc = 0.3   ! Guo et al. 2023 GRL
      parm_het_mortdoc = 0.5
      parm_het_mortpoc = 0.3   ! Guo et al. 2023 GRL
# endif 

#ifdef IODINE
      parm_iod1_cstflux = 5.0 * yps ! Constant excess Iodine flux in mmol m-2 -s
      parm_iod2_fz2 = 5.0 * yps  ! Iodine flux at z=z2 in mmol m-2 -s
      parm_iod2_z2 = 100.0 ! Reference depth for Iodine flux calculation
      parm_iod2_z1 = 50.0 ! Depth at which Iodine flux is maximum
      parm_iod2_b = -0.52 ! power law exponent
      parm_iod3_ric_ox = 10.0e-4 ! ratio of excess I produced to C oxidized in the sediments 
      parm_iod4_maxflux = 10.0 * yps  ! Max Iodine flux in mmol m-2 -s 
# endif
      parm_BSIbury           = 0.65 ! x1 default
      parm_Fe_scavenge_rate0 = 2.5  ! test initial scavenging rate 3 times higher
      parm_f_prod_sp_CaCO3   = 0.055 ! x1 default
      parm_POC_diss          = 88.0
      parm_SiO2_diss         = 250.0
      parm_CaCO3_diss        = 150.0

      parm_scalelen_z    = (/ 100.0, 250.0, 500.0, 700.0 /) ! DL: converted to m
      parm_scalelen_vals = (/     1.0,     2.9,     5.6,      5.7 /)

      kFe = (/ 0.025e-3, 0.05e-3, 0.025e-3 /) ! mmol Fe/m3
      kPO4 = (/ 0.0075, 0.06, 0.015 /)         ! mmol P/m3
      kDOP = (/ 0.22, 0.6, 0.05 /)          !
      kNO3 = (/ 0.2, 0.6, 2.0 /)            ! mmol N/m3   ! 0.11 0.45  8
      kNO2 = (/ 0.0300, 0.24, 0.8 /)        ! mmol N/m3  ! 0.0
      kNH4 = (/ 0.0025, 0.02, 0.0667 /)     ! mmol N/m3  ! 0.01 0.04  0.8
      kSiO3 = (/ 0.0, 0.6, 0.0 /)           ! mmol SiO3/m3
      Qp = (/ 0.00855, 0.00855, 0.00855 /)
      gQfe_0 = (/ 33.0e-6, 33.0e-6, 66.0e-6 /)
      gQfe_min = (/ 2.7e-6, 2.7e-6, 6.0e-6 /)
      alphaPI = (/ 0.4*dps, 0.31*dps, 0.33*dps /)  ! previous defaults, Chl. specific initial slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
      PCref = (/ 3.54*dps, 3.54*dps, 0.51*dps /) ! MF: sensitivity tests70%
      thetaN_max = (/ 2.5, 4.0, 2.5 /)
      loss_thres = (/ 0.02, 0.02, 0.02 /)
      loss_thres2 = (/ 0.0, 0.0, 0.001 /)
      temp_thres = (/ -10.0, -10.0, 16.0 /)
      mort = (/ 0.1*dps, 0.1*dps, 0.1*dps /)  ! non-grazing mortality (1/sec)
      mort2 = (/ 0.1*dps, 0.1*dps, 0.1*dps /)
      agg_rate_max = (/ 0.5, 0.5, 0.5 /)
      agg_rate_min = (/ 0.01, 0.02, 0.01 /)
      z_umax_0 = (/ 3.3*dps, 3.05*dps, 3.3*dps /) ! MF: sensitivity test 75%
      z_grz = (/ 1.2, 1.2, 1.2 /)      ! grazing coefficient (mmol C/m^3)
      graze_zoo = (/ 0.3, 0.25, 0.3 /)
      graze_poc = (/ 0.0, 0.4, 0.1 /)
      graze_doc = (/ 0.06, 0.06, 0.06 /)
      loss_poc = (/ 0.0, 0.0, 0.0 /)
      f_zoo_detr = (/ 0.2, 0.3, 0.1 /)


      auto_ind = sp_ind
      sname(auto_ind)         = 'sp'
      lname(auto_ind)         = 'Small Phyto'
      Nfixer(auto_ind)        = .false.
      exp_calcifier(auto_ind) = .false.
      grazee_ind(auto_ind)    = auto_ind
      C_ind(auto_ind)         = spC_ind_t
      Chl_ind(auto_ind)       = spChl_ind_t
      Fe_ind(auto_ind)        = spFe_ind_t
      Si_ind(auto_ind)        = 0
      imp_calcifier(auto_ind) = .true. ! this is over-written in case of BEC_COCCO
      CaCO3_ind(auto_ind)     = spCaCO3_ind_t

      ! More default parameters for diatoms:
      auto_ind = diat_ind
      sname(auto_ind)         = 'diat'
      lname(auto_ind)         = 'Diatom'
      Nfixer(auto_ind)        = .false.
      imp_calcifier(auto_ind) = .false.
      exp_calcifier(auto_ind) = .false.
      grazee_ind(auto_ind)    = auto_ind
      C_ind(auto_ind)         = diatC_ind_t
      Chl_ind(auto_ind)       = diatChl_ind_t
      Fe_ind(auto_ind)        = diatFe_ind_t
      Si_ind(auto_ind)        = diatSi_ind_t
      CaCO3_ind(auto_ind)     = 0

      ! More default parameters for diazotrophs:
      auto_ind = diaz_ind
      sname(auto_ind)         = 'diaz'
      lname(auto_ind)         = 'Diazotroph'
      Nfixer(auto_ind)        = .true.
      imp_calcifier(auto_ind) = .false.
      exp_calcifier(auto_ind) = .false.
      grazee_ind(auto_ind)    = auto_ind
      C_ind(auto_ind)         = diazC_ind_t
      Chl_ind(auto_ind)       = diazChl_ind_t
      Fe_ind(auto_ind)        = diazFe_ind_t
      Si_ind(auto_ind)        = 0
      CaCO3_ind(auto_ind)     = 0

#ifdef BEC_COCCO
      ! over-write sp information: no implicit calcifiers in this case    
      imp_calcifier(sp_ind) = .false.
      CaCO3_ind(sp_ind)     = 0

      ! Some more default parameters for coccolithophores:
      auto_ind = cocco_ind
      sname(auto_ind)         = 'cocco'
      lname(auto_ind)         = 'Coccolithophores'
      Nfixer(auto_ind)        = .false.
      imp_calcifier(auto_ind) = .false.
      exp_calcifier(auto_ind) = .true.
      grazee_ind(auto_ind)    = auto_ind
      C_ind(auto_ind)         = coccoC_ind_t
      Chl_ind(auto_ind)       = coccoChl_ind_t
      Fe_ind(auto_ind)        = coccoFe_ind_t
      Si_ind(auto_ind)        = 0
      CaCO3_ind(auto_ind)     = coccoCal_ind_t
#endif

#ifdef BIOPAR_NAMELIST
      !---------------------------------------------------------------------------
      !   read in namelist
      !---------------------------------------------------------------------------

      open(73,file=filename_biopar_bec_nml,status='old',iostat=status)
      ! Since this file is optional, it is okay to proceed if it does not exist or
      ! cannot be read.
      if (status .eq. 0) then
         read(73,NML=ecosys_parms_nml,iostat=status)
#ifdef MPI
         if (mynode.eq.0) then
#endif
            if (status .ne. 0) then
               print *, 'Could not read the biology parameter namelist!'
            else
               print *, 'Biology parameters were read from ',
     &              filename_biopar_bec_nml
            end if
#ifdef MPI
         end if  ! mynode.eq.0
#endif
         close(73)
      else  ! status .eq. 0
# ifdef MPI
         if (mynode.eq.0) then
# endif
            print *, 'Biology parameter namelist file was not found'
            print *, 'Default values will be used.'
# ifdef MPI
         end if
# endif
      end if  ! status .eq. 0
#ifndef BIOPAR_NOWRITE
# ifdef MPI
      ! Master process writes file with parameters used:
      if (mynode.eq.0) then
# endif
         open(74,file=filename_biopar_out_bec_nml,status='unknown',
     &           iostat=status)
         if (status .eq. 0) then
            write(74,NML=ecosys_parms_nml,iostat=status)
            if (status .ne. 0) then
               print *,'Could not write the biology parameter namelist!'
            end if
            close(74)
         else
            print *, 'Could not create output file ',
     &           filename_biopar_out_bec_nml
         end if
# ifdef MPI
      end if  ! (mynode.eq.0)
# endif
#endif /* BIOPAR_NOWRITE */
#endif /* BIOPAR_NAMELIST */


#ifdef BEC2_DIAG
      !---------------------------------------------------------------------------
      !   Initialize diagnostic variables:
      !---------------------------------------------------------------------------
# ifdef BEC2_DIAG_USER
      ! Allocate arrays for diagnostics:
      allocate(bec2_diag_2d(GLOBAL_2D_ARRAY,nr_bec2_diag_2d_user))
      allocate(bec2_diag_3d(GLOBAL_2D_ARRAY,N,nr_bec2_diag_3d_user))
#  ifdef AVERAGES
      allocate(bec2_diag_2d_avg(GLOBAL_2D_ARRAY,nr_bec2_diag_2d_user))
      allocate(bec2_diag_3d_avg(GLOBAL_2D_ARRAY,N,nr_bec2_diag_3d_user))
#   ifdef SLICE_AVG
      allocate(bec2_diag_2d_slavg(GLOBAL_2D_ARRAY,nr_bec2_diag_2d_user))
      allocate(bec2_diag_3d_slavg(GLOBAL_2D_ARRAY,nr_bec2_diag_3d_user))
#   endif
#  endif /* AVERAGES */
      MPI_master_only write (*,'(a,i3,a)') 'bec2_diag_2d allocated for ',
     &    nr_bec2_diag_2d_user, ' variables'
      MPI_master_only write (*,'(a,i3,a)') 'bec2_diag_3d allocated for ',
     &    nr_bec2_diag_3d_user, ' variables'
# endif /* BEC2_DIAG_USER */
      bec2_diag_2d = c0
      bec2_diag_3d = c0
#else /* BEC2_DIAG */
      ph_hist = c0
#endif /* BEC2_DIAG */

      !---------------------------------------------------------------------------
      !   Initialize ice fraction and atm. pressure field:
      !---------------------------------------------------------------------------
      ifrac = 0.0
      press = 1.

      lflux_gas_o2  = .TRUE.
# ifdef Ncycle_SY
      lflux_gas_n2o  = .TRUE.
      lflux_gas_n2  = .TRUE.
# endif
      lflux_gas_co2 = .TRUE.
      liron_flux = .TRUE.
#ifdef RIVER_LOAD_BIO
      lriver_load_n = .TRUE.
      lriver_load_p = .TRUE.
      lriver_load_alk = .TRUE.
      lriver_load_dic = .TRUE.
      lriver_load_si = .TRUE.
#endif
      lsource_sink  = .TRUE.  

      where(rmask==1)
         landmask=.TRUE.
      elsewhere 
         landmask=.FALSE.
      endwhere

#ifdef CCHEM_MOCSY
       optcon  = 'mol/m3'   ! input concentrations of DIC, alk, P, Si are in MOL/M3
       optt    = 'Tpot'     ! 'temp' is POTENTIAL temp [C]
       optp    = 'm'        ! input variable z_r is in METERS
       optb    = 'u74'      ! Uppstrom (1974) formulation for total boron
       optk1k2 = 'm95'      ! Millero (1995) formulations for k1 and k2
       optkf   = 'pf'       ! Perez & Fraga (1987) formulation for Kf
       optgas  = 'Pinsitu'  ! compute in situ values for fCO2 and pCO2 (i.e. accounting for pressure)
#endif


      end subroutine ecosys2_init

#else /* BIOLOGY_BEC2 */
      subroutine ecosys2_init()
      return 
      end
#endif /* BIOLOGY_BEC2 */
