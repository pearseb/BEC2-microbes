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
      parm_kao = 0.117 * dps ! max ammonium oxidation rate (mmol/m3/s)
      parm_kno = 0.05 * dps ! max nitrite oxidation rate (1/s)
      parm_ko2_ao = 0.333       ! Michaelis Menton O2 constant for ammonium oxidation (mmol m-3)
      parm_knh4_ao = 0.13    ! Michaelis Menton NH4 constant for ammonium oxidation (mmol m-3)
      parm_ko2_no = 0.778     ! Michaelis Menton O2 constant for nitrite oxidation (mmol m-3)
      parm_kno2_no = 1.0    ! Michaelis Menton NO2 constant for nitrite oxidation (mmol m-3)
      parm_kno3_den1 = 0.4      ! From optimization in 1D model, still needs some work
      parm_kno2_den2 = 0.05
      parm_kn2o_den3 = 0.2
      parm_ko2_oxic = 0.01
      parm_ko2_den1 = 5.0
      parm_ko2_den2 = 2.3
      parm_ko2_den3 = 0.11
      parm_koxic = 0.08 * dps
      parm_kden1 = 0.0215 * dps
      parm_kden2 = 0.008 * dps
      parm_kden3 = 0.0455 * dps
      parm_kax = 0.39 * dps
      parm_knh4_ax = 0.230
      parm_kno2_ax = 0.1
      parm_ko2_ax = 0.7
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
      kNH4 = (/ 0.0075, 0.06, 0.2 /)         ! mmol N/m3  ! 0.01 0.04  0.8
      kSiO3 = (/ 0.0, 0.6, 0.0 /)           ! mmol SiO3/m3
      Qp = (/ 0.008547, 0.008547, 0.008547 /)
      gQfe_0 = (/ 33.0e-6, 33.0e-6, 66.0e-6 /)
      gQfe_min = (/ 2.7e-6, 2.7e-6, 6.0e-6 /)
      alphaPI = (/ 0.4*dps, 0.31*dps, 0.33*dps /)  ! previous defaults, Chl. specific initial slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
      PCref = (/ 3.54*dps, 3.54*dps, 0.51*dps /) ! MF: sensitivity tests70%
      thetaN_max = (/ 2.5, 4.0, 2.5 /)
      loss_thres = (/ 0.02, 0.02, 0.02 /)
      loss_thres2 = (/ 0.0, 0.0, 0.001 /)
      temp_thres = (/ -10.0, -10.0, 16.0 /)
      mort = (/ 0.1*dps, 0.1*dps, 0.1*dps /)  ! non-grazing mortality (1/sec)
      mort2 = (/ 0.01*dps, 0.01*dps, 0.01*dps /)
      agg_rate_max = (/ 0.5, 0.5, 0.5 /)
      agg_rate_min = (/ 0.01, 0.02, 0.01 /)
      z_umax_0 = (/ 3.3*dps, 3.05*dps, 3.3*dps /) ! MF: sensitivity test 75%
      z_grz = (/ 1.2, 1.2, 1.2 /)      ! grazing coefficient (mmol C/m^3)
      graze_zoo = (/ 0.3, 0.25, 0.3 /)
      graze_poc = (/ 0.0, 0.4, 0.1 /)
      graze_doc = (/ 0.06, 0.06, 0.06 /)
      loss_poc = (/ 0.0, 0.0, 0.0 /)
      f_zoo_detr = (/ 0.1, 0.2, 0.1 /)


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
