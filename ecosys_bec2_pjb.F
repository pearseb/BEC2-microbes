#include "cppdefs.h"
#ifdef BIOLOGY_BEC2

       subroutine ecosys_bec2_tile(Istr,Iend,Jstr,Jend)

       implicit none

#include "param.h"  
#include "param_bec2.h"
#include "ecosys_bec2.h"
#include "ocean3d.h"
#include "forces.h"
#include "scalars.h"
#include "grid.h"

       integer i,j,k,m
       integer Istr,Iend,Jstr,Jend
       integer su_i1, su_i2, sv_j1, sv_j2

       real
     &    O2SAT_USED    ! used O2 saturation (mmol/m^3)
#if defined Ncycle_SY || defined N2O_NEV
     &    , N2OSAT_USED    ! used O2 saturation (mmol/m^3)
# endif
#ifdef Ncycle_SY
     &    ,N2SAT_USED    ! used N2 saturation (mmol/m^3)
# endif

       real,dimension(istr:iend,jstr:jend)::
     &    SSTT,   ! sea surface temperature (C)
     &    SSSS,    ! sea surface salinity (psu)
     &    SHF_QSW,             ! penetrative solar heat flux (W/m^2)
     &    QA_dust_def,       ! incoming deficit in the QA(dust) POC flux
     &    PAR_out,           ! photosynthetically available radiation (W/m^2)
     &    FICE_USED,   ! used ice fraction (non-dimensional)
     &    WS_USED,     ! used wind speed (m/s) ! HF: changed from (cm/s)
     &    XKW,         ! a * WS ** 2 (m/s) ! HF: changed from (cm/s)
     &    AP_USED,     ! used atmospheric pressure (atm)
     &    XKW_ICE,     ! common portion of piston vel., a*(1-fice)*u**2 (m/s) ! HF: changed from (cm/s)
     &    SCHMIDT_USED, ! used Schmidt number
     &    PV,           ! piston velocity (m/s) ! HF: changed from (cm/s)
     &    O2SAT_1atm,   ! O2 saturation @ 1 atm (mmol/m^3)
     &    FLUX,         ! tracer flux (mmol/m^2/s)
     &    sustr_rho,    ! wind stress in xi direction on rho grid (N/m^2)
     &    svstr_rho     ! wind stress in xi direction on rho grid (N/m^2)

       real STF(istr:iend,jstr:jend,ntrc_bio)

# ifdef CCHEM_MOCSY
       real ff, tk100, tk1002, co2star_air, mass_to_vol
       real,dimension(istr:iend)::
     &    fco2,         ! CO2 fugacity [uatm]
     &    hco3,         ! bicarbonate ion (HCO3-) concentration [mol/m3]
     &    co3,          ! carbonate ion (CO3-) concentration [mol/m3]
     &    omega_arag,   ! the aragonite saturation state
     &    omega_calc,   ! the calcite saturation state
     &    betad,        ! Revelle factor, i.e., dpCO2/pCO2 / dDIC/DIC
     &    rhosw,        ! in-situ density of seawater; rhoSW = f(s, t, p) in [kg/m3]
     &    p_ocean,      ! pressure [decibars]; p = f(depth, latitude)
     &    temp_is       ! in-situ temperature [degrees C]
#  ifdef CCHEM_TODEPTH
       real, dimension(istr:iend,jstr:jend,N)::  salt_depth, temp_depth
#  endif
# else  /* CCHEM_MOCSY */
       real,dimension(istr:iend)::
     &    PHLO          ! lower bound for ph in solver
     &    ,PHHI         ! upper bound for ph in solver
       real, PARAMETER :: 
     &    phlo_init = 7.0,  ! low bound for ph for no prev soln
     &    phhi_init = 9.0,  ! high bound for ph for no prev soln
     &    del_ph = 0.20     ! delta-ph for prev soln
# endif  /* CCHEM_MOCSY */

       real, PARAMETER :: 
     &    a = 8.6e-7       ! a = 0.31 cm/hr s^2/m^2 in (s/m)
       real,dimension(istr:iend)::
     &    PH_NEW,       ! computed PH from solver
     &    XCO2,         ! atmospheric co2 conc. (dry-air, 1 atm)
     &    CO2STAR_ROW,  ! CO2STAR from solver
     &    DCO2STAR_ROW, ! DCO2STAR from solver
     &    pco2oc_ROW,   ! pco2oc from solver
     &    DpCO2_ROW     ! DpCO2 from solver

# if defined Ncycle_SY || defined N2O_NEV
       real,dimension(istr:iend,jstr:jend)::
     &    N2OSAT_1atm   ! N2O saturation @ 1 atm (mmol/m^3)

       real,dimension(istr:iend)::
     &    XN2O         ! atmospheric n2o conc. (dry-air, 1 atm)
# endif
# if defined Ncycle_SY
       real,dimension(istr:iend,jstr:jend)::
     &    N2SAT_1atm   ! N2 saturation @ 1 atm (mmol/m^3)
# endif
!
! Interpolate/extrapolate wind stress to rho grid:
!
       do j=jstr,jend
          if (j == jstr) then
             sv_j1 = j+1
          else
             sv_j1 = j
          endif
          if (j .lt. jend) then
             sv_j2 = j+1
          else
             sv_j2 = j-1
          end if
          do i=istr,iend
             if (i == istr) then
                su_i1 = i+1
             else
                su_i1 = i
             end if
             if (i .lt. iend) then
                su_i2 = i+1
             else
                su_i2 = i-1
             end if
             sustr_rho(i,j) = 0.5*(sustr(su_i1,j)+sustr(su_i2,j))
             svstr_rho(i,j) = 0.5*(svstr(i,sv_j1)+svstr(i,sv_j2))
          end do
       end do

!
! Initialize local copies of tracer array:
!
       do m=1,ntrc_bio
          do k=1,N
             tracer(istr:iend,jstr:jend,k,m) =
     &              t(istr:iend,jstr:jend,k,nnew,1+ntrc_salt+ntrc_pas+m)             
          enddo 
       enddo
       sstt(istr:iend,jstr:jend) = t(istr:iend,jstr:jend,N,nnew,1)
       ssss(istr:iend,jstr:jend) = t(istr:iend,jstr:jend,N,nnew,2)
# if defined CCHEM_MOCSY && defined CCHEM_TODEPTH
       do k=1,N
          do j=jstr,jend
             do i=istr,iend
                temp_depth(i,j,k) = t(i,j,k,nnew,1)
                salt_depth(i,j,k) = t(i,j,k,nnew,2)
             end do
          end do
       end do
# endif /* CCHEM_MOCSY && CCHEM_TODEPTH */

!
! Get short wave radiation:
!
#ifdef DAILYPAR_BEC
       shf_qsw(istr:iend,jstr:jend) =
     &               srflx_dailyavg(istr:iend,jstr:jend)*rho0*Cp
#else
       shf_qsw(istr:iend,jstr:jend)=srflx(istr:iend,jstr:jend)*rho0*Cp
#endif


! ----------------- !
! 
!    Gas exchange   !
! ----------------- !
!                   !
!  i)   get ice concentration
!  ii)  get wind stress 
!  iii) N2O
!  iv)  N2
!  v)   O2
!  vi)  CO2


#ifdef ICEOBS
# define IFRAC sic /* Ice fraction has different name in this case */
#endif
       FICE_USED(istr:iend,jstr:jend) = IFRAC(istr:iend,jstr:jend)
       call WS(sustr_rho(istr:iend,jstr:jend),
     &    svstr_rho(istr:iend,jstr:jend),landmask(istr:iend,jstr:jend),
     &    ws_used(istr:iend,jstr:jend), istr,iend,jstr,jend)

       XKW(istr:iend,jstr:jend) = a * WS_USED(istr:iend,jstr:jend)*WS_USED(istr:iend,jstr:jend)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_2d_l(ws10m_idx_t))
     &     bec2_diag_2d(istr:iend,jstr:jend,bec2_diag_2d_idx(ws10m_idx_t)) = WS_USED(istr:iend,jstr:jend)
#  else
       bec2_diag_2d(istr:iend,jstr:jend,ws10m_idx_t) = WS_USED(istr:iend,jstr:jend)
#  endif
# endif

    !---------------------------------------------------------------------
    !   convert PRESS from Pascals to atm,no, it is atm already
    !----------------------------------------------------------------------
       AP_USED(istr:iend,jstr:jend) = PRESS(istr:iend,jstr:jend)

    !----------------------------------------------------------------------
    !   Set bad AP values to 1. This is necessary for runs restarting off
    !   a run in which the flux coupler didnt restart on AP correctly.
    !----------------------------------------------------------------------
       WHERE (AP_USED(istr:iend,jstr:jend) > 1.5 .OR. AP_USED(istr:iend,jstr:jend) < 0.5)
          AP_USED(istr:iend,jstr:jend) = c1
       END WHERE

# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_2d_l(xkw_idx_t))
     &     bec2_diag_2d(istr:iend,jstr:jend,bec2_diag_2d_idx(xkw_idx_t)) = XKW(istr:iend,jstr:jend)
       if (bec2_diag_2d_l(atmpress_idx_t))
     &     bec2_diag_2d(istr:iend,jstr:jend,bec2_diag_2d_idx(atmpress_idx_t)) = AP_USED(istr:iend,jstr:jend)
#  else
       bec2_diag_2d(istr:iend,jstr:jend,xkw_idx_t) = XKW(istr:iend,jstr:jend)
       bec2_diag_2d(istr:iend,jstr:jend,atmpress_idx_t) = AP_USED(istr:iend,jstr:jend)
#  endif
# endif

    !------------------------------------------------------------------------
    !   Compute XKW_ICE. XKW is zero over land, so XKW_ICE is too.
    !-------------------------------------------------------------------------

       IF (lflux_gas_o2 .OR. lflux_gas_co2) THEN
         XKW_ICE(istr:iend,jstr:jend) = XKW(istr:iend,jstr:jend)
         WHERE (FICE_USED(istr:iend,jstr:jend) > 0.2 .AND. FICE_USED(istr:iend,jstr:jend) < 0.9999)
            XKW_ICE(istr:iend,jstr:jend) = (c1 - FICE_USED(istr:iend,jstr:jend)) * XKW_ICE(istr:iend,jstr:jend)
         END WHERE
         WHERE (FICE_USED(istr:iend,jstr:jend) >= 0.9999)
           XKW_ICE(istr:iend,jstr:jend) = c0
         END WHERE
       END IF

#if defined Ncycle_SY || defined N2O_NEV
     !------------------------------------------------------------------------
     !   compute N2O flux
     !------------------------------------------------------------------------
        lflux_gas_n2o = .true.
       IF (lflux_gas_n2o) THEN
         call CSCHMIDT_N2O(SSTT(istr:iend,jstr:jend),
     &     landmask(istr:iend,jstr:jend),SCHMIDT_USED(istr:iend,jstr:jend),
     &     istr,iend,jstr,jend)


#if defined BEC2_DIAG && defined Ncycle_SY 
         bec2_diag_2d(istr:iend,jstr:jend,schmidt_n2o_idx_t) = SCHMIDT_USED(istr:iend,jstr:jend)
#endif

         call N2OSATU(SSTT(istr:iend,jstr:jend),
     &     SSSS(istr:iend,jstr:jend),landmask(istr:iend,jstr:jend),
     &     N2OSAT_1atm(istr:iend,jstr:jend),istr, iend,jstr,jend)
         do j=jstr,jend
            do i=istr,iend
               PV(i,j) = XKW_ICE(i,j) * SQRT(660.0 / SCHMIDT_USED(i,j))
# ifdef MASKING
     &                      * rmask(i,j)
# endif
        XN2O = 300.0 * 1e-9 !!! need to figure out namelists here
               N2OSAT_USED = AP_USED(i,j) * N2OSAT_1atm(i,j) * XN2O(i)
# ifdef MASKING
     &                      * rmask(i,j)
# endif

               !!! calculating gas exchange for all n2o tracer
               !!! we let n2o_atm trace the invasion of atmospheric n2o
               !!! and set the atmopspheric concentration of all other
               !!! tracer to 0 So that everything is additive.
# ifdef Ncycle_SY
               FLUX(i,j) = PV(i,j) * (N2OSAT_USED - tracer(i,j,N,n2o_ind_t))
               STF(i,j,n2o_ind_t) = FLUX(i,j)
# endif
# ifdef N2O_TRACER_DECOMP
               FLUX(i,j) = PV(i,j) * (c0 - tracer(i,j,N,n2o_ao1_ind_t))
               STF(i,j,n2o_ao1_ind_t) = FLUX(i,j)
               FLUX(i,j) = PV(i,j) * (c0 - tracer(i,j,N,n2o_soden_ind_t))
               STF(i,j,n2o_soden_ind_t) = FLUX(i,j)
               FLUX(i,j) = PV(i,j) * (c0 - tracer(i,j,N,n2o_siden_ind_t))
               STF(i,j,n2o_siden_ind_t) = FLUX(i,j)
               FLUX(i,j) = PV(i,j) * (N2OSAT_USED - tracer(i,j,N,n2o_atm_ind_t))
               STF(i,j,n2o_atm_ind_t) = FLUX(i,j)
# endif
# ifdef N2O_NEV 
               FLUX(i,j) = PV(i,j) * (N2OSAT_USED - tracer(i,j,N,n2o_nev_ind_t))
               STF(i,j,n2o_nev_ind_t) = FLUX(i,j)
#endif
               !!! Update tracer
# ifdef Ncycle_SY
               tracer(i,j,N,n2o_ind_t) = tracer(i,j,N,n2o_ind_t)+
     &                 stf(i,j,n2o_ind_t)*dt/Hz(i,j,N)
#endif

# ifdef N2O_TRACER_DECOMP
               tracer(i,j,N,n2o_ao1_ind_t) = tracer(i,j,N,n2o_ao1_ind_t)+
     &                 stf(i,j,n2o_ao1_ind_t)*dt/Hz(i,j,N)
               tracer(i,j,N,n2o_soden_ind_t) = tracer(i,j,N,n2o_soden_ind_t)+
     &                 stf(i,j,n2o_soden_ind_t)*dt/Hz(i,j,N)
               tracer(i,j,N,n2o_siden_ind_t) = tracer(i,j,N,n2o_siden_ind_t)+
     &                 stf(i,j,n2o_siden_ind_t)*dt/Hz(i,j,N)
               tracer(i,j,N,n2o_atm_ind_t) = tracer(i,j,N,n2o_atm_ind_t)+
     &                 stf(i,j,n2o_atm_ind_t)*dt/Hz(i,j,N)
# endif
# ifdef N2O_NEV
               tracer(i,j,N,n2o_nev_ind_t) = tracer(i,j,N,n2o_nev_ind_t)+
     &                 stf(i,j,n2o_nev_ind_t)*dt/Hz(i,j,N)
# endif

# if defined BEC2_DIAG && defined Ncycle_SY
               bec2_diag_2d(i,j,fgn2o_idx_t) = STF(i,j,n2o_ind_t)
               bec2_diag_2d(i,j,pvn2o_idx_t) = PV(i,j)
               bec2_diag_2d(i,j,n2osat_idx_t) = N2OSAT_USED
#ifdef N2O_TRACER_DECOMP
               bec2_diag_2d(i,j,fgn2o_ao1_idx_t) = STF(i,j,n2o_ao1_ind_t)
               bec2_diag_2d(i,j,fgn2o_siden_idx_t) = STF(i,j,n2o_siden_ind_t)
               bec2_diag_2d(i,j,fgn2o_soden_idx_t) = STF(i,j,n2o_soden_ind_t)
               bec2_diag_2d(i,j,fgn2o_atm_idx_t) = STF(i,j,n2o_atm_ind_t)
# endif
#endif

# ifdef BEC2_DIAG
#ifdef N2O_NEV
               bec2_diag_2d(i,j,fgn2o_nev_idx_t) = STF(i,j,n2o_nev_ind_t)
#endif
#endif
            end do  ! i
         end do  ! j
       END IF
# endif

# ifdef Ncycle_SY 
    !------------------------------------------------------------------------
    !   compute N2 flux
    !------------------------------------------------------------------------
        lflux_gas_n2 = .true.
       IF (lflux_gas_n2) THEN
         call CSCHMIDT_N2(SSTT(istr:iend,jstr:jend),
     &     landmask(istr:iend,jstr:jend),SCHMIDT_USED(istr:iend,jstr:jend),
     &     istr,iend,jstr,jend)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
         if (bec2_diag_2d_l(schmidt_n2_idx_t))
     &       bec2_diag_2d(istr:iend,jstr:jend,bec2_diag_2d_idx(schmidtn2_idx_t)) = SCHMIDT_USED(istr:iend,jstr:jend)
#  else
         bec2_diag_2d(istr:iend,jstr:jend,schmidt_n2_idx_t) = SCHMIDT_USED(istr:iend,jstr:jend)
#  endif
# endif
         call N2SATU(SSTT(istr:iend,jstr:jend),
     &     SSSS(istr:iend,jstr:jend),landmask(istr:iend,jstr:jend),
     &     N2SAT_1atm(istr:iend,jstr:jend),istr, iend,jstr,jend)
         do j=jstr,jend
            do i=istr,iend
               PV(i,j) = XKW_ICE(i,j) * SQRT(660.0 / SCHMIDT_USED(i,j))
# ifdef MASKING
     &                      * rmask(i,j)
# endif
               N2SAT_USED = AP_USED(i,j) * N2SAT_1atm(i,j) * 0.78084 !Molar ratio of N2 in air
# ifdef MASKING
     &                      * rmask(i,j)
# endif
!!!!!---------------Tracking only N2 excess for no2---------!!!!! to track full use:
!                FLUX(i,j) = PV(i,j) * (N2SAT_USED - tracer(i,j,N,n2_ind_t)) 
               FLUX(i,j) = PV(i,j) * (c0 - tracer(i,j,N,n2_ind_t)) 
               STF(i,j,n2_ind_t) = FLUX(i,j)
               tracer(i,j,N,n2_ind_t) = tracer(i,j,N,n2_ind_t)+
     &                 stf(i,j,n2_ind_t)*dt/Hz(i,j,N)
# ifdef BEC2_DIAG
               bec2_diag_2d(i,j,pvn2_idx_t) = PV(i,j)
               bec2_diag_2d(i,j,n2sat_idx_t) = N2SAT_USED
               bec2_diag_2d(i,j,fgn2_idx_t) = FLUX(i,j)
# endif
            end do  ! i
         end do  ! j
       END IF

# endif
    !------------------------------------------------------------------------
    !   compute O2 flux
    !------------------------------------------------------------------------

       IF (lflux_gas_o2) THEN
         call CSCHMIDT_O2(SSTT(istr:iend,jstr:jend),
     &     landmask(istr:iend,jstr:jend),SCHMIDT_USED(istr:iend,jstr:jend),
     &     istr,iend,jstr,jend)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
         if (bec2_diag_2d_l(schmidto2_idx_t))
     &       bec2_diag_2d(istr:iend,jstr:jend,bec2_diag_2d_idx(schmidto2_idx_t)) = SCHMIDT_USED(istr:iend,jstr:jend)
#  else
         bec2_diag_2d(istr:iend,jstr:jend,schmidto2_idx_t) = SCHMIDT_USED(istr:iend,jstr:jend)
#  endif
# endif
         call O2SATU(SSTT(istr:iend,jstr:jend),
     &     SSSS(istr:iend,jstr:jend),landmask(istr:iend,jstr:jend),
     &     O2SAT_1atm(istr:iend,jstr:jend),istr, iend,jstr,jend)
         do j=jstr,jend
            do i=istr,iend
               PV(i,j) = XKW_ICE(i,j) * SQRT(660.0 / SCHMIDT_USED(i,j))
# ifdef MASKING
     &                      * rmask(i,j)
# endif
               O2SAT_USED = AP_USED(i,j) * O2SAT_1atm(i,j)
# ifdef MASKING
     &                      * rmask(i,j)
# endif
               FLUX(i,j) = PV(i,j) * (O2SAT_USED - tracer(i,j,N,o2_ind_t))
               STF(i,j,o2_ind_t) = FLUX(i,j)
               tracer(i,j,N,o2_ind_t) = tracer(i,j,N,o2_ind_t)+
     &                 stf(i,j,o2_ind_t)*dt/Hz(i,j,N)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
               if (bec2_diag_2d_l(pvo2_idx_t))
     &             bec2_diag_2d(i,j,bec2_diag_2d_idx(pvo2_idx_t)) = PV(i,j)
               if (bec2_diag_2d_l(o2sat_idx_t))
     &             bec2_diag_2d(i,j,bec2_diag_2d_idx(o2sat_idx_t)) = O2SAT_USED
               if (bec2_diag_2d_l(fgo2_idx_t))
     &             bec2_diag_2d(i,j,bec2_diag_2d_idx(fgo2_idx_t)) = FLUX(i,j)
#  else
               bec2_diag_2d(i,j,pvo2_idx_t) = PV(i,j)
               bec2_diag_2d(i,j,o2sat_idx_t) = O2SAT_USED
               bec2_diag_2d(i,j,fgo2_idx_t) = FLUX(i,j)
#  endif
# endif
            end do  ! i
         end do  ! j
       END IF

    !------------------------------------------------------------------------
    !   compute CO2 flux, computing disequilibrium one row at a time
    !------------------------------------------------------------------------
       IF (lflux_gas_co2) THEN
          call CSCHMIDT_CO2(SSTT(istr:iend,jstr:jend),
     &       landmask(istr:iend,jstr:jend),
     &       SCHMIDT_USED(istr:iend,jstr:jend),
     &       istr,iend,jstr,jend)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
          if (bec2_diag_2d_l(schmidtco2_idx_t))
     &        bec2_diag_2d(istr:iend,jstr:jend,bec2_diag_2d_idx(schmidtco2_idx_t)) = SCHMIDT_USED(istr:iend,jstr:jend)
#  else
          bec2_diag_2d(istr:iend,jstr:jend,schmidtco2_idx_t) = SCHMIDT_USED(istr:iend,jstr:jend)
#  endif
# endif

          do j=jstr,jend
            do i=istr,iend
               PV(i,j) = XKW_ICE(i,j) * SQRT(660.0 / SCHMIDT_USED(i,j))
# ifdef MASKING
     &                      * rmask(i,j)
# endif
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
               if (bec2_diag_2d_l(pvco2_idx_t))
     &             bec2_diag_2d(i,j,bec2_diag_2d_idx(pvco2_idx_t)) = PV(i,j)
#  else
               bec2_diag_2d(i,j,pvco2_idx_t) = PV(i,j)
#  endif
# endif
            end do
          end do

! Set a default value for pCO2 in the air. Will be overwritten if
! VARIABLE_ANN_ATM_PCO2 and/or VARIABLE_ATM_PCO2 is defined:
          XCO2 = 370.0           
! OCMIP-2 mandated value (278 ppm, dry-air, 1 atm)
!HF           XCO2 = 278.0
# ifdef VARIABLE_ANN_ATM_PCO2
          ! Compute number of simulated days since the start of the simulation:
          tdays_pco2 = time*sec2day
          ! Compute where in which year we are:
          curr_year = start_year + get_frac_year(tdays_pco2)
          ! Determine atm pco2:
#  ifdef HINDCAST
          call get_atmco2(curr_year,atm_co2)
#  elif defined FUTURE_SCENARIO
          call get_atmco2(curr_year,futr_scen,atm_co2)
#  endif /* HINDCAST */
          XCO2 = atm_co2           
# endif
# ifdef VARIABLE_ATM_PCO2
      ! the factor deg2rad converts 360d into one full cycle of 2 pi
      ! the subtraction of 30d ensures that the maximum is reached
      ! at the end of April and the minimum at the end of October
          XCO2 = XCO2 + 2.9 * sin(deg2rad*(tdays-30.))
# endif /* VARIABLE_ATM_PCO2 */

          DO j = jstr,jend

! MF (added treatment of pCO2air forcing)
# ifdef PCO2AIR_FORCING
             XCO2(istr:iend) = pco2air(istr:iend,j)
# endif /* PCO2AIR_FORCING */

# ifdef CCHEM_MOCSY
             call vars(ph_new(istr:iend), pco2oc_ROW, fco2, co2star_row,  ! OUTPUT
     &          hco3, co3, omega_arag, omega_calc, betad, rhoSW, p_ocean, temp_is,           ! OUTPUT
     &          sstt(istr:iend,j), ssss(istr:iend,j), 0.001*tracer(istr:iend,j,N,alk_ind_t),       ! INPUT
     &          0.001*tracer(istr:iend,j,N,dic_ind_t), 0.001*tracer(istr:iend,j,N,sio3_ind_t),     ! INPUT
     &          0.001*tracer(istr:iend,j,N,po4_ind_t), AP_USED(istr:iend,j), -z_r(istr:iend,j,N),  ! INPUT
     &          latr(istr:iend,j), iend-istr+1, optCON, optT, optP, optB, optk1k2, optkf, optgas)  ! INPUT
             ! Unit conversion: mol/m3 -> mmol/m3
             co2star_row(istr:iend) = 1000 * co2star_row(istr:iend)
             ! Compute CO2star for air and DCO2star:
             do i = istr, iend
                ! For converting (mol/kg) -> (mmol/m^3). Note that rhoSW is kg/m^3.
                mass_to_vol = 1e3 * rhoSW(i)
                !---------------------------------------------------------------------
                !   f = k0(1-pH2O)*correction term for non-ideality
                !   Weiss  Price (1980, Mar. Chem., 8, 347-359;
                !                 Eq 13 with table 6 values)
                !---------------------------------------------------------------------
                tk100 = (T0_Kelvin + sstt(i,j)) * 1e-2
                tk1002 = tk100 * tk100
                ff = EXP(-162.8301 + 218.2968/tk100 + 
     &             90.9241*LOG(tk100) - 
     &             1.47696*tk1002 + ssss(i,j)*(.025695 - .025225*tk100 + 
     &             0.0049867*tk1002))
                co2star_air = (xco2(i)*1e-6) * ff * AP_USED(i,j) ! mol/kg
                co2star_air = co2star_air * mass_to_vol ! mmol/m^3
                DCO2STAR_ROW(i) = co2star_air - co2star_row(i) ! mmol/m^3
             end do
             ! Compute delta-pCO2:
             DpCO2_ROW = pco2oc_ROW - XCO2
#  ifdef BEC2_DIAG
#   ifdef CCHEM_TODEPTH
#    ifdef BEC2_DIAG_USER
             ! Multiply by 1000 to convert Mol/m3 -> mMol/m3:
             if (bec2_diag_3d_l(co3_idx_t))
     &           bec2_diag_3d(istr:iend,j,N,bec2_diag_3d_idx(co3_idx_t)) = 1000*co3
             if (bec2_diag_3d_l(hco3_idx_t))
     &           bec2_diag_3d(istr:iend,j,N,bec2_diag_3d_idx(hco3_idx_t)) = 1000*hco3
             if (bec2_diag_3d_l(co2star_idx_t))
     &           bec2_diag_3d(istr:iend,j,N,bec2_diag_3d_idx(co2star_idx_t)) = CO2STAR_ROW(istr:iend)
             if (bec2_diag_3d_l(pco2oc_idx_t))
     &           bec2_diag_3d(istr:iend,j,N,bec2_diag_3d_idx(pco2oc_idx_t)) = pco2oc_ROW(istr:iend)
             if (bec2_diag_3d_l(ph_idx_t))
     &           bec2_diag_3d(istr:iend,j,N,bec2_diag_3d_idx(ph_idx_t)) = ph_new(istr:iend)
#    else
             ! Multiply by 1000 to convert Mol/m3 -> mMol/m3:
             bec2_diag_3d(istr:iend,j,N,co3_idx_t) = 1000*co3
             bec2_diag_3d(istr:iend,j,N,hco3_idx_t) = 1000*hco3
             bec2_diag_3d(istr:iend,j,N,co2star_idx_t) = CO2STAR_ROW(istr:iend)
             bec2_diag_3d(istr:iend,j,N,pco2oc_idx_t) = pco2oc_ROW(istr:iend)
             bec2_diag_3d(istr:iend,j,N,ph_idx_t) = ph_new(istr:iend)
#    endif
#   else /* CCHEM_TODEPTH */
#    ifdef BEC2_DIAG_USER
             ! Multiply by 1000 to convert Mol/m3 -> mMol/m3:
             if (bec2_diag_2d_l(co3_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(co3_idx_t)) = 1000*co3
             if (bec2_diag_2d_l(hco3_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(hco3_idx_t)) = 1000*hco3
             if (bec2_diag_2d_l(co2star_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(co2star_idx_t)) = CO2STAR_ROW(istr:iend)
             if (bec2_diag_2d_l(pco2oc_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(pco2oc_idx_t)) = pco2oc_ROW(istr:iend)
             if (bec2_diag_2d_l(ph_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(ph_idx_t)) = ph_new(istr:iend)
#    else
             ! Multiply by 1000 to convert Mol/m3 -> mMol/m3:
             bec2_diag_2d(istr:iend,j,co3_idx_t) = 1000*co3
             bec2_diag_2d(istr:iend,j,hco3_idx_t) = 1000*hco3
             bec2_diag_2d(istr:iend,j,co2star_idx_t) = CO2STAR_ROW(istr:iend)
             bec2_diag_2d(istr:iend,j,pco2oc_idx_t) = pco2oc_ROW(istr:iend)
             bec2_diag_2d(istr:iend,j,ph_idx_t) = ph_new(istr:iend)
#    endif
#   endif /* CCHEM_TODEPTH */
#   ifdef BEC2_DIAG_USER
             if (bec2_diag_2d_l(pco2air_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(pco2air_idx_t)) = XCO2
             if (bec2_diag_2d_l(dco2star_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(dco2star_idx_t)) = DCO2STAR_ROW(istr:iend)
#   else
             bec2_diag_2d(istr:iend,j,pco2air_idx_t) = XCO2
             bec2_diag_2d(istr:iend,j,dco2star_idx_t) = DCO2STAR_ROW(istr:iend)
#   endif
#  endif /* BEC2_DIAG */
# else /* CCHEM_MOCSY */
#  ifdef BEC2_DIAG
#   ifdef BEC2_DIAG_USER
             ! pH is present in bec2_diag_2d, so bec2_diag_2d_idx(ph_idx_t)
             ! is defined:
             WHERE (BEC2_DIAG_2D(istr:iend,j,bec2_diag_2d_idx(ph_idx_t)) .ne. c0)
                PHLO(istr:iend) = BEC2_DIAG_2D(istr:iend,j,bec2_diag_2d_idx(ph_idx_t)) - del_ph
                PHHI(istr:iend) = BEC2_DIAG_2D(istr:iend,j,bec2_diag_2d_idx(ph_idx_t)) + del_ph
#   else
             WHERE (BEC2_DIAG_2D(istr:iend,j,ph_idx_t) .ne. c0)
                PHLO(istr:iend) = BEC2_DIAG_2D(istr:iend,j,ph_idx_t) - del_ph
                PHHI(istr:iend) = BEC2_DIAG_2D(istr:iend,j,ph_idx_t) + del_ph
#   endif
#  else /* BEC2_DIAG */
             WHERE (PH_HIST(istr:iend,j) .ne. c0)
                PHLO(istr:iend) = PH_HIST(istr:iend,j) - del_ph
                PHHI(istr:iend) = PH_HIST(istr:iend,j) + del_ph
#  endif /* BEC2_DIAG */ 
             ELSEWHERE
                PHLO = phlo_init
                PHHI = phhi_init
             END WHERE

             CALL co2calc_row(LANDMASK(istr:iend,j),
     &         SSTT(istr:iend,j), SSSS(istr:iend,j),
     &         tracer(istr:iend,j,N,dic_ind_t),
     &         tracer(istr:iend,j,N,alk_ind_t),
     &         tracer(istr:iend,j,N,po4_ind_t),
     &         tracer(istr:iend,j,N,sio3_ind_t),
     &         PHLO(istr:iend), PHHI(istr:iend), PH_NEW(istr:iend), XCO2(istr:iend),
     &         AP_USED(istr:iend,j), CO2STAR_ROW(istr:iend),
     &         DCO2STAR_ROW(istr:iend), pco2oc_ROW(istr:iend),
     &         DpCO2_ROW(istr:iend),istr,iend)
#  ifdef BEC2_DIAG
#   ifdef BEC2_DIAG_USER
             ! bec2_diag_2d_idx(ph_idx_t) is always defined:
             bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(ph_idx_t)) = PH_NEW(istr:iend)
             if (bec2_diag_2d_l(co2star_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(co2star_idx_t))  = CO2STAR_ROW(istr:iend)
             if (bec2_diag_2d_l(dco2star_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(dco2star_idx_t)) = DCO2STAR_ROW(istr:iend)
             if (bec2_diag_2d_l(pco2oc_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(pco2oc_idx_t))   = pco2oc_ROW(istr:iend)
             if (bec2_diag_2d_l(pco2air_idx_t))
     &           bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(pco2air_idx_t))  = XCO2
#   else
             bec2_diag_2d(istr:iend,j,ph_idx_t) = PH_NEW(istr:iend)
             bec2_diag_2d(istr:iend,j,co2star_idx_t)  = CO2STAR_ROW(istr:iend)
             bec2_diag_2d(istr:iend,j,dco2star_idx_t) = DCO2STAR_ROW(istr:iend)
             bec2_diag_2d(istr:iend,j,pco2oc_idx_t)   = pco2oc_ROW(istr:iend)
             bec2_diag_2d(istr:iend,j,pco2air_idx_t)  = XCO2
#   endif
#  else /* BEC2_DIAG */
             ph_hist(istr:iend,j) = PH_NEW(istr:iend)
             pCO2sw(istr:iend,j) = pco2oc_ROW(istr:iend)
#   ifdef PCO2AIR_FORCING
             pCO2air(istr:iend,j) = XCO2(istr:iend)
#   endif
#  endif /* BEC2_DIAG */
# endif /* CCHEM_MOCSY */

             FLUX(istr:iend,j) = PV(istr:iend,j) * DCO2STAR_ROW(istr:iend)
          END DO  ! j = jstr,jend

# if defined CCHEM_MOCSY && defined CCHEM_TODEPTH
          ! Compute carbon chemistry for k<N:
          do k=1,N-1
             do j = jstr,jend
                call vars(ph_new(istr:iend), pco2oc_ROW, fco2, co2star_row,                     ! OUTPUT
     &             hco3, co3, omega_arag, omega_calc, betad, rhoSW, p_ocean, temp_is,           ! OUTPUT
     &             temp_depth(istr:iend,j,k), salt_depth(istr:iend,j,k),                              ! INPUT
     &             0.001*tracer(istr:iend,j,k,alk_ind_t),                                             ! INPUT
     &             0.001*tracer(istr:iend,j,k,dic_ind_t), 0.001*tracer(istr:iend,j,k,sio3_ind_t),     ! INPUT
     &             0.001*tracer(istr:iend,j,k,po4_ind_t), AP_USED(istr:iend,j), -z_r(istr:iend,j,k),  ! INPUT
     &             latr(istr:iend,j), iend-istr+1, optCON, optT, optP, optB, optk1k2, optkf, optgas)  ! INPUT
#  ifdef BEC2_DIAG
#   ifdef BEC2_DIAG_USER
                ! Multiply by 1000 to convert Mol/m3 -> mMol/m3:
                if (bec2_diag_3d_l(co3_idx_t))
     &               bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(co3_idx_t)) = 1000*co3
                if (bec2_diag_3d_l(hco3_idx_t))
     &               bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(hco3_idx_t)) = 1000*hco3
                if (bec2_diag_3d_l(co2star_idx_t))
     &               bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(co2star_idx_t)) = 1000*CO2STAR_ROW(istr:iend)
                ! bec2_diag_3d_idx(pco2oc_idx_t) and bec2_diag_3d_idx(ph_idx_t)
                ! are always defined:
                bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(pco2oc_idx_t)) = pco2oc_ROW(istr:iend)
                bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(ph_idx_t)) = ph_new(istr:iend)
#   else
                ! Multiply by 1000 to convert Mol/m3 -> mMol/m3:
                bec2_diag_3d(istr:iend,j,k,co3_idx_t) = 1000*co3
                bec2_diag_3d(istr:iend,j,k,hco3_idx_t) = 1000*hco3
                bec2_diag_3d(istr:iend,j,k,co2star_idx_t) = 1000*CO2STAR_ROW(istr:iend)
                bec2_diag_3d(istr:iend,j,k,pco2oc_idx_t) = pco2oc_ROW(istr:iend)
                bec2_diag_3d(istr:iend,j,k,ph_idx_t) = ph_new(istr:iend)
#   endif /* BEC2_DIAG_USER */
#  endif /* BEC2_DIAG */
             end do  ! j = jstr,jend
          end do  ! k=1,N-1
# endif /* CCHEM_MOCSY && CCHEM_TODEPTH */

          STF(istr:iend,jstr:jend,dic_ind_t) = FLUX(istr:iend,jstr:jend)

          tracer(istr:iend,jstr:jend,N,dic_ind_t) =
     &       tracer(istr:iend,jstr:jend,N,dic_ind_t) +
     &       stf(istr:iend,jstr:jend,dic_ind_t)
     &       *dt/Hz(istr:iend,jstr:jend,N)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
          if (bec2_diag_2d_l(fgco2_idx_t))
     &        bec2_diag_2d(istr:iend,jstr:jend,bec2_diag_2d_idx(fgco2_idx_t)) = FLUX(istr:iend,jstr:jend)
#  else
          bec2_diag_2d(istr:iend,jstr:jend,fgco2_idx_t) = FLUX(istr:iend,jstr:jend)
#  endif
# endif
# ifdef WRITE_CO2FLX
          CO2flx(istr:iend,jstr:jend) = FLUX(istr:iend,jstr:jend)
# endif
       ENDIF  ! lflux_gas_co2

    !---------------------------------------------------------------------
    !   Surface fluxes:
    !----------------------------------------------------------------------

# ifdef RIVER_LOAD_N
       !CN: N loads from rivers
       if (lriver_load_n) then
          FLUX(istr:iend,jstr:jend)=din_river(istr:iend,jstr:jend)
       else
          FLUX = c0
       endif
      

       ! convert unit: from kg N -> mmol N  (with atomic mass of N:
       ! 14.00674 g/mol)
          FLUX(istr:iend,jstr:jend)=FLUX(istr:iend,jstr:jend)*71394.200220751
       ! convert from per day to per second!   
          FLUX(istr:iend,jstr:jend)=FLUX(istr:iend,jstr:jend)/86400
!       ! assumption: only 1% is bioavailable (very crude, just a test)
!          FLUX(istr:iend,jstr:jend)=FLUX(istr:iend,jstr:jend)*0.01

       ! add to tracer: NO3       
       tracer(istr:iend,jstr:jend,N,no3_ind_t) =
     &         tracer(istr:iend,jstr:jend,N,no3_ind_t) +
     &         FLUX(istr:iend,jstr:jend)
     &         *dt/Hz(istr:iend,jstr:jend,N)

!       ! add to tracer: DON       
!       tracer(istr:iend,jstr:jend,N,don_ind_t) =
!     &         tracer(istr:iend,jstr:jend,N,don_ind_t) +
!     &         0.5*FLUX(istr:iend,jstr:jend)
!     &         *dt/Hz(istr:iend,jstr:jend,N)

# endif /* RIVER_LOAD_N */

# ifdef RIVER_LOAD_P
       !CN: P loads from rivers
       if (lriver_load_p) then
          FLUX(istr:iend,jstr:jend)=dip_river(istr:iend,jstr:jend)
       else
          FLUX = c0
       endif
      
       ! onvert unit: from kg P -> mmol P  (with atomic mass of P:
       ! 30.97376 g/mol)
          FLUX(istr:iend,jstr:jend)=FLUX(istr:iend,jstr:jend)*32285.392538717
       ! convert from per day to per second!   
          FLUX(istr:iend,jstr:jend)=FLUX(istr:iend,jstr:jend)/86400
!       ! assumption: only 1% is bioavailable (very crude, just a test)
!          FLUX(istr:iend,jstr:jend)=FLUX(istr:iend,jstr:jend)*0.01

       ! add to tracer: PO4       
       tracer(istr:iend,jstr:jend,N,po4_ind_t) =
     &         tracer(istr:iend,jstr:jend,N,po4_ind_t) +
     &         FLUX(istr:iend,jstr:jend)
     &         *dt/Hz(istr:iend,jstr:jend,N)

!       ! add to tracer: DOP
!       tracer(istr:iend,jstr:jend,N,dop_ind_t) =
!     &         tracer(istr:iend,jstr:jend,N,dop_ind_t) +
!     &         0.25*FLUX(istr:iend,jstr:jend)
!     &         *dt/Hz(istr:iend,jstr:jend,N)
# endif /* RIVER_LOAD_P */

# ifdef RIVER_LOAD_SI_ALK_DIC   
       !CN: Si loads from rivers
       if (lriver_load_si) then
          FLUX(istr:iend,jstr:jend)=si_river(istr:iend,jstr:jend)
       else
          FLUX = c0
       endif
       ! convert from per day to per second!   
          FLUX(istr:iend,jstr:jend)=FLUX(istr:iend,jstr:jend)/86400
       ! add to tracer: SiO3
       tracer(istr:iend,jstr:jend,N,sio3_ind_t) =
     &         tracer(istr:iend,jstr:jend,N,sio3_ind_t) +
     &         FLUX(istr:iend,jstr:jend)
     &         *dt/Hz(istr:iend,jstr:jend,N)


       !CN: Alk loads from rivers
       if (lriver_load_alk) then
          FLUX(istr:iend,jstr:jend)=alk_river(istr:iend,jstr:jend)
       else
          FLUX = c0
       endif
       ! convert from per day to per second!   
          FLUX(istr:iend,jstr:jend)=FLUX(istr:iend,jstr:jend)/86400
       ! add to tracer: Alk
       tracer(istr:iend,jstr:jend,N,alk_ind_t) =
     &         tracer(istr:iend,jstr:jend,N,alk_ind_t) +
     &         FLUX(istr:iend,jstr:jend)
     &         *dt/Hz(istr:iend,jstr:jend,N)


       !CN: DIC loads from rivers
       if (lriver_load_dic) then
          FLUX(istr:iend,jstr:jend)=dic_river(istr:iend,jstr:jend)
       else
          FLUX = c0
       endif
       ! convert from per day to per second!   
          FLUX(istr:iend,jstr:jend)=FLUX(istr:iend,jstr:jend)/86400
       ! add to tracer: DIC
       tracer(istr:iend,jstr:jend,N,dic_ind_t) =
     &         tracer(istr:iend,jstr:jend,N,dic_ind_t) +
     &         FLUX(istr:iend,jstr:jend)
     &         *dt/Hz(istr:iend,jstr:jend,N)
# endif /* RIVER_LOAD_SI_ALK_DIC */


       !CN: iron fluxes in forcing file are in nmol Fe cm-2 s-1
       !     BEC needs them in mmol Fe m-2 s-1, convert here 
       !     (used to be done in get_iron.F)
       if (liron_flux) then
          FLUX(istr:iend,jstr:jend)=iron(istr:iend,jstr:jend)
       else
          FLUX = c0
       endif 
       FLUX = FLUX * parm_Fe_bioavail * nmol_cm2_to_mmol_m2 
       STF(istr:iend,jstr:jend,fe_ind_t) = FLUX(istr:iend,jstr:jend)
       tracer(istr:iend,jstr:jend,N,fe_ind_t) =
     &         tracer(istr:iend,jstr:jend,N,fe_ind_t) +
     &         stf(istr:iend,jstr:jend,fe_ind_t)
     &         *dt/Hz(istr:iend,jstr:jend,N)

# ifdef NOX_FORCING
       ! NOx: added to nitrate (convert units from kg N to mmol N)
       tracer(istr:iend,jstr:jend,N,no3_ind_t) =
     &         tracer(istr:iend,jstr:jend,N,no3_ind_t) +
     &         nox(istr:iend,jstr:jend)*71394.200220751
     &         *dt/Hz(istr:iend,jstr:jend,N)
# endif
# ifdef NHY_FORCING
       ! NHy: added to ammonium (convert units from kg N to mmol N)        
       tracer(istr:iend,jstr:jend,N,nh4_ind_t) =
     &         tracer(istr:iend,jstr:jend,N,nh4_ind_t) +
     &         nhy(istr:iend,jstr:jend)*71394.200220751
     &         *dt/Hz(istr:iend,jstr:jend,N)
# endif

# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_2d_l(ironflux_idx_t))
     &     bec2_diag_2d(istr:iend,jstr:jend,bec2_diag_2d_idx(ironflux_idx_t)) = FLUX(istr:iend,jstr:jend)
#  else
       bec2_diag_2d(istr:iend,jstr:jend,ironflux_idx_t) = FLUX(istr:iend,jstr:jend)
#  endif
# endif

!end of set_sflux part

    !---------------------------------------------------------------------
    ! Now we compute the tendencies for the biological tracers
    ! (stored in array DTRACER_MODULE).
    !
    ! NOTE: the k-loop needs to be in this order 
    ! because PAR etc. is initialized for k==N
    !---------------------------------------------------------------------

       do j = jstr, jend

# ifdef USE_EXPLICIT_VSINK
    !---------------------------------------------------------------------
    ! Compute vertical sinking:
    !---------------------------------------------------------------------
          call vert_sinking_bec2(istr,iend,j)
# endif
          do k = N,1,-1

          call ecosys_set_interior(k,t(istr:iend,j,k,nnew,1),
     &             SHF_QSW(istr:iend,j),  
     &             PAR_out(istr:iend,j),
     &             qa_dust_def(istr:iend,j),istr,iend,j,
     &             dt,dust(istr:iend,j))
          end do
       end do

    !---------------------------------------------------------------------
    ! Compute and write back new tracer values to ROMS:
    ! the explicit Euler method is used for computing the new
    ! biological tracer values.
    !---------------------------------------------------------------------

       do k=1,n
          do m=1,ntrc_bio
             t(istr:iend,jstr:jend,k,nnew,1+ntrc_salt+ntrc_pas+m) =
     &             tracer(istr:iend,jstr:jend,k,m)
     &             + dt * DTRACER_MODULE(istr:iend,jstr:jend,k,m)
          enddo 
       enddo

       end subroutine ecosys_bec2_tile


!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE ECOSYS_SET_INTERIOR   ********************
!
!-------------------------------------------------------------------------

       subroutine ecosys_set_interior(k,temp,SHF_QSW,  
     &  PAR_out,qa_dust_def,istr,iend,j,dt2,dust_flux)

        implicit none
#include "param.h"  
#include "param_bec2.h"
#include "ecosys_bec2.h"
#include "ocean3d.h"
#include "scalars.h"
#include "grid.h"

       integer istr,iend,i,j,k,n_loc,auto_ind,auto_ind2
       real dt2,tmp1,tmp2,tmp3,tmp4
       real temp(istr:iend), SHF_QSW(istr:iend),
     &    PAR_out(istr:iend), dust_flux(istr:iend)

       REAL, PARAMETER :: 
     &    epsC      = 1.00e-8,  ! small C concentration (mmol C/m^3)
     &    epsTinv   = 3.17e-8,  ! small inverse time scale (1/year) (1/sec)
     &    epsnondim = 1.00e-6,  ! small non-dimensional number (non-dim)
     &    cks = 9.0,         ! constant used in Fe quota modification
     &    cksi = 5.0         ! constant used in Si quota modification

       REAL, DIMENSION(istr:iend) :: 
     &    PO4_loc,         ! local copy of model PO4
     &    NO3_loc,         ! local copy of model NO3
     &    SiO3_loc,        ! local copy of model SiO3
     &    NH4_loc,         ! local copy of model NH4
     &    Fe_loc,          ! local copy of model Fe
     &    O2_loc,          ! local copy of model O2
     &    DOC_loc,         ! local copy of model DOC
     &    spC_loc,         ! local copy of model spC
     &    spChl_loc,       ! local copy of model spChl
     &    spCaCO3_loc,     ! local copy of model spCaCO3
     &    diatC_loc,       ! local copy of model diatC
     &    diatChl_loc,     ! local copy of model diatChl
     &    zooC_loc,        ! local copy of model zooC
     &    spFe_loc,        ! local copy of model spFe
     &    diatSi_loc,      ! local copy of model diatSi
     &    diatFe_loc,      ! local copy of model diatFe
     &    diazC_loc,       ! local copy of model diazC
     &    diazChl_loc,     ! local copy of model diazChl
     &    diazFe_loc,      ! local copy of model diazFe
     &    DON_loc,         ! local copy of model DON
     &    DOFe_loc,        ! local copy of model DOFe
     &    DOP_loc,         ! local copy of model DOP
     &    DOPr_loc,        ! local copy of model DOPr
     &    DONr_loc,        ! local copy of model DONr
# ifdef BEC_COCCO
     &    coccoC_loc,      ! local copy of model coccoC
     &    coccoChl_loc,    ! local copy of model coccoChl
     &    coccoCal_loc,    ! local copy of model coccoCal
     &    coccoFe_loc,     ! local copy of model coccoFe
     &    cal_loc,         ! local copy of model cal
# endif
# ifdef Ncycle_SY
     &    NO2_loc,         ! local copy of model no2
     &    N2O_atm_loc,     ! local copy of N2O from atmosphere
     &    N2O_ao1_loc,         ! local copy of model n2o prod. via ao1
     &    N2O_siden_loc,         ! local copy of model n2o sinked by den
     &    N2O_soden_loc,         ! local copy of model n2o sourced by den
     &    N2_loc,         ! local copy of model n2
     &    N2O_loc,         ! N2O_ao1_loc + N2O_soden_loc - N2O_siden_loc
     &    ammox,          ! ammonium oxidation (NH4 -> NO2) (mmol N/m^3/sec)
     &    nitrox,          ! nitrite oxidation (NO2 -> NO3) (mmolN/m^3/sec)
     &    anammox,          ! anaerobic ammonium oxidation (NH4 -> N2) (mmolN/m^3/sec)
     &    epsN,         ! small value for denitrification calculations
     &    R_oxic,        ! Pseudo oxic remineralization (NO3-> NO2, no units) 
     &    R_den1,        ! Pseudo denitrification 1 (NO3-> NO2, no units) 
     &    R_den2,        ! pseudo denitrification 2 (NO2-> N2O), no units) 
     &    R_den3,        ! Pseudo denitrification 3 (N2O-> N2, no units) 
     &    roxic,       ! portion of C remin done by denitrif that is done by denitrif 
     &    rden1,       ! portion of C remin done by denitrif that is done by denitrif 
     &    rden2,       ! portion of C remin done by denitrif that is done by denitrif 
     &    rden3,       ! portion of C remin done by denitrif that is done by denitrif 
     &    Rden,       ! portion of C remin done by denitrif 
     &    DENITRIF1,       ! NO3-->NO2 consumed during denitrification
     &    DENITRIF2,       ! NO2-->N2O consumed during denitrification
     &    DENITRIF3,       ! N2O-->N2 consumed during denitrification
     &    n2o_nh4_yield1,   ! N2O yield during ammonium oxidation (0 --> 1)  
     &    no2_nh4_yield1,   ! NO2 yield during ammonium oxidation (0 --> 1)
# endif
# ifdef N2O_NEV
     &    N2O_nev_loc,         ! local copy of model n2o using Nevisson
     &    n2o_prod_nev,
     &    n2o_cons_nev,
# endif
# ifdef EXPLICIT_MICROBES
     &    AOA_loc,         ! local copy of model AOA
     &    NOB_loc,         ! local copy of model NOB
     &    AOX_loc,         ! local copy of model AOX
     &    AER_loc,         ! local copy of model AER
     &    NAR_loc,         ! local copy of model NAR
     &    NAI_loc,         ! local copy of model NAI
     &    NIR_loc,         ! local copy of model NIR
     &    NIO_loc,         ! local copy of model NIO
     &    NOS_loc,         ! local copy of model NOS
     &    NAO_loc,         ! local copy of model NAO
     &    DOCr_loc,        ! local copy of model DOCr
     &    szooC_loc,       ! local copy of model small zooC
     &    DOC_tot,         ! total DOC
     &    AOAprime,        ! AOA biomass minus minimum biomass threshold 
     &    NOBprime,        ! NOB biomass minus minimum biomass threshold
     &    AOXprime,        ! AOX biomass minus minimum biomass threshold
     &    AERprime,        ! AER biomass minus minimum biomass threshold 
     &    NARprime,        ! NAR biomass minus minimum biomass threshold 
     &    NAIprime,        ! NAI biomass minus minimum biomass threshold 
     &    NIRprime,        ! NIR biomass minus minimum biomass threshold
     &    NIOprime,        ! NIO biomass minus minimum biomass threshold
     &    NOSprime,        ! NOS biomass minus minimum biomass threshold
     &    NAOprime,        ! NAO biomass minus minimum biomass threshold 
     &    SZOOprime,       ! SZOO biomass minus minimum biomass threshold
     &    zuptoxy,         ! uptake rate of oxygen by microbe
     &    zuptnh4,         ! uptake rate of ammonium by microbe
     &    zuptno2,         ! uptake rate of nitrite by microbe
     &    zuptdoc,         ! uptake rate of dissolved organic matter by microbe
     &    zuptno3,         ! uptake rate of nitrate by microbe
     &    zuptn2o,         ! uptake rate of nitrous oxide by microbe
     &    zlimaer_doc,     ! limitation of AER growth by DOC
     &    zlimaer_docr,    ! limitation of AER growth by DOCr
     &    zlimnar_doc,     ! limitation of NAR growth by DOC
     &    zlimnar_docr,    ! limitation of NAR growth by DOCr
     &    zlimnai_doc,     ! limitation of NAI growth by DOC
     &    zlimnai_docr,    ! limitation of NAI growth by DOCr
     &    zlimnir_doc,     ! limitation of NIR growth by DOC
     &    zlimnir_docr,    ! limitation of NIR growth by DOCr
     &    zlimnio_doc,     ! limitation of NIO growth by DOC
     &    zlimnio_docr,    ! limitation of NIO growth by DOCr
     &    zlimnos_doc,     ! limitation of NOS growth by DOC
     &    zlimnos_docr,    ! limitation of NOS growth by DOCr
     &    zlimnao_doc,     ! limitation of NAO growth by DOC
     &    zlimnao_docr,    ! limitation of NAO growth by DOCr
     &    zmuaoa,          ! realised growth rate of AOA
     &    zmunob,          ! realised growth rate of NOB
     &    zmunob_ana,      ! realised growth rate of NOB (anaerobic)
     &    zmuaox,          ! realised growth rate of AOX
     &    zmuaer,          ! realised growth rate of AER
     &    zmunar_aer,      ! realised aerobic growth rate of NAR
     &    zmunar_ana,      ! realised anaerobic growth rate of NAR
     &    zmunar_fac,      ! facultative growth? (1 == yes)
     &    zmunar,          ! realised growth rate of NAR
     &    zmunai_aer,      ! realised aerobic growth rate of NAI
     &    zmunai_ana,      ! realised anaerobic growth rate of NAI
     &    zmunai_fac,      ! facultative growth? (1 == yes)
     &    zmunai,          ! realised growth rate of NAI
     &    zmunir_aer,      ! realised aerobic growth rate of NIR
     &    zmunir_ana,      ! realised anaerobic growth rate of NIR
     &    zmunir_fac,      ! facultative growth? (1 == yes)
     &    zmunir,          ! realised growth rate of NIR
     &    zmunio_aer,      ! realised aerobic growth rate of NIO
     &    zmunio_ana,      ! realised anaerobic growth rate of NIO
     &    zmunio_fac,      ! facultative growth? (1 == yes)
     &    zmunio,          ! realised growth rate of NIO
     &    zmunos_aer,      ! realised aerobic growth rate of NOS
     &    zmunos_ana,      ! realised anaerobic growth rate of NOS
     &    zmunos_fac,      ! facultative growth? (1 == yes)
     &    zmunos,          ! realised growth rate of NOS
     &    zmunao_aer,      ! realised aerobic growth rate of NAO
     &    zmunao_ana,      ! realised anaerobic growth rate of NAO
     &    zmunao_fac,      ! facultative growth? (1 == yes)
     &    zmunao,          ! realised growth rate of NAO
     &    zbioaoa,         ! biomass accumulation rate of AOA
     &    zbionob,         ! biomass accumulation rate of NOB
     &    zbionob_ana,     ! biomass accumulation rate of NOB (anaerobic)
     &    zbioaox,         ! biomass accumulation rate of AOX
     &    zbioaer,         ! biomass accumulation rate of AER
     &    zbionar,         ! biomass accumulation rate of NAR
     &    zbionai,         ! biomass accumulation rate of NAI
     &    zbionir,         ! biomass accumulation rate of NIR
     &    zbionio,         ! biomass accumulation rate of NIO
     &    zbionos,         ! biomass accumulation rate of NOS
     &    zbionao,         ! biomass accumulation rate of NAO
     &    zmufnob,         ! f-value for NOB growth (1 = limited by NO2, 0 = limited by O2)
     &    yaoa_n2o,        ! N2O production yield during ammonia oxidation
     &    zmortaoa,        ! mortality loss rate of AOA
     &    zmortaoa_doc,    ! mortality loss rate of AOA routed to DOC
     &    zmortaoa_dic,    ! mortality loss rate of AOA routed to DIC
     &    zmortaoa_poc,    ! mortality loss rate of AOA routed to POC
     &    zmortaoa_don,    ! mortality loss rate of AOA routed to DON
     &    zmortaoa_dop,    ! mortality loss rate of AOA routed to DOP
     &    zmortnob2,       ! temporary mortality loss rate of NOB (remove later!)
     &    zmortnob,        ! mortality loss rate of NOB
     &    zmortnob_doc,    ! mortality loss rate of NOB routed to DOC
     &    zmortnob_dic,    ! mortality loss rate of NOB routed to DIC
     &    zmortnob_poc,    ! mortality loss rate of NOB routed to POC
     &    zmortnob_don,    ! mortality loss rate of NOB routed to DON
     &    zmortnob_dop,    ! mortality loss rate of NOB routed to DOP
     &    zmortaox,        ! mortality loss rate of AOX
     &    zmortaox_doc,    ! mortality loss rate of AOX routed to DOC
     &    zmortaox_dic,    ! mortality loss rate of AOX routed to DIC
     &    zmortaox_poc,    ! mortality loss rate of AOX routed to POC
     &    zmortaox_don,    ! mortality loss rate of AOX routed to DON
     &    zmortaox_dop,    ! mortality loss rate of AOX routed to DOP
     &    zmortaer,        ! mortality loss rate of AER
     &    zmortaer_doc,    ! mortality loss rate of AER routed to DOC
     &    zmortaer_dic,    ! mortality loss rate of AER routed to DIC
     &    zmortaer_poc,    ! mortality loss rate of AER routed to POC
     &    zmortaer_don,    ! mortality loss rate of AER routed to DON
     &    zmortaer_dop,    ! mortality loss rate of AER routed to DOP
     &    zmortnar,        ! mortality loss rate of NAR
     &    zmortnar_doc,    ! mortality loss rate of NAR routed to DOC
     &    zmortnar_dic,    ! mortality loss rate of NAR routed to DIC
     &    zmortnar_poc,    ! mortality loss rate of NAR routed to POC
     &    zmortnar_don,    ! mortality loss rate of NAR routed to DON
     &    zmortnar_dop,    ! mortality loss rate of NAR routed to DOP
     &    zmortnai,        ! mortality loss rate of NAI
     &    zmortnai_doc,    ! mortality loss rate of NAI routed to DOC
     &    zmortnai_dic,    ! mortality loss rate of NAI routed to DIC
     &    zmortnai_poc,    ! mortality loss rate of NAI routed to POC
     &    zmortnai_don,    ! mortality loss rate of NAI routed to DON
     &    zmortnai_dop,    ! mortality loss rate of NAI routed to DOP
     &    zmortnir,        ! mortality loss rate of NIR
     &    zmortnir_doc,    ! mortality loss rate of NIR routed to DOC
     &    zmortnir_dic,    ! mortality loss rate of NIR routed to DIC
     &    zmortnir_poc,    ! mortality loss rate of NIR routed to POC
     &    zmortnir_don,    ! mortality loss rate of NIR routed to DON
     &    zmortnir_dop,    ! mortality loss rate of NIR routed to DOP
     &    zmortnio,        ! mortality loss rate of NIO
     &    zmortnio_doc,    ! mortality loss rate of NIO routed to DOC
     &    zmortnio_dic,    ! mortality loss rate of NIO routed to DIC
     &    zmortnio_poc,    ! mortality loss rate of NIO routed to POC
     &    zmortnio_don,    ! mortality loss rate of NIO routed to DON
     &    zmortnio_dop,    ! mortality loss rate of NIO routed to DOP
     &    zmortnos,        ! mortality loss rate of NOS
     &    zmortnos_doc,    ! mortality loss rate of NOS routed to DOC
     &    zmortnos_dic,    ! mortality loss rate of NOS routed to DIC
     &    zmortnos_poc,    ! mortality loss rate of NOS routed to POC
     &    zmortnos_don,    ! mortality loss rate of NOS routed to DON
     &    zmortnos_dop,    ! mortality loss rate of NOS routed to DOP
     &    zmortnao,        ! mortality loss rate of NAO
     &    zmortnao_doc,    ! mortality loss rate of NAO routed to DOC
     &    zmortnao_dic,    ! mortality loss rate of NAO routed to DIC
     &    zmortnao_poc,    ! mortality loss rate of NAO routed to POC
     &    zmortnao_don,    ! mortality loss rate of NAO routed to DON
     &    zmortnao_dop,    ! mortality loss rate of NAO routed to DOP
     &    zzooaoa,         ! grazing mortality loss rate of AOA
     &    zzooaoa_zoo,     ! zzooaoa routed to zoo (mmol C/m^3/sec)
     &    zzooaoa_poc,     ! zzooaoa routed to poc (mmol C/m^3/sec)
     &    zzooaoa_doc,     ! zzooaoa routed to doc (mmol C/m^3/sec)
     &    zzooaoa_dic,     ! zzooaoa routed to dic (mmol C/m^3/sec)
     &    zzooaoa_excN,    ! zzooaoa N routed to zoo but excreted (mmol N/m^3/sec)
     &    zzooaoa_excP,    ! zzooaoa P routed to zoo but excreted (mmol P/m^3/sec)
     &    zzooaoa_excFe,   ! zzooaoa Fe routed to zoo but becomes dFe (mmol Fe/m^3/sec)
     &    zzoonob,         ! grazing mortality loss rate of NOB
     &    zzoonob_zoo,     ! zzoonob routed to zoo (mmol C/m^3/sec)
     &    zzoonob_poc,     ! zzoonob routed to poc (mmol C/m^3/sec)
     &    zzoonob_doc,     ! zzoonob routed to doc (mmol C/m^3/sec)
     &    zzoonob_dic,     ! zzoonob routed to dic (mmol C/m^3/sec)
     &    zzoonob_excN,    ! zzoonob N routed to zoo but excreted (mmol N/m^3/sec)
     &    zzoonob_excP,    ! zzoonob P routed to zoo but excreted (mmol P/m^3/sec)
     &    zzoonob_excFe,   ! zzoonob Fe routed to zoo but becomes dFe (mmol Fe/m^3/sec)
     &    zzooaox,         ! grazing mortality loss rate of AOX
     &    zzooaox_zoo,     ! zzooaox routed to zoo (mmol C/m^3/sec)
     &    zzooaox_poc,     ! zzooaox routed to poc (mmol C/m^3/sec)
     &    zzooaox_doc,     ! zzooaox routed to doc (mmol C/m^3/sec)
     &    zzooaox_dic,     ! zzooaox routed to dic (mmol C/m^3/sec)
     &    zzooaox_excN,    ! zzooaox N routed to zoo but excreted (mmol N/m^3/sec)
     &    zzooaox_excP,    ! zzooaox P routed to zoo but excreted (mmol P/m^3/sec)
     &    zzooaox_excFe,   ! zzooaox Fe routed to zoo but becomes dFe (mmol Fe/m^3/sec)
     &    zszooaer,         ! grazing mortality loss rate of AER
     &    zszooaer_szoo,     ! zszooaer routed to zoo (mmol C/m^3/sec)
     &    zszooaer_poc,     ! zszooaer routed to poc (mmol C/m^3/sec)
     &    zszooaer_doc,     ! zszooaer routed to doc (mmol C/m^3/sec)
     &    zszooaer_dic,     ! zszooaer routed to dic (mmol C/m^3/sec)
     &    zszooaer_excN,    ! zszooaer N routed to zoo but excreted (mmol N/m^3/sec)
     &    zszooaer_excP,    ! zszooaer P routed to zoo but excreted (mmol P/m^3/sec)
     &    zszooaer_excFe,   ! zszooaer Fe routed to zoo but becomes dFe (mmol Fe/m^3/sec)
     &    zszoonar,         ! grazing mortality loss rate of NAR
     &    zszoonar_szoo,     ! zszoonar routed to zoo (mmol C/m^3/sec)
     &    zszoonar_poc,     ! zszoonar routed to poc (mmol C/m^3/sec)
     &    zszoonar_doc,     ! zszoonar routed to doc (mmol C/m^3/sec)
     &    zszoonar_dic,     ! zszoonar routed to dic (mmol C/m^3/sec)
     &    zszoonar_excN,    ! zszoonar N routed to zoo but excreted (mmol N/m^3/sec)
     &    zszoonar_excP,    ! zszoonar P routed to zoo but excreted (mmol P/m^3/sec)
     &    zszoonar_excFe,   ! zszoonar Fe routed to zoo but becomes dFe (mmol Fe/m^3/sec)
     &    zszoonai,         ! grazing mortality loss rate of NAI
     &    zszoonai_szoo,     ! zszoonai routed to zoo (mmol C/m^3/sec)
     &    zszoonai_poc,     ! zszoonai routed to poc (mmol C/m^3/sec)
     &    zszoonai_doc,     ! zszoonai routed to doc (mmol C/m^3/sec)
     &    zszoonai_dic,     ! zszoonai routed to dic (mmol C/m^3/sec)
     &    zszoonai_excN,    ! zszoonai N routed to zoo but excreted (mmol N/m^3/sec)
     &    zszoonai_excP,    ! zszoonai P routed to zoo but excreted (mmol P/m^3/sec)
     &    zszoonai_excFe,   ! zszoonai Fe routed to zoo but becomes dFe (mmol Fe/m^3/sec)
     &    zszoonir,         ! grazing mortality loss rate of NIR
     &    zszoonir_szoo,     ! zszoonir routed to zoo (mmol C/m^3/sec)
     &    zszoonir_poc,     ! zszoonir routed to poc (mmol C/m^3/sec)
     &    zszoonir_doc,     ! zszoonir routed to doc (mmol C/m^3/sec)
     &    zszoonir_dic,     ! zszoonir routed to dic (mmol C/m^3/sec)
     &    zszoonir_excN,    ! zszoonir N routed to zoo but excreted (mmol N/m^3/sec)
     &    zszoonir_excP,    ! zszoonir P routed to zoo but excreted (mmol P/m^3/sec)
     &    zszoonir_excFe,   ! zszoonir Fe routed to zoo but becomes dFe (mmol Fe/m^3/sec)
     &    zszoonio,         ! grazing mortality loss rate of NIO
     &    zszoonio_szoo,     ! zszoonio routed to zoo (mmol C/m^3/sec)
     &    zszoonio_poc,     ! zszoonio routed to poc (mmol C/m^3/sec)
     &    zszoonio_doc,     ! zszoonio routed to doc (mmol C/m^3/sec)
     &    zszoonio_dic,     ! zszoonio routed to dic (mmol C/m^3/sec)
     &    zszoonio_excN,    ! zszoonio N routed to zoo but excreted (mmol N/m^3/sec)
     &    zszoonio_excP,    ! zszoonio P routed to zoo but excreted (mmol P/m^3/sec)
     &    zszoonio_excFe,   ! zszoonio Fe routed to zoo but becomes dFe (mmol Fe/m^3/sec)
     &    zszoonos,         ! grazing mortality loss rate of NOS
     &    zszoonos_szoo,     ! zszoonos routed to zoo (mmol C/m^3/sec)
     &    zszoonos_poc,     ! zszoonos routed to poc (mmol C/m^3/sec)
     &    zszoonos_doc,     ! zszoonos routed to doc (mmol C/m^3/sec)
     &    zszoonos_dic,     ! zszoonos routed to dic (mmol C/m^3/sec)
     &    zszoonos_excN,    ! zszoonos N routed to zoo but excreted (mmol N/m^3/sec)
     &    zszoonos_excP,    ! zszoonos P routed to zoo but excreted (mmol P/m^3/sec)
     &    zszoonos_excFe,   ! zszoonos Fe routed to zoo but becomes dFe (mmol Fe/m^3/sec)
     &    zszoonao,         ! grazing mortality loss rate of NAO
     &    zszoonao_szoo,     ! zszoonao routed to zoo (mmol C/m^3/sec)
     &    zszoonao_poc,     ! zszoonao routed to poc (mmol C/m^3/sec)
     &    zszoonao_doc,     ! zszoonao routed to doc (mmol C/m^3/sec)
     &    zszoonao_dic,     ! zszoonao routed to dic (mmol C/m^3/sec)
     &    zszoonao_excN,    ! zszoonao N routed to zoo but excreted (mmol N/m^3/sec)
     &    zszoonao_excP,    ! zszoonao P routed to zoo but excreted (mmol P/m^3/sec)
     &    zszoonao_excFe,   ! zszoonao Fe routed to zoo but becomes dFe (mmol Fe/m^3/sec)
     &    qDON,            ! DON to DOC ratio
     &    qDOP,            ! DOP to DOC ratio
     &    qDONr,           ! DONr to DOCr ratio
     &    qDOPr,           ! DOPr to DOCr ratio
     &    DIN_release,     ! release of DIN by heterotrophic bacteria
     &    DIP_release,     ! release of DIP by heterotrophic bacteria
     &    DIN_uptake,      ! supplementary uptake of DIN by heterotrophic bacteria
     &    DIP_uptake,      ! supplementary uptake of DIP by heterotrophic bacteria
     &    dFe_require,     ! required uptake of dFe by heterotrophic bacteria
     &    dFe_uptake,      ! supplementary uptake of dFe by heterotrophic bacteria
     &    szoo_loss,       ! mortality & higher trophic grazing on small zooplankton (mmol C/m^3/sec)
     &    szoo_loss_doc,   ! szoo_loss routed to doc (mmol C/m^3/sec)
     &    szoo_loss_dic,   ! szoo_loss routed to dic (mmol C/m^3/sec)
     &    szoo_loss_poc,   ! szoo_loss routed to poc (mmol C/m^3/sec)
     &    DENITRIF4,       ! NO3-->N2O consumed during denitrification
     &    DENITRIF5,       ! NO2-->N2 consumed during denitrification
     &    DENITRIF6,       ! NO3-->N2 consumed during denitrification
     &    nitrox_ana,      ! anaerobic nitrite oxidation (NO2 -> NO3) (mmolN/m^3/sec)
# endif /* EXPLICIT_MICROBES */
     &    QA_dust_def,     ! incoming deficit in the QA(dust) POC flux
     &    NITRIF,          ! nitrification (NH4 -> NO3) (mmol N/m^3/sec)
     &    DENITRIF,        ! water column nitrification (NO3 -> N2) (mmol N/m^3/sec)
     &    SED_DENITRIF,    ! sedimentary denitrification (mmol N/m^3/sec)
     &    OTHER_REMIN,     ! organic C remin not due oxic or denitrif (mmolC/m^3/sec)
     &    flux_oxidated,   ! carbon oxidation rate in sediment
     &    O2_CONSUMPTION,  ! O2 consumption
     &    O2_PRODUCTION    ! O2 production
# ifdef IODINE
     &    ,iod_sed_flux1
     &    ,iod_sed_flux2
     &    ,iod_sed_flux3
     &    ,iod_sed_flux4
# endif

       REAL, DIMENSION(istr:iend) ::
     &    Tfunc,           ! temp response function GD98 (non-dim)
     &    Tfunc_diat,      ! temp response function GD98 (non-dim) for diatoms
# ifdef BEC_COCCO
     &    Tfunc_cocco,     ! temp response function GD98 (non-dim) for coccos
     &    light_lim_calc,  ! light limitation factor for calc of coccos
# endif
     &    PAR_in,          ! photosynthetically available radiation (W/m^2)
     &    KPARdz,          ! PAR adsorption coefficient (non-dim)
     &    PAR_lay,         ! average PAR over layer depth (W/m^2)
     &    f_loss_thres,    ! fraction of grazing loss reduction at depth
     &    C_loss_thres,    ! bio-C threshold at which losses go to zero (mmol C/m^3)
     &    VFe,             ! C-specific Fe uptake (non-dim)
     &    VSiO3,           ! C-specific SiO3 uptake (non-dim)
     &    f_nut,           ! nut limitation factor, modifies C fixation (non-dim)
     &    PCmax,           ! max value of PCphoto at temperature TEMP (1/sec)
     &    light_lim,       ! light limitation factor
     &    PCphoto,         ! C-specific rate of photosynth. (1/sec)
     &    pChl,            ! Chl synth. regulation term (mg Chl/mmol N)
     &    VNC,             ! C-specific N uptake rate (mmol N/mmol C/sec)
     &    z_umax,          ! max. zoo growth rate at local T (1/sec)
     &    remaining_P,     ! used in routing P from autotrophs w/ Qp different from Qp_zoo_pom
     &    Zprime,          ! used to limit zoo mort at low biomass (mmol C/m^3)
     &    zoo_loss,        ! mortality & higher trophic grazing on zooplankton (mmol C/m^3/sec)
     &    zoo_loss_doc,    ! zoo_loss routed to doc (mmol C/m^3/sec)
     &    zoo_loss_dic,    ! zoo_loss routed to dic (mmol C/m^3/sec)
     &    zoo_loss_poc,    ! zoo_loss routed to poc (mmol C/m^3/sec)
     &    f_zoo_detr_loc,  ! frac of zoo losses into large detrital pool (non-dim)
     &    DOC_prod,        ! production of DOC (mmol C/m^3/sec)
     &    DON_prod,        ! production of dissolved organic N
     &    DOFe_prod,       ! production of dissolved organic Fe
     &    DOP_prod,        ! production of dissolved organic P
     &    DOC_remin,       ! remineralization of DOC (mmol C/m^3/sec)
     &    DON_remin,       ! portion of DON remineralized
     &    DOFe_remin,      ! portion of DOFe remineralized
     &    DONr_remin,      ! poirtion of refractory DON remineralized
# ifdef EXPLICIT_MICROBES
     &    DOCr_remin,      ! portion of refractory DOC remineralized
     &    z_o2lim,         ! degree of oxygen limitation of microzooplankton grazing
     &    aoa_parlim,      ! degree of light limitation of ammonia oxidation
# endif   
     &    POC_prod,        ! POC production (mmol C/m^3/sec)
     &    P_CaCO3_prod,    ! large detritus CaCO3, base units = mmol CaCO3
     &    P_SiO2_prod,     ! prod. of large detritus SiO2
     &    Fe_scavenge,     ! loss of dissolved iron, scavenging (mmol Fe/m^3/sec)
     &    P_iron_prod,     ! P iron production (base units = mmol Fe)
     &    Fe_scavenge_rate, ! annual scavenging rate of iron as % of ambient
     &    WORK1,WORK2,WORK3      ! temporary arrays

       REAL, DIMENSION(istr:iend,autotroph_cnt) ::
     &    gQfe,            ! fe/C for growth
     &    Qfe,             ! init fe/C ratio (mmolFe/mmolC)
     &    Qsi,             ! Diatom initial Si/C ratio (mmol Si/mmol C)
     &    gQsi,            ! diatom Si/C ratio for growth (new biomass)
     &    thetaC,          ! local Chl/C ratio (mg Chl/mmol C)
     &    Pprime,          ! used to limit autotroph mort at low biomass (mmol C/m^3)
     &    C_loc,           ! amount of carbon in each autotroph species
     &    Chl_loc,         ! amount of chlorophyll in each autotroph species
     &    autoFe_loc,      ! amount of iron in each autotroph species
     &    Si_loc,          ! amount of Si in each autotroph species
     &    CaCO3_loc,       ! amount of CaCO3 in each autotroph species
     &    QCaCO3,          ! CaCO3/C ratio (mmol CaCO3/mmol C)
# ifdef Ncycle_SY
     &    VNOX,            ! NO3 + NO2 uptake rate (non-dim)
     &    VNO3,            ! NO3 uptake rate (non-dim)
     &    VNO2,            ! NO2 uptake rate (non-dim)
# else 
     &    VNO3,            ! NO3 uptake rate (non-dim)
# endif 
     &    VNH4,            ! NH4 uptake rate (non-dim)
     &    VNtot,           ! total N uptake rate (non-dim)
     &    VPO4,            ! C-specific PO4 uptake (non-dim)
     &    VDOP,            ! C-specific DOP uptake rate (non-dim)
     &    VPtot,           ! total P uptake rate/limitation (non-dim)
     &    photoC,          ! C-fixation (mmol C/m^3/sec)
     &    photoFe,         ! iron uptake
     &    photoSi,         ! silicon uptake (mmol Si/m^3/sec)
     &    NO3_V,           ! nitrate uptake (mmol NO3/m^3/sec)
#ifdef Ncycle_SY
     &    NO2_V,           ! nitrate uptake (mmol NO3/m^3/sec)
# endif
     &    NH4_V,           ! ammonium uptake (mmol NH4/m^3/sec)
     &    PO4_V,           ! PO4 uptake (mmol PO4/m^3/sec)
     &    DOP_V,           ! DOP uptake (mmol DOP/m^3/sec)
     &    CaCO3_PROD,      ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
     &    photoacc,        ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
     &    auto_loss,       ! autotroph non-grazing mort (mmol C/m^3/sec)
     &    auto_loss_poc,   ! auto_loss routed to poc (mmol C/m^3/sec)
     &    auto_loss_doc,   ! auto_loss routed to doc (mmol C/m^3/sec)
     &    auto_loss_dic,   ! auto_loss routed to dic (mmol C/m^3/sec)
     &    auto_agg,        ! autotroph aggregation (mmol C/m^3/sec)
     &    auto_graze,      ! autotroph grazing rate (mmol C/m^3/sec)
     &    auto_graze_zoo,  ! auto_graze routed to zoo (mmol C/m^3/sec)
     &    auto_graze_poc,  ! auto_graze routed to poc (mmol C/m^3/sec)
     &    auto_graze_doc,  ! auto_graze routed to doc (mmol C/m^3/sec)
     &    auto_graze_dic,  ! auto_graze routed to dic (mmol C/m^3/sec)
     &    remaining_P_dop, ! remaining_P from mort routed to DOP pool
     &    remaining_P_dip, ! remaining_P from mort routed to remin
     &    Nfix,            ! total Nitrogen fixation (mmol N/m^3/sec)
     &    Nexcrete         ! fixed N excretion


    !---------------------------------------------------------------------------
    !   exit immediately if computations are not to be performed
    !---------------------------------------------------------------------------

       IF (.NOT. lsource_sink) RETURN

    !---------------------------------------------------------------------------
    !   create local copies of model tracers
    !   treat negative values as zero
    !---------------------------------------------------------------------------
       PO4_loc  = MAX(c0, TRACER(istr:iend,j,k,po4_ind_t))
       NO3_loc  = MAX(c0, TRACER(istr:iend,j,k,no3_ind_t))
       SiO3_loc = MAX(c0, TRACER(istr:iend,j,k,sio3_ind_t))
       NH4_loc  = MAX(c0, TRACER(istr:iend,j,k,nh4_ind_t))
       Fe_loc   = MAX(c0, TRACER(istr:iend,j,k,fe_ind_t))
       O2_loc   = MAX(c0, TRACER(istr:iend,j,k,o2_ind_t))
       DOC_loc  = MAX(c0, TRACER(istr:iend,j,k,doc_ind_t))
       spC_loc  = MAX(c0, TRACER(istr:iend,j,k,spC_ind_t))
       spChl_loc= MAX(c0, TRACER(istr:iend,j,k,spChl_ind_t))
       spCaCO3_loc  = MAX(c0,TRACER(istr:iend,j,k,spCaCO3_ind_t))
       diatC_loc= MAX(c0, TRACER(istr:iend,j,k,diatC_ind_t))
       diatChl_loc  = MAX(c0, TRACER(istr:iend,j,k,diatChl_ind_t))
       zooC_loc  = MAX(c0, TRACER(istr:iend,j,k,zooC_ind_t))
       spFe_loc  = MAX(c0, TRACER(istr:iend,j,k,spFe_ind_t))
       diatSi_loc = MAX(c0, TRACER(istr:iend,j,k,diatSi_ind_t))
       diatFe_loc  = MAX(c0, TRACER(istr:iend,j,k,diatFe_ind_t))
       diazC_loc    = MAX(c0, TRACER(istr:iend,j,k,diazC_ind_t))
       diazChl_loc  = MAX(c0, TRACER(istr:iend,j,k,diazChl_ind_t))
       diazFe_loc   = MAX(c0, TRACER(istr:iend,j,k,diazFe_ind_t))
       DON_loc  = MAX(c0, TRACER(istr:iend,j,k,don_ind_t))
       DOFe_loc = MAX(c0, TRACER(istr:iend,j,k,dofe_ind_t))
       DOP_loc  = MAX(c0, TRACER(istr:iend,j,k,dop_ind_t))
       DOPr_loc = MAX(c0, TRACER(istr:iend,j,k,dopr_ind_t))
       DONr_loc = MAX(c0, TRACER(istr:iend,j,k,donr_ind_t))
# ifdef BEC_COCCO
       coccoC_loc = MAX(c0, TRACER(istr:iend,j,k,coccoc_ind_t))
       coccoChl_loc = MAX(c0, TRACER(istr:iend,j,k,coccochl_ind_t))
       coccoCal_loc = MAX(c0, TRACER(istr:iend,j,k,coccocal_ind_t))
       coccoFe_loc = MAX(c0, TRACER(istr:iend,j,k,coccofe_ind_t))
       cal_loc = MAX(c0, TRACER(istr:iend,j,k,cal_ind_t))
# endif
# ifdef Ncycle_SY
       NO2_loc  = MAX(c0, TRACER(istr:iend,j,k,no2_ind_t))
# ifdef N2O_TRACER_DECOMP
       N2O_atm_loc = MAX(c0, TRACER(istr:iend,j,k,n2o_atm_ind_t))
       N2O_ao1_loc  = MAX(c0, TRACER(istr:iend,j,k,n2o_ao1_ind_t))
       N2O_siden_loc  = MAX(c0, TRACER(istr:iend,j,k,n2o_siden_ind_t))
       N2O_soden_loc  = MAX(c0, TRACER(istr:iend,j,k,n2o_soden_ind_t))
# endif
       N2_loc  = MAX(c0, TRACER(istr:iend,j,k,n2_ind_t))
       N2O_loc = MAX(c0, TRACER(istr:iend,j,k,n2o_ind_t))
# endif
# ifdef N2O_NEV
       N2O_nev_loc  = MAX(c0, TRACER(istr:iend,j,k,n2o_nev_ind_t))
# endif
# ifdef EXPLICIT_MICROBES
       AOA_loc  = MAX(c0, TRACER(istr:iend,j,k,aoa_ind_t))
       NOB_loc  = MAX(c0, TRACER(istr:iend,j,k,nob_ind_t))
       AOX_loc  = MAX(c0, TRACER(istr:iend,j,k,aox_ind_t))
       AER_loc  = MAX(c0, TRACER(istr:iend,j,k,aer_ind_t))
       NAR_loc  = MAX(c0, TRACER(istr:iend,j,k,nar_ind_t))
       NAI_loc  = MAX(c0, TRACER(istr:iend,j,k,nai_ind_t))
       NIR_loc  = MAX(c0, TRACER(istr:iend,j,k,nir_ind_t))
       NIO_loc  = MAX(c0, TRACER(istr:iend,j,k,nio_ind_t))
       NOS_loc  = MAX(c0, TRACER(istr:iend,j,k,nos_ind_t))
       NAO_loc  = MAX(c0, TRACER(istr:iend,j,k,nao_ind_t))
       DOCr_loc  = MAX(c0, TRACER(istr:iend,j,k,docr_ind_t))
       szooC_loc  = MAX(c0, TRACER(istr:iend,j,k,szooc_ind_t))
# endif /* EXPLICIT_MICROBES */

#ifdef BEC2_DIAG
# ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(totchl_idx_t))
     &     bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(totchl_idx_t)) = tracer(istr:iend,j,k,spchl_ind_t)
     &         + tracer(istr:iend,j,k,diatchl_ind_t) + tracer(istr:iend,j,k,diazchl_ind_t)
       if (bec2_diag_3d_l(totphytoc_idx_t))
     &     bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(totphytoc_idx_t)) = tracer(istr:iend,j,k,spc_ind_t)
     &         + tracer(istr:iend,j,k,diatc_ind_t) + tracer(istr:iend,j,k,diazc_ind_t)
       if (bec2_diag_3d_l(spczero_idx_t))
     &     bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(spczero_idx_t)) = max(c0,-TRACER(istr:iend,j,k,spC_ind_t))/dt
       if (bec2_diag_3d_l(diatczero_idx_t))
     &     bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatczero_idx_t)) = max(c0,-TRACER(istr:iend,j,k,diatC_ind_t))/dt
       if (bec2_diag_3d_l(diazczero_idx_t))
     &     bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diazczero_idx_t)) = max(c0,-TRACER(istr:iend,j,k,diazC_ind_t))/dt
       if (bec2_diag_3d_l(doczero_idx_t))
     &     bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(doczero_idx_t)) = max(c0,-TRACER(istr:iend,j,k,doc_ind_t))/dt
       if (bec2_diag_3d_l(zooczero_idx_t))
     &     bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(zooczero_idx_t)) = max(c0,-TRACER(istr:iend,j,k,zooC_ind_t))/dt
       if (bec2_diag_3d_l(spcaco3zero_idx_t))
     &     bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(spcaco3zero_idx_t)) = max(c0,-TRACER(istr:iend,j,k,spCaCO3_ind_t))/dt
# else
       bec2_diag_3d(istr:iend,j,k,totchl_idx_t) = tracer(istr:iend,j,k,spchl_ind_t)
     &     + tracer(istr:iend,j,k,diatchl_ind_t) + tracer(istr:iend,j,k,diazchl_ind_t)
       bec2_diag_3d(istr:iend,j,k,totphytoc_idx_t) = tracer(istr:iend,j,k,spc_ind_t)
     &     + tracer(istr:iend,j,k,diatc_ind_t) + tracer(istr:iend,j,k,diazc_ind_t)
       bec2_diag_3d(istr:iend,j,k,spczero_idx_t) = max(c0,-TRACER(istr:iend,j,k,spC_ind_t))/dt
       bec2_diag_3d(istr:iend,j,k,diatczero_idx_t) = max(c0,-TRACER(istr:iend,j,k,diatC_ind_t))/dt
       bec2_diag_3d(istr:iend,j,k,diazczero_idx_t) = max(c0,-TRACER(istr:iend,j,k,diazC_ind_t))/dt
       bec2_diag_3d(istr:iend,j,k,doczero_idx_t) = max(c0,-TRACER(istr:iend,j,k,doc_ind_t))/dt
       bec2_diag_3d(istr:iend,j,k,zooczero_idx_t) = max(c0,-TRACER(istr:iend,j,k,zooC_ind_t))/dt
       bec2_diag_3d(istr:iend,j,k,spcaco3zero_idx_t) = max(c0,-TRACER(istr:iend,j,k,spCaCO3_ind_t))/dt
# endif
#endif

    !---------------------------------------------------------------------------
    !   If any phyto box are zero, set others to zeros.
    !---------------------------------------------------------------------------

       WHERE (spC_loc == c0 .OR. spChl_loc == c0 .OR. spFe_loc == c0)
          spC_loc = c0
          spChl_loc = c0
          spCaCO3_loc = c0
          spFe_loc = c0
       END WHERE

       WHERE (diatC_loc == c0 .OR. diatChl_loc == c0 .OR. 
     &     diatFe_loc == c0 
     &    .OR. diatSi_loc == c0)
          diatC_loc = c0
          diatChl_loc = c0
          diatFe_loc = c0
          diatSi_loc = c0
       END WHERE
    
       WHERE (diazC_loc == c0 .OR. diazChl_loc == c0 .OR. 
     &   diazFe_loc == c0)
          diazC_loc = c0
          diazChl_loc = c0
          diazFe_loc = c0
       END WHERE

# ifdef BEC_COCCO
       WHERE (coccoC_loc == c0 .OR. coccoChl_loc == c0 .OR. coccoFe_loc == c0)
          coccoC_loc = c0
          coccoChl_loc = c0
          coccoCal_loc = c0
          coccoFe_loc = c0
       END WHERE
# endif

    !---------------------------------------------------------------------------
    ! Initialize C_loc, Chl_loc, CaCO3_loc, autoFe_loc:
    !---------------------------------------------------------------------------

       C_loc(istr:iend,sp_ind) = spC_loc
       C_loc(istr:iend,diat_ind) = diatC_loc
       C_loc(istr:iend,diaz_ind) = diazC_loc
       Chl_loc(istr:iend,sp_ind) = spChl_loc
       Chl_loc(istr:iend,diat_ind) = diatChl_loc
       Chl_loc(istr:iend,diaz_ind) = diazChl_loc
       CaCO3_loc(istr:iend,sp_ind) = spCaCO3_loc
       CaCO3_loc(istr:iend,diat_ind) = c0
       CaCO3_loc(istr:iend,diaz_ind) = c0
       Si_loc(istr:iend,sp_ind) = c0
       Si_loc(istr:iend,diat_ind) = diatSi_loc
       Si_loc(istr:iend,diaz_ind) = c0
       autoFe_loc(istr:iend,sp_ind) = spFe_loc
       autoFe_loc(istr:iend,diat_ind) = diatFe_loc
       autoFe_loc(istr:iend,diaz_ind) = diazFe_loc
# ifdef BEC_COCCO
       C_loc(istr:iend,cocco_ind) = coccoC_loc
       Chl_loc(istr:iend,cocco_ind) = coccoChl_loc
       CaCO3_loc(istr:iend,cocco_ind) = coccoCal_loc
       Si_loc(istr:iend,cocco_ind) = c0
       autoFe_loc(istr:iend,cocco_ind) = coccoFe_loc
# endif

    !---------------------------------------------------------!
    !---------------------------------------------------------!
    !---------------------------------------------------------!
    !   (1.) Stoichiometry of Phytoplankton (Fe, Si, CaCO3)   !
    !---------------------------------------------------------!
    !---------------------------------------------------------!
    !---------------------------------------------------------!

    ! ---------------------------------------------------------- !
    !             Begin loop over phytoplankton                  !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
       do auto_ind = 1, autotroph_cnt
          thetaC(istr:iend,auto_ind) = Chl_loc(istr:iend,auto_ind) / (C_loc(istr:iend,auto_ind) + epsC)
          Qfe(istr:iend,auto_ind) = autoFe_loc(istr:iend,auto_ind) / (C_loc(istr:iend,auto_ind) + epsC)
          if (kSiO3(auto_ind) > c0) then
             Qsi(istr:iend,auto_ind) = min(Si_loc(istr:iend,auto_ind) / (C_loc(istr:iend,auto_ind) + epsC), gQsi_max)
          end if
       end do
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    !                End loop over phytoplankton                 !
    ! ---------------------------------------------------------- !

    !---------------------------------------------------------------------------
    !   DETERMINE NEW ELEMENTAL RATIOS FOR GROWTH (NEW BIOMASS)
    ! Modify these initial ratios under low ambient iron conditions
    ! Modify the initial si/C ratio under low ambient Si conditions
    !---------------------------------------------------------------------------

    ! ---------------------------------------------------------- !
    !             Begin loop over phytoplankton                  !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
       do auto_ind = 1, autotroph_cnt
          gQfe(istr:iend,auto_ind) = gQfe_0(auto_ind)
          where (Fe_loc < cks * kFe(auto_ind))
             gQfe(istr:iend,auto_ind) =
     &           max((gQfe(istr:iend,auto_ind) * Fe_loc / (cks * kFe(auto_ind))),
     &           gQfe_min(auto_ind))
          end where

          if (kSiO3(auto_ind) > c0) then
             gQsi(istr:iend,auto_ind) = gQsi_0
             where ((Fe_loc < cksi * kFe(auto_ind)) .and. (Fe_loc > c0) .and.
     &           (SiO3_loc > (cksi * kSiO3(auto_ind))))
                gQsi(istr:iend,auto_ind) = min((gQsi(istr:iend,auto_ind) * cksi * kFe(auto_ind) / Fe_loc), gQsi_max)
             end where

             where (Fe_loc == c0)
                gQsi(istr:iend,auto_ind) = gQsi_max
             end where

             where (SiO3_loc < (cksi * kSiO3(auto_ind)))
                gQsi(istr:iend,auto_ind) = max((gQsi(istr:iend,auto_ind) * SiO3_loc / (cksi * kSiO3(auto_ind))),
     &                                  gQsi_min)
             end where
          endif

    !-----------------------------------------------------------------------
    !  QCaCO3 is the percentage of sp organic matter which is associated
    !  with coccolithophores
    !-----------------------------------------------------------------------

          if (CaCO3_ind(auto_ind) > 0) then
#  ifdef BEC_COCCO
             if (exp_calcifier(auto_ind)) then
                QCaCO3(istr:iend,auto_ind) = 1.00  ! 0.433
             else 
                QCaCO3(istr:iend,auto_ind) = 
     &              min( CaCO3_loc(istr:iend,auto_ind) 
     &                 / (C_loc(istr:iend,auto_ind) + epsC), QCaCO3_max )
             endif
#  else /* case not BEC_COCCO  */
             QCaCO3(istr:iend,auto_ind) = min(CaCO3_loc(istr:iend,auto_ind) 
     &                      / (C_loc(istr:iend,auto_ind) + epsC), QCaCO3_max)
#  endif /* BEC_COCCO */
          end if
       end do  ! auto_ind = 1, autotroph_cnt
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    !                End loop over phytoplankton                 !
    ! ---------------------------------------------------------- !

# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(spqcaco3_idx_t))
     &     bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(spqcaco3_idx_t)) = QCaCO3(istr:iend,sp_ind)
#  else
       bec2_diag_3d(istr:iend,j,k,spqcaco3_idx_t) = QCaCO3(istr:iend,sp_ind)
#  endif
# endif



    !------------------------------------------!
    !------------------------------------------!
    !------------------------------------------!
    !   (2.) PAR and phytoplankton shading     !
    !------------------------------------------!
    !------------------------------------------!
    !------------------------------------------!
    
    !---------------------------------------------------------------------------
    !   various k==N initializations
    !
    !   f_qsw_par   fraction of incoming SW -> PAR (non-dim)
    !---------------------------------------------------------------------------

       if (k == N) then
          PAR_out = MAX(c0, f_qsw_par * SHF_QSW)
# ifdef MASKING
     &            * rmask(istr:iend,j)
# endif

# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
          if (bec2_diag_2d_l(parinc_idx_t))
     &        bec2_diag_2d(istr:iend,j,bec2_diag_2d_idx(parinc_idx_t)) = PAR_out(istr:iend)
#  else
          bec2_diag_2d(istr:iend,j,parinc_idx_t) = PAR_out(istr:iend)
#  endif
# else
          PARinc(istr:iend,j) = PAR_out(istr:iend)
# endif

          CALL init_particulate_terms(QA_dust_def,istr,iend,j,dust_flux)

       end if   


    !-----------------------------------------------------------------------
    !  compute PAR related quantities
    !  Morel, Maritorena, JGR, Vol 106, No. C4, pp 7163--7180, 2001
    !-----------------------------------------------------------------------

       do i = istr, iend
          PAR_in(i) = PAR_out(i)
# ifdef BEC_COCCO
          KPARdz(i) = max(spChl_loc(i)+diatChl_loc(i)+diazChl_loc(i)
     &        +coccoChl_loc(i), 0.02)
# else
          KPARdz(i) = max(spChl_loc(i)+diatChl_loc(i)+diazChl_loc(i), 0.02)
# endif
          ! DL: in the following computation of KPARdz, it was taken into 
          ! account that Hz is in m
          if (KPARdz(i) < 0.13224) then
             KPARdz(i) = 0.0919*(KPARdz(i)**0.3536) * Hz(i,j,k)
          else
             KPARdz(i) = 0.1131*(KPARdz(i)**0.4562) * Hz(i,j,k)
          end if
          PAR_out(i) = PAR_in(i) * exp(-KPARdz(i))
          PAR_lay(i) = PAR_in(i) * (c1 - exp(-KPARdz(i))) / KPARdz(i)

# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
          if (bec2_diag_3d_l(par_idx_t))
     &        bec2_diag_3d(i,j,k,bec2_diag_3d_idx(par_idx_t)) = PAR_lay(i)
#  else
          bec2_diag_3d(i,j,k,par_idx_t) = PAR_lay(i)
#  endif
# else
          PAR(i,j,k) = PAR_lay(i)
# endif

       end do



    !------------------------------------------!
    !------------------------------------------!
    !------------------------------------------!
    !   (3.) Temperature-dependent growth      !
    !------------------------------------------!
    !------------------------------------------!
    !------------------------------------------!

    !-----------------------------------------------------------------------
    !  Tref = 30.0 reference temperature (deg. C)
    !
    !  Using q10 formulation with Q10 value of 2.0 (Doney et al., 1996).
    !  growth, mort and grazing rates scaled by Tfunc where they are computed
    !-----------------------------------------------------------------------

       Tfunc = Q_10**(((temp + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)
# ifdef BEC_COCCO
       Tfunc_cocco = 1.45**(((temp + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)
# endif /* BEC_COCCO */
       Tfunc_diat  = 1.55**(((temp + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)



    !--------------------------------------------------!
    !--------------------------------------------------!
    !--------------------------------------------------!
    !   (4.) Depth-dependent mortality thresholds      !
    !--------------------------------------------------!
    !--------------------------------------------------!
    !--------------------------------------------------!

    !-----------------------------------------------------------------------
    !  calculate the loss threshold interpolation factor
    !-----------------------------------------------------------------------

       do i = istr, iend
          if (-z_r(i,j,k) > thres_z1) then
             if (-z_r(i,j,k) < thres_z2) then
                !               = ( 150 - depth ) / ( 150 - 100 )
                f_loss_thres(i) = (thres_z2 + z_r(i,j,k))/(thres_z2 - thres_z1)
             else
                f_loss_thres(i) = c0
             endif
          else
             f_loss_thres(i) = c1
          endif
       end do
       ! f_loss_thres = 1 above 100 metres (thres_z1)
       ! f_loss_thres = 0 below 150 metres (thres_z2)
       ! f_loss_thres > 0 but < 1 between 100 and 150 metres (thres_z1 < z < thres_z2)

    !-----------------------------------------------------------------------
    !  Compute Pprime for all autotrophs, used for loss terms
    !-----------------------------------------------------------------------

       do auto_ind = 1, autotroph_cnt
          where (temp < temp_thres(auto_ind))
             C_loss_thres = f_loss_thres * loss_thres2(auto_ind)
          elsewhere
             C_loss_thres = f_loss_thres * loss_thres(auto_ind)
          end where
          ! C_loss_thres = (/ 0.02, 0.02, 0.02 /) when z < 100m and where T > (/ -10, -10, 16 /)
          ! C_loss_thres < (/ 0.02, 0.02, 0.02 /) btw 100-150m and where T > (/ -10, -10, 16 /)
          ! C_loss_thres = (/ 0.,   0.,   0.   /) when z > 150m and where T > (/ -10, -10, 16 /)
          Pprime(istr:iend,auto_ind) = max(C_loc(istr:iend,auto_ind) - C_loss_thres, c0)
          ! Pprime is zero when biomass is at or below threshold
          ! When biomass is high, Pprime ~= C_loc
          ! Because C_loss_thres is equal to 0 below 150 metres, phytoplankton that are present
          ! beneath 150 metres can be grazed, die and aggregate to extinction. However, above 
          ! 150 metres, phytoplankton losses are zero when biomass is <= C_loss_thres, meaning that
          ! they cannot go to extinction.
       end do



    !--------------------------------------------!
    !--------------------------------------------!
    !--------------------------------------------!
    !   (5.) Ecosystem growth and mortality      !
    !--------------------------------------------!
    !--------------------------------------------!
    !--------------------------------------------!

      !-------------------------------------!
      !   (5a) Nutrient limitation terms    !
      !-------------------------------------!

    !-----------------------------------------------------------------------
    !  Get relative nutrient uptake rates for autotrophs,
    !  min. relative uptake rate modifies C fixation in the manner
    !  that the min. cell quota does in GD98.
    !-----------------------------------------------------------------------

          !--------------------------------------------------------------
          !                        *PJB* 
          !--------------------------------------------------------------
          ! These are limitation terms for phytoplankton growth  
          ! 
          !    VNO3 = [0,1]
          !    VNO2 = [0,1]
          !    VNH4 = [0,1]
          !    VNH4 + VNO2 + VNO3 = VNtot = [0,1]
          ! 
          !    VFe = [0,1]
          ! 
          !    VPO4 = [0,1]
          !    VDOP = [0,1]
          !    VPO4 + VDOP = VPtot = [0,1]
          ! 
          !    VSiO3 = [0,1]
          ! 
          !    f_nut = MIN(VNtot, VPtot, VFe, VSiO3) = [0,1]
          ! 
          !--------------------------------------------------------------
          ! kNH4 = (/ 0.0025, 0.02, 0.0667 /)     ! mmol N/m3
          ! kNO2 = (/ 0.0300, 0.24, 0.8 /)        ! mmol N/m3
          ! kNO3 = (/ 0.2,    0.6,  2.0 /)        ! mmol N/m3
          ! kFe = (/ 0.025e-3, 0.05e-3, 0.025e-3 /) ! mmol Fe/m3
          ! kPO4 = (/ 0.0075, 0.06, 0.015 /)         ! mmol P/m3
          ! kDOP = (/ 0.22, 0.6, 0.05 /)
          ! kSiO3 = (/ 0.0, 0.6, 0.0 /)           ! mmol SiO3/m3



    ! ---------------------------------------------------------- !
    !             Begin loop over phytoplankton                  !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
       do auto_ind = 1, autotroph_cnt
# ifdef Ncycle_SY
          VNO3(istr:iend,auto_ind) = (NO3_loc / kNO3(auto_ind)) /
     &        (c1 + (NO3_loc / kNO3(auto_ind)) + (NO2_loc / kNO2(auto_ind)) + (NH4_loc / kNH4(auto_ind)))
          VNO2(istr:iend,auto_ind) = (NO2_loc / kNO2(auto_ind)) /
     &        (c1 + (NO3_loc / kNO3(auto_ind)) + (NO2_loc / kNO2(auto_ind)) + (NH4_loc / kNH4(auto_ind)))
          VNH4(istr:iend,auto_ind) = (NH4_loc / kNH4(auto_ind)) /
     &        (c1 + (NO3_loc / kNO3(auto_ind)) + (NO2_loc / kNO2(auto_ind)) + (NH4_loc / kNH4(auto_ind)))
          VNtot(istr:iend,auto_ind) = VNO3(istr:iend,auto_ind) + VNO2(istr:iend,auto_ind) + VNH4(istr:iend,auto_ind)
          if (Nfixer(auto_ind)) VNtot(istr:iend,auto_ind) = c1
# else
          VNO3(istr:iend,auto_ind) = (NO3_loc / kNO3(auto_ind)) /
     &        (c1 + (NO3_loc / kNO3(auto_ind)) + (NH4_loc / kNH4(auto_ind)))
          VNH4(istr:iend,auto_ind) = (NH4_loc / kNH4(auto_ind)) /
     &        (c1 + (NO3_loc / kNO3(auto_ind)) + (NH4_loc / kNH4(auto_ind)))
          VNtot(istr:iend,auto_ind) = VNO3(istr:iend,auto_ind) + VNH4(istr:iend,auto_ind)
          if (Nfixer(auto_ind)) VNtot(istr:iend,auto_ind) = c1
# endif

          ! Iron
          VFe = Fe_loc / (Fe_loc + kFe(auto_ind))
          f_nut = min(VNtot(istr:iend,auto_ind), VFe)

          ! Phosphorus
          VPO4(istr:iend,auto_ind) = (PO4_loc / kPO4(auto_ind)) /
     &           (c1 + (PO4_loc / kPO4(auto_ind)) + (DOP_loc / kDOP(auto_ind)))
          VDOP(istr:iend,auto_ind) = (DOP_loc / kDOP(auto_ind)) /
     &           (c1 + (PO4_loc / kPO4(auto_ind)) + (DOP_loc / kDOP(auto_ind)))
          VPtot(istr:iend,auto_ind) = VPO4(istr:iend,auto_ind) + VDOP(istr:iend,auto_ind)
          f_nut = min(f_nut, VPtot(istr:iend,auto_ind))

          ! silicic acid !
          if (kSiO3(auto_ind) > c0) then
             VSiO3 = SiO3_loc / (SiO3_loc + kSiO3(auto_ind))
             f_nut = min(f_nut, VSiO3)
          endif


# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
          if (auto_ind == sp_ind) then
             if (bec2_diag_3d_l(spfeuptake_idx_t))
     &           bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(spfeuptake_idx_t)) = VFe
             if (bec2_diag_3d_l(sppo4uptake_idx_t))
     &           bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(sppo4uptake_idx_t)) = VPO4
          else if (auto_ind == diat_ind) then
             if (bec2_diag_3d_l(diatfeuptake_idx_t))
     &           bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatfeuptake_idx_t)) = VFe
             if (bec2_diag_3d_l(diatpo4uptake_idx_t))
     &           bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatpo4uptake_idx_t)) = VPO4
             if (bec2_diag_3d_l(diatsio3uptake_idx_t))
     &           bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatsio3uptake_idx_t)) = VSiO3
          else
             if (bec2_diag_3d_l(diazfeuptake_idx_t))
     &           bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diazfeuptake_idx_t)) = VFe
             if (bec2_diag_3d_l(diazpo4uptake_idx_t))
     &           bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diazpo4uptake_idx_t)) = VPO4
          end if
#  else
          if (auto_ind == sp_ind) then
             bec2_diag_3d(istr:iend,j,k,spfeuptake_idx_t) = VFe
             bec2_diag_3d(istr:iend,j,k,sppo4uptake_idx_t) =VPO4(istr:iend,sp_ind)
          else if (auto_ind == diat_ind) then
             bec2_diag_3d(istr:iend,j,k,diatfeuptake_idx_t) = VFe
             bec2_diag_3d(istr:iend,j,k,diatpo4uptake_idx_t) =VPO4(istr:iend,diat_ind)
             bec2_diag_3d(istr:iend,j,k,diatsio3uptake_idx_t) = VSiO3
#   ifdef BEC_COCCO
          else if (auto_ind == cocco_ind) then
             bec2_diag_3d(istr:iend,j,k,coccofeuptake_idx_t) = VFe
             bec2_diag_3d(istr:iend,j,k,coccopo4uptake_idx_t) =VPO4(istr:iend,cocco_ind)
#   endif
          else
             bec2_diag_3d(istr:iend,j,k,diazfeuptake_idx_t) = VFe
             bec2_diag_3d(istr:iend,j,k,diazpo4uptake_idx_t) =VPO4(istr:iend,diaz_ind)
          end if
#  endif
# endif


      !---------------------------------------!
      !   (5b) Apply nutrient limitation      !
      !---------------------------------------!

          if (auto_ind==diat_ind) then
             PCmax = PCref(auto_ind) * f_nut * Tfunc_diat
# ifdef BEC_COCCO
          else if (auto_ind==cocco_ind) then
             PCmax = PCref(auto_ind) * f_nut * Tfunc_cocco
# endif /* BEC_COCCO */
          else
             PCmax = PCref(auto_ind) * f_nut * Tfunc
          endif
          where (temp < temp_thres(auto_ind)) PCmax = c0


      !------------------------------------!
      !   (5c) Apply light limitation      !
      !------------------------------------!

          light_lim = c1 - 
     &                exp((-c1 * alphaPI(auto_ind) * thetaC(istr:iend,auto_ind) * PAR_lay) /
     &                      (PCmax + epsTinv))
          PCphoto = PCmax * light_lim
# ifdef BEC_COCCO
          if (auto_ind==cocco_ind) then
             where (PAR_lay < 10) light_lim = c0
             PCphoto = PCmax * light_lim
!   reduce cocco growth at low temperatures
             where (temp < CaCO3_temp_thres1)
     &          PCphoto =PCphoto * max((temp-CaCO3_temp_thres2), c0) /
     &                    (CaCO3_temp_thres1-CaCO3_temp_thres2)
          end if
# endif /* BEC_COCCO */


# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
          if (auto_ind == sp_ind .and. bec2_diag_3d_l(splightlim_idx_t)) then
             bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(splightlim_idx_t)) = light_lim
          else if (auto_ind == diat_ind .and. bec2_diag_3d_l(diatlightlim_idx_t)) then
             bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatlightlim_idx_t)) = light_lim
          else if (bec2_diag_3d_l(diazlightlim_idx_t)) then
             ! Diazotrophs:
             bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diazlightlim_idx_t)) = light_lim
          end if
#  else
          if (auto_ind == sp_ind) then
             bec2_diag_3d(istr:iend,j,k,splightlim_idx_t) = light_lim
          else if (auto_ind == diat_ind) then
             bec2_diag_3d(istr:iend,j,k,diatlightlim_idx_t) = light_lim
#   ifdef BEC_COCCO 
          else if (auto_ind == cocco_ind) then
             bec2_diag_3d(istr:iend,j,k,coccolightlim_idx_t) = light_lim
#   endif
          else
             ! Diazotrophs:
             bec2_diag_3d(istr:iend,j,k,diazlightlim_idx_t) = light_lim
          end if
#  endif
# endif /* BEC2_DIAG */



      !------------------------------------!
      !   (5d) Compute biomass growth      !
      !------------------------------------!

          photoC(istr:iend,auto_ind) = PCphoto * C_loc(istr:iend,auto_ind)


      !-----------------------------------!
      !   (5e) Nutrient uptake rates      !
      !-----------------------------------!

    !-----------------------------------------------------------------------
    !  Get nutrient uptakes by autotrophs based on calculated C fixation
    !  total N uptake VNC is used in photoadaption
    !-----------------------------------------------------------------------

       ! Nitrogen !
# ifdef Ncycle_SY
       epsN = 1.0e-30 ! small value to avoid dividing by 0
          where (VNtot(istr:iend,auto_ind) > c0)
             where (NO3_loc .le. c0)
                NO3_V(istr:iend,auto_ind) = c0
             elsewhere
                NO3_V(istr:iend,auto_ind) = (VNO3(istr:iend,auto_ind) / VNtot(istr:iend,auto_ind)) 
     &                                        * photoC(istr:iend,auto_ind) * Q 
             end where
             where (NO2_loc .le. c0)
                NO2_V(istr:iend,auto_ind) = c0
             elsewhere
                NO2_V(istr:iend,auto_ind) = (VNO2(istr:iend,auto_ind) / VNtot(istr:iend,auto_ind)) 
     &                                        * photoC(istr:iend,auto_ind) * Q 
             end where
             NH4_V(istr:iend,auto_ind) = (VNH4(istr:iend,auto_ind) / VNtot(istr:iend,auto_ind))
     &                    * photoC(istr:iend,auto_ind) * Q

             VNC = PCphoto * Q
          elsewhere
             NO3_V(istr:iend,auto_ind) = c0
             NO2_V(istr:iend,auto_ind) = c0
             NH4_V(istr:iend,auto_ind) = c0
             VNC = c0
          end where
# else
          where (VNtot(istr:iend,auto_ind) > c0)
             NO3_V(istr:iend,auto_ind) = (VNO3(istr:iend,auto_ind) / VNtot(istr:iend,auto_ind))
     &                   * photoC(istr:iend,auto_ind) * Q
             NH4_V(istr:iend,auto_ind) = (VNH4(istr:iend,auto_ind) / VNtot(istr:iend,auto_ind))
     &                   * photoC(istr:iend,auto_ind) * Q
             VNC = PCphoto * Q
          elsewhere
             NO3_V(istr:iend,auto_ind) = c0
             NH4_V(istr:iend,auto_ind) = c0
             VNC = c0
          end where
# endif

       ! Phosphorus !
          where (VPtot(istr:iend,auto_ind) > c0)
             PO4_V(istr:iend,auto_ind) = (VPO4(istr:iend,auto_ind) / VPtot(istr:iend,auto_ind)) 
     &                    * photoC(istr:iend,auto_ind) * Qp(auto_ind)
             DOP_V(istr:iend,auto_ind) = (VDOP(istr:iend,auto_ind) / VPtot(istr:iend,auto_ind)) 
     &                    * photoC(istr:iend,auto_ind) * Qp(auto_ind)
          elsewhere
             PO4_V(istr:iend,auto_ind) = c0
             DOP_V(istr:iend,auto_ind) = c0
          end where

          ! Iron
          photoFe(istr:iend,auto_ind) = photoC(istr:iend,auto_ind) * gQfe(istr:iend,auto_ind)

          ! Silicate !
          if (kSiO3(auto_ind) > c0) then
             photoSi(istr:iend,auto_ind) = photoC(istr:iend,auto_ind) * gQsi(istr:iend,auto_ind)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
             if (auto_ind == diat_ind .and. bec2_diag_3d_l(diatsiuptake_idx_t)) then
                bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatsiuptake_idx_t)) = photoSi(istr:iend,auto_ind)
             end if
#  else
             if (auto_ind == diat_ind) then
                bec2_diag_3d(istr:iend,j,k,diatsiuptake_idx_t) = photoSi(istr:iend,auto_ind)
             end if
#  endif
# endif
          endif

# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
          if (bec2_diag_3d_l(ironuptakesp_idx_t+auto_ind-1))
     &        bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(ironuptakesp_idx_t+auto_ind-1)) = photoFe(istr:iend,auto_ind)
#  else
#   ifdef BEC_COCCO
          if (auto_ind == sp_ind) then
            bec2_diag_3d(istr:iend,j,k,ironuptakesp_idx_t) = photoFe(istr:iend,sp_ind)
          else if (auto_ind == diat_ind) then
            bec2_diag_3d(istr:iend,j,k,ironuptakediat_idx_t) = photoFe(istr:iend,diat_ind)
          else if (auto_ind == cocco_ind) then
            bec2_diag_3d(istr:iend,j,k,ironuptakecocco_idx_t) = photoFe(istr:iend,cocco_ind)
          else
            bec2_diag_3d(istr:iend,j,k,ironuptakediaz_idx_t) = photoFe(istr:iend,diaz_ind)
          endif
#   else /* case not BEC_COCCO */
          bec2_diag_3d(istr:iend,j,k,ironuptakesp_idx_t+auto_ind-1) = photoFe(istr:iend,auto_ind)
#   endif
#  endif
# endif



      !-----------------------------!
      !   (5f) Photoadaptation      !
      !-----------------------------!

    !-----------------------------------------------------------------------
    !  calculate pChl, (used in photoadapt., GD98)
    !  2.3   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
    !  GD 98 Chl. synth. term
    !-----------------------------------------------------------------------

          WORK1 = alphaPI(auto_ind) * thetaC(istr:iend,auto_ind) * PAR_lay

    ! CN (April 2015): changed from c0 to value slightly larger than c0 (parameter
    ! PAR_thres_pChl, BEC blew up in SO setup otherwise):
          where (WORK1 > PAR_thres_pChl)
             pChl = thetaN_max(auto_ind) * PCphoto / WORK1
             photoacc(istr:iend,auto_ind) = (pChl * VNC / thetaC(istr:iend,auto_ind))
     &                       * Chl_loc(istr:iend,auto_ind)
          elsewhere
             photoacc(istr:iend,auto_ind) = c0
          end where


# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
          if (bec2_diag_3d_l(spphotoacc_idx_t))
     &        bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(spphotoacc_idx_t)) = photoacc(istr:iend,sp_ind)
          if (bec2_diag_3d_l(diatphotoacc_idx_t))
     &        bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatphotoacc_idx_t)) = photoacc(istr:iend,diat_ind)
          if (bec2_diag_3d_l(diazphotoacc_idx_t))
     &        bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diazphotoacc_idx_t)) = photoacc(istr:iend,diaz_ind)
#  else
          bec2_diag_3d(istr:iend,j,k,spphotoacc_idx_t) = photoacc(istr:iend,sp_ind)
          bec2_diag_3d(istr:iend,j,k,diatphotoacc_idx_t) = photoacc(istr:iend,diat_ind)
          bec2_diag_3d(istr:iend,j,k,diazphotoacc_idx_t) = photoacc(istr:iend,diaz_ind)
#   ifdef BEC_COCCO
          bec2_diag_3d(istr:iend,j,k,coccophotoacc_idx_t) = photoacc(istr:iend,cocco_ind)
#   endif
#  endif
# endif



      !---------------------------------------!
      !   (5f) Implicit CaCO3 production      !
      !---------------------------------------!


    !-----------------------------------------------------------------------
    !  implicit CaCO3 Production, parameterized as function of small phyto production:
    !    - decrease CaCO3 as function of nutrient limitation
    !    - decrease CaCO3 prod at low temperatures
    !    - increase CaCO3 prod under bloom conditions
    !  Maximum calcification rate is 40% of primary production.
    !
    !  explicit CaCO3 Production, parametrized as function of photoC of PFT
    !    think about: temp dependence, f_nut dependence, CaCO3:Corg
    !    ratio (right now set to one)
    !-----------------------------------------------------------------------

          if (imp_calcifier(auto_ind)) then
             CaCO3_PROD(istr:iend,auto_ind) = parm_f_prod_sp_CaCO3 * photoC(istr:iend,auto_ind)
             CaCO3_PROD(istr:iend,auto_ind) = CaCO3_PROD(istr:iend,auto_ind) * f_nut

             where (temp < CaCO3_temp_thres1)
     &          CaCO3_PROD(istr:iend,auto_ind) = CaCO3_PROD(istr:iend,auto_ind) * max((temp-CaCO3_temp_thres2), c0) /
     &                    (CaCO3_temp_thres1-CaCO3_temp_thres2)

             where (C_loc(istr:iend,auto_ind) > CaCO3_sp_thres)
     &          CaCO3_PROD(istr:iend,auto_ind) = min((CaCO3_PROD(istr:iend,auto_ind)
     &                     * C_loc(istr:iend,auto_ind) / CaCO3_sp_thres),
     &                    (f_photosp_CaCO3 * photoC(istr:iend,auto_ind)))
          elseif (exp_calcifier(auto_ind)) then
!!  calcification: use different light_lim than for photosynthesis (see
!!  Zondervan 2007)
!             QCaCO3(istr:iend,auto_ind) = QCaCO3(istr:iend,auto_ind)*light_lim_calc
!             CaCO3_PROD(istr:iend,auto_ind) = QCaCO3(istr:iend,auto_ind)*photoC(istr:iend,auto_ind)*(1-light_lim)             
             CaCO3_PROD(istr:iend,auto_ind) = QCaCO3(istr:iend,auto_ind)*photoC(istr:iend,auto_ind)
          !   if (k==42) then
          !   if (mynode==222) then
          !     print*,'surface, mynode==222:'
          !     print*,'QCaCO3',QCaCO3(istr:iend,auto_ind)
          !     print*,'photoC',photoC(istr:iend,auto_ind)
          !     print*,'photoC',photoC(istr:iend,cocco_ind)
          !     print*,'CaCO3_PROD',CaCO3_PROD(istr:iend,auto_ind)
          !   endif
          !   endif
          else
             CaCO3_PROD(istr:iend,auto_ind) = c0
          endif
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
          if (auto_ind == sp_ind .and. bec2_diag_3d_l(caco3prod_idx_t)) then
             bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(caco3prod_idx_t)) = CaCO3_PROD(istr:iend,auto_ind)
          end if
#  else
          if (auto_ind == sp_ind) then
             bec2_diag_3d(istr:iend,j,k,caco3prod_idx_t) = CaCO3_PROD(istr:iend,auto_ind)
          end if
#   ifdef BEC_COCCO
          if (auto_ind == cocco_ind) then
             bec2_diag_3d(istr:iend,j,k,caco3prodcocco_idx_t) = CaCO3_PROD(istr:iend,cocco_ind)
          end if
#   endif
#  endif
# endif

      !-------------------------------------------------------!
      !   (5g) Autotroph mortality and aggregation losses     !
      !-------------------------------------------------------!

          !--------------------------------------------------------------
          !                        *PJB* 
          !--------------------------------------------------------------
          ! Mortality of the phytoplankton groups takes the form of:
          !    mortality = mort * B(phyto)
          ! 
          !   where...
          ! mort = (/ 0.1*dps, 0.1*dps, 0.1*dps /) * Tfunc 
          ! Tfunc is [0.2 --> 1.0] from 0 to 30 degC
          !
          ! For diatoms, Tfunc is higher at lower temperatures:
          ! Tfunc is [0.27 --> 1.0] from 0 to 30 degC
          !
          !--------------------------------------------------------------
       if (auto_ind == diat_ind) then
          auto_loss(istr:iend,auto_ind) =  mort(auto_ind) * Pprime(istr:iend,auto_ind) * Tfunc_diat
# ifdef BEC_COCCO
       elseif (auto_ind == cocco_ind) then
          auto_loss(istr:iend,auto_ind) =  mort(auto_ind) * Pprime(istr:iend,auto_ind) * Tfunc_cocco
# endif /* BEC_COCCO */
       else
          auto_loss(istr:iend,auto_ind) =  mort(auto_ind) * Pprime(istr:iend,auto_ind) * Tfunc
       endif


          !--------------------------------------------------------------
          !                        *PJB* 
          !--------------------------------------------------------------
          ! Aggregation of the phytoplankton groups takes the form of:
          !    aggregation = mort2 * B(phyto)**2
          ! 
          !   where...
          ! mort2 = (/ 0.01*dps, 0.01*dps, 0.01*dps /)   Quadratic mort
          !
          ! However, this aggregation rate cannot exceed a maximum rate or
          ! a minimum rate, set as equal to:
          !    max_aggregation = agg_max * B(phyto)
          !    min_aggregation = agg_min * B(phyto)
          !
          !   where...
          ! agg_max = (/ 0.5,  0.5,  0.5  /)
          ! agg_min = (/ 0.01, 0.02, 0.01 /)
          !
          ! This means that losses through aggregation cannot exceed half
          ! of biomass, and losses through aggregation must always occur 
          ! at 1% of biomass when biomass is > than the minimum threshold
          ! set earlier by C_loss_thres
          !
          ! However, it is very unlikely that aggregation losses will
          ! approach 50 % of biomass, because for this to occur the 
          ! biomass of phytoplankton would need to be 50 mmol C / m3 
          !
          !--------------------------------------------------------------
          auto_agg(istr:iend,auto_ind) = min(agg_rate_max(auto_ind) * dps * Pprime(istr:iend,auto_ind), 
     &                               mort2(auto_ind) * Pprime(istr:iend,auto_ind) * Pprime(istr:iend,auto_ind))
          auto_agg(istr:iend,auto_ind) = max((agg_rate_min(auto_ind) * dps) * Pprime(istr:iend,auto_ind),
     &                               auto_agg(istr:iend,auto_ind))

      !--------------------------------------!
      !   (5h) Autotroph grazing losses      !
      !--------------------------------------!

          ! ---------------------------------------------------------------------------------
          !                        *PJB* 
          ! ---------------------------------------------------------------------------------
          ! Grazing by zooplankton takes the simple form of:
          !    total_grazing = mu(zoo) * B(zoo) *  B(prey) / ( B(prey) + K )
          !
          !   where...
          ! mu = parm_zoo_mumax = (/ 3.3*dps, 3.05*dps, 3.3*dps /) * Tfunc 
          ! K  = parm_kzoo_phy    = (/ 1.2, 1.2, 1.2 /)      ! grazing coefficient (mmol C/m^3)
          ! Tfunc is [0.2 --> 1.0] from 0 to 30 degC
          !
          ! Grazing losses to a specific phyto is calculated as a proportion of their biomass
          ! relative to the total biomass of phytoplankton, such that:
          !    diatom_grazing = (B(diat) / B(phyto)) * total_grazing
          !
          ! ---------------------------------------------------------------------------------

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Meso-zooplankton eats phytoplankton and chemoautotrophs !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !!! 1. RETRIEVE THE BIOMASS OF ALL PREY !!!
          WORK1 = c0
# ifdef EXPLICIT_MICROBES
          AOAprime = MAX(AOA_loc - parm_che_bmin, c0)
          NOBprime = MAX(NOB_loc - parm_che_bmin, c0)
          AOXprime = MAX(AOX_loc - parm_che_bmin, c0)
# else 
          AOAprime = c0
          NOBprime = c0
          AOXprime = c0
# endif /* EXPLICIT_MICROBES */
          WORK1 = WORK1 + sum(Pprime(istr:iend,:), dim=2) + AOAprime + NOBprime + AOXprime
          ! O2 limitation term
          z_o2lim = min(1.0, max(0.0, (exp(-O2_loc/10.0)) ))
          z_umax = z_umax_0(auto_ind) * Tfunc * (c1-z_o2lim)
          where (WORK1 > c0)
             auto_graze(istr:iend,auto_ind) = z_umax * zooC_loc *
     &                                        (Pprime(istr:iend,auto_ind) / (WORK1 + z_grz(auto_ind)))
             ! Grazing on chemoautotrophs is half as fast as grazing on phytoplankton
             zzooaoa = z_umax_0(sp_ind)*0.5 * Tfunc * (c1-z_o2lim) * (AOAprime(istr:iend) / (WORK1 + parm_kzoo_aoa))
             zzoonob = z_umax_0(sp_ind)*0.5 * Tfunc * (c1-z_o2lim) * (NOBprime(istr:iend) / (WORK1 + parm_kzoo_nob))
             zzooaox = z_umax_0(sp_ind)*0.5 * Tfunc * (c1-z_o2lim) * (AOXprime(istr:iend) / (WORK1 + parm_kzoo_aox))
          elsewhere
             auto_graze(istr:iend,auto_ind) = c0
             zzooaoa = c0
             zzoonob = c0
             zzooaox = c0
          end where


# ifdef EXPLICIT_MICROBES 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          ! Nano zooplankton eats bacterial heterotrophs ! 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          AERprime = MAX(AER_loc - parm_het_bmin, c0)
          NARprime = MAX(NAR_loc - parm_het_bmin, c0)
          NAIprime = MAX(NAI_loc - parm_het_bmin, c0)
          NIRprime = MAX(NIR_loc - parm_het_bmin, c0)
          NIOprime = MAX(NIO_loc - parm_het_bmin, c0)
          NOSprime = MAX(NOS_loc - parm_het_bmin, c0)
          NAOprime = MAX(NAO_loc - parm_het_bmin, c0)

          WORK1 = AERprime + NARprime + NAIprime + NIRprime + NIOprime + NOSprime + NAOprime
          z_umax = parm_szoo_mumax * Tfunc

          where (WORK1 > c0)
             zszooaer = z_umax * szooC_loc * (AERprime / (WORK1 + parm_kszoo_aer))
             zszoonar = z_umax * szooC_loc * (NARprime / (WORK1 + parm_kszoo_nar))
             zszoonai = z_umax * szooC_loc * (NAIprime / (WORK1 + parm_kszoo_nai))
             zszoonir = z_umax * szooC_loc * (NIRprime / (WORK1 + parm_kszoo_nir))
             zszoonio = z_umax * szooC_loc * (NIOprime / (WORK1 + parm_kszoo_nio))
             zszoonos = z_umax * szooC_loc * (NOSprime / (WORK1 + parm_kszoo_nos))
             zszoonao = z_umax * szooC_loc * (NAOprime / (WORK1 + parm_kszoo_nao))
          elsewhere
             zszooaer = c0
             zszoonar = c0
             zszoonai = c0
             zszoonir = c0
             zszoonio = c0
             zszoonos = c0
             zszoonao = c0
          end where
# endif /* EXPLICIT_MICROBES */ 

      !---------------------------------------------!
      !   (5i) Nitrogen fixation by diazotrophs     !
      !---------------------------------------------!

          ! ----------------------------------------------------------------
          !                        *PJB* 
          ! ----------------------------------------------------------------
          !  N-fixation of diazotrophs takes the form:
          !    Nfix = (photoC * Q * r_Nfix_photo)
          !
          !   where...
          ! photoC       = ...     ! C-fixation (not N-limited)
          ! Q            = 0.137   ! N/C ratio (mmol/mmol) of phyto & zoo 
          ! r_Nfix_photo = 1.25    ! N fix relative to C fix (non-dim)
          !
          ! 
          ! This total Nfix term is then altered by accounting for the N
          ! that is taken up from the environment (NH4 + NO2 + NO3), versus
          ! that which is fixed internally via N2.
          !  
          ! Hence, actual N2 --> Norg is: 
          !    Nfix = (photoC * Q * r_Nfix_photo)
          !           - NO3_V - NO2_V - NH4_V
          !
          !   where...
          ! NO3_V  = [0-1]*Nfix   ! Fraction of N-fixed biomass via NO3 uptake
          ! NO2_V  = [0-1]*Nfix   ! Fraction of N-fixed biomass via NO2 uptake 
          ! NH4_V  = [0-1]*Nfix   ! Fraction of N-fixed biomass via NH4 uptake
          !
          ! The sum of NO3_V + NO2_V + NH4_V is only near 1 when N is high,
          ! and even then, there is always a small fraction of N2 --> Norg. 
          ! If there was no lower temperature threshold for diazotroph growth,
          ! Then this would mean that nitrogen fixation would happen at high
          ! latitudes, albeit at low rates.
          !
          ! N-excrete is always 0.25 of total N-fixation 
          !   [set by "r_Nfix_photo"]
          ! ----------------------------------------------------------------

# ifdef Ncycle_SY
          if (Nfixer(auto_ind)) then
             WORK1 = photoC(istr:iend,auto_ind) * Q
             Nfix(istr:iend,auto_ind)     = (WORK1 * r_Nfix_photo) - NO3_V(istr:iend,auto_ind) - NO2_V(istr:iend,auto_ind)
     &                 - NH4_V(istr:iend,auto_ind)
             ! This Nexcrete calculation seems silly to me, since it works out to 
             !  always be 25% of the total N2 --> Norg term (PJB)
             Nexcrete(istr:iend,auto_ind) = Nfix(istr:iend,auto_ind) + NO3_V(istr:iend,auto_ind) + NO2_V(istr:iend,auto_ind)
     &                 + NH4_V(istr:iend,auto_ind) - WORK1
          endif
# else
          if (Nfixer(auto_ind)) then
             WORK1 = photoC(istr:iend,auto_ind) * Q
             Nfix(istr:iend,auto_ind)     = (WORK1 * r_Nfix_photo) - NO3_V(istr:iend,auto_ind)
     &                 - NH4_V(istr:iend,auto_ind)
             Nexcrete(istr:iend,auto_ind) = Nfix(istr:iend,auto_ind) + NO3_V(istr:iend,auto_ind)
     &                 + NH4_V(istr:iend,auto_ind) - WORK1
          endif
#endif


# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
          if (auto_ind == diaz_ind .and. bec2_diag_3d_l(diaznfix_idx_t)) then
             bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diaznfix_idx_t)) = Nfix(istr:iend,auto_ind)
          end if
#  else
          if (auto_ind == diaz_ind) then
             bec2_diag_3d(istr:iend,j,k,diaznfix_idx_t) = Nfix(istr:iend,auto_ind)
          end if
#  endif
# endif


      !---------------------------------------------------------!
      !   (5j) Route grazing terms to Zoo, POC, DIC and DOC     !
      !---------------------------------------------------------!
      !---------------------------------------------------------!
      !  NOTES:
      !   - all aggregation goes to POC
      !   - currently assumes that 33% of grazed caco3 is remineralized
      !   - if autotrophs(sp_ind)%graze_zoo ever changes, coefficients on routing grazed sp must change!
      !   - min.%C routed to POC from grazing for ballast requirements = 0.4 * Qcaco3
      !   - min.%C routed from sp_loss = 0.59 * QCaCO3, or P_CaCO3%rho
      !   - if autotrophs(diat_ind)%graze_zoo is changed, coeff.s for poc,doc and dic must change!
      !-----------------------------------------------------------------------
      
          ! ----------------------------------------------------------------
          !                        *PJB* 
          ! ----------------------------------------------------------------
          !  
          !    Zooplankton grazing loss of phytoplankton is routed to 4 pools:
          !      - Zoo
          !      - POC
          !      - DOC
          !      - DIC
          !
          !   where...
          ! graze_zoo = (/ 0.30, 0.25, 0.30 /)
          ! graze_doc = (/ 0.06, 0.06, 0.06 /)
          ! graze_poc = (/ 0.00, 0.40, 0.10 /)
          !                  |  
          !                  |
          !                  V
          !                 This is actually not zero when we  
          !                 consider implicit calcifiers
          !  
          ! The remaining fraction of grazing is routed to DIC
          !                SP    Diat  Diaz
          !                64%   29%   54%
          !
          ! Temporary routing of bacterial/archaeal biomass to ZOO,
          ! POC, DOC and DIC the same as for small phytoplankton. This
          ! assumes that bacteria of similar size to Prochlorococcus (SP)
          !
          ! ----------------------------------------------------------------

          ! 1. Route grazing of phytoplankton to Mesozooplpankton
          auto_graze_zoo(istr:iend,auto_ind) = graze_zoo(auto_ind) * auto_graze(istr:iend,auto_ind)
          ! 2. Route grazing of phytoplankton to POC
          if (imp_calcifier(auto_ind)) then
             auto_graze_poc(istr:iend,auto_ind) = auto_graze(istr:iend,auto_ind)
     &       * max((caco3_poc_min * QCaCO3(istr:iend,auto_ind)),
     &             min(spc_poc_fac * max(1.0,Pprime(istr:iend,auto_ind)), f_graze_sp_poc_lim))
          else     
             auto_graze_poc(istr:iend,auto_ind) = graze_poc(auto_ind) * auto_graze(istr:iend,auto_ind)
          endif
          ! 3. Route grazing of phytoplankton to DOC
          auto_graze_doc(istr:iend,auto_ind) = graze_doc(auto_ind) * auto_graze(istr:iend,auto_ind)
          ! 4. Route grazing of phytoplankton to DIC
          auto_graze_dic(istr:iend,auto_ind) = auto_graze(istr:iend,auto_ind)
     &        - (auto_graze_zoo(istr:iend,auto_ind) + auto_graze_poc(istr:iend,auto_ind)
     &        + auto_graze_doc(istr:iend,auto_ind)) 

     
     
          ! 1. Route loss to POC
          if (imp_calcifier(auto_ind)) then
             auto_loss_poc(istr:iend,auto_ind) = QCaCO3(istr:iend,auto_ind) * auto_loss(istr:iend,auto_ind)
          else
             auto_loss_poc(istr:iend,auto_ind) = loss_poc(auto_ind) * auto_loss(istr:iend,auto_ind)
          endif
          ! 2. Route loss to DOC
          auto_loss_doc(istr:iend,auto_ind) = (c1 - parm_labile_ratio) * (auto_loss(istr:iend,auto_ind)
     &         - auto_loss_poc(istr:iend,auto_ind)) 
          ! 3. Route loss to DIC
          auto_loss_dic(istr:iend,auto_ind) = parm_labile_ratio * (auto_loss(istr:iend,auto_ind)
     &        - auto_loss_poc(istr:iend,auto_ind))

    !-----------------------------------------------------------------------
    ! P from some autotrophs w/ Qp different from Qp_zoo_pom must be routed differently than other
    ! elements to ensure that sinking detritus and zooplankton pools get their fixed P/C ratios.
    ! The remaining P is split evenly between DOP and PO4.
    !-----------------------------------------------------------------------

          if (Qp(auto_ind) /= Qp_zoo_pom) then
             remaining_P = ((auto_graze(istr:iend,auto_ind) + auto_loss(istr:iend,auto_ind) 
     &                  + auto_agg(istr:iend,auto_ind)) * Qp(auto_ind))
     &                  - ((auto_graze_zoo(istr:iend,auto_ind)) * Qp_zoo_pom)
     &                  - ((auto_graze_poc(istr:iend,auto_ind) + auto_loss_poc(istr:iend,auto_ind) 
     &                  + auto_agg(istr:iend,auto_ind)) * Qp_zoo_pom)

             remaining_P_dop(istr:iend,auto_ind) = (c1 - parm_labile_ratio) * remaining_P
             remaining_P_dip(istr:iend,auto_ind) = parm_labile_ratio * remaining_P
          endif


       end do  ! auto_ind = 1, autotroph_cnt
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    !               End loop over phytoplankton                  !
    ! ---------------------------------------------------------- !


# ifdef EXPLICIT_MICROBES
      !-------------------------------------!
      !   (6) Route grazing of microbes to: !
      !        - Grazers (ZOO and SZOO)     !
      !        - POC, DOC and DIC           !
      !        - all Fe goes to dFe or POFe !
      !-------------------------------------!

       ! 1. Route grazing to Zooplpankton (0.3)
       zzooaoa_zoo = graze_zoo(sp_ind) * zzooaoa
       zzoonob_zoo = graze_zoo(sp_ind) * zzoonob
       zzooaox_zoo = graze_zoo(sp_ind) * zzooaox
       zszooaer_szoo = graze_zoo(sp_ind) * zszooaer
       zszoonar_szoo = graze_zoo(sp_ind) * zszoonar
       zszoonai_szoo = graze_zoo(sp_ind) * zszoonai
       zszoonir_szoo = graze_zoo(sp_ind) * zszoonir
       zszoonio_szoo = graze_zoo(sp_ind) * zszoonio
       zszoonos_szoo = graze_zoo(sp_ind) * zszoonos
       zszoonao_szoo = graze_zoo(sp_ind) * zszoonao

       ! 2. Route grazing to POC (0.0) THIS SHOULD BE ZERO
       zzooaoa_poc = graze_poc(sp_ind) * zzooaoa
       zzoonob_poc = graze_poc(sp_ind) * zzoonob
       zzooaox_poc = graze_poc(sp_ind) * zzooaox
       zszooaer_poc = graze_poc(sp_ind) * zszooaer
       zszoonar_poc = graze_poc(sp_ind) * zszoonar
       zszoonai_poc = graze_poc(sp_ind) * zszoonai
       zszoonir_poc = graze_poc(sp_ind) * zszoonir
       zszoonio_poc = graze_poc(sp_ind) * zszoonio
       zszoonos_poc = graze_poc(sp_ind) * zszoonos
       zszoonao_poc = graze_poc(sp_ind) * zszoonao

       ! 3. Route grazing to DOC, DON, DOP (not DOFe) (0.06)
       zzooaoa_doc = graze_doc(sp_ind) * zzooaoa
       zzoonob_doc = graze_doc(sp_ind) * zzoonob
       zzooaox_doc = graze_doc(sp_ind) * zzooaox
       zszooaer_doc = graze_doc(sp_ind) * zszooaer
       zszoonar_doc = graze_doc(sp_ind) * zszoonar
       zszoonai_doc = graze_doc(sp_ind) * zszoonai
       zszoonir_doc = graze_doc(sp_ind) * zszoonir
       zszoonio_doc = graze_doc(sp_ind) * zszoonio
       zszoonos_doc = graze_doc(sp_ind) * zszoonos
       zszoonao_doc = graze_doc(sp_ind) * zszoonao

       ! 4. Route grazing to DIC, NH4, PO4, dFe (1 - 0.3 - 0.0 - 0.06 = 0.64)
       zzooaoa_dic = zzooaoa - zzooaoa_zoo - zzooaoa_poc - zzooaoa_doc
       zzoonob_dic = zzoonob - zzoonob_zoo - zzoonob_poc - zzoonob_doc
       zzooaox_dic = zzooaox - zzooaox_zoo - zzooaox_poc - zzooaox_doc
       zszooaer_dic = zszooaer - zszooaer_szoo - zszooaer_poc - zszooaer_doc
       zszoonar_dic = zszoonar - zszoonar_szoo - zszoonar_poc - zszoonar_doc
       zszoonai_dic = zszoonai - zszoonai_szoo - zszoonai_poc - zszoonai_doc
       zszoonir_dic = zszoonir - zszoonir_szoo - zszoonir_poc - zszoonir_doc
       zszoonio_dic = zszoonio - zszoonio_szoo - zszoonio_poc - zszoonio_doc
       zszoonos_dic = zszoonos - zszoonos_szoo - zszoonos_poc - zszoonos_doc
       zszoonao_dic = zszoonao - zszoonao_szoo - zszoonao_poc - zszoonao_doc

       ! Because zooplankton and microbes can have different elemental stoichiometries, 
       ! we need to account for the spill over effect of zooplankton consuming N-rich, 
       ! P-rich or Fe-rich foods. If consuming rich foods (high quality), then the 
       ! zooplankton will excrete the excess N, P and Fe --> DOM. If consuming poor quality 
       ! food (C:N of prey > C:N of Zoo), then these numbers will be negative, and the 
       ! zooplankton will need to consume extra prey. 
       ! However, supplementary feeding is yet to be coded, so:
       !    C:N, C:P and C:Fe of Microbes <= C:N, C:P and C:Fe of Zoo
       zzooaoa_excN = zzooaoa_zoo/parm_aoa_CN - zzooaoa_zoo*Q +
     &                zzooaoa_poc/parm_aoa_CN - zzooaoa_poc*Q
       zzoonob_excN = zzoonob_zoo/parm_nob_CN - zzoonob_zoo*Q +
     &                zzoonob_poc/parm_nob_CN - zzoonob_poc*Q
       zzooaox_excN = zzooaox_zoo/parm_aox_CN - zzooaox_zoo*Q +
     &                zzooaox_poc/parm_aox_CN - zzooaox_poc*Q
       zzooaoa_excP = zzooaoa_zoo/parm_aoa_CP - zzooaoa_zoo*Qp_zoo_pom +
     &                zzooaoa_poc/parm_aoa_CP - zzooaoa_poc*Qp_zoo_pom
       zzoonob_excP = zzoonob_zoo/parm_nob_CP - zzoonob_zoo*Qp_zoo_pom +
     &                zzoonob_poc/parm_nob_CP - zzoonob_poc*Qp_zoo_pom
       zzooaox_excP = zzooaox_zoo/parm_aox_CP - zzooaox_zoo*Qp_zoo_pom +
     &                zzooaox_poc/parm_aox_CP - zzooaox_poc*Qp_zoo_pom
       zzooaoa_excFe = zzooaoa_zoo/parm_aoa_CFe - zzooaoa_zoo*Qfe_zoo  ! extra iron excreted during zoo feeding,
       zzoonob_excFe = zzoonob_zoo/parm_nob_CFe - zzoonob_zoo*Qfe_zoo  ! but particles remain Fe-rich  
       zzooaox_excFe = zzooaox_zoo/parm_aox_CFe - zzooaox_zoo*Qfe_zoo

       zszooaer_excN = zszooaer_szoo/parm_aer_CN - zszooaer_szoo*Q +
     &                 zszooaer_poc/parm_aer_CN - zszooaer_poc*Q
       zszoonar_excN = zszoonar_szoo/parm_nar_CN - zszoonar_szoo*Q +
     &                 zszoonar_poc/parm_nar_CN - zszoonar_poc*Q
       zszoonai_excN = zszoonai_szoo/parm_nai_CN - zszoonai_szoo*Q +
     &                 zszoonai_poc/parm_nai_CN - zszoonai_poc*Q
       zszoonir_excN = zszoonir_szoo/parm_nir_CN - zszoonir_szoo*Q +
     &                 zszoonir_poc/parm_nir_CN - zszoonir_poc*Q
       zszoonio_excN = zszoonio_szoo/parm_nio_CN - zszoonio_szoo*Q +
     &                 zszoonio_poc/parm_nio_CN - zszoonio_poc*Q
       zszoonos_excN = zszoonos_szoo/parm_nos_CN - zszoonos_szoo*Q +
     &                 zszoonos_poc/parm_nos_CN - zszoonos_poc*Q
       zszoonao_excN = zszoonao_szoo/parm_nao_CN - zszoonao_szoo*Q +
     &                 zszoonao_poc/parm_nao_CN - zszoonao_poc*Q
       zszooaer_excP = zszooaer_szoo/parm_aer_CP - zszooaer_szoo*Qp_zoo_pom +
     &                 zszooaer_poc/parm_aer_CP - zszooaer_poc*Qp_zoo_pom
       zszoonar_excP = zszoonar_szoo/parm_nar_CP - zszoonar_szoo*Qp_zoo_pom +
     &                 zszoonar_poc/parm_nar_CP - zszoonar_poc*Qp_zoo_pom
       zszoonai_excP = zszoonai_szoo/parm_nai_CP - zszoonai_szoo*Qp_zoo_pom +
     &                 zszoonai_poc/parm_nai_CP - zszoonai_poc*Qp_zoo_pom
       zszoonir_excP = zszoonir_szoo/parm_nir_CP - zszoonir_szoo*Qp_zoo_pom +
     &                 zszoonir_poc/parm_nir_CP - zszoonir_poc*Qp_zoo_pom
       zszoonio_excP = zszoonio_szoo/parm_nio_CP - zszoonio_szoo*Qp_zoo_pom +
     &                 zszoonio_poc/parm_nio_CP - zszoonio_poc*Qp_zoo_pom
       zszoonos_excP = zszoonos_szoo/parm_nos_CP - zszoonos_szoo*Qp_zoo_pom +
     &                 zszoonos_poc/parm_nos_CP - zszoonos_poc*Qp_zoo_pom
       zszoonao_excP = zszoonao_szoo/parm_nao_CP - zszoonao_szoo*Qp_zoo_pom +
     &                 zszoonao_poc/parm_nao_CP - zszoonao_poc*Qp_zoo_pom
       zszooaer_excFe = zszooaer_szoo/parm_aer_CFe - zszooaer_szoo*Qfe_zoo  ! extra iron excreted during zoo feeding,
       zszoonar_excFe = zszoonar_szoo/parm_nar_CFe - zszoonar_szoo*Qfe_zoo  ! extra iron excreted during zoo feeding,
       zszoonai_excFe = zszoonai_szoo/parm_nai_CFe - zszoonai_szoo*Qfe_zoo  ! extra iron excreted during zoo feeding,
       zszoonir_excFe = zszoonir_szoo/parm_nir_CFe - zszoonir_szoo*Qfe_zoo  ! but particles remain Fe-rich
       zszoonio_excFe = zszoonio_szoo/parm_nio_CFe - zszoonio_szoo*Qfe_zoo  ! but particles remain Fe-rich
       zszoonos_excFe = zszoonos_szoo/parm_nos_CFe - zszoonos_szoo*Qfe_zoo
       zszoonao_excFe = zszoonao_szoo/parm_nao_CFe - zszoonao_szoo*Qfe_zoo  ! extra iron excreted during zoo feeding,

       ! remove the extra Fe from the zooplankton  
# endif /* EXPLICIT_MICROBES */


      !------------------------------------!
      !   (7) Losses of zooplankton        !
      !------------------------------------!

    !-----------------------------------------------------------------------
    !  get fractional factor for routing of zoo losses, based on food supply
    !  more material is routed to large detrital pool when diatoms eaten
    !-----------------------------------------------------------------------
          !--------------------------------------------------------------
          !                        *PJB* 
          !--------------------------------------------------------------
          ! Losses of zooplankton to detritus (i.e., POC) are scaled by 
          ! what they eat. If they eat more diatoms, then more of their
          ! losses (mortality) is routed to POC
          !    
          !   f_zoo_detr_loc = fraction detritus from total zooplankton grazing
          ! 
          ! where...
          !    f_zoo_detr = (/ 0.1, 0.2, 0.1 /)
          !    epsC = 1.00e-8,  ! small C concentration (mmol C/m^3)
          !    epsTinv = 3.17e-8,  ! small inverse time scale (1/year) (1/sec)
          ! 
          !--------------------------------------------------------------

       WORK1 = c0
       WORK2 = c0
       do auto_ind = 1, autotroph_cnt
          WORK1 = WORK1 + f_zoo_detr(auto_ind) * (auto_graze(istr:iend,auto_ind) + epsC * epsTinv)
          WORK2 = WORK2 + (auto_graze(istr:iend,auto_ind) + epsC * epsTinv)
       end do
       ! [ 0.1 when all grazing is small phytos or diazotrophs --> 0.2 when all grazing is diatoms]
       f_zoo_detr_loc = WORK1 / WORK2

          !--------------------------------------------------------------
          !                        *PJB* 
          !--------------------------------------------------------------
          ! Now calculate losses of zooplankton via quadratic closure and
          ! how these losses are partitioned to DOC, DIC, and POC
          ! 
          ! where...
          !    f_loss_thres = 1 above 100 metres (thres_z1)
          !    f_loss_thres = 0 below 150 metres (thres_z2)
          !    f_loss_thres > 0 but < 1 between 100 and 150 metres (thres_z1 < z < thres_z2)
          !    loss_thres_zoo = 0.06,  ! concentration where losses go to zero
          !    parm_z_mort_0  = 0.1 * dps,  ! 0.1 per day  |  Linear mortality
          !    parm_z_mort2_0 = 0.4 * dps,  ! 0.4 per day  |  Quadratic mortality
          !    parm_labile_ratio = 0.85     ! 85% of non-detrital zoo losses go to DIC
          ! 
          !--------------------------------------------------------------

       !!!!!!!!!!!!!!!!!!!!!
       ! Large zooplankton !
       !!!!!!!!!!!!!!!!!!!!!

       ! Calculate the concentration beneath which losses do not occur 
       !   [0 at depth --> 0.06 in upp 100 metres]
       C_loss_thres = f_loss_thres * loss_thres_zoo
       Zprime = max(zooC_loc - C_loss_thres, c0)

       ! Calculate loss of zooplankton (both quadratic and linear mortality)
       zoo_loss = (parm_z_mort2_0 * Zprime**1.5 + parm_z_mort_0 * Zprime) * Tfunc

       ! Route zooplankton losses to DOC, DIC and POC
       zoo_loss_doc = (c1 - parm_labile_ratio) * (c1 - f_zoo_detr_loc) * zoo_loss
       zoo_loss_dic = parm_labile_ratio * (c1 - f_zoo_detr_loc) * zoo_loss
       zoo_loss_poc = f_zoo_detr_loc * zoo_loss 

# ifdef EXPLICIT_MICROBES
       !!!!!!!!!!!!!!!!!!!!!
       ! Small zooplankton !
       !!!!!!!!!!!!!!!!!!!!!

       ! Calculate the concentration beneath which losses do not occur 
       !   [0 at depth --> 0.01 in upp 100 metres]
       C_loss_thres = f_loss_thres * parm_szoo_losst
       SZOOprime = max(szooC_loc - C_loss_thres, c0)

       ! Calculate loss of zooplankton (both quadratic and linear mortality)
       szoo_loss = (parm_szoo_qmort*(SZOOprime*szooC_loc) + parm_szoo_lmort*SZOOprime) * Tfunc

       ! Route zooplankton losses to DOC, DIC and POC
         ! parm_szoo_fdetr = 0.15
         ! parm_szoo_labil = 0.85
       !   Losses to DOC = 0.85 * 0.15 = 0.1275 
       !   Losses to DIC = 0.85 * 0.85 = 0.7225
       !   Losses to POC =               0.15  
       szoo_loss_doc = szoo_loss * (c1 - parm_szoo_fdetr) * (c1 - parm_szoo_labil)
       szoo_loss_dic = szoo_loss * (c1 - parm_szoo_fdetr) * parm_szoo_labil
       szoo_loss_poc = szoo_loss * parm_szoo_fdetr
# endif /* EXPLICIT_MICROBES */


# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(grazesp_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(grazesp_idx_t)) = auto_graze(istr:iend,sp_ind)
       if (bec2_diag_3d_l(grazediat_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(grazediat_idx_t)) = auto_graze(istr:iend,diat_ind)
       if (bec2_diag_3d_l(grazediaz_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(grazediaz_idx_t)) = auto_graze(istr:iend,diaz_ind)
       if (bec2_diag_3d_l(sploss_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(sploss_idx_t)) = auto_loss(istr:iend,sp_ind)
       if (bec2_diag_3d_l(diatloss_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatloss_idx_t)) = auto_loss(istr:iend,diat_ind)
       if (bec2_diag_3d_l(diazloss_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diazloss_idx_t)) = auto_loss(istr:iend,diaz_ind)
       if (bec2_diag_3d_l(spagg_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(spagg_idx_t)) = auto_agg(istr:iend,sp_ind)
       if (bec2_diag_3d_l(diatagg_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatagg_idx_t)) = auto_agg(istr:iend,diat_ind)
       if (bec2_diag_3d_l(diazagg_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diazagg_idx_t)) = auto_agg(istr:iend,diaz_ind)
       if (bec2_diag_3d_l(photocsp_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(photocsp_idx_t)) = photoC(istr:iend,sp_ind)
       if (bec2_diag_3d_l(photocdiat_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(photocdiat_idx_t)) = photoC(istr:iend,diat_ind)
       if (bec2_diag_3d_l(photocdiaz_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(photocdiaz_idx_t)) = photoC(istr:iend,diaz_ind)
       if (bec2_diag_3d_l(totprod_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(totprod_idx_t)) = sum(photoC,dim=2)
       if (bec2_diag_3d_l(spnlim_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(spnlim_idx_t)) = VNtot(istr:iend,sp_ind)
       if (bec2_diag_3d_l(diatnlim_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatnlim_idx_t)) = VNtot(istr:iend,diat_ind)
       if (bec2_diag_3d_l(zooloss_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(zooloss_idx_t)) = zoo_loss
       if (bec2_diag_3d_l(zoolossdic_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(zoolossdic_idx_t)) = zoo_loss_dic
# ifdef Ncycle_SY
       if (bec2_diag_3d_l(spno2uptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(spno2uptake_idx_t)) = NO2_V(istr:iend,sp_ind)
       if (bec2_diag_3d_l(diatno2uptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatno2uptake_idx_t)) = NO2_V(istr:iend,diat_ind)
       if (bec2_diag_3d_l(diazno2uptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diazno2uptake_idx_t)) = NO2_V(istr:iend,diaz_ind)
# endif
# ifdef EXPLICIT_MICROBES
       if (bec2_diag_3d_l(aoa2zoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aoa2zoo_idx_t)) = zzooaoa(istr:iend)
       if (bec2_diag_3d_l(nob2zoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nob2zoo_idx_t)) = zzoonob(istr:iend)
       if (bec2_diag_3d_l(aox2zoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aox2zoo_idx_t)) = zzooaox(istr:iend)
       if (bec2_diag_3d_l(aoa2zoo_excN_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aoa2zoo_excN_idx_t)) = zzooaoa_excN(istr:iend)
       if (bec2_diag_3d_l(nob2zoo_excN_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nob2zoo_excN_idx_t)) = zzoonob_excN(istr:iend)
       if (bec2_diag_3d_l(aox2zoo_excN_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aox2zoo_excN_idx_t)) = zzooaox_excN(istr:iend)
       if (bec2_diag_3d_l(aoa2zoo_excP_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aoa2zoo_excP_idx_t)) = zzooaoa_excP(istr:iend)
       if (bec2_diag_3d_l(nob2zoo_excP_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nob2zoo_excP_idx_t)) = zzoonob_excP(istr:iend)
       if (bec2_diag_3d_l(aox2zoo_excP_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aox2zoo_excP_idx_t)) = zzooaox_excP(istr:iend)
       if (bec2_diag_3d_l(aoa2zoo_excFe_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aoa2zoo_excFe_idx_t)) = zzooaoa_excFe(istr:iend)
       if (bec2_diag_3d_l(nob2zoo_excFe_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nob2zoo_excFe_idx_t)) = zzoonob_excFe(istr:iend)
       if (bec2_diag_3d_l(aox2zoo_excFe_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aox2zoo_excFe_idx_t)) = zzooaox_excFe(istr:iend)
       if (bec2_diag_3d_l(aer2szoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aer2szoo_idx_t)) = zszooaer(istr:iend)
       if (bec2_diag_3d_l(nar2szoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nar2szoo_idx_t)) = zszoonar(istr:iend)
       if (bec2_diag_3d_l(nai2szoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nai2szoo_idx_t)) = zszoonai(istr:iend)
       if (bec2_diag_3d_l(nir2szoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nir2szoo_idx_t)) = zszoonir(istr:iend)
       if (bec2_diag_3d_l(nio2szoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nio2szoo_idx_t)) = zszoonio(istr:iend)
       if (bec2_diag_3d_l(nos2szoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nos2szoo_idx_t)) = zszoonos(istr:iend)
       if (bec2_diag_3d_l(nao2szoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nao2szoo_idx_t)) = zszoonao(istr:iend)
       if (bec2_diag_3d_l(aer2szoo_excN_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aer2szoo_excN_idx_t)) = zszooaer_excN(istr:iend)
       if (bec2_diag_3d_l(nar2szoo_excN_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nar2szoo_excN_idx_t)) = zszoonar_excN(istr:iend)
       if (bec2_diag_3d_l(nai2szoo_excN_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nai2szoo_excN_idx_t)) = zszoonai_excN(istr:iend)
       if (bec2_diag_3d_l(nir2szoo_excN_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nir2szoo_excN_idx_t)) = zszoonir_excN(istr:iend)
       if (bec2_diag_3d_l(nio2szoo_excN_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nio2szoo_excN_idx_t)) = zszoonio_excN(istr:iend)
       if (bec2_diag_3d_l(nos2szoo_excN_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nos2szoo_excN_idx_t)) = zszoonos_excN(istr:iend)
       if (bec2_diag_3d_l(nao2szoo_excN_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nao2szoo_excN_idx_t)) = zszoonao_excN(istr:iend)
       if (bec2_diag_3d_l(aer2szoo_excP_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aer2szoo_excP_idx_t)) = zszooaer_excP(istr:iend)
       if (bec2_diag_3d_l(nar2szoo_excP_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nar2szoo_excP_idx_t)) = zszoonar_excP(istr:iend)
       if (bec2_diag_3d_l(nai2szoo_excP_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nai2szoo_excP_idx_t)) = zszoonai_excP(istr:iend)
       if (bec2_diag_3d_l(nir2szoo_excP_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nir2szoo_excP_idx_t)) = zszoonir_excP(istr:iend)
       if (bec2_diag_3d_l(nio2szoo_excP_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nio2szoo_excP_idx_t)) = zszoonio_excP(istr:iend)
       if (bec2_diag_3d_l(nos2szoo_excP_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nos2szoo_excP_idx_t)) = zszoonos_excP(istr:iend)
       if (bec2_diag_3d_l(nao2szoo_excP_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nao2szoo_excP_idx_t)) = zszoonao_excP(istr:iend)
       if (bec2_diag_3d_l(aer2szoo_excFe_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aer2szoo_excFe_idx_t)) = zszooaer_excFe(istr:iend)
       if (bec2_diag_3d_l(nar2szoo_excFe_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nar2szoo_excFe_idx_t)) = zszoonar_excFe(istr:iend)
       if (bec2_diag_3d_l(nai2szoo_excFe_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nai2szoo_excFe_idx_t)) = zszoonai_excFe(istr:iend)
       if (bec2_diag_3d_l(nir2szoo_excFe_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nir2szoo_excFe_idx_t)) = zszoonir_excFe(istr:iend)
       if (bec2_diag_3d_l(nio2szoo_excFe_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nio2szoo_excFe_idx_t)) = zszoonio_excFe(istr:iend)
       if (bec2_diag_3d_l(nos2szoo_excFe_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nos2szoo_excFe_idx_t)) = zszoonos_excFe(istr:iend)
       if (bec2_diag_3d_l(nao2szoo_excFe_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nao2szoo_excFe_idx_t)) = zszoonao_excFe(istr:iend)
       if (bec2_diag_3d_l(szooloss_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(szooloss_idx_t)) = szoo_loss
# endif /* EXPLICIT_MICROBES */
       if (bec2_diag_3d_l(spno3uptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(spno3uptake_idx_t)) = NO3_V(istr:iend,sp_ind)
       if (bec2_diag_3d_l(diatno3uptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatno3uptake_idx_t)) = NO3_V(istr:iend,diat_ind)
       if (bec2_diag_3d_l(diazno3uptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diazno3uptake_idx_t)) = NO3_V(istr:iend,diaz_ind)
       if (bec2_diag_3d_l(spnh4uptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(spnh4uptake_idx_t)) = NH4_V(istr:iend,sp_ind)
       if (bec2_diag_3d_l(diatnh4uptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diatnh4uptake_idx_t)) = NH4_V(istr:iend,diat_ind)
       if (bec2_diag_3d_l(diaznh4uptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diaznh4uptake_idx_t)) = NH4_V(istr:iend,diaz_ind)
       if (bec2_diag_3d_l(grazedicsp_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(grazedicsp_idx_t)) = auto_graze_dic(istr:iend,sp_ind)
       if (bec2_diag_3d_l(grazedicdiat_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(grazedicdiat_idx_t)) = auto_graze_dic(istr:iend,diat_ind)
       if (bec2_diag_3d_l(grazedicdiaz_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(grazedicdiaz_idx_t)) = auto_graze_dic(istr:iend,diaz_ind)
       if (bec2_diag_3d_l(lossdicsp_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(lossdicsp_idx_t)) = auto_loss_dic(istr:iend,sp_ind)
       if (bec2_diag_3d_l(lossdicdiat_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(lossdicdiat_idx_t)) = auto_loss_dic(istr:iend,diat_ind)
       if (bec2_diag_3d_l(lossdicdiaz_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(lossdicdiaz_idx_t)) = auto_loss_dic(istr:iend,diaz_ind)
       if (bec2_diag_3d_l(grazespzoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(grazespzoo_idx_t)) = auto_graze_zoo(istr:iend,sp_ind)
       if (bec2_diag_3d_l(grazediatzoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(grazediatzoo_idx_t)) = auto_graze_zoo(istr:iend,diat_ind)
       if (bec2_diag_3d_l(grazediazzoo_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(grazediazzoo_idx_t)) = auto_graze_zoo(istr:iend,diaz_ind)
#  else
       bec2_diag_3d(istr:iend,j,k,grazesp_idx_t) = auto_graze(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,grazediat_idx_t) = auto_graze(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,grazediaz_idx_t) = auto_graze(istr:iend,diaz_ind)
       bec2_diag_3d(istr:iend,j,k,sploss_idx_t) = auto_loss(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,diatloss_idx_t) = auto_loss(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,diazloss_idx_t) = auto_loss(istr:iend,diaz_ind)
       bec2_diag_3d(istr:iend,j,k,spagg_idx_t) = auto_agg(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,diatagg_idx_t) = auto_agg(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,diazagg_idx_t) = auto_agg(istr:iend,diaz_ind)
       bec2_diag_3d(istr:iend,j,k,photocsp_idx_t) = photoC(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,photocdiat_idx_t) = photoC(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,photocdiaz_idx_t) = photoC(istr:iend,diaz_ind)
       bec2_diag_3d(istr:iend,j,k,totprod_idx_t) = sum(photoC,dim=2)
       bec2_diag_3d(istr:iend,j,k,spnlim_idx_t) = VNtot(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,diatnlim_idx_t) = VNtot(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,spplim_idx_t) = VPtot(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,diatplim_idx_t) = VPtot(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,diazplim_idx_t) = VPtot(istr:iend,diaz_ind)
       bec2_diag_3d(istr:iend,j,k,zooloss_idx_t) = zoo_loss
       bec2_diag_3d(istr:iend,j,k,zoolossdic_idx_t) = zoo_loss_dic
# ifdef Ncycle_SY
       bec2_diag_3d(istr:iend,j,k,spno2uptake_idx_t) = NO2_V(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,diatno2uptake_idx_t) = NO2_V(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,diazno2uptake_idx_t) = NO2_V(istr:iend,diaz_ind)
# endif
# ifdef EXPLICIT_MICROBES
       bec2_diag_3d(istr:iend,j,k,aoa2zoo_idx_t) = zzooaoa(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nob2zoo_idx_t) = zzoonob(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aox2zoo_idx_t) = zzooaox(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aoa2zoo_excN_idx_t) = zzooaoa_excN(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nob2zoo_excN_idx_t) = zzoonob_excN(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aox2zoo_excN_idx_t) = zzooaox_excN(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aoa2zoo_excP_idx_t) = zzooaoa_excP(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nob2zoo_excP_idx_t) = zzoonob_excP(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aox2zoo_excP_idx_t) = zzooaox_excP(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aoa2zoo_excFe_idx_t) = zzooaoa_excFe(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nob2zoo_excFe_idx_t) = zzoonob_excFe(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aox2zoo_excFe_idx_t) = zzooaox_excFe(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aer2szoo_idx_t) = zszooaer(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nar2szoo_idx_t) = zszoonar(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nai2szoo_idx_t) = zszoonai(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nir2szoo_idx_t) = zszoonir(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nio2szoo_idx_t) = zszoonio(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nos2szoo_idx_t) = zszoonos(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nao2szoo_idx_t) = zszoonao(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aer2szoo_excN_idx_t) = zszooaer_excN(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nar2szoo_excN_idx_t) = zszoonar_excN(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nai2szoo_excN_idx_t) = zszoonai_excN(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nir2szoo_excN_idx_t) = zszoonir_excN(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nio2szoo_excN_idx_t) = zszoonio_excN(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nos2szoo_excN_idx_t) = zszoonos_excN(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nao2szoo_excN_idx_t) = zszoonao_excN(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aer2szoo_excP_idx_t) = zszooaer_excP(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nar2szoo_excP_idx_t) = zszoonar_excP(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nai2szoo_excP_idx_t) = zszoonai_excP(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nir2szoo_excP_idx_t) = zszoonir_excP(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nio2szoo_excP_idx_t) = zszoonio_excP(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nos2szoo_excP_idx_t) = zszoonos_excP(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nao2szoo_excP_idx_t) = zszoonao_excP(istr:iend)
       bec2_diag_3d(istr:iend,j,k,aer2szoo_excFe_idx_t) = zszooaer_excFe(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nar2szoo_excFe_idx_t) = zszoonar_excFe(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nai2szoo_excFe_idx_t) = zszoonai_excFe(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nir2szoo_excFe_idx_t) = zszoonir_excFe(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nio2szoo_excFe_idx_t) = zszoonio_excFe(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nos2szoo_excFe_idx_t) = zszoonos_excFe(istr:iend)
       bec2_diag_3d(istr:iend,j,k,nao2szoo_excFe_idx_t) = zszoonao_excFe(istr:iend)
       bec2_diag_3d(istr:iend,j,k,szooloss_idx_t) = szoo_loss
# endif /* EXPLICIT_MICROBES */
       bec2_diag_3d(istr:iend,j,k,spno3uptake_idx_t) = NO3_V(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,diatno3uptake_idx_t) = NO3_V(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,diazno3uptake_idx_t) = NO3_V(istr:iend,diaz_ind)
       bec2_diag_3d(istr:iend,j,k,spnh4uptake_idx_t) = NH4_V(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,diatnh4uptake_idx_t) = NH4_V(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,diaznh4uptake_idx_t) = NH4_V(istr:iend,diaz_ind)
       bec2_diag_3d(istr:iend,j,k,grazedicsp_idx_t) = auto_graze_dic(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,grazedicdiat_idx_t) = auto_graze_dic(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,grazedicdiaz_idx_t) = auto_graze_dic(istr:iend,diaz_ind)
       bec2_diag_3d(istr:iend,j,k,lossdicsp_idx_t) = auto_loss_dic(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,lossdicdiat_idx_t) = auto_loss_dic(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,lossdicdiaz_idx_t) = auto_loss_dic(istr:iend,diaz_ind)
       bec2_diag_3d(istr:iend,j,k,grazespzoo_idx_t) = auto_graze_zoo(istr:iend,sp_ind)
       bec2_diag_3d(istr:iend,j,k,grazediatzoo_idx_t) = auto_graze_zoo(istr:iend,diat_ind)
       bec2_diag_3d(istr:iend,j,k,grazediazzoo_idx_t) = auto_graze_zoo(istr:iend,diaz_ind)
#   ifdef BEC_COCCO
       bec2_diag_3d(istr:iend,j,k,grazecocco_idx_t) = auto_graze(istr:iend,cocco_ind)
       bec2_diag_3d(istr:iend,j,k,coccoloss_idx_t) = auto_loss(istr:iend,cocco_ind)
       bec2_diag_3d(istr:iend,j,k,coccoagg_idx_t) = auto_agg(istr:iend,cocco_ind)
       bec2_diag_3d(istr:iend,j,k,photoccocco_idx_t) = photoC(istr:iend,cocco_ind)
       bec2_diag_3d(istr:iend,j,k,cocconlim_idx_t) = VNtot(istr:iend,cocco_ind)
       bec2_diag_3d(istr:iend,j,k,coccoplim_idx_t) = VPtot(istr:iend,cocco_ind)
       bec2_diag_3d(istr:iend,j,k,coccono3uptake_idx_t) = NO3_V(istr:iend,cocco_ind)
       bec2_diag_3d(istr:iend,j,k,cocconh4uptake_idx_t) = NH4_V(istr:iend,cocco_ind)
       bec2_diag_3d(istr:iend,j,k,coccograzedic_idx_t) = auto_graze_dic(istr:iend,cocco_ind)
       bec2_diag_3d(istr:iend,j,k,coccolossdic_idx_t) = auto_loss_dic(istr:iend,cocco_ind)
       bec2_diag_3d(istr:iend,j,k,grazecoccozoo_idx_t) = auto_graze_zoo(istr:iend,cocco_ind)
#   endif /* BEC_COCCO */
#  endif
# endif /* BEC2_DIAG */


    !------------------------------------------------------------------------------------!
    !   (8) Linear and Quadratic mortality terms for chemoautrophs and heterotrophs      !
    !------------------------------------------------------------------------------------!

# ifdef EXPLICIT_MICROBES
        AOAprime = MAX(AOA_loc - parm_che_bmin, c0)
        NOBprime = MAX(NOB_loc - parm_che_bmin, c0)
        AOXprime = MAX(AOX_loc - parm_che_bmin, c0)
        AERprime = MAX(AER_loc - parm_het_bmin, c0)
        NARprime = MAX(NAR_loc - parm_het_bmin, c0)
        NAIprime = MAX(NAI_loc - parm_het_bmin, c0)
        NIRprime = MAX(NIR_loc - parm_het_bmin, c0)
        NIOprime = MAX(NIO_loc - parm_het_bmin, c0)
        NOSprime = MAX(NOS_loc - parm_het_bmin, c0)
        NAOprime = MAX(NAO_loc - parm_het_bmin, c0)
        zmortaoa = ( AOAprime * parm_aoa_lmort + (AOAprime * AOA_loc) * parm_aoa_qmort )* Tfunc
        zmortnob = ( NOBprime * parm_nob_lmort + (NOBprime * NOB_loc) * parm_nob_qmort )* Tfunc
        zmortaox = ( AOXprime * parm_aox_lmort + (AOXprime * AOX_loc) * parm_aox_qmort )* Tfunc
        zmortaer = ( AERprime * parm_aer_lmort + (AERprime * AER_loc) * parm_aer_qmort )* Tfunc
        zmortnar = ( NARprime * parm_nar_lmort + (NARprime * NAR_loc) * parm_nar_qmort )* Tfunc
        zmortnai = ( NAIprime * parm_nai_lmort + (NAIprime * NAI_loc) * parm_nai_qmort )* Tfunc
        zmortnir = ( NIRprime * parm_nir_lmort + (NIRprime * NIR_loc) * parm_nir_qmort )* Tfunc
        zmortnio = ( NIOprime * parm_nio_lmort + (NIOprime * NIO_loc) * parm_nio_qmort )* Tfunc
        zmortnos = ( NOSprime * parm_nos_lmort + (NOSprime * NOS_loc) * parm_nos_qmort )* Tfunc
        zmortnao = ( NAOprime * parm_nao_lmort + (NAOprime * NAO_loc) * parm_nao_qmort )* Tfunc

!        ! Add extra strong mortality to NOB below 1 µM O2
!        where (O2_loc .lt. c1 )
!            !zmortnob2 = NOB_loc * parm_munob
!            zmortnob = zmortnob + NOB_loc * parm_munob
!        endwhere

        ! Route mortality to different pools, conserving constant C:N:P of implicit POM, excess --> DOM
        zmortaoa_doc = zmortaoa * parm_che_mortdoc * (c1-parm_che_mortpoc)
        zmortaoa_dic = zmortaoa * (c1-parm_che_mortdoc) * (c1-parm_che_mortpoc)
        zmortaoa_poc = zmortaoa * parm_che_mortpoc
        zmortaoa_don = zmortaoa/parm_aoa_CN - (zmortaoa_dic/parm_aoa_CN + zmortaoa_poc*Q)
        zmortaoa_dop = zmortaoa/parm_aoa_CP - (zmortaoa_dic/parm_aoa_CP + zmortaoa_poc*Qp_zoo_pom)
        zmortnob_doc = zmortnob * parm_che_mortdoc * (c1-parm_che_mortpoc)
        zmortnob_dic = zmortnob * (c1-parm_che_mortdoc) * (c1-parm_che_mortpoc)
        zmortnob_poc = zmortnob * parm_che_mortpoc
        zmortnob_don = zmortnob/parm_nob_CN - (zmortnob_dic/parm_nob_CN + zmortnob_poc*Q)
        zmortnob_dop = zmortnob/parm_nob_CP - (zmortnob_dic/parm_nob_CP + zmortnob_poc*Qp_zoo_pom)
        zmortaox_doc = zmortaox * parm_che_mortdoc * (c1-parm_che_mortpoc)
        zmortaox_dic = zmortaox * (c1-parm_che_mortdoc) * (c1-parm_che_mortpoc)
        zmortaox_poc = zmortaox * parm_che_mortpoc
        zmortaox_don = zmortaox/parm_aox_CN - (zmortaox_dic/parm_aox_CN + zmortaox_poc*Q)
        zmortaox_dop = zmortaox/parm_aox_CP - (zmortaox_dic/parm_aox_CP + zmortaox_poc*Qp_zoo_pom)
        zmortaer_doc = zmortaer * parm_het_mortdoc * (c1-parm_het_mortpoc)
        zmortaer_dic = zmortaer * (c1-parm_het_mortdoc) * (c1-parm_het_mortpoc)
        zmortaer_poc = zmortaer * parm_het_mortpoc
        zmortaer_don = zmortaer/parm_aer_CN - (zmortaer_dic/parm_aer_CN + zmortaer_poc*Q)
        zmortaer_dop = zmortaer/parm_aer_CP - (zmortaer_dic/parm_aer_CP + zmortaer_poc*Qp_zoo_pom)
        zmortnar_doc = zmortnar * parm_het_mortdoc * (c1-parm_het_mortpoc)
        zmortnar_dic = zmortnar * (c1-parm_het_mortdoc) * (c1-parm_het_mortpoc)
        zmortnar_poc = zmortnar * parm_het_mortpoc
        zmortnar_don = zmortnar/parm_nar_CN - (zmortnar_dic/parm_nar_CN + zmortnar_poc*Q)
        zmortnar_dop = zmortnar/parm_nar_CP - (zmortnar_dic/parm_nar_CP + zmortnar_poc*Qp_zoo_pom)
        zmortnai_doc = zmortnai * parm_het_mortdoc * (c1-parm_het_mortpoc)
        zmortnai_dic = zmortnai * (c1-parm_het_mortdoc) * (c1-parm_het_mortpoc)
        zmortnai_poc = zmortnai * parm_het_mortpoc
        zmortnai_don = zmortnai/parm_nai_CN - (zmortnai_dic/parm_nai_CN + zmortnai_poc*Q)
        zmortnai_dop = zmortnai/parm_nai_CP - (zmortnai_dic/parm_nai_CP + zmortnai_poc*Qp_zoo_pom)
        zmortnir_doc = zmortnir * parm_het_mortdoc * (c1-parm_het_mortpoc)
        zmortnir_dic = zmortnir * (c1-parm_het_mortdoc) * (c1-parm_het_mortpoc)
        zmortnir_poc = zmortnir * parm_het_mortpoc
        zmortnir_don = zmortnir/parm_nir_CN - (zmortnir_dic/parm_nir_CN + zmortnir_poc*Q)
        zmortnir_dop = zmortnir/parm_nir_CP - (zmortnir_dic/parm_nir_CP + zmortnir_poc*Qp_zoo_pom)
        zmortnio_doc = zmortnio * parm_het_mortdoc * (c1-parm_het_mortpoc)
        zmortnio_dic = zmortnio * (c1-parm_het_mortdoc) * (c1-parm_het_mortpoc)
        zmortnio_poc = zmortnio * parm_het_mortpoc
        zmortnio_don = zmortnio/parm_nio_CN - (zmortnio_dic/parm_nio_CN + zmortnio_poc*Q)
        zmortnio_dop = zmortnio/parm_nio_CP - (zmortnio_dic/parm_nio_CP + zmortnio_poc*Qp_zoo_pom)
        zmortnos_doc = zmortnos * parm_het_mortdoc * (c1-parm_het_mortpoc)
        zmortnos_dic = zmortnos * (c1-parm_het_mortdoc) * (c1-parm_het_mortpoc)
        zmortnos_poc = zmortnos * parm_het_mortpoc
        zmortnos_don = zmortnos/parm_nos_CN - (zmortnos_dic/parm_nos_CN + zmortnos_poc*Q)
        zmortnos_dop = zmortnos/parm_nos_CP - (zmortnos_dic/parm_nos_CP + zmortnos_poc*Qp_zoo_pom)
        zmortnao_doc = zmortnao * parm_het_mortdoc * (c1-parm_het_mortpoc)
        zmortnao_dic = zmortnao * (c1-parm_het_mortdoc) * (c1-parm_het_mortpoc)
        zmortnao_poc = zmortnao * parm_het_mortpoc
        zmortnao_don = zmortnao/parm_nao_CN - (zmortnao_dic/parm_nao_CN + zmortnao_poc*Q)
        zmortnao_dop = zmortnao/parm_nao_CP - (zmortnao_dic/parm_nao_CP + zmortnao_poc*Qp_zoo_pom)
        ! NOTE:
        !  You may notice that the routing to DOFe is not calculated here, unlike DON and DOP.
        !  This is on purpose. The C:Fe ratio of sinking particulate matter is not conserved like 
        !  for C:N:P. This means that the Fe-rich character of the microbes can be routed to POM.
        !  Also, all biomass that goes to DOM should go to dFe for Fe.
# endif /* EXPLICIT_MICROBES */



    !------------------------------------!
    !   (9)  Compute terms for DOM       !
    !------------------------------------!

       ! DOC !
       DOC_prod = zoo_loss_doc + sum(auto_loss_doc, dim=2) + sum(auto_graze_doc, dim=2)
# ifdef EXPLICIT_MICROBES
     &            + zzooaoa_doc + zzoonob_doc + zzooaox_doc
     &            + zszooaer_doc + zszoonar_doc + zszoonai_doc + zszoonir_doc 
     &            + zszoonio_doc + zszoonos_doc + zszoonao_doc
     &            + zmortaoa_doc + zmortnob_doc + zmortaox_doc
     &            + zmortaer_doc + zmortnar_doc + zmortnai_doc + zmortnir_doc 
     &            + zmortnio_doc + zmortnos_doc + zmortnao_doc
     &            + POC_remin(istr:iend,j)
     &            + szoo_loss_doc
# endif /* EXPLICIT_MICROBES */

       ! DON !
       DON_prod = Q * (zoo_loss_doc + sum(auto_loss_doc, dim=2) + sum(auto_graze_doc, dim=2))
# ifdef EXPLICIT_MICROBES
     &            + zzooaoa_doc/parm_aoa_CN + zmortaoa_don
     &            + zzoonob_doc/parm_nob_CN + zmortnob_don
     &            + zzooaox_doc/parm_aox_CN + zmortaox_don
     &            + zszooaer_doc/parm_aer_CN + zmortaer_don
     &            + zszoonar_doc/parm_nar_CN + zmortnar_don
     &            + zszoonai_doc/parm_nai_CN + zmortnai_don
     &            + zszoonir_doc/parm_nir_CN + zmortnir_don
     &            + zszoonio_doc/parm_nio_CN + zmortnio_don
     &            + zszoonos_doc/parm_nos_CN + zmortnos_don
     &            + zszoonao_doc/parm_nao_CN + zmortnao_don
     &            + POC_remin(istr:iend,j) * Q
     &            + szoo_loss_doc * Q
# endif /* EXPLICIT_MICROBES */
  
       ! DOP !
       DOP_prod = Qp_zoo_pom * zoo_loss_doc
# ifdef EXPLICIT_MICROBES
     &            + zzooaoa_doc/parm_aoa_CP + zmortaoa_dop
     &            + zzoonob_doc/parm_nob_CP + zmortnob_dop
     &            + zzooaox_doc/parm_aox_CP + zmortaox_dop
     &            + zszooaer_doc/parm_aer_CP + zmortaer_dop
     &            + zszoonar_doc/parm_nar_CP + zmortnar_dop
     &            + zszoonai_doc/parm_nai_CP + zmortnai_dop
     &            + zszoonir_doc/parm_nir_CP + zmortnir_dop
     &            + zszoonio_doc/parm_nio_CP + zmortnio_dop
     &            + zszoonos_doc/parm_nos_CP + zmortnos_dop
     &            + zszoonao_doc/parm_nao_CP + zmortnao_dop
     &            + POC_remin(istr:iend,j) * Qp_zoo_pom
     &            + szoo_loss_doc * Qp_zoo_pom
# endif /* EXPLICIT_MICROBES */
       do auto_ind = 1, autotroph_cnt
          if (Qp(auto_ind) == Qp_zoo_pom) then
             DOP_prod = DOP_prod + Qp(auto_ind) * (auto_loss_doc(istr:iend,auto_ind)
     &                  + auto_graze_doc(istr:iend,auto_ind))
          else
             DOP_prod = DOP_prod + remaining_P_dop(istr:iend,auto_ind)
          endif
       end do

       ! DOFe !
       ! PJB - really, all Fe from biomass should be going to the dFe pool... There is
       !       no way to differentiate DOFe and Fe' when measuring dFe (< 0.2 µm filter) 
       ! DOFe_prod is therefore routed to dFe, not to DOFe, but we collect it here anyway
       DOFe_prod = Qfe_zoo * zoo_loss_doc
# ifdef EXPLICIT_MICROBES
     &            + (zzooaoa_doc + zmortaoa_doc) / parm_aoa_CFe
     &            + (zzoonob_doc + zmortnob_doc) / parm_nob_CFe
     &            + (zzooaox_doc + zmortaox_doc) / parm_aox_CFe
     &            + (zszooaer_doc + zmortaer_doc) / parm_aer_CFe
     &            + (zszoonar_doc + zmortnar_doc) / parm_nar_CFe
     &            + (zszoonai_doc + zmortnai_doc) / parm_nai_CFe
     &            + (zszoonir_doc + zmortnir_doc) / parm_nir_CFe
     &            + (zszoonio_doc + zmortnio_doc) / parm_nio_CFe
     &            + (zszoonos_doc + zmortnos_doc) / parm_nos_CFe
     &            + (zszoonao_doc + zmortnao_doc) / parm_nao_CFe
     &            + szoo_loss_doc * Qfe_zoo
# endif /* EXPLICIT_MICROBES */
       ! Add phytoplankton mortality --> DOFe
       do auto_ind = 1, autotroph_cnt
          DOFe_prod = DOFe_prod + Qfe(istr:iend,auto_ind) * 
     &    (auto_loss_doc(istr:iend,auto_ind) + auto_graze_doc(istr:iend,auto_ind))
       end do


# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(docprod_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(docprod_idx_t)) = DOC_prod
       if (bec2_diag_3d_l(donprod_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(donprod_idx_t)) = DON_prod
       if (bec2_diag_3d_l(dopprod_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(dopprod_idx_t)) = DOP_prod
       if (bec2_diag_3d_l(dofeprod_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(dofeprod_idx_t)) = DOFe_prod
#  else
       bec2_diag_3d(istr:iend,j,k,docprod_idx_t) = DOC_prod
       bec2_diag_3d(istr:iend,j,k,donprod_idx_t) = DON_prod
       bec2_diag_3d(istr:iend,j,k,dopprod_idx_t) = DOP_prod
       bec2_diag_3d(istr:iend,j,k,dofeprod_idx_t) = DOFe_prod
#  endif
# endif


# ifndef EXPLICIT_MICROBES

       DOC_remin  = DOC_loc  * DOC_reminR * 0.0685
       DON_remin  = DON_loc  * DON_reminR * 0.1
       DOFe_remin = DOFe_loc * DOFe_reminR * 0.05
       DOP_remin(istr:iend,j)  = DOP_loc  * DOP_reminR * 0.05

    !-----------------------------------------------------------------------
    !  Refractory remin rate due to photochemistry
    !  below euphotic zone remin rate sharply decrease
    !-----------------------------------------------------------------------

       where (PAR_lay > 1.0)
          DONr_remin = DONr_loc * DONr_reminR
          DOPr_remin(istr:iend,j) = DOPr_loc * DOPr_reminR
       elsewhere
          DONr_remin = DONr_loc * (c1/(365.0*670.0)) * dps ! 1/670 yrs
          DOPr_remin(istr:iend,j) = DOPr_loc * (c1/(365.0*460.0)) * dps ! 1/460 yrs
       end where

# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(docremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(docremin_idx_t)) = DOC_remin
       if (bec2_diag_3d_l(donremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(donremin_idx_t)) = DON_remin
       if (bec2_diag_3d_l(doferemin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(doferemin_idx_t)) = DOFe_remin
       if (bec2_diag_3d_l(dopremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(dopremin_idx_t)) = DOP_remin(istr:iend,j)
       if (bec2_diag_3d_l(donrremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(donrremin_idx_t)) = DONr_remin
#  else
       bec2_diag_3d(istr:iend,j,k,docremin_idx_t) = DOC_remin
       bec2_diag_3d(istr:iend,j,k,donremin_idx_t) = DON_remin
       bec2_diag_3d(istr:iend,j,k,doferemin_idx_t) = DOFe_remin
       bec2_diag_3d(istr:iend,j,k,dopremin_idx_t) = DOP_remin(istr:iend,j)
       bec2_diag_3d(istr:iend,j,k,donrremin_idx_t) = DONr_remin
#  endif
# endif

# endif /* not EXPLICIT_MICROBES */


    !---------------------------! 
    !   (10)  terms for POM     !
    !---------------------------!

       ! POC !
       POC_prod = zoo_loss_poc + sum(auto_graze_poc, dim=2)
     &            + sum(auto_agg, dim=2) + sum(auto_loss_poc, dim=2)
# ifdef EXPLICIT_MICROBES
     &            + zmortaoa_poc + zmortnob_poc + zmortaox_poc
     &            + zmortaer_poc + zmortnar_poc + zmortnai_poc + zmortnir_poc 
     &            + zmortnio_poc + zmortnos_poc + zmortnao_poc
     &            + zzooaoa_poc + zzoonob_poc + zzooaox_poc
     &            + zszooaer_poc + zszoonar_poc + zszoonai_poc + zszoonir_poc 
     &            + zszoonio_poc + zszoonos_poc + zszoonao_poc
     &            + szoo_loss_poc
# endif /* EXPLICIT_MICROBES */

       ! CaCO3   ( 33% of CaCO3 is remin when phyto are grazed )
       P_CaCO3_prod = ((c1 - f_graze_CaCO3_REMIN) * auto_graze(istr:iend,sp_ind) +
     &                auto_loss(istr:iend,sp_ind) + auto_agg(istr:iend,sp_ind)) * QCaCO3(istr:iend,sp_ind)
# ifdef BEC_COCCO
     &                + ((c1 - f_graze_CaCO3_REMIN) * auto_graze(istr:iend,cocco_ind) +
     &                auto_loss(istr:iend,cocco_ind) + auto_agg(istr:iend,cocco_ind)) * QCaCO3(istr:iend,cocco_ind)
# endif

       ! Biogenic Silicate   ( 60% of grazed diatom SiO2 is remineralised )
       P_SiO2_prod = c0
       do auto_ind = 1, autotroph_cnt
          if (kSiO3(auto_ind) > c0) then
             P_SiO2_prod = P_SiO2_prod + Qsi(istr:iend,auto_ind)
     &          * ((c1 - f_graze_si_remin) * auto_graze(istr:iend,auto_ind) + auto_agg(istr:iend,auto_ind)
     &          + loss_poc(auto_ind) * auto_loss(istr:iend,auto_ind))
          endif
       end do



# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(pocprod_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(pocprod_idx_t)) = POC_prod(istr:iend)
       if (bec2_diag_3d_l(sio2prod_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(sio2prod_idx_t)) = P_SiO2_prod(istr:iend)
       if (bec2_diag_3d_l(pcaco3prod_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(pcaco3prod_idx_t)) = P_CaCO3_prod(istr:iend)
#  else
       bec2_diag_3d(istr:iend,j,k,pocprod_idx_t) = POC_prod(istr:iend)
       bec2_diag_3d(istr:iend,j,k,sio2prod_idx_t) = P_SiO2_prod(istr:iend)
       bec2_diag_3d(istr:iend,j,k,pcaco3prod_idx_t) = P_CaCO3_prod(istr:iend)
#  endif
# endif

    !-----------------------------------------------------------------------!
    !  (11) Iron scavenging 
    !       a) compute in terms of loss per year per unit iron (/year)
    !       b) scale by sinking POMx10 + Dust + bSi + CaCO3 flux
    !       c) increase scavenging rate at higher iron (>0.6nM)
    !       d) convert to net loss per second    
    !-----------------------------------------------------------------------!

       Fe_scavenge_rate = parm_Fe_scavenge_rate0

       Fe_scavenge_rate = Fe_scavenge_rate *
     &     ((POC_sflux_out(istr:iend,j) + POC_hflux_out(istr:iend,j)) * POC_mass * 10.0 +
     &      (P_CaCO3_sflux_out(istr:iend,j) + P_CaCO3_hflux_out(istr:iend,j)) * P_CaCO3_mass +
     &      (P_SiO2_sflux_out(istr:iend,j) + P_SiO2_hflux_out(istr:iend,j)) * P_SiO2_mass +
     &      (dust_sflux_out(istr:iend,j) + dust_hflux_out(istr:iend,j)) * dust_fescav_scale)

       where (Fe_loc > Fe_scavenge_thres1)
          Fe_scavenge_rate = Fe_scavenge_rate + (Fe_loc - Fe_scavenge_thres1) * fe_max_scale2
       end where

       Fe_scavenge = yps * Fe_loc * Fe_scavenge_rate

       P_iron_prod = zoo_loss_poc * Qfe_zoo + Fe_scavenge
# ifdef EXPLICIT_MICROBES
     &               + zmortaoa_poc/parm_aoa_CFe + zmortnob_poc/parm_nob_CFe + zmortaox_poc/parm_aox_CFe
     &               + zmortaer_poc/parm_aer_CFe 
     &               + zmortnar_poc/parm_nar_CFe + zmortnai_poc/parm_nai_CFe + zmortnir_poc/parm_nir_CFe 
     &               + zmortnio_poc/parm_nio_CFe + zmortnos_poc/parm_nos_CFe + zmortnao_poc/parm_nao_CFe
     &               + zzooaoa_poc/parm_aoa_CFe + zzoonob_poc/parm_nob_CFe + zzooaox_poc/parm_aox_CFe
     &               + zszooaer_poc/parm_aer_CFe 
     &               + zszoonar_poc/parm_nar_CFe + zszoonai_poc/parm_nai_CFe + zszoonir_poc/parm_nir_CFe 
     &               + zszoonio_poc/parm_nio_CFe + zszoonos_poc/parm_nos_CFe + zszoonao_poc/parm_nao_CFe
     &               + szoo_loss_poc * Qfe_zoo
# endif /* EXPLICIT_MICROBES */
       ! Add phytoplankton mortality and grazing terms
       do auto_ind = 1, autotroph_cnt
          P_iron_prod = P_iron_prod + Qfe(istr:iend,auto_ind) 
     &       * (auto_agg(istr:iend,auto_ind) + auto_graze_poc(istr:iend,auto_ind) + auto_loss_poc(istr:iend,auto_ind))
       end do


# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(pironprod_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(pironprod_idx_t)) = P_iron_prod
       if (bec2_diag_3d_l(fescavenge_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(fescavenge_idx_t)) = Fe_scavenge
       if (bec2_diag_3d_l(fescavengerate_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(fescavengerate_idx_t)) = Fe_scavenge_rate
#  else
       bec2_diag_3d(istr:iend,j,k,pironprod_idx_t) = P_iron_prod
       bec2_diag_3d(istr:iend,j,k,fescavenge_idx_t) = Fe_scavenge
       bec2_diag_3d(istr:iend,j,k,fescavengerate_idx_t) = Fe_scavenge_rate
#  endif
# endif


    !-----------------------------------------------------! 
    !   (12)  sinking of POC and its remineralisation     !
    !-----------------------------------------------------!

# ifdef IODINE
       call compute_particulate_terms(k,QA_dust_def, TEMP, O2_loc, NO3_loc,
     &        POC_prod,P_iron_prod,P_CaCO3_prod,P_SiO2_prod,SED_DENITRIF,
     &        OTHER_REMIN,iod_sed_flux1,iod_sed_flux2,iod_sed_flux3,
     &        iod_sed_flux4, istr,iend,j)
# else
       call compute_particulate_terms(k,QA_dust_def, TEMP, O2_loc, NO3_loc,POC_prod,
     &        P_iron_prod,P_CaCO3_prod,P_SiO2_prod,SED_DENITRIF, OTHER_REMIN,istr,iend,j)
# endif


    !----------------------------------------------------!
    !   (13) Nitrification and anaerobic metabolisms     !
    !----------------------------------------------------!

          ! ---------------------------------------------!
          !                        *PJB* 
          ! ---------------------------------------------!
          !  
          !
          ! ---------------------------------------------!


# ifdef Ncycle_SY

        ! initialise arrays
        ammox     = c0
        nitrox    = c0
        anammox   = c0
        DENITRIF1 = c0
        DENITRIF2 = c0
        DENITRIF3 = c0

        ! ------------------------------------------!
        !     Implicit nitrification and anammox
        ! ------------------------------------------!
# ifndef EXPLICIT_MICROBES

        ! ammonium oxidation 
        where (NH4_loc .gt. c0 .and. O2_loc .gt. c0 .and. PAR_out < parm_nitrif_par_lim)
            ammox = parm_kao * (O2_loc / (parm_ko2_ao + O2_loc)) 
     &      * (NH4_loc/(parm_knh4_ao + NH4_loc))
                where (PAR_in > parm_nitrif_par_lim)
                        ammox = ammox * log(PAR_out / parm_nitrif_par_lim) / (-KPARdz)
                end where
        elsewhere
            ammox = c0
        end where    
        ! find N2O yield associated with ammonia oxidation
        where (O2_loc .gt. 1e-5)
                n2o_nh4_yield1 = c1/(2.0 + 2.0 / ((parm_n2o_ji_a / O2_loc + parm_n2o_ji_b)/100.0));
                no2_nh4_yield1 = c1/(1.0 * ((parm_n2o_ji_a / O2_loc +parm_n2o_ji_b)/100.0) + 1.0);
        elsewhere
                n2o_nh4_yield1 = 0.5
                no2_nh4_yield1 = c0
        end where

        ! nitrite oxidation
        where (NO2_loc .gt. c0 .and. O2_loc .gt. c0 .and. PAR_out < parm_nitrif_par_lim)
               nitrox = parm_kno * (O2_loc / (parm_ko2_no + O2_loc))
     &         * (NO2_loc/(parm_kno2_no + NO2_loc)) 
               where (PAR_in > parm_nitrif_par_lim)
                        nitrox = nitrox * log(PAR_out / parm_nitrif_par_lim) / (-KPARdz)
               end where
        elsewhere
            nitrox = c0
        end where

        ! anammox
        where (NH4_loc.gt.c0 .and. NO2_loc.gt.c0)
             anammox = parm_kax * (NH4_loc / (parm_knh4_ax + NH4_loc)) *
     &                (NO2_loc / (parm_kno2_ax + NO2_loc)) * exp(-max(O2_loc,c0)/parm_ko2_ax)
        elsewhere
             anammox = c0
        end where

          ! --------------------------------!
          !     Implicit denitrification 
          ! --------------------------------!

       ! split denitrification into three components i.e., one done by no3, one by no2 and one by n2o
       ! Total organic matter denitrified in C units 
       R_oxic = parm_koxic * (O2_loc / (parm_ko2_oxic + O2_loc)) + epsN
       R_den1 = parm_kden1 * (NO3_loc / (parm_kno3_den1 + NO3_loc)) * exp(-max(O2_loc,c0)/parm_ko2_den1) + epsN
       R_den2 = parm_kden2 * (NO2_loc / (parm_kno2_den2 + NO2_loc)) * exp(-max(O2_loc,c0)/parm_ko2_den2) + epsN
       R_den3 = parm_kden3 * (N2O_loc / (parm_kn2o_den3 + N2O_loc)) * exp(-max(O2_loc,c0)/parm_ko2_den3) + epsN
       roxic = R_oxic / (R_den1 + R_den2 + R_den3 + R_oxic)
       rden1 = R_den1 / (R_den1 + R_den2 + R_den3 + R_oxic)
       rden2 = R_den2 / (R_den1 + R_den2 + R_den3 + R_oxic)
       rden3 = R_den3 / (R_den1 + R_den2 + R_den3 + R_oxic)

       Rden =  (DOC_remin + POC_remin(istr:iend,j)
     &                - OTHER_REMIN - SED_DENITRIF * denitrif_C_N)

       ! DENITRIF1 NO3 --> NO2
       where (NO3_loc + (nitrox - rden1 * Rden * denitrif_NO3_C -
     &         SED_DENITRIF - sum(NO3_V, dim=2)) * dt .ge. c0)
                DENITRIF1 = rden1 * Rden * denitrif_NO3_C
       elsewhere (NO3_loc .le. c0)
                DENITRIF1 =  c0
       elsewhere
                DENITRIF1 = NO3_loc / dt + nitrox - SED_DENITRIF -
     &          sum(NO3_V, dim=2)
       end where

       ! DENITRIF2 NO2 -->0.5  N2O
       where (NO2_loc  +  (ammox * no2_nh4_yield1 - nitrox +
     &   DENITRIF1 - rden2 * Rden * denitrif_NO2_C - anammox -
     &   sum(NO2_V, dim=2)) * dt .ge. c0)
            DENITRIF2 = rden2 * Rden * denitrif_NO2_C
       elsewhere (NO2_loc .le. c0)
            DENITRIF2 = c0
       elsewhere
            DENITRIF2 = NO2_loc / dt + ammox * no2_nh4_yield1 - nitrox +
     &      DENITRIF1
       end where

       ! DENITRIF3 N2O --> N2
       where (N2O_loc + (ammox * n2o_nh4_yield1 + 0.5 * DENITRIF2
     &   - rden3 * Rden * denitrif_N2O_C) * dt  .ge. c0)
            DENITRIF3 = rden3 * Rden * denitrif_N2O_C
        elsewhere (N2O_loc .le. c0)
            DENITRIF3 = c0
        elsewhere
            DENITRIF3 = N2O_loc / dt + ammox * n2o_nh4_yield1 + 0.5 * DENITRIF2
        end where
# endif /* not EXPLICIT_MICROBES */


# ifdef EXPLICIT_MICROBES
          ! -----------------------------------------------------------!
          !     Explicit nitrification, anammox and denitrification 
          ! -----------------------------------------------------------!

        ! initialise arrays (chemoautotrophs)
        zmuaoa  = c0
        zbioaoa = c0
        zmunob  = c0
        zmunob_ana = c0
        zbionob = c0
        zbionob_ana = c0
        zmuaox  = c0
        zbioaox = c0
        zmufnob    = c0

        ! determine N2O yield from ammonia oxidation using empirical relationship with O2 (Frey et al. 2023)
        yaoa_n2o = min(3.0, (0.2/(O2_loc + epsN) + 0.5))/100.0

        ! apply light limitation function to ammonia oxidation 
        !aoa_parlim = c1 - (PAR_out/PAR_in) / (PAR_out/PAR_in + parm_kaoa_par)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ammonia oxidising archaea metabolism !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        zuptnh4 = c0
        zuptoxy = c0
        where (NH4_loc.gt.c0 .and. O2_loc.gt.c0)
           zuptoxy = parm_aoa_po2 * O2_loc     !!! diffusive limit to O2 uptake
           !zuptoxy = parm_muaoa * parm_yaoa_oxy * ( O2_loc / ( O2_loc + parm_kaoa_oxy ) )
           zuptnh4 = parm_muaoa * parm_yaoa_nh4 * ( NH4_loc / ( NH4_loc + parm_kaoa_nh4 ) )
        endwhere
        zmuaoa = max(c0, min( (zuptoxy/parm_yaoa_oxy), (zuptnh4/parm_yaoa_nh4) ) ) * Tfunc !* aoa_parlim
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! nitrite oxidising bacteria metabolism !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        zuptno2 = c0
        zuptoxy = c0
        where (NH4_loc.gt.c0 .and. NO2_loc.gt.c0 .and. O2_loc.gt.c0)
            zuptoxy = parm_nob_po2 * O2_loc     !!! diffusive limit to O2 uptake
            !zuptoxy = (parm_munob * parm_ynob_oxy) * ( O2_loc / ( O2_loc + parm_knob_oxy ) )
            zuptno2 = (parm_munob * parm_ynob_no2) * ( NO2_loc / ( NO2_loc + parm_knob_no2 ) )
        endwhere
        zmunob = max(c0, min( (zuptoxy/parm_ynob_oxy), (zuptno2/parm_ynob_no2) ) ) * Tfunc
        where ((zuptoxy/parm_ynob_oxy) .gt. (zuptno2/parm_ynob_no2)) zmufnob = c1
        ! try simulating some anaerobic nitrite oxidation (only consumes NO2)
        where (NH4_loc.gt.c0 .and. NO2_loc.gt.c0 .and. O2_loc.eq.c0)
            zuptno2 = (parm_munob * parm_ynob_no2 * 0.05) * ( NO2_loc / ( NO2_loc + parm_knob_no2 ) )
            zmunob_ana = (zuptno2/parm_ynob_no2) * Tfunc * c0
        endwhere
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! anaerobic ammonium oxidising archaea metabolism !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        zuptnh4 = c0
        zuptno2 = c0
        where (NH4_loc.gt.c0 .and. NO2_loc.gt.c0)
            zuptnh4 = ( parm_muaox * parm_yaox_nh4) * ( NH4_loc / ( NH4_loc + parm_kaox_nh4 ) )
            zuptno2 = ( parm_muaox * parm_yaox_no2) * ( NO2_loc / ( NO2_loc + parm_kaox_no2 ) )
        endwhere
        zmuaox = max(c0, min( (zuptnh4/parm_yaox_nh4), (zuptno2/parm_yaox_no2) ) ) * Tfunc

        ! CHEMOAUTOS : growth rate, accumulation of C biomass and rates
        zbioaoa = zmuaoa * AOA_loc
        zbionob = zmunob * NOB_loc
        zbionob_ana = zmunob_ana * NOB_loc
        zbioaox = zmuaox * AOX_loc
        ammox   = zbioaoa * parm_yaoa_nh4
        nitrox  = zbionob * parm_ynob_no2 
        nitrox_ana  = zbionob_ana * parm_ynob_no2 
        anammox = zbioaox * parm_yaox_nh4   ! measured by NH4 use


        ! initialise arrays (heterotrophs)
        zmuaer = c0
        zbioaer = c0
        zmunar_aer = c0
        zmunar_ana = c0
        zmunar_fac = c0
        zmunar  = c0
        zbionar = c0
        zmunai_aer = c0
        zmunai_ana = c0
        zmunai_fac = c0
        zmunai  = c0
        zbionai = c0
        zmunir_aer = c0
        zmunir_ana = c0
        zmunir_fac = c0
        zmunir  = c0
        zbionir = c0
        zmunio_aer = c0
        zmunio_ana = c0
        zmunio_fac = c0
        zmunio  = c0
        zbionio = c0
        zmunos_aer = c0
        zmunos_ana = c0
        zmunos_fac = c0
        zmunos  = c0
        zbionos = c0
        zmunao_aer = c0
        zmunao_ana = c0
        zmunao_fac = c0
        zmunao  = c0
        zbionao = c0

        DOC_tot = DOC_loc + DOCr_loc

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Obligate aerobic heterotrophs !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! limitation term for DOC uptake (DOC and DOCr)
        zlimaer_doc = (DOC_loc/parm_kaer_doc) / (1.0 + (DOC_loc/parm_kaer_doc) + (DOCr_loc/parm_kaer_docr))
        zlimaer_docr= (DOCr_loc/parm_kaer_docr)/(1.0 + (DOC_loc/parm_kaer_doc) + (DOCr_loc/parm_kaer_docr))
        zuptdoc = c0
        zuptoxy = c0
        where (DOC_tot.gt.c0 .and. O2_loc.gt.c0)
           zuptoxy = parm_aer_po2 * O2_loc     !!! diffusive limit to O2 uptake
           zuptdoc = parm_Vmax_doc * ( zlimaer_doc + zlimaer_docr )
        endwhere
        zmuaer = max(c0, min( (zuptoxy/parm_yaer_oxy), (zuptdoc/parm_yaer_doc) ) ) * Tfunc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Nitrate reducers (NO3-->NO2) !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! limitation term for DOC uptake (DOC and DOCr)
        zlimnar_doc = (DOC_loc/parm_knar_doc) / (1.0 + (DOC_loc/parm_knar_doc) + (DOCr_loc/parm_knar_docr))
        zlimnar_docr= (DOCr_loc/parm_knar_docr)/(1.0 + (DOC_loc/parm_knar_doc) + (DOCr_loc/parm_knar_docr))
        zuptdoc = c0
        zuptoxy = c0
        where (DOC_tot.gt.c0 .and. O2_loc.gt.c0)
           zuptoxy = parm_nar_po2 * O2_loc     !!! diffusive limit to O2 uptake
           zuptdoc = parm_Vmax_doc * ( zlimnar_doc + zlimnar_docr )
        endwhere
        zmunar_aer = max(c0, min( (zuptoxy/parm_ynar_oxy), (zuptdoc/parm_ynar_aer) ) ) * Tfunc
        zuptdoc = c0
        zuptno3 = c0
        where (DOC_tot.gt.c0 .and. NO3_loc.gt.c0)
           zuptno3 = parm_Vmax_no3 * ( NO3_loc / ( NO3_loc + parm_knar_no3 ) )
           zuptdoc = parm_Vmax_doc * ( zlimnar_doc + zlimnar_docr )
        endwhere
        zmunar_ana = max(c0, min( (zuptno3/parm_ynar_no3), (zuptdoc/parm_ynar_ana) ) ) * Tfunc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Nitrite reducers (NO2-->N2O) !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! limitation term for DOC uptake (DOC and DOCr)
        zlimnir_doc = (DOC_loc/parm_knir_doc) / (1.0 + (DOC_loc/parm_knir_doc) + (DOCr_loc/parm_knir_docr))
        zlimnir_docr= (DOCr_loc/parm_knir_docr)/(1.0 + (DOC_loc/parm_knir_doc) + (DOCr_loc/parm_knir_docr))
        zuptdoc = c0
        zuptoxy = c0
        where (DOC_tot.gt.c0 .and. O2_loc.gt.c0)
           zuptoxy = parm_nir_po2 * O2_loc     !!! diffusive limit to O2 uptake
           zuptdoc = parm_Vmax_doc * ( zlimnir_doc + zlimnir_docr )
        endwhere
        zmunir_aer = max(c0, min( (zuptoxy/parm_ynir_oxy), (zuptdoc/parm_ynir_aer) ) ) * Tfunc
        zuptdoc = c0
        zuptno2 = c0
        where (DOC_tot.gt.c0 .and. NO2_loc.gt.c0)
           zuptno2 = parm_Vmax_no2 * ( NO2_loc / ( NO2_loc + parm_knir_no2 ) )
           zuptdoc = parm_Vmax_doc * ( zlimnir_doc + zlimnir_docr )
        endwhere
        zmunir_ana = max(c0, min( (zuptno2/parm_ynir_no2), (zuptdoc/parm_ynir_ana) ) ) * Tfunc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Nitrous oxide reducers (N2O-->N2) !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! limitation term for DOC uptake (DOC and DOCr)
        zlimnos_doc = (DOC_loc/parm_knos_doc) / (1.0 + (DOC_loc/parm_knos_doc) + (DOCr_loc/parm_knos_docr))
        zlimnos_docr= (DOCr_loc/parm_knos_docr)/(1.0 + (DOC_loc/parm_knos_doc) + (DOCr_loc/parm_knos_docr))
        zuptdoc = c0
        zuptoxy = c0
        where (DOC_tot.gt.c0 .and. O2_loc.gt.c0)
           zuptoxy = parm_nos_po2 * O2_loc     !!! diffusive limit to O2 uptake
           zuptdoc = parm_Vmax_doc * ( zlimnos_doc + zlimnos_docr )
        endwhere
        zmunos_aer = max(c0, min( (zuptoxy/parm_ynos_oxy), (zuptdoc/parm_ynos_aer) ) ) * Tfunc
        zuptdoc = c0
        zuptn2o = c0
        where (DOC_tot.gt.c0 .and. N2O_loc.gt.c0)
           zuptn2o = parm_nos_pn2o * N2O_loc   !!! diffusive limit to N2O uptake
           zuptdoc = parm_Vmax_doc * ( zlimnos_doc + zlimnos_docr )
        endwhere
        zmunos_ana = max(c0, min( (zuptn2o/parm_ynos_n2o), (zuptdoc/parm_ynos_ana) ) ) * Tfunc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Nitrate reducers (NO3-->N2O) !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! limitation term for DOC uptake (DOC and DOCr)
        zlimnai_doc = (DOC_loc/parm_knai_doc) / (1.0 + (DOC_loc/parm_knai_doc) + (DOCr_loc/parm_knai_docr))
        zlimnai_docr= (DOCr_loc/parm_knai_docr)/(1.0 + (DOC_loc/parm_knai_doc) + (DOCr_loc/parm_knai_docr))
        zuptdoc = c0
        zuptoxy = c0
        where (DOC_tot.gt.c0 .and. O2_loc.gt.c0)
           zuptoxy = parm_nai_po2 * O2_loc     !!! diffusive limit to O2 uptake
           zuptdoc = parm_Vmax_doc * ( zlimnai_doc + zlimnai_docr )
        endwhere
        zmunai_aer = max(c0, min( (zuptoxy/parm_ynai_oxy), (zuptdoc/parm_ynai_aer) ) ) * Tfunc
        zuptdoc = c0
        zuptno3 = c0
        where (DOC_tot.gt.c0 .and. NO3_loc.gt.c0)
           zuptno3 = parm_Vmax_no3 * ( NO3_loc / ( NO3_loc + parm_knai_no3 ) )
           zuptdoc = parm_Vmax_doc * ( zlimnai_doc + zlimnai_docr )
        endwhere
        zmunai_ana = max(c0, min( (zuptno3/parm_ynai_no3), (zuptdoc/parm_ynai_ana) ) ) * Tfunc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Nitrite reducers (NO2-->N2O) !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! limitation term for DOC uptake (DOC and DOCr)
        zlimnio_doc = (DOC_loc/parm_knio_doc) / (1.0 + (DOC_loc/parm_knio_doc) + (DOCr_loc/parm_knio_docr))
        zlimnio_docr= (DOCr_loc/parm_knio_docr)/(1.0 + (DOC_loc/parm_knio_doc) + (DOCr_loc/parm_knio_docr))
        zuptdoc = c0
        zuptoxy = c0
        where (DOC_tot.gt.c0 .and. O2_loc.gt.c0)
           zuptoxy = parm_nio_po2 * O2_loc     !!! diffusive limit to O2 uptake
           zuptdoc = parm_Vmax_doc * ( zlimnio_doc + zlimnio_docr )
        endwhere
        zmunio_aer = max(c0, min( (zuptoxy/parm_ynio_oxy), (zuptdoc/parm_ynio_aer) ) ) * Tfunc
        zuptdoc = c0
        zuptno2 = c0
        where (DOC_tot.gt.c0 .and. NO2_loc.gt.c0)
           zuptno2 = parm_Vmax_no2 * ( NO2_loc / ( NO2_loc + parm_knio_no2 ) )
           zuptdoc = parm_Vmax_doc * ( zlimnio_doc + zlimnio_docr )
        endwhere
        zmunio_ana = max(c0, min( (zuptno2/parm_ynio_no2), (zuptdoc/parm_ynio_ana) ) ) * Tfunc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Nitrate reducers (NO3-->N2) !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! limitation term for DOC uptake (DOC and DOCr)
        zlimnao_doc = (DOC_loc/parm_knao_doc) / (1.0 + (DOC_loc/parm_knao_doc) + (DOCr_loc/parm_knao_docr))
        zlimnao_docr= (DOCr_loc/parm_knao_docr)/(1.0 + (DOC_loc/parm_knao_doc) + (DOCr_loc/parm_knao_docr))
        zuptdoc = c0
        zuptoxy = c0
        where (DOC_tot.gt.c0 .and. O2_loc.gt.c0)
           zuptoxy = parm_nao_po2 * O2_loc     !!! diffusive limit to O2 uptake
           zuptdoc = parm_Vmax_doc * ( zlimnao_doc + zlimnao_docr )
        endwhere
        zmunao_aer = max(c0, min( (zuptoxy/parm_ynao_oxy), (zuptdoc/parm_ynao_aer) ) ) * Tfunc
        zuptdoc = c0
        zuptno3 = c0
        where (DOC_tot.gt.c0 .and. NO3_loc.gt.c0)
           zuptno3 = parm_Vmax_no3 * ( NO3_loc / ( NO3_loc + parm_knao_no3 ) )
           zuptdoc = parm_Vmax_doc * ( zlimnao_doc + zlimnao_docr )
        endwhere
        zmunao_ana = max(c0, min( (zuptno3/parm_ynao_no3), (zuptdoc/parm_ynao_ana) ) ) * Tfunc
      
         
# ifdef FACULTATIVE_MICROBES
        ! Save occurance of facultative anaerobic heterotrophy to array (1 if anaer)
        where (zmunar_aer.lt.zmunar_ana) zmunar_fac = c1
        where (zmunai_aer.lt.zmunai_ana) zmunai_fac = c1
        where (zmunir_aer.lt.zmunir_ana) zmunir_fac = c1
        where (zmunio_aer.lt.zmunio_ana) zmunio_fac = c1
        where (zmunos_aer.lt.zmunos_ana) zmunos_fac = c1
        where (zmunao_aer.lt.zmunao_ana) zmunao_fac = c1
        ! Determine whether growth is faster on oxygen or nitrogen
        zmunar = max(zmunar_aer, zmunar_ana)
        zmunai = max(zmunai_aer, zmunai_ana)
        zmunir = max(zmunir_aer, zmunir_ana)
        zmunio = max(zmunio_aer, zmunio_ana)
        zmunos = max(zmunos_aer, zmunos_ana)
        zmunao = max(zmunao_aer, zmunao_ana)
# else
        ! Make growth anaerobic
        zmunar_fac = c1
        zmunai_fac = c1
        zmunir_fac = c1
        zmunio_fac = c1
        zmunos_fac = c1
        zmunao_fac = c1
        ! Make growth anaerobic
        zmunar = zmunar_ana
        zmunai = zmunai_ana
        zmunir = zmunir_ana
        zmunio = zmunio_ana
        zmunos = zmunos_ana
        zmunao = zmunao_ana
# endif
        ! accumulation of C biomass and rates
        zbioaer = zmuaer * AER_loc
        zbionar = zmunar * NAR_loc
        zbionai = zmunai * NAI_loc
        zbionir = zmunir * NIR_loc
        zbionio = zmunio * NIO_loc
        zbionos = zmunos * NOS_loc
        zbionao = zmunao * NAO_loc
        DENITRIF1 = zbionar * zmunar_fac * parm_ynar_no3
        DENITRIF2 = zbionir * zmunir_fac * parm_ynir_no2
        DENITRIF3 = zbionos * zmunos_fac * parm_ynos_n2o
        DENITRIF4 = zbionai * zmunai_fac * parm_ynai_no3
        DENITRIF5 = zbionio * zmunio_fac * parm_ynio_no2
        DENITRIF6 = zbionao * zmunao_fac * parm_ynao_no3


      ! Compute remineralisation (consumption) of DOM --> Inorganics
      ! NOTE:
      !  We treat the remineralisation of C,N,P differently from Fe.
      !  N and P are remineralised in their ratio with C, since heterotrophs get N and P from DOM.
      !  Fe, however, is assimilated as dFe by bacteria, not consumed as part of DOM. No DOFe pool.

        ! For N and P, first get the elemental ratios of DOM in the environment.
        qDON = ( DON_loc + epsN ) / ( DOC_loc + epsN )
        qDOP = ( DOP_loc + epsN ) / ( DOC_loc + epsN )
        qDONr = ( DONr_loc + epsN ) / ( DOCr_loc + epsN )
        qDOPr = ( DOPr_loc + epsN ) / ( DOCr_loc + epsN )

        ! Then calculate DOC consumption first based on the yields of the bacteria * biomass growth
        DOC_remin = zbioaer * parm_yaer_doc * ((zlimaer_doc + epsN) / ( zlimaer_doc + zlimaer_docr + epsN)) +
     &              ( zbionar * (parm_ynar_aer * (c1-zmunar_fac)) + zbionar * (parm_ynar_ana * zmunar_fac))
     &              * ( (zlimnar_doc + epsN) / ( zlimnar_doc + zlimnar_docr + epsN)) +
     &              ( zbionai * (parm_ynai_aer * (c1-zmunai_fac)) + zbionai * (parm_ynai_ana * zmunai_fac))
     &              * ( (zlimnai_doc + epsN) /( zlimnai_doc + zlimnai_docr + epsN)) +
     &              ( zbionir * (parm_ynir_aer * (c1-zmunir_fac)) + zbionir * (parm_ynir_ana * zmunir_fac))
     &              * ( (zlimnir_doc + epsN) /( zlimnir_doc + zlimnir_docr + epsN)) +
     &              ( zbionio * (parm_ynio_aer * (c1-zmunio_fac)) + zbionio * (parm_ynio_ana * zmunio_fac))
     &              * ( (zlimnio_doc + epsN) /( zlimnio_doc + zlimnio_docr + epsN)) +
     &              ( zbionos * (parm_ynos_aer * (c1-zmunos_fac)) + zbionos * (parm_ynos_ana * zmunos_fac))
     &              * ( (zlimnos_doc + epsN) /( zlimnos_doc + zlimnos_docr + epsN)) +
     &              ( zbionao * (parm_ynao_aer * (c1-zmunao_fac)) + zbionao * (parm_ynao_ana * zmunao_fac))
     &              * ( (zlimnao_doc + epsN) /( zlimnao_doc + zlimnao_docr + epsN))
        DOCr_remin= zbioaer * parm_yaer_doc * ((zlimaer_docr + epsN) / ( zlimaer_doc + zlimaer_docr + epsN)) +
     &              ( zbionar * (parm_ynar_aer * (c1-zmunar_fac)) + zbionar * (parm_ynar_ana * zmunar_fac))
     &              * ( (zlimnar_docr + epsN)/( zlimnar_doc + zlimnar_docr + epsN)) +
     &              ( zbionai * (parm_ynai_aer * (c1-zmunai_fac)) + zbionai * (parm_ynai_ana * zmunai_fac))
     &              * ( (zlimnai_docr + epsN)/( zlimnai_doc + zlimnai_docr + epsN)) +
     &              ( zbionir * (parm_ynir_aer * (c1-zmunir_fac)) + zbionir * (parm_ynir_ana * zmunir_fac))
     &              * ( (zlimnir_docr + epsN)/( zlimnir_doc + zlimnir_docr + epsN)) +
     &              ( zbionio * (parm_ynio_aer * (c1-zmunio_fac)) + zbionio * (parm_ynio_ana * zmunio_fac))
     &              * ( (zlimnio_docr + epsN)/( zlimnio_doc + zlimnio_docr + epsN)) +
     &              ( zbionos * (parm_ynos_aer * (c1-zmunos_fac)) + zbionos * (parm_ynos_ana * zmunos_fac))
     &              * ( (zlimnos_docr + epsN)/( zlimnos_doc + zlimnos_docr + epsN)) +
     &              ( zbionao * (parm_ynao_aer * (c1-zmunao_fac)) + zbionao * (parm_ynao_ana * zmunao_fac))
     &              * ( (zlimnao_docr + epsN)/( zlimnao_doc + zlimnao_docr + epsN))

        ! Apply the elemental ratios of DOM to the DOC consumption
        !   - This assumes that microbes indiscriminately munch on DOM that is available
        !     i.e., that they consume whatever DOM is available in the stoich that is available
        !   - This approach should also ensure that DON and DOP never become negative
        DON_remin = DOC_remin * qDON
        DOP_remin(istr:iend,j) = DOC_remin * qDOP
        DONr_remin = DOCr_remin * qDONr
        DOPr_remin(istr:iend,j) = DOCr_remin * qDOPr

      ! However, because bacteria munch indiscriminantly on available DOM, we must also account 
      ! for how variations in the stoichiometry of the DOM meets their nutritional needs.
      ! If DOM is poor in N & P, then less NH4 & PO4 will be released.
      ! If DOM is rich in N & P, then more NH4 & PO4 will be released.
      ! NOTE: DOM has to be extremely poor in nutritional quality for bacteria to take up inorganic nutrients
        DIN_release = max(c0, (DON_remin+DONr_remin) - 
     &                (zbioaer/parm_aer_CN +
     &                 zbionar/parm_nar_CN + zbionai/parm_nai_CN + zbionir/parm_nir_CN 
     &               + zbionio/parm_nio_CN + zbionos/parm_nos_CN + zbionao/parm_nao_CN) )
        DIP_release = max(c0, (DOP_remin(istr:iend,j)+DOPr_remin(istr:iend,j)) -
     &                (zbioaer/parm_aer_CP +
     &                 zbionar/parm_nar_CP + zbionai/parm_nai_CP + zbionir/parm_nir_CP 
     &               + zbionio/parm_nio_CP + zbionos/parm_nos_CP + zbionao/parm_nao_CP) )
        DIN_uptake = min(c0, (DON_remin+DONr_remin) -
     &                (zbioaer/parm_aer_CN +
     &                 zbionar/parm_nar_CN + zbionai/parm_nai_CN + zbionir/parm_nir_CN 
     &               + zbionio/parm_nio_CN + zbionos/parm_nos_CN + zbionao/parm_nao_CN) )
        DIP_uptake = min(c0, (DOP_remin(istr:iend,j)+DOPr_remin(istr:iend,j)) -
     &                (zbioaer/parm_aer_CP +
     &                 zbionar/parm_nar_CP + zbionai/parm_nai_CP + zbionir/parm_nir_CP 
     &               + zbionio/parm_nio_CP + zbionos/parm_nos_CP + zbionao/parm_nao_CP) )

        ! For Fe, bacteria assimilate dFe from their environment, competing with phytos.
        !  1. Determine how much Fe is needed.
        dFe_require = (zbioaer/parm_aer_CFe
     &               + zbionar/parm_nar_CFe + zbionai/parm_nai_CFe + zbionir/parm_nir_CFe 
     &               + zbionio/parm_nio_CFe + zbionos/parm_nos_CFe + zbionao/parm_nao_CFe)
        !  2. Determine if enouch dFe is available
        dFe_uptake  = dFe_require * c1 

        ! TEST THIS LATER !
        !where (dFe_require.gt.c0)
        !   dFe_uptake  = dFe_require * min(1.0, Fe_loc/(dFe_require * dt))
        !   DOFe_require = (c1 - dFe_uptake/dFe_require) * dFe_require
        !elsewhere
        !   dFe_uptake  = c0
        !   DOFe_require = c0
        !endwhere
        !!  3. If not, they take it from the particulate pool about to release to dFe
        !where (DOFe_require.gt.c0)
        !   P_iron_remin(istr:iend,j) = P_iron_remin(istr:iend,j) - DOFe_require
        !elsewhere
        !   DOFe_uptake  = c0
        !endwhere


# endif /* EXPLICIT_MICROBES */



! Save chemoautotrophic and heterotrophic metabolisms output
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(ammox_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(ammox_idx_t)) = ammox
       if (bec2_diag_3d_l(nitrox_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nitrox_idx_t)) = nitrox
       if (bec2_diag_3d_l(anammox_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(anammox_idx_t)) = anammox
       if (bec2_diag_3d_l(denitrif1_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(denitrif1_idx_t)) = DENITRIF1
       if (bec2_diag_3d_l(denitrif2_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(denitrif2_idx_t)) = DENITRIF2
       if (bec2_diag_3d_l(denitrif3_idx_t))
     &   bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(denitrif3_idx_t)) = DENITRIF3
       if (bec2_diag_3d_l(docremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(docremin_idx_t)) = DOC_remin
       if (bec2_diag_3d_l(donremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(donremin_idx_t)) = DON_remin
       if (bec2_diag_3d_l(dopremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(dopremin_idx_t)) = DOP_remin(istr:iend,j)
       if (bec2_diag_3d_l(donrremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(donrremin_idx_t)) = DONr_remin
#   ifdef EXPLICIT_MICROBES
       if (bec2_diag_3d_l(denitrif4_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(denitrif4_idx_t)) = DENITRIF4
       if (bec2_diag_3d_l(denitrif5_idx_t))
     &   bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(denitrif5_idx_t)) = DENITRIF5
       if (bec2_diag_3d_l(denitrif6_idx_t))
     &   bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(denitrif6_idx_t)) = DENITRIF6
       if (bec2_diag_3d_l(docrremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(docrremin_idx_t)) = DOCr_remin
       if (bec2_diag_3d_l(doprremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(doprremin_idx_t)) = DOPr_remin(istr:iend,j)
       if (bec2_diag_3d_l(aoamu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aoamu_idx_t)) = zmuaoa
       if (bec2_diag_3d_l(nobmu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nobmu_idx_t)) = zmunob
       if (bec2_diag_3d_l(nobmuana_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nobmuana_idx_t)) = zmunob_ana
       if (bec2_diag_3d_l(nobfmu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nobfmu_idx_t)) = zmufnob
       if (bec2_diag_3d_l(aoxmu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aoxmu_idx_t)) = zmuaox
       if (bec2_diag_3d_l(aoabio_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aoabio_idx_t)) = zbioaoa
       if (bec2_diag_3d_l(nobbio_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nobbio_idx_t)) = zbionob
       if (bec2_diag_3d_l(nobbioana_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nobbioana_idx_t)) = zbionob_ana
       if (bec2_diag_3d_l(aoxbio_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aoxbio_idx_t)) = zbioaox
       if (bec2_diag_3d_l(aermu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aermu_idx_t)) = zmuaer
       if (bec2_diag_3d_l(narmu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(narmu_idx_t)) = zmunar
       if (bec2_diag_3d_l(naimu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(naimu_idx_t)) = zmunai
       if (bec2_diag_3d_l(nirmu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nirmu_idx_t)) = zmunir
       if (bec2_diag_3d_l(niomu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(niomu_idx_t)) = zmunio
       if (bec2_diag_3d_l(nosmu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nosmu_idx_t)) = zmunos
       if (bec2_diag_3d_l(naomu_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(naomu_idx_t)) = zmunao
       if (bec2_diag_3d_l(aerbio_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aerbio_idx_t)) = zbioaer
       if (bec2_diag_3d_l(narbio_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(narbio_idx_t)) = zbionar
       if (bec2_diag_3d_l(naibio_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(naibio_idx_t)) = zbionai
       if (bec2_diag_3d_l(nirbio_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nirbio_idx_t)) = zbionir
       if (bec2_diag_3d_l(niobio_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(niobio_idx_t)) = zbionio
       if (bec2_diag_3d_l(nosbio_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nosbio_idx_t)) = zbionos
       if (bec2_diag_3d_l(naobio_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(naobio_idx_t)) = zbionao
       if (bec2_diag_3d_l(narana_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(narana_idx_t)) = zmunar_fac
       if (bec2_diag_3d_l(naiana_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(naiana_idx_t)) = zmunai_fac
       if (bec2_diag_3d_l(nirana_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nirana_idx_t)) = zmunir_fac
       if (bec2_diag_3d_l(nioana_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nioana_idx_t)) = zmunio_fac
       if (bec2_diag_3d_l(nosana_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nosana_idx_t)) = zmunos_fac
       if (bec2_diag_3d_l(naoana_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(naoana_idx_t)) = zmunao_fac
       if (bec2_diag_3d_l(dinrelease_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(dinrelease_idx_t)) = DIN_release
       if (bec2_diag_3d_l(diprelease_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(diprelease_idx_t)) = DIP_release
       if (bec2_diag_3d_l(dinuptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(dinuptake_idx_t)) = DIN_uptake
       if (bec2_diag_3d_l(dipuptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(dipuptake_idx_t)) = DIP_uptake
       if (bec2_diag_3d_l(dfeuptake_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(dfeuptake_idx_t)) = dFe_uptake
       if (bec2_diag_3d_l(aoamort_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aoamort_idx_t)) = zmortaoa
       if (bec2_diag_3d_l(nobmort_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nobmort_idx_t)) = zmortnob
       if (bec2_diag_3d_l(aoxmort_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aoxmort_idx_t)) = zmortaox
       if (bec2_diag_3d_l(aermort_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(aermort_idx_t)) = zmortaer
       if (bec2_diag_3d_l(narmort_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(narmort_idx_t)) = zmortnar
       if (bec2_diag_3d_l(naimort_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(naimort_idx_t)) = zmortnai
       if (bec2_diag_3d_l(nirmort_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nirmort_idx_t)) = zmortnir
       if (bec2_diag_3d_l(niomort_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(niomort_idx_t)) = zmortnio
       if (bec2_diag_3d_l(nosmort_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nosmort_idx_t)) = zmortnos
       if (bec2_diag_3d_l(naomort_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(naomort_idx_t)) = zmortnao
#   endif /* EXPLICIT_MICROBES */
#  else     
      bec2_diag_3d(istr:iend,j,k,ammox_idx_t) = ammox
      bec2_diag_3d(istr:iend,j,k,nitrox_idx_t) = nitrox
      bec2_diag_3d(istr:iend,j,k,anammox_idx_t) = anammox
      bec2_diag_3d(istr:iend,j,k,denitrif1_idx_t) = DENITRIF1
      bec2_diag_3d(istr:iend,j,k,denitrif2_idx_t) = DENITRIF2
      bec2_diag_3d(istr:iend,j,k,denitrif3_idx_t) = DENITRIF3
      bec2_diag_3d(istr:iend,j,k,docremin_idx_t) = DOC_remin
      bec2_diag_3d(istr:iend,j,k,donremin_idx_t) = DON_remin
      bec2_diag_3d(istr:iend,j,k,dopremin_idx_t) = DOP_remin(istr:iend,j)
      bec2_diag_3d(istr:iend,j,k,donrremin_idx_t) = DONr_remin
#    ifdef EXPLICIT_MICROBES
       bec2_diag_3d(istr:iend,j,k,denitrif4_idx_t) = DENITRIF4
       bec2_diag_3d(istr:iend,j,k,denitrif5_idx_t) = DENITRIF5
       bec2_diag_3d(istr:iend,j,k,denitrif6_idx_t) = DENITRIF6
       bec2_diag_3d(istr:iend,j,k,nitroxana_idx_t) = nitrox_ana
       bec2_diag_3d(istr:iend,j,k,docrremin_idx_t) = DOCr_remin
       bec2_diag_3d(istr:iend,j,k,doprremin_idx_t) = DOPr_remin(istr:iend,j)
       bec2_diag_3d(istr:iend,j,k,aoamu_idx_t) = zmuaoa
       bec2_diag_3d(istr:iend,j,k,nobmuana_idx_t) = zmunob_ana
       bec2_diag_3d(istr:iend,j,k,nobfmu_idx_t) = zmufnob
       bec2_diag_3d(istr:iend,j,k,aoxmu_idx_t) = zmuaox
       bec2_diag_3d(istr:iend,j,k,aoabio_idx_t) = zbioaoa
       bec2_diag_3d(istr:iend,j,k,nobbioana_idx_t) = zbionob_ana
       bec2_diag_3d(istr:iend,j,k,aoxbio_idx_t) = zbioaox
       bec2_diag_3d(istr:iend,j,k,aermu_idx_t) = zmuaer
       bec2_diag_3d(istr:iend,j,k,narmu_idx_t) = zmunar
       bec2_diag_3d(istr:iend,j,k,naimu_idx_t) = zmunai
       bec2_diag_3d(istr:iend,j,k,nirmu_idx_t) = zmunir
       bec2_diag_3d(istr:iend,j,k,niomu_idx_t) = zmunio
       bec2_diag_3d(istr:iend,j,k,nosmu_idx_t) = zmunos
       bec2_diag_3d(istr:iend,j,k,naomu_idx_t) = zmunao
       bec2_diag_3d(istr:iend,j,k,aerbio_idx_t) = zbioaer
       bec2_diag_3d(istr:iend,j,k,narbio_idx_t) = zbionar
       bec2_diag_3d(istr:iend,j,k,naibio_idx_t) = zbionai
       bec2_diag_3d(istr:iend,j,k,nirbio_idx_t) = zbionir
       bec2_diag_3d(istr:iend,j,k,niobio_idx_t) = zbionio
       bec2_diag_3d(istr:iend,j,k,nosbio_idx_t) = zbionos
       bec2_diag_3d(istr:iend,j,k,naobio_idx_t) = zbionao
       bec2_diag_3d(istr:iend,j,k,narana_idx_t) = zmunar_fac
       bec2_diag_3d(istr:iend,j,k,naiana_idx_t) = zmunai_fac
       bec2_diag_3d(istr:iend,j,k,nirana_idx_t) = zmunir_fac
       bec2_diag_3d(istr:iend,j,k,nioana_idx_t) = zmunio_fac
       bec2_diag_3d(istr:iend,j,k,nosana_idx_t) = zmunos_fac
       bec2_diag_3d(istr:iend,j,k,naoana_idx_t) = zmunao_fac
       bec2_diag_3d(istr:iend,j,k,dinrelease_idx_t) = DIN_release
       bec2_diag_3d(istr:iend,j,k,diprelease_idx_t) = DIP_release
       bec2_diag_3d(istr:iend,j,k,dinuptake_idx_t) = DIN_uptake
       bec2_diag_3d(istr:iend,j,k,dipuptake_idx_t) = DIP_uptake
       bec2_diag_3d(istr:iend,j,k,dfeuptake_idx_t) = dFe_uptake
       bec2_diag_3d(istr:iend,j,k,aoamort_idx_t) = zmortaoa
       bec2_diag_3d(istr:iend,j,k,nobmort_idx_t) = zmortnob
       bec2_diag_3d(istr:iend,j,k,aoxmort_idx_t) = zmortaox
       bec2_diag_3d(istr:iend,j,k,aermort_idx_t) = zmortaer
       bec2_diag_3d(istr:iend,j,k,narmort_idx_t) = zmortnar
       bec2_diag_3d(istr:iend,j,k,naimort_idx_t) = zmortnai
       bec2_diag_3d(istr:iend,j,k,nirmort_idx_t) = zmortnir
       bec2_diag_3d(istr:iend,j,k,niomort_idx_t) = zmortnio
       bec2_diag_3d(istr:iend,j,k,nosmort_idx_t) = zmortnos
       bec2_diag_3d(istr:iend,j,k,naomort_idx_t) = zmortnao
#    endif /* EXPLICIT_MICROBES */
#  endif  
# endif     

# else /* Ncycle_SY */

       WORK1 = (O2_loc - parm_o2_min) / parm_o2_min_delta
       WORK1 = min(max(WORK1,c0),c1)
       where (PAR_out < parm_nitrif_par_lim)
          NITRIF = WORK1 * parm_kappa_nitrif * NH4_loc
          where (PAR_in > parm_nitrif_par_lim)
             NITRIF = NITRIF * log(PAR_out / parm_nitrif_par_lim) / (-KPARdz)
          end where
       elsewhere
          NITRIF = c0
       end where

       WORK1 = ((parm_o2_min + parm_o2_min_delta) - O2_loc) / parm_o2_min_delta
       WORK1 = min(max(WORK1,c0),c1)

       WORK1 = merge(c0, WORK1, NO3_loc == c0)
       DENITRIF = WORK1 * ((DOC_remin
     &                  + POC_remin(istr:iend,j)
     &                  - OTHER_REMIN) / denitrif_C_N - SED_DENITRIF)

       where (NO3_loc < parm_denitrif_NO3_limit)
            DENITRIF = 1/parm_denitrif_NO3_limit * NO3_loc * DENITRIF
       end where


# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(nitrif_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(nitrif_idx_t)) = NITRIF
       if (bec2_diag_3d_l(denitrif_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(denitrif_idx_t)) = DENITRIF
#  else     
       bec2_diag_3d(istr:iend,j,k,nitrif_idx_t) = NITRIF
       bec2_diag_3d(istr:iend,j,k,denitrif_idx_t) = DENITRIF
#  endif  
# endif   

# endif /* Ncycle_SY */


# ifdef NO3_RESTORE
       if (k .eq. N)
             restore_no3 = c1/2592000 * (no3_restore(istr:iend,j,k)-NO3_loc)
       else
             restore_no3 = c0
# endif



    !-------------------------------!
    !-------------------------------!
    !-------------------------------!
    !   (14) SOURCES AND SINKS      !
    !-------------------------------!
    !-------------------------------!
    !-------------------------------!

# ifdef Ncycle_SY


     !!!!!!!!!!!!!!!!!
     ! Nitrate (NO3) !
     !!!!!!!!!!!!!!!!!
       DTRACER_MODULE(istr:iend,j,k,no3_ind_t) = nitrox + nitrox_ana
     &                                           - DENITRIF1 - DENITRIF4 - DENITRIF6
     &                                           - SED_DENITRIF
     &                                           - sum(NO3_V, dim=2)
# ifdef NO3_RESTORE
     &                                           + restore_no3
# endif
# ifdef EXPLICIT_MICROBES
     &                                           + ( zbioaox * parm_paox_no3 )
# endif /* EXPLICIT_MICROBES */


     !!!!!!!!!!!!!!!!!
     ! Nitrite (NO2) !
     !!!!!!!!!!!!!!!!!
       DTRACER_MODULE(istr:iend,j,k,no2_ind_t) = DENITRIF1 - DENITRIF2 - DENITRIF5
     &                                           - nitrox - nitrox_ana
     &                                           - sum(NO2_V, dim=2)
# ifdef EXPLICIT_MICROBES
     &                                           - (zbioaox * parm_yaox_no2)
     &                                           + ( zbioaoa * (parm_yaoa_nh4 - c1/parm_aoa_CN) * (c1-yaoa_n2o) )
# else
     &                                           - r_no2tonh4_ax * anammox
     &                                           + ammox * no2_nh4_yield1
# endif /* EXPLICIT_MICROBES */


     !!!!!!!!!!!!!!!!!!
     ! Ammonium (NH4) !
     !!!!!!!!!!!!!!!!!!
       DTRACER_MODULE(istr:iend,j,k,nh4_ind_t) = zoo_loss_dic * Q
     &                                           + sum(auto_loss_dic, dim=2) * Q
     &                                           + sum(auto_graze_dic, dim=2) * Q
     &                                           - sum(NH4_V, dim=2)
     &                                           - ammox - anammox
# ifdef EXPLICIT_MICROBES
     &                                           + DIN_release - DIN_uptake
     &                                           - ( zbionob / parm_nob_CN )
     &                                           - ( zbionob_ana / parm_nob_CN )
     &                                           + zmortaoa_dic/parm_aoa_CN
     &                                           + zmortnob_dic/parm_nob_CN
     &                                           + zmortaox_dic/parm_aox_CN
     &                                           + zmortaer_dic/parm_aer_CN
     &                                           + zmortnar_dic/parm_nar_CN
     &                                           + zmortnai_dic/parm_nai_CN
     &                                           + zmortnir_dic/parm_nir_CN
     &                                           + zmortnio_dic/parm_nio_CN
     &                                           + zmortnos_dic/parm_nos_CN
     &                                           + zmortnao_dic/parm_nao_CN
     &                                           + zzooaoa_dic/parm_aoa_CN
     &                                           + zzoonob_dic/parm_nob_CN
     &                                           + zzooaox_dic/parm_aox_CN
     &                                           + (zzooaoa_excN + zzoonob_excN + zzooaox_excN)
     &                                           + zszooaer_dic/parm_aer_CN
     &                                           + zszoonar_dic/parm_nar_CN
     &                                           + zszoonai_dic/parm_nai_CN
     &                                           + zszoonir_dic/parm_nir_CN
     &                                           + zszoonio_dic/parm_nio_CN
     &                                           + zszoonos_dic/parm_nos_CN
     &                                           + zszoonao_dic/parm_nao_CN
     &                                           + (zszoonar_excN
     &                                           +  zszoonar_excN + zszoonai_excN + zszoonir_excN 
     &                                           +  zszoonio_excN + zszoonos_excN + zszoonao_excN)
     &                                           + Q * szoo_loss_dic
# else
     &                                           + DON_remin + DONr_remin
     &                                           + POC_remin(istr:iend,j) * (c1 - DONrefract) * Q
# endif /* EXPLICIT_MICROBES */


     !!!!!!!!!!!!!!!!!!!!!!!
     ! Nitrous Oxide (N2O) !
     !!!!!!!!!!!!!!!!!!!!!!!
       DTRACER_MODULE(istr:iend,j,k,n2o_ind_t) = 0.5*DENITRIF2 + 0.5*DENITRIF4 - DENITRIF3
# ifdef EXPLICIT_MICROBES 
     &                                           + ( zbioaoa * (parm_yaoa_nh4 - c1/parm_aoa_CN) * yaoa_n2o )
# else 
     &                                           + ammox * n2o_nh4_yield1
# endif /* EXPLICIT_MICROBES */


     !!!!!!!!!!!!!!!!!!!
     ! Dinitrogen (N2) !
     !!!!!!!!!!!!!!!!!!!
       DTRACER_MODULE(istr:iend,j,k,n2_ind_t) = DENITRIF3 + 0.5*DENITRIF5 + 0.5*DENITRIF6
     &                                          + 0.5*SED_DENITRIF + anammox


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Nitrous Oxide components !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ifdef N2O_TRACER_DECOMP
       WORK1 = N2O_ao1_loc + N2O_siden_loc + N2O_atm_loc + N2O_soden_loc + 4 * epsN
       DTRACER_MODULE(istr:iend,j,k,n2o_ao1_ind_t) = ammox *
# ifdef EXPLICIT_MICROBES
     &                                           yaoa_n2o
# else
     &                                           n2o_nh4_yield1
# endif /* EXPLICIT_MICROBES */
     &      - (N2O_ao1_loc + epsN) / WORK1 * DENITRIF3
       DTRACER_MODULE(istr:iend,j,k,n2o_soden_ind_t) = 0.5 * DENITRIF2 -
     &      (N2O_soden_loc + epsN) / WORK1 * DENITRIF3
       DTRACER_MODULE(istr:iend,j,k,n2o_siden_ind_t) =
     &      - (N2O_siden_loc + epsN) / WORK1 * DENITRIF3
       DTRACER_MODULE(istr:iend,j,k,n2o_atm_ind_t) =
     &      - (N2O_atm_loc + epsN) / WORK1 * DENITRIF3 
 
       DTRACER_MODULE(istr:iend,j,k,n2_sed_ind_t) = 
     &      0.5 * SED_DENITRIF 
 
# ifdef BEC2_DIAG 
      bec2_diag_3d(istr:iend,j,k,n2oammox_idx_t) = ammox * 
# ifdef EXPLICIT_MICROBES 
     &                                           yaoa_n2o 
# else 
     &                                           n2o_nh4_yield1 
# endif /* EXPLICIT_MICROBES */ 
      WORK1 = N2O_ao1_loc + N2O_siden_loc + N2O_atm_loc + 
     & N2O_soden_loc + 4 * epsN 
      bec2_diag_3d(istr:iend,j,k,n2oao1_cons_idx_t) = 
     &     (N2O_ao1_loc + epsN) / WORK1 * DENITRIF3 
      bec2_diag_3d(istr:iend,j,k,n2osoden_cons_idx_t) = 
     &     (N2O_soden_loc + epsN) / WORK1 * DENITRIF3 
      bec2_diag_3d(istr:iend,j,k,n2osiden_cons_idx_t) = 
     &     (N2O_siden_loc + epsN) / WORK1 * DENITRIF3 
      bec2_diag_3d(istr:iend,j,k,n2oatm_cons_idx_t) = 
     &     (N2O_atm_loc + epsN) / WORK1 * DENITRIF3 
 
# endif 
# endif 


# else /* Ncycle_SY */
       DTRACER_MODULE(istr:iend,j,k,no3_ind_t) = NITRIF - DENITRIF - 
     &                                 SED_DENITRIF - sum(NO3_V, dim=2)
# ifdef NO3_RESTORE
     &                                 + restore_no3
# endif
       DTRACER_MODULE(istr:iend,j,k,nh4_ind_t) = -sum(NH4_V, dim=2) - NITRIF + DON_remin + DONr_remin
     &     + Q * (zoo_loss_dic + sum(auto_loss_dic, dim=2) + sum(auto_graze_dic, dim=2)
     &     + POC_remin(istr:iend,j) * (c1 - DONrefract))
# endif /* Ncycle_SY */


       do auto_ind = 1, autotroph_cnt
          if (Nfixer(auto_ind))
     &        DTRACER_MODULE(istr:iend,j,k,nh4_ind_t) = DTRACER_MODULE(istr:iend,j,k,nh4_ind_t) 
     &                                          + Nexcrete(istr:iend,auto_ind)
       end do


    !-----------------------------------------------------------------------
    ! Excess Iodine 
    !-----------------------------------------------------------------------
# ifdef IODINE
       DTRACER_MODULE(istr:iend,j,k,iodxs1_ind_t) = iod_sed_flux1
       DTRACER_MODULE(istr:iend,j,k,iodxs2_ind_t) = iod_sed_flux2
       DTRACER_MODULE(istr:iend,j,k,iodxs3_ind_t) = iod_sed_flux3
       DTRACER_MODULE(istr:iend,j,k,iodxs4_ind_t) = iod_sed_flux4
# endif


    !-----------------------------------------------------------------------
    !  dissolved iron
    !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,fe_ind_t) = P_iron_remin(istr:iend,j)
     &                                          + (Qfe_zoo * zoo_loss_dic) 
     &                                          + DOFe_remin 
     &                                          - sum(photoFe, dim=2) - Fe_scavenge
# ifdef EXPLICIT_MICROBES
     &                                          - dFe_uptake + DOFe_prod
     &                                          - zbioaoa/parm_aoa_CFe
     &                                          - zbionob/parm_nob_CFe
     &                                          - zbionob_ana/parm_nob_CFe
     &                                          - zbioaox/parm_aox_CFe
     &                                          + zmortaoa_dic/parm_aoa_CFe
     &                                          + zmortnob_dic/parm_nob_CFe
     &                                          + zmortaox_dic/parm_aox_CFe
     &                                          + zmortaer_dic/parm_aer_CFe
     &                                          + zmortnar_dic/parm_nar_CFe
     &                                          + zmortnai_dic/parm_nai_CFe
     &                                          + zmortnir_dic/parm_nir_CFe
     &                                          + zmortnio_dic/parm_nio_CFe
     &                                          + zmortnos_dic/parm_nos_CFe
     &                                          + zmortnao_dic/parm_nao_CFe
     &                                          + zzooaoa_dic/parm_aoa_CFe
     &                                          + zzoonob_dic/parm_nob_CFe
     &                                          + zzooaox_dic/parm_aox_CFe
     &                                          + zszooaer_dic/parm_aer_CFe
     &                                          + zszoonar_dic/parm_nar_CFe
     &                                          + zszoonai_dic/parm_nai_CFe
     &                                          + zszoonir_dic/parm_nir_CFe
     &                                          + zszoonio_dic/parm_nio_CFe
     &                                          + zszoonos_dic/parm_nos_CFe
     &                                          + zszoonao_dic/parm_nao_CFe
     &                                          + szoo_loss_dic * Qfe_zoo
     &                                          + (zzooaoa_excFe + zzoonob_excFe + zzooaox_excFe)
     &                                          + (zszooaer_excFe 
     &                                          +  zszoonar_excFe + zszoonai_excFe + zszoonir_excFe 
     &                                           + zszoonio_excFe + zszoonos_excFe + zszoonao_excFe)
# else
     &                                          + DOFe_remin
# endif /* EXPLICIT_MICROBES */

       do auto_ind = 1, autotroph_cnt
          DTRACER_MODULE(istr:iend,j,k,fe_ind_t) = DTRACER_MODULE(istr:iend,j,k,fe_ind_t)
     &      + (Qfe(istr:iend,auto_ind) * (auto_loss_dic(istr:iend,auto_ind) + auto_graze_dic(istr:iend,auto_ind)))
     &      + auto_graze_zoo(istr:iend,auto_ind) * (Qfe(istr:iend,auto_ind)-Qfe_zoo)
       end do

    !-----------------------------------------------------------------------
    !  dissolved SiO3
    !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,sio3_ind_t) = P_SiO2_remin(istr:iend,j)

       do auto_ind = 1, autotroph_cnt
          if (kSiO3(auto_ind) > 0) then
             DTRACER_MODULE(istr:iend,j,k,sio3_ind_t) = DTRACER_MODULE(istr:iend,j,k,sio3_ind_t) 
     &                 - photoSi(istr:iend,auto_ind) + Qsi(istr:iend,auto_ind)
     &                 * (f_graze_si_remin * auto_graze(istr:iend,auto_ind)
     &                 + (c1 - loss_poc(auto_ind)) * auto_loss(istr:iend,auto_ind))
          endif
       end do

    !-----------------------------------------------------------------------
    !  phosphate
    !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,po4_ind_t) = Qp_zoo_pom * zoo_loss_dic - sum(PO4_V, dim=2)
# ifdef EXPLICIT_MICROBES
     &                                           + DIP_release - DIP_uptake
     &                                           - zbioaoa/parm_aoa_CP
     &                                           - zbionob/parm_nob_CP
     &                                           - zbionob_ana/parm_nob_CP
     &                                           - zbioaox/parm_aox_CP
     &                                           + zmortaoa_dic/parm_aoa_CP
     &                                           + zmortnob_dic/parm_nob_CP
     &                                           + zmortaox_dic/parm_aox_CP
     &                                           + zmortaer_dic/parm_aer_CP
     &                                           + zmortnar_dic/parm_nar_CP
     &                                           + zmortnai_dic/parm_nai_CP
     &                                           + zmortnir_dic/parm_nir_CP
     &                                           + zmortnio_dic/parm_nio_CP
     &                                           + zmortnos_dic/parm_nos_CP
     &                                           + zmortnao_dic/parm_nao_CP
     &                                           + zzooaoa_dic/parm_aoa_CP
     &                                           + zzoonob_dic/parm_nob_CP
     &                                           + zzooaox_dic/parm_aox_CP
     &                                           + zszooaer_dic/parm_aer_CP
     &                                           + zszoonar_dic/parm_nar_CP
     &                                           + zszoonai_dic/parm_nai_CP
     &                                           + zszoonir_dic/parm_nir_CP
     &                                           + zszoonio_dic/parm_nio_CP
     &                                           + zszoonos_dic/parm_nos_CP
     &                                           + zszoonao_dic/parm_nao_CP
     &                                           + szoo_loss_dic * Qp_zoo_pom
     &                                           + (zszooaer_excP
     &                                           +  zszoonar_excP + zszoonai_excP + zszoonir_excP 
     &                                           +  zszoonio_excP + zszoonos_excP + zszoonao_excP)
     &                                           + (zzooaoa_excP + zzoonob_excP + zzooaox_excP)
# else
     &                                           + DOP_remin(istr:iend,j) + DOPr_remin(istr:iend,j)
     &                                           + Qp_zoo_pom * (c1 - DOPrefract) * POC_remin(istr:iend,j)
# endif /* EXPLICIT_MICROBES */

       do auto_ind = 1, autotroph_cnt
          if (Qp(auto_ind) == Qp_zoo_pom) then
             DTRACER_MODULE(istr:iend,j,k,po4_ind_t) = DTRACER_MODULE(istr:iend,j,k,po4_ind_t)
     &       + Qp(auto_ind) * (auto_loss_dic(istr:iend,auto_ind) + auto_graze_dic(istr:iend,auto_ind))
          else
             DTRACER_MODULE(istr:iend,j,k,po4_ind_t) = DTRACER_MODULE(istr:iend,j,k,po4_ind_t) 
     &                  + remaining_P_dip(istr:iend,auto_ind)
          endif
       end do

    !-----------------------------------------------------------------------
    !  autotroph Carbon
    !  autotroph Chlorophyll
    !  autotroph Fe
    !  autotroph Si
    !  autotroph CaCO3
    !-----------------------------------------------------------------------

       do auto_ind = 1, autotroph_cnt
          WORK1 = auto_graze(istr:iend,auto_ind) + auto_loss(istr:iend,auto_ind) + auto_agg(istr:iend,auto_ind)

          n_loc = C_ind(auto_ind)
          DTRACER_MODULE(istr:iend,j,k,n_loc) = photoC(istr:iend,auto_ind) - WORK1

          n_loc = Chl_ind(auto_ind)
          DTRACER_MODULE(istr:iend,j,k,n_loc) = photoacc(istr:iend,auto_ind) - thetaC(istr:iend,auto_ind) * WORK1

          n_loc = Fe_ind(auto_ind)
          DTRACER_MODULE(istr:iend,j,k,n_loc) =  photoFe(istr:iend,auto_ind) - Qfe(istr:iend,auto_ind) * WORK1

          n_loc = Si_ind(auto_ind)
          if (n_loc > 0) then
             DTRACER_MODULE(istr:iend,j,k,n_loc) =  photoSi(istr:iend,auto_ind) - Qsi(istr:iend,auto_ind) * WORK1
          endif

          n_loc = CaCO3_ind(auto_ind)
          if (n_loc > 0) then
             DTRACER_MODULE(istr:iend,j,k,n_loc) = CaCO3_PROD(istr:iend,auto_ind) 
     &              - QCaCO3(istr:iend,auto_ind) * WORK1
          endif
       end do

    !-----------------------------------------------------------------------
    !  heterotrophs and chemoautotrophs
    !-----------------------------------------------------------------------

# ifdef EXPLICIT_MICROBES
       DTRACER_MODULE(istr:iend,j,k,aoa_ind_t) = zbioaoa - zzooaoa - zmortaoa
       DTRACER_MODULE(istr:iend,j,k,nob_ind_t) = zbionob + zbionob_ana - zzoonob - zmortnob! - zmortnob2
       DTRACER_MODULE(istr:iend,j,k,aox_ind_t) = zbioaox - zzooaox - zmortaox
       DTRACER_MODULE(istr:iend,j,k,aer_ind_t) = zbioaer - zszooaer - zmortaer
       DTRACER_MODULE(istr:iend,j,k,nar_ind_t) = zbionar - zszoonar - zmortnar
       DTRACER_MODULE(istr:iend,j,k,nai_ind_t) = zbionai - zszoonai - zmortnai
       DTRACER_MODULE(istr:iend,j,k,nir_ind_t) = zbionir - zszoonir - zmortnir
       DTRACER_MODULE(istr:iend,j,k,nio_ind_t) = zbionio - zszoonio - zmortnio
       DTRACER_MODULE(istr:iend,j,k,nos_ind_t) = zbionos - zszoonos - zmortnos
       DTRACER_MODULE(istr:iend,j,k,nao_ind_t) = zbionao - zszoonao - zmortnao
# endif /* EXPLICIT_MICROBES */

    !-----------------------------------------------------------------------
    !  zoo Carbon
    !-----------------------------------------------------------------------

# ifdef EXPLICIT_MICROBES
       DTRACER_MODULE(istr:iend,j,k,zooC_ind_t) = sum(auto_graze_zoo, dim=2) - zoo_loss 
     &                                          + zzooaoa_zoo + zzoonob_zoo + zzooaox_zoo
       DTRACER_MODULE(istr:iend,j,k,szooc_ind_t) = zszooaer_szoo 
     &                                           + zszoonar_szoo + zszoonai_szoo + zszoonir_szoo 
     &                                           + zszoonio_szoo + zszoonos_szoo + zszoonao_szoo
     &                                           - szoo_loss
# else
       DTRACER_MODULE(istr:iend,j,k,zooC_ind_t) = sum(auto_graze_zoo, dim=2) - zoo_loss
# endif /* EXPLICIT_MICROBES */

    !-----------------------------------------------------------------------
    !  dissolved organic matter
    !  from sinking remin small fraction to refractory pool
    !-----------------------------------------------------------------------

# ifdef EXPLICIT_MICROBES
       ! labile DOC !
       DTRACER_MODULE(istr:iend,j,k,doc_ind_t) = DOC_prod * (c1-DOCrefract) - DOC_remin
       ! labile DON !
       DTRACER_MODULE(istr:iend,j,k,don_ind_t) = DON_prod * (c1-DONrefract) - DON_remin
       ! labile DOP !
       DTRACER_MODULE(istr:iend,j,k,dop_ind_t) = DOP_prod * (c1-DOPrefract) - DOP_remin(istr:iend,j) - sum(DOP_V, dim=2)
       ! labile DOFe !
       DTRACER_MODULE(istr:iend,j,k,dofe_ind_t) = c0  ! There is no such distinction between DOFe and dFe
       ! refractory DOC
       DTRACER_MODULE(istr:iend,j,k,docr_ind_t) = DOC_prod * DOCrefract - DOCr_remin
       ! refractory DON
       DTRACER_MODULE(istr:iend,j,k,donr_ind_t) = DON_prod * DONrefract - DONr_remin
       ! refractory DOP
       DTRACER_MODULE(istr:iend,j,k,dopr_ind_t) = DOP_prod * DOPrefract - DOPr_remin(istr:iend,j)
# else
       DTRACER_MODULE(istr:iend,j,k,doc_ind_t) = DOC_prod - DOC_remin
       DTRACER_MODULE(istr:iend,j,k,don_ind_t) = (DON_prod * (c1 - DONrefract)) - DON_remin
       DTRACER_MODULE(istr:iend,j,k,dop_ind_t) = (DOP_prod * (c1 - DOPrefract)) - DOP_remin(istr:iend,j)
       DTRACER_MODULE(istr:iend,j,k,dofe_ind_t) = DOFe_prod - DOFe_remin
       DTRACER_MODULE(istr:iend,j,k,donr_ind_t) = (DON_prod * DONrefract) - DONr_remin
     &                                            + ( POC_remin(istr:iend,j) * DONrefract * Q)
       DTRACER_MODULE(istr:iend,j,k,dopr_ind_t) = (DOP_prod * DOPrefract) - DOPr_remin(istr:iend,j)
     &                                            + ( POC_remin(istr:iend,j) * DOPrefract * Qp_zoo_pom)
# endif /* EXPLICIT_MICROBES */


    !-----------------------------------------------------------------------
    !   dissolved inorganic Carbon
    !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,dic_ind_t) =
     &       sum(auto_loss_dic, dim=2) + sum(auto_graze_dic, dim=2) - sum(photoC, dim=2)
     &       + P_CaCO3_remin(istr:iend,j)
     &       + zoo_loss_dic
# ifdef EXPLICIT_MICROBES
     &       - zbioaoa + zzooaoa_dic + zmortaoa_dic
     &       - zbionob - zbionob_ana + zzoonob_dic + zmortnob_dic
     &       - zbioaox + zzooaox_dic + zmortaox_dic
     &       + (zbioaer*(parm_yaer_doc-c1))
     &       + (zbionar*(parm_ynar_aer-c1)*(c1-zmunar_fac)) + (zbionar*(parm_ynar_ana-c1)*zmunar_fac)
     &       + (zbionai*(parm_ynai_aer-c1)*(c1-zmunai_fac)) + (zbionai*(parm_ynai_ana-c1)*zmunai_fac)
     &       + (zbionir*(parm_ynir_aer-c1)*(c1-zmunir_fac)) + (zbionir*(parm_ynir_ana-c1)*zmunir_fac)
     &       + (zbionio*(parm_ynio_aer-c1)*(c1-zmunio_fac)) + (zbionio*(parm_ynio_ana-c1)*zmunio_fac)
     &       + (zbionos*(parm_ynos_aer-c1)*(c1-zmunos_fac)) + (zbionos*(parm_ynos_ana-c1)*zmunos_fac)
     &       + (zbionao*(parm_ynao_aer-c1)*(c1-zmunao_fac)) + (zbionao*(parm_ynao_ana-c1)*zmunao_fac)
     &       + zmortaer_dic + zszooaer_dic
     &       + zmortnar_dic + zszoonar_dic
     &       + zmortnai_dic + zszoonai_dic
     &       + zmortnir_dic + zszoonir_dic
     &       + zmortnio_dic + zszoonio_dic
     &       + zmortnos_dic + zszoonos_dic
     &       + zmortnao_dic + zszoonao_dic
     &       + szoo_loss_dic
# else 
     &       + DOC_remin
     &       + POC_remin(istr:iend,j)
# endif /* EXPLICIT_MICROBES */

       do auto_ind = 1, autotroph_cnt
          if (CaCO3_ind(auto_ind) > 0)
     &        DTRACER_MODULE(istr:iend,j,k,dic_ind_t) = DTRACER_MODULE(istr:iend,j,k,dic_ind_t)
     &        + f_graze_CaCO3_REMIN * auto_graze(istr:iend,auto_ind) * QCaCO3(istr:iend,auto_ind)
     &        - CaCO3_PROD(istr:iend,auto_ind)
       end do

    !-----------------------------------------------------------------------
    !  alkalinity
    !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,alk_ind_t) = -DTRACER_MODULE(istr:iend,j,k,no3_ind_t)
     &                                           - DTRACER_MODULE(istr:iend,j,k,no2_ind_t)
     &                                           + DTRACER_MODULE(istr:iend,j,k,nh4_ind_t)
     &                                           + c2 * P_CaCO3_remin(istr:iend,j)

       do auto_ind = 1, autotroph_cnt
          if (CaCO3_ind(auto_ind) > 0)
     &        DTRACER_MODULE(istr:iend,j,k,alk_ind_t) = DTRACER_MODULE(istr:iend,j,k,alk_ind_t)
     &        + c2 * (f_graze_CaCO3_REMIN * auto_graze(istr:iend,auto_ind) * QCaCO3(istr:iend,auto_ind)
     &        - CaCO3_PROD(istr:iend,auto_ind))
       end do

    !-----------------------------------------------------------------------
    !  oxygen
    !-----------------------------------------------------------------------

       O2_PRODUCTION = c0
# ifdef Ncycle_SY
       do auto_ind = 1, autotroph_cnt
          if (.not. Nfixer(auto_ind)) then
             where (photoC(istr:iend,auto_ind) > c0)
                O2_PRODUCTION = O2_PRODUCTION + photoC(istr:iend,auto_ind) *
     &          ((NO3_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) 
     &          + NO2_V(istr:iend,auto_ind) + NH4_V(istr:iend,auto_ind))) / parm_Red_D_C_O2
     &          + (NH4_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) 
     &          + NO2_V(istr:iend,auto_ind) + NH4_V(istr:iend,auto_ind))) / parm_Remin_D_C_O2
     &          + (NO2_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) 
     &          + NO2_V(istr:iend,auto_ind) + NH4_V(istr:iend,auto_ind))) / parm_Red_D_C_O2_NO2V)
             end where
          else
             where (photoC(istr:iend,auto_ind) > c0)
                O2_PRODUCTION = O2_PRODUCTION + photoC(istr:iend,auto_ind) *
     &          ((NO3_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) 
     &             + NO2_V(istr:iend,auto_ind) + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Red_D_C_O2
     &           + (NO2_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) + NO2_V(istr:iend,auto_ind) 
     &           + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Red_D_C_O2_NO2V
     &           + (NH4_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) + NO2_V(istr:iend,auto_ind) 
     &           + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Remin_D_C_O2
     &           + (Nfix(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) + NO2_V(istr:iend,auto_ind)
     &           + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Red_D_C_O2_diaz)
             end where
          endif
       end do
#else
       do auto_ind = 1, autotroph_cnt
          if (.not. Nfixer(auto_ind)) then
             where (photoC(istr:iend,auto_ind) > c0)
                O2_PRODUCTION = O2_PRODUCTION + photoC(istr:iend,auto_ind) *
     &          ((NO3_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) 
     &          + NH4_V(istr:iend,auto_ind))) / parm_Red_D_C_O2
     &          + (NH4_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) 
     &          + NH4_V(istr:iend,auto_ind))) / parm_Remin_D_C_O2)
             end where
          else
             where (photoC(istr:iend,auto_ind) > c0)
                O2_PRODUCTION = O2_PRODUCTION + photoC(istr:iend,auto_ind) *
     &          ((NO3_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) 
     &             + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Red_D_C_O2
     &           + (NH4_V(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) 
     &           + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Remin_D_C_O2
     &           + (Nfix(istr:iend,auto_ind) / (NO3_V(istr:iend,auto_ind) 
     &           + NH4_V(istr:iend,auto_ind) + Nfix(istr:iend,auto_ind))) / parm_Red_D_C_O2_diaz)
             end where
          endif
       end do
# endif

# ifdef Ncycle_SY

       where (O2_loc .le. c0)
          WORK1 = c0
       elsewhere
          WORK1 = roxic
       endwhere

       O2_CONSUMPTION = WORK1 * (c0 - (SED_DENITRIF*denitrif_C_N) - OTHER_REMIN) / parm_Remin_D_C_O2
     &                + (zoo_loss_dic + sum(auto_loss_dic, dim=2) + sum(auto_graze_dic, dim=2)) / parm_Remin_D_C_O2
# ifdef EXPLICIT_MICROBES
     &   + zbioaoa*parm_yaoa_oxy
     &   + zbionob*parm_ynob_oxy
     &   + ( zzooaoa_dic + zzoonob_dic + zzooaox_dic ) / parm_Remin_D_C_O2
     &   + zbioaer*parm_yaer_oxy
     &   + zbionar*parm_ynar_oxy*(c1-zmunar_fac)
     &   + zbionai*parm_ynai_oxy*(c1-zmunai_fac)
     &   + zbionir*parm_ynir_oxy*(c1-zmunir_fac)
     &   + zbionio*parm_ynio_oxy*(c1-zmunio_fac)
     &   + zbionos*parm_ynos_oxy*(c1-zmunos_fac)
     &   + zbionao*parm_ynao_oxy*(c1-zmunao_fac)
     &   + ( zszooaer_dic + zszoonar_dic + zszoonai_dic + zszoonir_dic 
     &     + zszoonio_dic + zszoonos_dic + zszoonao_dic ) / parm_Remin_D_C_O2
     &   +  szoo_loss_dic / parm_Remin_D_C_O2
# else
     &   + 1.5 * ammox + 0.5 * nitrox
     &   + DOC_remin / parm_Remin_D_C_O2
     &   + POC_remin(istr:iend,j) / parm_Remin_D_C_O2
# endif /* EXPLICIT_MICROBES */

# else

       WORK1 = (O2_loc - parm_o2_min) / parm_o2_min_delta
       WORK1 = min(max(WORK1,c0),c1)

       O2_CONSUMPTION = WORK1 *(( POC_remin(istr:iend,j) + DOC_remin - 
     &                            SED_DENITRIF*denitrif_C_N - OTHER_REMIN + 
     &                            zoo_loss_dic + sum(auto_loss_dic, dim=2) + sum(auto_graze_dic, dim=2) ) / parm_Remin_D_C_O2 
     &                            + (c2 * NITRIF))

# endif

       DTRACER_MODULE(istr:iend,j,k,o2_ind_t) = O2_PRODUCTION - O2_CONSUMPTION

# ifdef BEC2_DIAG
      bec2_diag_3d(istr:iend,j,k,o2prod_idx_t) = O2_PRODUCTION
      bec2_diag_3d(istr:iend,j,k,o2cons_idx_t) = O2_CONSUMPTION
# endif        



# ifdef N2O_NEV
        ! Nevison Param: N2O production 
        where (O2_loc .gt. O2_crit_nev .and. PAR_out < parm_nitrif_par_lim)
                n2o_prod_nev = (parm_Red_D_C_O2/parm_Red_D_C_N) * (parm_n2o_nev_a1/O2_loc + parm_n2o_nev_a2)
     &          * exp(-z_r(i,j,k)/z_scale_nev) * (O2_CONSUMPTION-O2_PRODUCTION)
               where (PAR_in > parm_nitrif_par_lim)
                        n2o_prod_nev = n2o_prod_nev * log(PAR_out/parm_nitrif_par_lim) / (-KPARdz)
               endwhere
        elsewhere 
                n2o_prod_nev = c0        
        end where
        ! First order consumption
        where (O2_loc .le. O2_crit_nev .and. O2_loc .gt. c0)
                n2o_cons_nev = N2O_cons_tau_nev * N2O_nev_loc 
        elsewhere
                n2o_cons_nev = c0
        endwhere

       DTRACER_MODULE(istr:iend,j,k,n2o_nev_ind_t) = n2o_prod_nev - n2o_cons_nev
# ifdef BEC2_DIAG
      bec2_diag_3d(istr:iend,j,k,n2oprodnev_idx_t) = n2o_prod_nev
      bec2_diag_3d(istr:iend,j,k,n2oconsnev_idx_t) = n2o_cons_nev
# endif
# endif

# ifdef BEC_COCCO
    !-----------------------------------------------------------------------
    ! CN, June/July 2015: Implementation of coccolilthophores 
    !-----------------------------------------------------------------------


    !--------------------------------
    ! 
    ! detached CaCO3: from explicit calcifiers only
    !
    !--------------------------------

       do auto_ind = 1, autotroph_cnt
          if (exp_calcifier(auto_ind)) then
            DTRACER_MODULE(istr:iend,j,k,cal_ind_t) = c0
!            DTRACER_MODULE(istr:iend,j,k,cal_ind_t) = P_CaCO3_prod(istr:iend) - P_CaCO3_remin(istr:iend,j)
!     &          - P_CaCO3_sed_loss(istr:iend,j)
          endif
       end do
# endif /* BEC_COCCO */

# ifdef USE_EXPLICIT_VSINK
    !-----------------------------------------------------------------------
    !  Particulate pools
    !-----------------------------------------------------------------------

       DTRACER_MODULE(istr:iend,j,k,dusthard_ind_t) = -dusthard_remin(istr:iend,j)
       DTRACER_MODULE(istr:iend,j,k,pochard_ind_t) = POC_gamma*POC_prod - POChard_remin(istr:iend,j)
       DTRACER_MODULE(istr:iend,j,k,pcaco3hard_ind_t) = P_CaCO3_gamma*P_CaCO3_prod - P_CaCO3hard_remin(istr:iend,j)
       DTRACER_MODULE(istr:iend,j,k,psio2hard_ind_t) = P_SiO2_gamma*P_SiO2_prod - P_SiO2hard_remin(istr:iend,j)
       DTRACER_MODULE(istr:iend,j,k,pironhard_ind_t) = P_iron_gamma*P_iron_prod - P_ironhard_remin(istr:iend,j)
       DTRACER_MODULE(istr:iend,j,k,dustsoft_ind_t) = -dustsoft_remin(istr:iend,j)
       DTRACER_MODULE(istr:iend,j,k,pocsoft_ind_t) = (c1-POC_gamma)*POC_prod - POCsoft_remin(istr:iend,j)
       DTRACER_MODULE(istr:iend,j,k,pcaco3soft_ind_t) = (c1-P_CaCO3_gamma)*P_CaCO3_prod - 
     &                                                  P_CaCO3soft_remin(istr:iend,j)
       DTRACER_MODULE(istr:iend,j,k,psio2soft_ind_t) = (c1-P_SiO2_gamma)*P_SiO2_prod - P_SiO2soft_remin(istr:iend,j)
       DTRACER_MODULE(istr:iend,j,k,pironsoft_ind_t) = (c1-P_iron_gamma)*P_iron_prod - P_ironsoft_remin(istr:iend,j)
# endif /* USE_EXPLICIT_VSINK */

       return
       end subroutine ecosys_set_interior


!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE INIT_PARTICULATE_TERMS   ****************
!
!-------------------------------------------------------------------------

        SUBROUTINE init_particulate_terms(QA_dust_def,istr,
     &      iend,j,net_dust_in)

        implicit none
#include "param.h" 
#include "param_bec2.h"
#include "ecosys_bec2.h"
#include "scalars.h"
#include "ocean3d.h"

    !---------------------------------------------------------------------------
    !   Set incoming fluxes (put into outgoing flux for first level usage).
    !   Set dissolution length, production fraction and mass terms.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------
        integer istr, iend, j
        REAL, DIMENSION(istr:iend), INTENT(OUT) :: 
     &          QA_dust_def     ! incoming deficit in the QA(dust) POC flux

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    !DL: the incoming dust flux has the same unit as in the ROMS forcing file.
    ! It is assumed to be kg/m2/s.
    !
        REAL, DIMENSION(istr:iend) :: 
     &         net_dust_in        ! net incoming dust flux

    !-----------------------------------------------------------------------
    !  parameters, from Armstrong et al. 2000
    !
    !  July 2002, length scale for excess POC and bSI modified by temperature
    !  Value given here is at Tref of 30 deg. C, JKM
    ! 
    !   diss       dissolution length for soft subclass
    !   gamma      fraction of production -> hard subclass
    !   mass       mass of 1e6 base units in kg    (WAS: 1e9 base units in g)
    !   rho        QA mass ratio of POC to this particle class
    !
    !   Base units:
    !     POC:        mmol C      (WAS: nmol C)
    !     P_CaCO3:    mmol CaCO3  (WAS: nmol CaCO3)
    !     P_SiO2:     mmol SiO2   (WAS: nmol SiO2)
    !     dust:       kg dust     (WAS: g dust)
    !     P_iron:     mmol Fe     (WAS: nmol Fe)
    !
    !  Units of fluxes:
    !     sflux_in:    incoming flux of soft subclass (base units/m^2/sec)
    !     hflux_in:    incoming flux of hard subclass (base units/m^2/sec)
    !     prod:        production term (base units/m^3/sec)
    !     sflux_out:   outgoing flux of soft subclass (base units/m^2/sec)
    !     hflux_out:   outgoing flux of hard subclass (base units/m^2/sec)
    !     remin:       remineralization term (base units/m^3/sec)
    !    NOTE: Area/volume units were cm^2 and cm^3!
    !-----------------------------------------------------------------------

        POC_diss      = parm_POC_diss   ! diss. length (m), modified by TEMP
        POC_mass      = 12.01           ! molecular weight of POC
        POC_gamma     = c0

        P_CaCO3_diss  = parm_CaCO3_diss ! diss. length (m)
        P_CaCO3_gamma = 0.30            ! prod frac -> hard subclass
        P_CaCO3_mass  = 100.09          ! molecular weight of CaCO3
        P_CaCO3_rho   = 0.05 * P_CaCO3_mass / POC_mass ! QA mass ratio for CaCO3

        P_SiO2_diss   = parm_SiO2_diss  ! diss. length (m), modified by TEMP
        P_SiO2_gamma  = 0.030           ! prod frac -> hard subclass
        P_SiO2_mass   = 60.08           ! molecular weight of SiO2
        P_SiO2_rho    = 0.05 * P_SiO2_mass / POC_mass ! QA mass ratio for SiO2

        dust_diss     = 200.0           ! diss. length (m) (DL: changed from cm)
        dust_gamma    = 0.97            ! prod frac -> hard subclass
        dust_mass     = 1.0e6           ! base units are already kg
        dust_rho      = 0.05 * dust_mass / POC_mass ! QA mass ratio for dust

        P_iron_gamma  = c0              ! prod frac -> hard subclass

    !-----------------------------------------------------------------------
    !  Set incoming fluxes
    !-----------------------------------------------------------------------

        P_CaCO3_sflux_out = c0
        P_CaCO3_hflux_out = c0
        P_SiO2_sflux_out  = c0
        P_SiO2_hflux_out  = c0

        dust_sflux_out(istr:iend,j) = (c1 - dust_gamma) * net_dust_in
        dust_hflux_out(istr:iend,j) = dust_gamma * net_dust_in

        P_iron_sflux_out = c0
        P_iron_hflux_out = c0

    !-----------------------------------------------------------------------
    !  Hard POC is QA flux and soft POC is excess POC.
    !
    !  Note (MF):
    !  These names are convenient given the mineral associated soft and
    !  hard components, but rather confusing when reading about the
    !  particulate/ballast models in Armstrong et al 2002, and Lima et
    !  al 2014 (there is no such thing as POC hard, POC_gamma
    !  doesnt exist here or is set to c0 in CESM BEC). Soft/Hard actually
    !  applies only to the mineral associated components. Think about 
    !  POC_hflux as being the QA component of POC and POC_sflux the excess POC.
    !-----------------------------------------------------------------------

        POC_sflux_out = c0
        POC_hflux_out = c0

    !-----------------------------------------------------------------------
    !  Compute initial QA(dust) POC flux deficit.
    !-----------------------------------------------------------------------

        QA_dust_def = dust_rho * (dust_sflux_out(istr:iend,j) + dust_hflux_out(istr:iend,j))

    !-----------------------------------------------------------------------
    !  Apply dust forcing to surface layer (only if explicitly tracking
    !  dust as a tracer)
    !-----------------------------------------------------------------------
# ifdef USE_EXPLICIT_VSINK
        tracer(istr:iend,j,N,dusthard_ind_t) = tracer(istr:iend,j,N,dusthard_ind_t)
     &                       + dt*dust_gamma*net_dust_in/Hz(istr:iend,j,N)
        tracer(istr:iend,j,N,dustsoft_ind_t) = tracer(istr:iend,j,N,dustsoft_ind_t)
     &                       + dt*(c1-dust_gamma)*net_dust_in/Hz(istr:iend,j,N)
# endif

        return
        end subroutine init_particulate_terms


!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE COMPUTE_PARTICULATE_TERMS   ****************
!
!-------------------------------------------------------------------------
#ifdef IODINE
       subroutine compute_particulate_terms(k,QA_dust_def, TEMP, O2_loc, NO3_loc, 
     &        POC_prod,P_iron_prod,P_CaCO3_prod,P_SiO2_prod,SED_DENITRIF,
     &        OTHER_REMIN,iod_sed_flux1,iod_sed_flux2,iod_sed_flux3,
     &        iod_sed_flux4, istr,iend,j)
#else
       subroutine compute_particulate_terms(k,QA_dust_def, TEMP, O2_loc, NO3_loc, 
     &        POC_prod,P_iron_prod,P_CaCO3_prod,P_SiO2_prod,SED_DENITRIF,
     &        OTHER_REMIN,istr,iend,j)
#endif

       implicit none
#include "param.h"
#include "ecosys_bec2.h"
#include "param_bec2.h"
#include "ocean3d.h"
#include "grid.h"

       integer istr, iend, i, j, k, n_loc, iB
       real, dimension(istr:iend) ::
     &   temp, O2_loc, NO3_loc
     &   ,POC_prod,P_CaCO3_prod,P_SiO2_prod,P_iron_prod
     &   ,QA_dust_def   ! incoming deficit in the QA (dust) POC flux
     &   ,dust_remin    ! remineralization of dust (base units/m^3/sec)
     &   ,SED_DENITRIF  ! sedimentary denitrification (mmolN/m^3/s)
     &   ,OTHER_REMIN   ! sedimentary remin not due to oxic or denitrification
     &   ,flux_oxidated ! carbon oxidation rate in sediment
     &   ,FESEDFLUX     ! sedimentary Fe inputs
#ifdef IODINE
     &   ,iod_sed_flux1     ! sedimentary I input -- param 1
     &   ,iod_sed_flux2     ! sedimentary I input -- param 2
     &   ,iod_sed_flux3     ! sedimentary I input -- param 3
     &   ,iod_sed_flux4     ! sedimentary I input -- param 4
#endif
     &   ,scalelength       
     &   ,TfuncS         ! temperature scaling from soft POM remin (right now just applied to ballast)
#ifdef TDEP_REMIN
     &   ,Tfunc_soft     ! temperature scaling from soft POM remin (Laufkoetter 2017)
#endif
# ifndef USE_EXPLICIT_VSINK
     &   ,DECAY_Hard       ! scaling factor for dissolution of Hard Ballast
     &   ,DECAY_HardDust   ! scaling factor for dissolution of Hard dust
# endif
       real
# ifndef USE_EXPLICIT_VSINK
     &   poc_diss_loc,      ! diss. length used (m)
     &   sio2_diss_loc,     ! diss. length varies spatially with O2 (m)
     &   caco3_diss_loc,
     &   dust_diss_loc,
     &   decay_POC_E,       ! scaling factor for dissolution of excess POC
     &   decay_SiO2,        ! scaling factor for dissolution of SiO2
     &   decay_CaCO3,       ! scaling factor for dissolution of CaCO3
     &   decay_dust,        ! scaling factor for dissolution of dust
# else
     &   pocsoft_remin_loc, pochard_remin_loc,  ! Local remin rates subject to adjustments
     &   pcaco3soft_remin_loc, pcaco3hard_remin_loc,
     &   psio2soft_remin_loc, psio2hard_remin_loc,
     &   pironsoft_remin_loc, pironhard_remin_loc,
     &   dustsoft_remin_loc, dusthard_remin_loc,
# endif
     &   POC_PROD_avail,    ! POC production available for excess POC flux (mmol/m^3/s)
     &   new_QA_dust_def,   ! outgoing deficit in the QA(dust) POC flux
     &   flux, flux_alt,    ! temp variables used to update sinking flux
     &   dz_loc, dzr_loc    ! Hz and its inverse at a particular i,j,k location

       logical
     &   poc_error      ! POC error flag

    !-----------------------------------------------------------------------
    !  incoming fluxes are outgoing fluxes from previous level
    !-----------------------------------------------------------------------

       P_CaCO3_sflux_in(istr:iend,j) = P_CaCO3_sflux_out(istr:iend,j)
       P_CaCO3_hflux_in(istr:iend,j) = P_CaCO3_hflux_out(istr:iend,j)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(caco3fluxin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(caco3fluxin_idx_t)) = P_CaCO3_sflux_in(istr:iend,j)
     &          + P_CaCO3_hflux_in(istr:iend,j)
#  else
       bec2_diag_3d(istr:iend,j,k,caco3fluxin_idx_t) = P_CaCO3_sflux_in(istr:iend,j)
     &          + P_CaCO3_hflux_in(istr:iend,j)
#  endif
# endif

       P_SiO2_sflux_in(istr:iend,j) = P_SiO2_sflux_out(istr:iend,j)
       P_SiO2_hflux_in(istr:iend,j) = P_SiO2_hflux_out(istr:iend,j)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(sio2fluxin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(sio2fluxin_idx_t)) = P_SiO2_sflux_in(istr:iend,j)
     &          + P_SiO2_hflux_in(istr:iend,j)
#  else
       bec2_diag_3d(istr:iend,j,k,sio2fluxin_idx_t) = P_SiO2_sflux_in(istr:iend,j)
     &          + P_SiO2_hflux_in(istr:iend,j)
#  endif
# endif

       dust_sflux_in(istr:iend,j) = dust_sflux_out(istr:iend,j)
       dust_hflux_in(istr:iend,j) = dust_hflux_out(istr:iend,j)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(dustfluxin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(dustfluxin_idx_t)) = dust_sflux_in(istr:iend,j)
     &          + dust_hflux_in(istr:iend,j)
#  else
       bec2_diag_3d(istr:iend,j,k,dustfluxin_idx_t) = dust_sflux_in(istr:iend,j)
     &          + dust_hflux_in(istr:iend,j)
#  endif
# endif

       POC_sflux_in(istr:iend,j) = POC_sflux_out(istr:iend,j)
       POC_hflux_in(istr:iend,j) = POC_hflux_out(istr:iend,j)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(pocfluxin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(pocfluxin_idx_t)) = POC_sflux_in(istr:iend,j)
     &          + POC_hflux_in(istr:iend,j)
#  else
       bec2_diag_3d(istr:iend,j,k,pocfluxin_idx_t) = POC_sflux_in(istr:iend,j)
     &          + POC_hflux_in(istr:iend,j)
#  endif
# endif

       P_iron_sflux_in(istr:iend,j) = P_iron_sflux_out(istr:iend,j)
       P_iron_hflux_in(istr:iend,j) = P_iron_hflux_out(istr:iend,j)
# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(pironfluxin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(pironfluxin_idx_t)) = P_iron_sflux_in(istr:iend,j)
     &          + P_iron_hflux_in(istr:iend,j)
#  else
       bec2_diag_3d(istr:iend,j,k,pironfluxin_idx_t) = P_iron_sflux_in(istr:iend,j)
     &          + P_iron_hflux_in(istr:iend,j)
#  endif
# endif

    !-----------------------------------------------------------------------
    !  initialize loss to sediments = 0 and local copy of percent sed
    !-----------------------------------------------------------------------

       P_iron_sed_loss(istr:iend,j) = c0
       POC_sed_loss(istr:iend,j) = c0
       P_CaCO3_sed_loss(istr:iend,j) = c0
       P_SiO2_sed_loss(istr:iend,j) = c0
       dust_sed_loss(istr:iend,j) = c0
       SED_DENITRIF = c0
       OTHER_REMIN = c0
# ifdef USE_EXPLICIT_VSINK
       FESEDFLUX(istr:iend) = c0
# endif
#ifdef IODINE
       iod_sed_flux1(istr:iend) = c0
       iod_sed_flux2(istr:iend) = c0
       iod_sed_flux3(istr:iend) = c0
       iod_sed_flux4(istr:iend) = c0
#endif
    !-----------------------------------------------------------------------
    !  compute scalelength and decay factors
    !-----------------------------------------------------------------------

       do i = istr, iend
          if (-z_w(i,j,k) < parm_scalelen_z(1)) then
             scalelength(i) = parm_scalelen_vals(1)
          else if (-z_w(i,j,k) >= parm_scalelen_z(size(parm_scalelen_z))) then
             scalelength(i) = parm_scalelen_vals(size(parm_scalelen_z))
          else
             do n_loc = 2, size(parm_scalelen_z)
                if (-z_w(i,j,k) < parm_scalelen_z(n_loc)) then
                   scalelength(i) = parm_scalelen_vals(n_loc-1)
     &                 + (parm_scalelen_vals(n_loc) - parm_scalelen_vals(n_loc-1))
     &                 * (-z_w(i,j,k) - parm_scalelen_z(n_loc-1))/
     &                 (parm_scalelen_z(n_loc) - parm_scalelen_z(n_loc-1))
                   exit
                endif
             end do
          endif
# ifndef USE_EXPLICIT_VSINK      
          DECAY_Hard(i)     = exp(-Hz(i,j,k) / 4.0e4)  ! DL: denominator converted to m
          DECAY_HardDust(i) = exp(-Hz(i,j,k) / 1.2e5)  ! DL: denominator converted to m
# endif      
       end do

    !----------------------------------------------------------------------
    !   Tref = 30.0 reference temperature (deg. C)
    !-----------------------------------------------------------------------

       TfuncS = Q_10**(((TEMP + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)
#ifdef TDEP_REMIN
       Tfunc_soft = exp(parm_ktfunc_soft*(TEMP+T0_Kelvin)) 
#endif

       poc_error = .false.

       do i = istr, iend
          dz_loc = Hz(i,j,k)
         ! dzr_loc = c1/dz_loc ! moved this line out of the if-loop
          if (LANDMASK(i,j)) then
             dzr_loc = c1 / dz_loc

# ifdef USE_EXPLICIT_VSINK             
             ! MF: Create local copies of remineralization lengths
             ! to adjust later on. Only soft components get adjusted
             ! to environmental conditions later! 
             pochard_remin_loc = param_pochard_remin
             pocsoft_remin_loc = param_pocsoft_remin
             psio2hard_remin_loc = param_psio2hard_remin
             psio2soft_remin_loc = param_psio2soft_remin
             pcaco3hard_remin_loc = param_pcaco3hard_remin
             pcaco3soft_remin_loc = param_pcaco3soft_remin
             pironhard_remin_loc = param_pironhard_remin
             pironsoft_remin_loc = param_pironsoft_remin
             dusthard_remin_loc = param_dusthard_remin
             dustsoft_remin_loc = param_dustsoft_remin
# else
             poc_diss_loc = POC_diss
             sio2_diss_loc = P_SiO2_diss
             caco3_diss_loc = P_CaCO3_diss
             dust_diss_loc = dust_diss
# endif

    !-----------------------------------------------------------------------
    !  increase POC diss length scale where O2 concentrations are low
    !-----------------------------------------------------------------------

             if ((O2_loc(i) >= 5.0) .and. (O2_loc(i) < 40.0)) then
# ifdef USE_EXPLICIT_VSINK             
                pocsoft_remin_loc = pocsoft_remin_loc/(c1+(parm_lowo2_remin_factor-c1)*(40.0 - O2_loc(i))/35.0)
# else
                poc_diss_loc = POC_diss*(c1+(parm_lowo2_remin_factor-c1)*(40.0 - O2_loc(i))/35.0)
# endif
             else if (O2_loc(i) < 5.0) then
# ifdef USE_EXPLICIT_VSINK             
                pocsoft_remin_loc = pocsoft_remin_loc / parm_lowo2_remin_factor
# else
                poc_diss_loc = POC_diss * parm_lowo2_remin_factor
# endif                
             endif

    !-----------------------------------------------------------------------
    !  apply scalelength factor to length scales
    !-----------------------------------------------------------------------

# ifdef USE_EXPLICIT_VSINK             
             pocsoft_remin_loc = pocsoft_remin_loc / scalelength(i)
             pcaco3soft_remin_loc = pcaco3soft_remin_loc / scalelength(i)
             psio2soft_remin_loc = psio2soft_remin_loc / scalelength(i)
             pironsoft_remin_loc = pironsoft_remin_loc / scalelength(i)
             dustsoft_remin_loc = dustsoft_remin_loc / scalelength(i)
# elif defined TDEP_REMIN    
             poc_diss_loc = scalelength(i) * poc_diss_loc / Tfunc_soft(i)
             sio2_diss_loc = scalelength(i) * sio2_diss_loc / Tfunc_soft(i)
             caco3_diss_loc = scalelength(i) * caco3_diss_loc / Tfunc_soft(i)
             dust_diss_loc = scalelength(i) * dust_diss_loc / Tfunc_soft(i)
# else
             poc_diss_loc = scalelength(i) * poc_diss_loc 
             sio2_diss_loc = scalelength(i) * sio2_diss_loc 
             caco3_diss_loc = scalelength(i) * caco3_diss_loc 
             dust_diss_loc = scalelength(i) * dust_diss_loc 
# endif

    !-----------------------------------------------------------------------
    !  apply temperature dependence to sio2_diss length scale
    !-----------------------------------------------------------------------
             ! MF: Not sure whether this should still be done? In CESM
             ! BEC, they dont do this anymore, but only apply the
             ! scaling factors, which in the end also represent some
             ! kind of temperature dependence!
# ifdef USE_EXPLICIT_VSINK             
             psio2soft_remin_loc = psio2soft_remin_loc * TfuncS(i)
# else    
             sio2_diss_loc = sio2_diss_loc / TfuncS(i)
# endif
    !-----------------------------------------------------------------------
    !  decay_POC_E and decay_SiO2 set locally, modified by O2
    !-----------------------------------------------------------------------

# ifndef USE_EXPLICIT_VSINK             
             decay_POC_E = exp(-dz_loc / poc_diss_loc)
             decay_SiO2  = exp(-dz_loc / sio2_diss_loc)
             decay_CaCO3 = exp(-dz_loc / caco3_diss_loc)
             decay_dust  = exp(-dz_loc / dust_diss_loc)

    !-----------------------------------------------------------------------
    !  Set outgoing fluxes for non-iron pools.
    !  The outoing fluxes for ballast materials are from the
    !  solution of the coresponding continuous ODE across the model
    !  level. The ODE has a constant source term and linear decay.
    !  It is assumed that there is no sub-surface dust production.
    !-----------------------------------------------------------------------

             P_CaCO3_sflux_out(i,j) = P_CaCO3_sflux_in(i,j) * decay_CaCO3 +
     &           P_CaCO3_prod(i) * ((c1 - P_CaCO3_gamma) * (c1 - decay_CaCO3)
     &             * caco3_diss_loc)

             P_CaCO3_hflux_out(i,j) = P_CaCO3_hflux_in(i,j) * DECAY_Hard(i) +
     &           P_CaCO3_prod(i) * (P_CaCO3_gamma * dz_loc)

             P_SiO2_sflux_out(i,j) = P_SiO2_sflux_in(i,j) * decay_SiO2 +
     &           P_SiO2_prod(i) * ((c1 - P_SiO2_gamma) * (c1 - decay_SiO2)
     &             * sio2_diss_loc)

             P_SiO2_hflux_out(i,j) = P_SiO2_hflux_in(i,j) * DECAY_Hard(i) +
     &           P_SiO2_prod(i) * (P_SiO2_gamma * dz_loc)

             dust_sflux_out(i,j) = dust_sflux_in(i,j) * decay_dust

             dust_hflux_out(i,j) = dust_hflux_in(i,j) * DECAY_HardDust(i)

# else

    !-----------------------------------------------------------------------
    !  Set outgoing fluxes.
    !  In the case of explicit sinking, these are actually computed in
    !  vsink_BEC2.F and stored in VSinkFlux(:,:,N+1,nsink) at w levels.
    !  Incoming fluxes for level n at VSinkFlux(:,:,n,sink_tracer) and
    !  outgoing fluxes for level n at VSinkFlux(:,:,n-1,sink_tracer).
    !  It is assumed that there is no sub-surface dust production.
    !  Downward fluxes are negative in VSinkFlux!
    !-----------------------------------------------------------------------

             P_CaCO3_sflux_out(i,j) = -VSinkFlux(i,j,k-1,iPCACO3SOFT_VSink)
             P_CaCO3_hflux_out(i,j) = -VSinkFlux(i,j,k-1,iPCACO3HARD_VSink)

             P_SiO2_sflux_out(i,j) = -VSinkFlux(i,j,k-1,iPSIO2SOFT_VSink)
             P_SiO2_hflux_out(i,j) = -VSinkFlux(i,j,k-1,iPSIO2HARD_VSink)

             P_iron_sflux_out(i,j) = -VSinkFlux(i,j,k-1,iPIRONSOFT_VSink)
             P_iron_hflux_out(i,j) = -VSinkFlux(i,j,k-1,iPIRONHARD_VSink)

             dust_sflux_out(i,j) = -VSinkFlux(i,j,k-1,iDUSTSOFT_VSink)
             dust_hflux_out(i,j) = -VSinkFlux(i,j,k-1,iDUSTHARD_VSink)

# endif /* USE_EXPLICIT_VSINK */

    !-----------------------------------------------------------------------
    !  Compute how much POC_PROD is available for deficit reduction
    !  and excess POC flux after subtracting off fraction of non-dust
    !  ballast production from net POC_PROD.
    !-----------------------------------------------------------------------
 
             POC_PROD_avail = POC_prod(i) -
     &          P_CaCO3_rho * P_CaCO3_prod(i) -
     &          P_SiO2_rho * P_SiO2_prod(i)

    !-----------------------------------------------------------------------
    !  Check for POC production bounds violations
    !-----------------------------------------------------------------------

             if (POC_PROD_avail < c0) then
                print *, 'subroutine compute_particulate_terms:',
     &            ' non_dust ballast production exceeds POC production'
                print*,'k',k
                print *, 'POC_prod_avail: ', POC_prod_avail
                print*,'POC_prod(i)',POC_prod(i)
                print*,'P_CaCO3_rho*P_CaCO3_prod(i)',P_CaCO3_rho * P_CaCO3_prod(i)
                print*,'P_SiO2_rho * P_SiO2_prod(i)',P_SiO2_rho * P_SiO2_prod(i)
                poc_error = .true.
             endif

    !-----------------------------------------------------------------------
    !  Compute 1st approximation to new QA_dust_def, the QA_dust
    !  deficit leaving the cell (implicit). In the case of explicit 
    !  sinking, new_QA_dust_def equals dust flux leaving the layer.
    !  Ignore POC_PROD_avail at this stage.
    !-----------------------------------------------------------------------

             if (QA_dust_def(i) > 0) then
# ifdef USE_EXPLICIT_VSINK
                new_QA_dust_def = dust_rho * (dust_sflux_out(i,j) +
     &                            dust_hflux_out(i,j))
# else
                new_QA_dust_def = QA_dust_def(i) *
     &             (dust_sflux_out(i,j) + dust_hflux_out(i,j)) / 
     &             (dust_sflux_in(i,j) + dust_hflux_in(i,j))
# endif
             else
                new_QA_dust_def = c0
             endif

    !-----------------------------------------------------------------------
    !  Use POC_PROD_avail to reduce new_QA_dust_def.
    !-----------------------------------------------------------------------

        ! MF: Discuss the whole QA_dust_def treatment for the explicit
        ! case !!!

             if (new_QA_dust_def > c0) then
                new_QA_dust_def = new_QA_dust_def - POC_PROD_avail * dz_loc
                if (new_QA_dust_def < c0) then
                   POC_PROD_avail = -new_QA_dust_def * dzr_loc
                   new_QA_dust_def = c0
                else
                   POC_PROD_avail = c0
                endif
             endif

             QA_dust_def(i) = new_QA_dust_def

    !-----------------------------------------------------------------------
    !  Compute outgoing POC fluxes. QA POC flux is computing using
    !  ballast fluxes and new_QA_dust_def. If no QA POC flux came in
    !  and no production occured, then no QA POC flux goes out. This
    !  shortcut is present to avoid roundoff cancellation errors from
    !  the dust_rho * dust_flux_out - QA_dust_def computation.
    !  Any POC_PROD_avail still remaining goes into excess POC flux.
    !-----------------------------------------------------------------------

        ! HOW DO I DEAL WITH THIS???? I am not sure how this section
        ! would be dealt with in the explicit case. Since
        ! POChard is not actually produced (POC_gamma=0) it is
        ! assembeled from the different rho fractions of ballast
        ! material, so I assume this calculation is 'OK'?!? How
        ! about POCsoft, can I just read it from VSinkFlux, and/or
        ! whats the contribution of POC_PROD_avail???

             if (POC_hflux_in(i,j) == c0 .and. POC_prod(i) == c0) then
                POC_hflux_out(i,j) = c0
             else
                POC_hflux_out(i,j) = P_CaCO3_rho *
     &             (P_CaCO3_sflux_out(i,j) + P_CaCO3_hflux_out(i,j)) +
     &             P_SiO2_rho *
     &             (P_SiO2_sflux_out(i,j) + P_SiO2_hflux_out(i,j)) +
     &             dust_rho *
     &             (dust_sflux_out(i,j) + dust_hflux_out(i,j)) -
     &             new_QA_dust_def
                POC_hflux_out(i,j) = max(POC_hflux_out(i,j), c0)
             endif

# ifdef USE_EXPLICIT_VSINK
             POC_sflux_out(i,j) = -VSinkFlux(i,j,k-1,iPOCSOFT_VSink)
# else             
             POC_sflux_out(i,j) = POC_sflux_in(i,j) * decay_POC_E +
     &          POC_PROD_avail *((c1 - decay_POC_E) * poc_diss)
# endif

    !-----------------------------------------------------------------------
    !  Compute remineralization terms. It is assumed that there is no
    !  sub-surface dust production.
    !-----------------------------------------------------------------------

# ifdef USE_EXPLICIT_VSINK
             P_CaCO3hard_remin(i,j) = pcaco3hard_remin_loc *
     &                      dps * tracer(i,j,k,pcaco3hard_ind_t)
             P_CaCO3soft_remin(i,j) = pcaco3soft_remin_loc *
     &                      dps * tracer(i,j,k,pcaco3soft_ind_t)
             P_CaCO3_remin(i,j) = P_CaCO3hard_remin(i,j)+P_CaCO3soft_remin(i,j)

             POChard_remin(i,j) = pochard_remin_loc *
     &                      dps * tracer(i,j,k,pochard_ind_t)
             POCsoft_remin(i,j) = pocsoft_remin_loc *
     &                      dps * tracer(i,j,k,pocsoft_ind_t)
             POC_remin(i,j) = POChard_remin(i,j)+POCsoft_remin(i,j)

             P_SiO2hard_remin(i,j) = psio2hard_remin_loc *
     &                      dps * tracer(i,j,k,psio2hard_ind_t)
             P_SiO2soft_remin(i,j) = psio2soft_remin_loc *
     &                      dps * tracer(i,j,k,psio2soft_ind_t)
             P_SiO2_remin(i,j) = P_SiO2hard_remin(i,j)+P_SiO2soft_remin(i,j)

             P_ironhard_remin(i,j) = pironhard_remin_loc *
     &                      dps * tracer(i,j,k,pironhard_ind_t)
             P_ironsoft_remin(i,j) = pironsoft_remin_loc *
     &                      dps * tracer(i,j,k,pironsoft_ind_t)
             P_iron_remin(i,j) = P_ironhard_remin(i,j)+P_ironsoft_remin(i,j)

             dusthard_remin(i,j) = dusthard_remin_loc *
     &                      dps * tracer(i,j,k,dusthard_ind_t)
             dustsoft_remin(i,j) = dustsoft_remin_loc *
     &                      dps * tracer(i,j,k,dustsoft_ind_t)
             dust_remin(i) = dusthard_remin(i,j)+dustsoft_remin(i,j)
# else    
             P_CaCO3_remin(i,j) = P_CaCO3_prod(i) +
     &          ((P_CaCO3_sflux_in(i,j) - P_CaCO3_sflux_out(i,j)) +
     &          (P_CaCO3_hflux_in(i,j) - P_CaCO3_hflux_out(i,j))) * dzr_loc

             P_SiO2_remin(i,j) = P_SiO2_prod(i) +
     &          ((P_SiO2_sflux_in(i,j) - P_SiO2_sflux_out(i,j)) +
     &          (P_SiO2_hflux_in(i,j) - P_SiO2_hflux_out(i,j))) * dzr_loc

             POC_remin(i,j) = POC_prod(i) +
     &          ((POC_sflux_in(i,j) - POC_sflux_out(i,j)) +
     &          (POC_hflux_in(i,j) - POC_hflux_out(i,j))) * dzr_loc

             dust_remin(i) =
     &          ((dust_sflux_in(i,j) - dust_sflux_out(i,j)) +
     &          (dust_hflux_in(i,j) - dust_hflux_out(i,j))) * dzr_loc
# endif /* USE_EXPLICIT_VSINK */


    !-----------------------------------------------------------------------
    !  Compute iron remineralization and flux out.
    !-----------------------------------------------------------------------
        
        ! MF: Discuss Fe remin in the case of explicit sinking

# ifdef USE_EXPLICIT_VSINK
             P_iron_remin(i,j) = P_iron_remin(i,j) + POC_remin(i,j)*parm_Red_Fe_C
# else
             if ((POC_sflux_in(i,j)+POC_hflux_in(i,j)) == c0) then
                P_Iron_remin(i,j) = (POC_remin(i,j) * parm_Red_Fe_C)
             else
                P_Iron_remin(i,j) = (POC_remin(i,j) *
     &             (P_iron_sflux_in(i,j) + P_iron_hflux_in(i,j)) /
     &             (POC_sflux_in(i,j) + POC_hflux_in(i,j)))
             endif
             P_Iron_remin(i,j) = P_Iron_remin(i,j)
     &             + P_iron_sflux_in(i,j) * 1.5e-5
             ! MF: What does the hard-coded 1.5e-5 stand for? Applied to
             ! a flux of mmol Fe m-2 s-1? CESM BEC says Fe desorption
             ! from particulate pool. If this is the case, shouldnt it
             ! be treate like the other mineral pools with a dissolution
             ! length and an associated decay?

             P_iron_sflux_out(i,j) = P_iron_sflux_in(i,j) + dz_loc *
     &           ((c1 - P_iron_gamma) * P_iron_prod(i) - P_Iron_remin(i,j))

             ! MF: The following piece does not really pertain to explicit
             ! sinking of particulate Fe?!?
             if (P_iron_sflux_out(i,j) < c0) then
                P_iron_sflux_out(i,j) = c0
                P_Iron_remin(i,j) = P_iron_sflux_in(i,j) * dzr_loc +
     &            (c1 - P_iron_gamma) * P_iron_prod(i)
             endif
# endif 

    !-----------------------------------------------------------------------
    !  Compute iron release from dust remin/dissolution
    !
    !  dust remin gDust = 0.035 / 55.847 * 1.0e9 = 626712.0 nmolFe
    !                      gFe     molFe     nmolFe
    !  Also add in Fe source from sediments if applicable to this cell.
    !-----------------------------------------------------------------------


             P_Iron_remin(i,j) = P_Iron_remin(i,j)
     &           + dust_remin(i) * dust_to_Fe

# ifdef BEC2_DIAG
             bec2_diag_3d(i,j,k,pironremin_idx_t) = P_Iron_remin(i,j)
# endif

             P_iron_hflux_out(i,j) = P_iron_hflux_in(i,j)

          else  ! LANDMASK(i,j)
             P_CaCO3_sflux_out(i,j) = c0
             P_CaCO3_hflux_out(i,j) = c0
             P_CaCO3_remin(i,j) = c0

             P_SiO2_sflux_out(i,j) = c0
             P_SiO2_hflux_out(i,j) = c0
             P_SiO2_remin(i,j) = c0

             dust_sflux_out(i,j) = c0
             dust_hflux_out(i,j) = c0
             dust_remin(i) = c0

             POC_sflux_out(i,j) = c0
             POC_hflux_out(i,j) = c0
             POC_remin(i,j) = c0

             P_iron_sflux_out(i,j) = c0
             P_iron_hflux_out(i,j) = c0
             P_Iron_remin(i,j) = c0
          endif  ! LANDMASK(i,j)

    !-----------------------------------------------------------------------
    !  Bottom Sediments Cell?
    !  If so compute sedimentary burial and denitrification N losses.
    !  Using empirical relations from Bohlen et al., 2012 (doi:10.1029/2011GB004198) for Sed Denitrification
    !  OTHER_REMIN estimates organic matter remineralized in the sediments
    !      by the processes other than oxic remin and denitrification (SO4 and CO2,
    !      etc..)
    !      based on Soetaert et al., 1996, varies between 10% and 50%
    !      0.4_r8 is a coefficient with units mmolC/cm2/yr sinking flux,
    !      OTHER_REMIN is 50% above this high flux value,
    !      In special case where bottom O2 has been depleted to < 1.0 uM,
    !               all sedimentary remin is due to DENITRIFICATION + OTHER_REMIN
    !  POC burial from Dunne et al. 2007 (doi:10.1029/2006GB002907), maximum of 80% burial efficiency imposed
    !  Bsi preservation in sediments based on
    !     Ragueneau et al. 2000 (doi:10.1016/S0921-8181(00)00052-7)
    !  Calcite is preserved in sediments above the lysocline, dissolves below.
    !       Here a constant depth is used for lysocline.
    !-----------------------------------------------------------------------

          if (LANDMASK(i,j) .and. (k == 1)) then

             flux = POC_sflux_out(i,j)+POC_hflux_out(i,j) ! mmol C/m^2/s

             if (flux > c0) then
               flux_alt = flux*spd ! convert to mmol C/m^2/day

               POC_sed_loss(i,j) = flux * min(0.8, parm_POMbury
     &             * (0.013 + 0.53 * flux_alt*flux_alt / (7.0 + flux_alt)**2))

            ! CN: set denitr to zero if NO3<=0 in bottom water (to not
            ! further deplete already not existing NO3)
            ! if this works: actually, this might still take out too
            ! much NO3 if NO3_BW is just slightly above zero. Include
            ! something like "if SED_DENITRIF > NO3_loc,
            ! reduce SED_DENITRIF to available NO3"
            ! check: is some NO3 threshold implemented for water col
            ! denitrif? could also cause neg NO3...
            ! FOR NOW: NO3min=5 (see Simon's draft!! he uses 2, but said
            ! it's kind of arbitrary...), below this
            ! sed_denitr is reduced! same for water column denitrf!

            ! CN (Oct 2015): constrain SED_DENITRIF by what is left to
            ! be remineralized after POC_sed_loss

            ! SY: Changed to Michaelis Menton limitation by NO3
               SED_DENITRIF(i) = dzr_loc * min(flux * (0.06 + 0.19 *
     &            0.99**(O2_loc(i)-NO3_loc(i))),
     &            flux-POC_sed_loss(i,j)*denitrif_C_N) * 
     &            (NO3_loc(i)/(NO3_loc(i)+parm_sed_denitrif_NO3_limit))

#ifdef BEC2_DIAG
# ifdef BEC2_DIAG_USER
               if (bec2_diag_2d_l(seddenitrif_idx_t))
     &            bec2_diag_2d(i,j,bec2_diag_2d_idx(seddenitrif_idx_t)) = SED_DENITRIF(i)
# else
               bec2_diag_2d(i,j,seddenitrif_idx_t) = SED_DENITRIF(i)
               bec2_diag_2d(i,j,fluxtosed_idx_t) = flux
               bec2_diag_2d(i,j,pocsedloss_idx_t) = POC_sed_loss(i,j)
# endif
#endif

               flux_alt = flux*1e-4*spd*365.0 ! convert to mmol C/cm^2/year
               OTHER_REMIN(i) = dzr_loc
     &             * min(min(0.1 + flux_alt,0.5) * (flux - POC_sed_loss(i,j)),
     &                   (flux - POC_sed_loss(i,j) - (SED_DENITRIF(i)*dz_loc*denitrif_C_N)))

    !----------------------------------------------------------------------------------
    !              if bottom water O2 is depleted, assume all remin is denitrif + other               
    !----------------------------------------------------------------------------------

               if (O2_loc(i) < c1) then
                 OTHER_REMIN(i) = dzr_loc *
     &                (flux - POC_sed_loss(i,j) - (SED_DENITRIF(i)*dz_loc*denitrif_C_N))
              endif

#ifdef BEC2_DIAG
               bec2_diag_2d(i,j,otherremin_idx_t) = OTHER_REMIN(i)
#endif

# ifdef CONSTANT_FE_SEDFLUX
             ! Add diffusive iron flux if depth < 1100.0m, based on
             ! Johnson et al.1999, value of 5.0 umolFe/m2/day.
             !
             ! 2.0 umolFe/m2/day,  2.3148e-6  nmolFe/cm2/sec
             ! 2.0 umolFe/m2/day,  2.3148e-3  mmolFe/m2/sec
             !-----------------------------------------------------------------
             ! CN: take value from Moore 2004 (in mmolFe/m2/sec,see old BEC code)
              if (z_r(i,j,k) < 1100.) then
                FESEDFLUX(i) = 2.0*0.001/86400
              endif
#  else
            ! CN (Sep 2015): FESEDFLUX parametrization from Dale 2015
            ! carbon oxidation rate in mmol m-2 s-1 
            flux_oxidated(i) = flux - POC_sed_loss(i,j)

            ! convert from mmol m-2 s-1 to mmol m-2 d-1
            flux_oxidated(i) = flux_oxidated(i) * spd  ! day2sec=86400

            ! constant 170 in umol m-2 d-1, flux_oxidated in mmol m-2
            ! day-1, o2 in muM
            if (flux_oxidated(i) == 0) then
                FESEDFLUX(i) = 0 
            elseif (flux_oxidated(i) > 0) then
                FESEDFLUX(i) = 170.0 * tanh(flux_oxidated(i)/O2_loc(i))
            elseif (flux_oxidated(i) > 0 .and. O2_loc(i)==0) then
                FESEDFLUX(i) = 170.0 

            endif
            ! convert from umol m-2 d-1 to mmol m-2 s-1
            if (latr(i,j) > 5.0 .and. latr(i,j) < 20.0) then
                FESEDFLUX(i) = FESEDFLUX(i)*0.001/86400*1.0     !sec2day=1/86400
            else
                FESEDFLUX(i) = FESEDFLUX(i)*0.001/86400*1.0
            endif
# endif

# ifdef IODINE
            ! -----------------------------------
            ! Excess Iodine sedimentary emissions
            ! -----------------------------------
            
            ! Parameters (defined in ecosys_bec2_init.F)
            ! -----------------------------------
            !     parm_iod1_cstflux = 5.0 * yps ! Constant excess Iodine flux in mmol m-2 -s
            !     parm_iod2_fz2 = 5.0 * yps  ! Iodine flux at z=z2 in mmol m-2 -s
            !     parm_iod2_z2 = 100.0 ! Reference depth for Iodine flux calculation
            !     parm_iod2_z1 = 50.0 ! Depth at which Iodine flux is maximum
            !     parm_iod2_b = -0.52 ! power law exponent
            !     parm_iod3_ric_ox = 10.0e-4 ! ratio of excess I produced to C oxidized in the sediments
            !     parm_iod4_maxflux = 10.0 * yps  ! Max Iodine flux in mmol m-2 -s
            ! -----------------------------------


            ! Parameterization 1: constant source
            ! -----------------------------------
            ! *assume constant emmision rate of 5 mmol m-2 y-1
            iod_sed_flux1(i) = parm_iod1_cstflux * dzr_loc

            ! Parameterization 2: depth dependant power law 
            ! -----------------------------------
            ! *assume emmision rate of 5 mmol m-2 y-1
            ! at z = 100 m. Assume saturation at 50m.
            if (-z_w(i,j,k) .le. parm_iod2_z1) then 
                iod_sed_flux2(i) = parm_iod2_fz2 * dzr_loc*(parm_iod2_z1
     &          /parm_iod2_z2) ** parm_iod2_b
            elseif (-z_w(i,j,k) > parm_iod2_z1) then 
                iod_sed_flux2(i) = parm_iod2_fz2*dzr_loc* (-z_w(i,j,k)
     &          /parm_iod2_z2) ** parm_iod2_b
            endif

            ! Parameterization 3: proportional to oxidated C flux 
            ! -----------------------------------
            ! units of mmol m-2 s-1
            iod_sed_flux3(i) = parm_iod3_ric_ox * dzr_loc * 
     &                         (flux -POC_sed_loss(i,j))


            ! Parameterization 4: f(oxidated C flux, O2) 
            ! -----------------------------------
            ! Following Fe sed flux parameterization:
            !     *flux_oxidated should be in mmol m-2 d-1
            !     and o2 in mmol m-3
            ! Assume max emmision rate is 10 mmol m-2 y-1

            ! First convert oxidated flux from mmol m-2 s-1 to mmol m-2 d-1
            flux_oxidated(i) = (flux -POC_sed_loss(i,j))*spd

            if (flux_oxidated(i) == c0) then
                iod_sed_flux4(i) = c0  
            elseif (flux_oxidated(i) > c0 .and. O2_loc(i)>c0) then
                iod_sed_flux4(i) = parm_iod4_maxflux * tanh(flux_oxidated(i)/O2_loc(i))
            elseif (flux_oxidated(i) > c0 .and. O2_loc(i)==c0) then 
                iod_sed_flux4(i) = parm_iod4_maxflux
            endif
# endif

# ifdef BEC2_DIAG
             bec2_diag_2d(i,j,fesedflux_idx_t) = FESEDFLUX(i)
# endif
            endif  ! flux > c0

            flux = P_SiO2_sflux_out(i,j)+P_SiO2_hflux_out(i,j)
            flux_alt = flux*spd ! convert to mmol/m^2/day
            ! first compute burial efficiency, then compute loss to sediments
            if (flux_alt > c1) then
               P_SiO2_sed_loss(i,j) = 0.3 * flux_alt - 0.06
            else
               P_SiO2_sed_loss(i,j) = 0.04
            endif
            P_SiO2_sed_loss(i,j) = flux * parm_BSIbury * P_SiO2_sed_loss(i,j)

# ifdef BEC2_DIAG
            bec2_diag_2d(i,j,sio2fluxtosed_idx_t) = flux
            bec2_diag_2d(i,j,sio2sedloss_idx_t) = P_SiO2_sed_loss(i,j)
# endif

            flux = c0  ! set flux back to zero
            if (-z_w(i,j,k) < 3300.0) then
               flux = P_CaCO3_sflux_out(i,j) + P_CaCO3_hflux_out(i,j)
               P_CaCO3_sed_loss(i,j) = flux
            endif

# ifdef BEC2_DIAG
            bec2_diag_2d(i,j,caco3fluxtosed_idx_t) = flux
            bec2_diag_2d(i,j,caco3sedloss_idx_t) = P_CaCO3_sed_loss(i,j)
# endif

    !----------------------------------------------------------------------------------
    !  Update sinking fluxes and remin fluxes, accounting for sediments.
    !  flux used to hold sinking fluxes before update.
    !----------------------------------------------------------------------------------

            flux = P_CaCO3_sflux_out(i,j) + P_CaCO3_hflux_out(i,j)
            if (flux > c0) then
               P_CaCO3_remin(i,j) = P_CaCO3_remin(i,j)
     &             + ((flux - P_CaCO3_sed_loss(i,j)) * dzr_loc)
            endif

            flux = P_SiO2_sflux_out(i,j) + P_SiO2_hflux_out(i,j)
            if (flux > c0) then
               P_SiO2_remin(i,j) = P_SiO2_remin(i,j)
     &             + ((flux - P_SiO2_sed_loss(i,j)) * dzr_loc)
            endif

            flux = POC_sflux_out(i,j) + POC_hflux_out(i,j)
            if (flux > c0) then
               POC_remin(i,j) = POC_remin(i,j)
     &             + ((flux - POC_sed_loss(i,j)) * dzr_loc)
            endif

    !-----------------------------------------------------------------------
    !   Remove all Piron and dust that hits bottom, sedimentary Fe source 
    !        accounted for by FESEDFLUX elsewhere.
    !-----------------------------------------------------------------------


            flux = c0  ! set flux back to zero
            flux = (P_iron_sflux_out(i,j) + P_iron_hflux_out(i,j))
            if (flux > c0) then
               P_iron_sed_loss(i,j) = flux
            endif

            dust_sed_loss(i,j) = dust_sflux_out(i,j) + dust_hflux_out(i,j)

# ifdef BEC2_DIAG
             bec2_diag_2d(i,j,pironfluxtosed_idx_t) = flux
             bec2_diag_2d(i,j,dustfluxtosed_idx_t) = dust_sed_loss(i,j)
# endif

    !-----------------------------------------------------------------------
    !   Bottom layer: set all outgoing fluxes to 0.0
    !-----------------------------------------------------------------------

            if (k == 1) then
               P_CaCO3_sflux_out(i,j) = c0
               P_CaCO3_hflux_out(i,j) = c0

               P_SiO2_sflux_out(i,j) = c0
               P_SiO2_hflux_out(i,j) = c0

               dust_sflux_out(i,j) = c0
               dust_hflux_out(i,j) = c0

               POC_sflux_out(i,j) = c0
               POC_hflux_out(i,j) = c0

               P_iron_sflux_out(i,j) = c0
               P_iron_hflux_out(i,j) = c0
            endif

             P_Iron_remin(i,j) = P_Iron_remin(i,j)
     &           + (FESEDFLUX(i) * dzr_loc)
         end if  ! LANDMASK(i,j) .and. (k == 1)

       end do  ! i = istr, iend

# ifdef BEC2_DIAG
#  ifdef BEC2_DIAG_USER
       if (bec2_diag_3d_l(pironremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(pironremin_idx_t)) = P_Iron_remin(istr:iend,j)
       if (bec2_diag_3d_l(caco3remin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(caco3remin_idx_t)) = P_CaCO3_remin(istr:iend,j)
       if (bec2_diag_3d_l(sio2remin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(sio2remin_idx_t)) = P_SiO2_remin(istr:iend,j)
       if (bec2_diag_3d_l(pocremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(pocremin_idx_t)) = POC_remin(istr:iend,j)
       if (bec2_diag_3d_l(dustremin_idx_t))
     &    bec2_diag_3d(istr:iend,j,k,bec2_diag_3d_idx(dustremin_idx_t)) = dust_remin(istr:iend)
#  else
       bec2_diag_3d(istr:iend,j,k,pironremin_idx_t) = P_Iron_remin(istr:iend,j)
       bec2_diag_3d(istr:iend,j,k,caco3remin_idx_t) = P_CaCO3_remin(istr:iend,j)
       bec2_diag_3d(istr:iend,j,k,sio2remin_idx_t) = P_SiO2_remin(istr:iend,j)
       bec2_diag_3d(istr:iend,j,k,pocremin_idx_t) = POC_remin(istr:iend,j)
       bec2_diag_3d(istr:iend,j,k,dustremin_idx_t) = dust_remin(istr:iend)
#  endif
# endif

       return
       end subroutine compute_particulate_terms


!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE WS            ****************************
!
!-------------------------------------------------------------------------

       subroutine WS(SMFTX, SMFTY,landmask,work,istr,iend,jstr,jend)
    !-------------------------------------------------------------------------
    !   result declaration
    !-------------------------------------------------------------------------
        implicit none

! needed for scalars.h
#include "param.h"  
! needed for rho0
#include "scalars.h"       
        integer::istr,iend,jstr,jend
        REAL, DIMENSION(istr:iend,jstr:jend), INTENT(OUT)  :: WORK ! 10m windspeed (m/s)

    !-------------------------------------------------------------------------
    !   Compute 10m windspeed from surface velocity fluxes where LAND_MASK is
    !   true. Give zero where LAND_MASK is false.
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !   Derivation
    !
    !   stress (N/m^2) is a function of surface velocity fluxes
    !   (taux,tauy) = (SMFTX,SMFTY) / momentum_factor
    !
    !   turbulent velocity scale (m/s) is a function of stress
    !   ustar^2 = sqrt(taux^2 + tauy^2) / rho_air
    !           = sqrt(SMFTX^2 + SMFTY^2) / (momentum_factor rho_air)
    !
    !   turbulent velocity scale (m/s) is a function of 10m windspeed (m/s)
    !   ustar^2 = coef_1 ws + coef_2 WS^2 + coef_3 WS^3
    !   ustar^2 = c_d WS^2
    !
    !   Compute ustar^2 from velocity fluxes and solve cubic equation
    !   for 10m windspeed. Initial guess for iterative solver comes
    !   from assuming a constant drag coefficient c_d.
    !
    !-----------------------------------------------------------------------
        real,parameter::
     &    rho_air   = 1.2            ! ambient air density (kg/m^3)
!HF not needed, use rho0     &   momentum_factor  = 10.0
!HF not needed     &   cmperm        = 100.        ! cm per meter

    !-----------------------------------------------------------------------
    !   arguments
    !-------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: 
     &    SMFTX,  ! 'zonal' surface velocity flux (m^2/s^2)
     &    SMFTY    ! 'meridional' surface velocity flux (m^2/s^2)
! HF: both were in (cm^2/s^s)

    !------------------------------------------------------------------------
    !   local variables
    !-------------------------------------------------------------------------

        REAL, PARAMETER :: 
     &    coef_1  = 0.0027,     ! (m/s)
     &    coef_2  = 0.000142,   ! (non-dimensional)
     &    coef_3  = 0.0000764,  ! (s/m)
     &    c_d     = 1.7e-3      ! drag coefficient for initial u

        REAL, DIMENSION(istr:iend,jstr:jend) :: 
     &    ustar_squared  ! square of turbulent velocity scale (m^2/s^2)
        LOGICAL,DIMENSION(istr:iend,jstr:jend) :: landmask
    !------------------------------------------------------------------------
    !   solve for WS in (m/s)
    !-------------------------------------------------------------------------

       WHERE (landmask)
         ustar_squared = SQRT(SMFTX**2 + SMFTY**2) * rho0 / 
     &         rho_air

         WORK = SQRT(ustar_squared / c_d)

       !----------------------------------------------------------------------
       !   perform 3 Newton iterations
       !   loop constructs are not allowed inside WHERE constructs
       !----------------------------------------------------------------------

          WORK = WORK - (WORK*(coef_1 + WORK*(coef_2 + WORK*coef_3)) 
     &         - ustar_squared) / 
     &         (coef_1 + WORK*(2*coef_2 + WORK*3*coef_3))

          WORK = WORK - (WORK*(coef_1 + WORK*(coef_2 + WORK*coef_3)) 
     &          - ustar_squared) / 
     &         (coef_1 + WORK*(2*coef_2 + WORK*3*coef_3))

          WORK = WORK - (WORK*(coef_1 + WORK*(coef_2 + WORK*coef_3))
     &           - ustar_squared) / 
     &       (coef_1 + WORK*(2*coef_2 + WORK*3*coef_3))

        ELSEWHERE
           WORK = 0.0
        END WHERE

    !-------------------------------------------------------------------------
    !   HF: DO NOT convert to (cm/s)
    !-------------------------------------------------------------------------

        RETURN
        END subroutine WS

# if defined Ncycle_SY || defined N2O_NEV
!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE CSCHMIDT_N2O   ****************************
!
!-------------------------------------------------------------------------

         subroutine CSCHMIDT_N2O(SSTT,landmask,SCHMIDT_N2O,istr,
     &          iend,jstr,jend)
        implicit none
!! needed for c0
!#include "param_bec2.h"
        
       integer::istr,iend,jstr,jend
         LOGICAL,DIMENSION(istr:iend,jstr:jend) :: landmask

     !---------------------------------------------------------------------------
     !   Compute Schmidt number of N2O in seawater as function of SST
     !   where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   ref : Sarmiento and Gruber 2006 book (OBCD, table 3.3.1 page 85 
     !---------------------------------------------------------------------------
     !
     !DL: over land, do NOT assign zero to the result, but use very large value instead
     ! in order to prevent division by zero.

     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

           REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: SSTT

     !---------------------------------------------------------------------------
     !   result declaration
     !---------------------------------------------------------------------------

           REAL, DIMENSION(istr:iend,jstr:jend) :: SCHMIDT_N2O

     !---------------------------------------------------------------------------
     !   coefficients in expansion
     !---------------------------------------------------------------------------

          REAL, PARAMETER :: 
     &      a = 2301.1, 
     &      b = 151.1, 
     &      c = 4.7364, 
     &      d = 0.059431

         WHERE (LANDMASK)
           SCHMIDT_N2O = a + SSTT* (-b + SSTT * (c + SSTT * (-d)))
         ELSEWHERE
            SCHMIDT_N2O = 1e30
          END WHERE
          return
          END subroutine CSCHMIDT_N2O

          subroutine N2OSATU(SSTT, SSSS,landmask, N2OSAT,istr,
     &         iend,jstr,jend)
        implicit none
! needed for c0
#include "param_bec2.h"
       integer::istr,iend,jstr,jend
       LOGICAL,DIMENSION(istr:iend,jstr:jend) :: landmask

     !---------------------------------------------------------------------------
     !
     !   Computes n2o saturation concentration at 1 atm n2o pressure
     !   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
     !   in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   FROM WEISS AND PRICE (1980), Marine chemistry.
     !
     !---------------------------------------------------------------------------

     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: 
     &    SSTT,  ! sea surface temperature (C)
     &    SSSS    ! sea surface salinity (psu)

     !---------------------------------------------------------------------------
     !   result declaration
     !---------------------------------------------------------------------------

         REAL, DIMENSION(istr:iend,jstr:jend) :: N2OSAT

     !---------------------------------------------------------------------------
     !   local variables
     !---------------------------------------------------------------------------

         REAL, DIMENSION(istr:iend,jstr:jend) :: TS

     !---------------------------------------------------------------------------
     !   coefficients in expansion
     !---------------------------------------------------------------------------

        REAL, PARAMETER :: 
     &    a_1 = -165.8802, 
     &    a_2 = 222.8743, 
     &    a_3 = 92.0792, 
     &    a_4 = -1.48425, 
     &    b_1 = -0.056235,
     &    b_2 = 0.031619, 
     &    b_3 = -0.0048472 

        WHERE (LANDMASK)
           TS = T0_Kelvin + SSTT

         N2OSAT =  exp(a_1 + a_2 * (100.0/TS) + a_3 * log(TS/100.0) +
     &          a_4 * (TS/100.0)**2 + SSSS * (b_1 +b_2 * (TS/100.0) + b_3 *
     &          (TS/100.0)**2)) 
        ELSEWHERE
          N2OSAT = c0
       END WHERE
     !---------------------------------------------------------------------------
     !   Convert from mol/l to mmol/m^3
     !---------------------------------------------------------------------------

        N2OSAT = N2OSAT * 1000 * 1000

        return
        END subroutine N2OSATU

     !---------------------------------------------------------------------------
     !---------------------------------------------------------------------------
     !   CSCHMIDT_N2
     !---------------------------------------------------------------------------
     !---------------------------------------------------------------------------


         subroutine CSCHMIDT_N2(SSTT,landmask,SCHMIDT_N2,istr,
     &          iend,jstr,jend)
        implicit none
        
       integer::istr,iend,jstr,jend
         LOGICAL,DIMENSION(istr:iend,jstr:jend) :: landmask

    !---------------------------------------------------------------------------
    !   Compute Schmidt number of N2 in seawater as function of SST
    !   where LAND_MASK is true. Give zero where LAND_MASK is false.
    !
    !   ref : Sarmiento and Gruber 2006 (Table 3.3.1)
    !         
    !---------------------------------------------------------------------------
    !
    !DL: over land, do NOT assign zero to the result, but use very large value instead
    ! in order to prevent division by zero.

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: SSTT

    !---------------------------------------------------------------------------
    !   result declaration
    !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend) :: SCHMIDT_N2

    !---------------------------------------------------------------------------
    !   coefficients in expansion
    !---------------------------------------------------------------------------

          REAL, PARAMETER :: 
     &      a = 2206.1,   
     &      b = 144.86, 
     &      c = 4.5413, 
     &      d = 0.056988

         WHERE (LANDMASK)
           SCHMIDT_N2 = a + SSTT* (-b + SSTT * (c + SSTT * (-d)))
         ELSEWHERE
            SCHMIDT_N2 = 1e30
          END WHERE
          return
          END subroutine CSCHMIDT_N2

  !*****************************************************************************


          subroutine N2SATU(SSTT, SSSS, landmask, N2SAT,istr,
     &         iend,jstr,jend)
        implicit none
! needed for c0
#include "param_bec2.h"
       integer::istr,iend,jstr,jend
       LOGICAL,DIMENSION(istr:iend,jstr:jend) :: landmask

     !---------------------------------------------------------------------------
     !
     !   Computes dinitrogen saturation concentration at 1 atm total pressure
     !   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
     !   in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   FROM Sarmiento and Gruber 2006 
     !   --> Solubility : table 3.2.2
     !   --> apply correction for moist air following Panel 3.2.1
     !
     !   N2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
     !   Simon Yang, Feb 2018
     !---------------------------------------------------------------------------
 
     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

       REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: 
     &    SSTT,  ! sea surface temperature (C)
     &    SSSS    ! sea surface salinity (psu)

    !---------------------------------------------------------------------------
    !   result declaration
    !---------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend) :: N2SAT

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend) :: TS, BUNSEN, PH2O_to_P

    !---------------------------------------------------------------------------
    !   coefficients in expansion
    !---------------------------------------------------------------------------

        REAL, PARAMETER :: 
     &    a_1 = -59.6274, 
     &    a_2 = 85.7761, 
     &    a_3 = 24.3696, 
     &    b_1 = -0.051580, 
     &    b_2 = 0.026329, 
     &    b_3 = -0.0037252,
     &    vbarn2o = 22.4136 ! volar molume of N2 at standard temperature and pressure

        WHERE (LANDMASK)
           TS = T0_Kelvin + SSTT

          BUNSEN=  exp(a_1 + a_2 * (100.0/TS) + a_3 * log(TS/100.0) +
     &          SSSS * (b_1 +b_2 * (TS/100.0) + b_3 *
     &          (TS/100.0)**2))
        
!!!! Calculate Ph2o/P to correct for moist air 
        PH2O_to_P = exp(24.4543-67.4509*(100.0/TS)-4.8489*log(TS/100.0)
     &              -0.000544*SSSS)

        N2SAT = ( BUNSEN / vbarn2o ) * (c1 - PH2O_to_P)

        ELSEWHERE
          N2SAT = c0
       END WHERE
    !---------------------------------------------------------------------------
    !   Convert from ml/l to mmol/m^3
    !---------------------------------------------------------------------------

        N2SAT = N2SAT * 1e6

        return
        END subroutine N2SATU

  !*****************************************************************************
# endif
     !---------------------------------------------------------------------------
     !   Subroutine CSCHMIDT_O2
     !---------------------------------------------------------------------------

         subroutine CSCHMIDT_O2(SSTT,landmask,SCHMIDT_O2,istr,
     &          iend,jstr,jend)
        implicit none
!! needed for c0
!#include "param_bec2.h"
        
       integer::istr,iend,jstr,jend
         LOGICAL,DIMENSION(istr:iend,jstr:jend) :: landmask

    !---------------------------------------------------------------------------
    !   Compute Schmidt number of O2 in seawater as function of SST
    !   where LAND_MASK is true. Give zero where LAND_MASK is false.
    !
    !   ref : Keeling et al, Global Biogeochem. Cycles, Vol. 12,
    !         No. 1, pp. 141-163, March 1998
    !---------------------------------------------------------------------------
    !
    !DL: over land, do NOT assign zero to the result, but use very large value instead
    ! in order to prevent division by zero.

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: SSTT

    !---------------------------------------------------------------------------
    !   result declaration
    !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend) :: SCHMIDT_O2

    !---------------------------------------------------------------------------
    !   coefficients in expansion
    !---------------------------------------------------------------------------

          REAL, PARAMETER :: 
     &      a = 1638.0, 
     &      b = 81.83, 
     &      c = 1.483, 
     &      d = 0.008004

         WHERE (LANDMASK)
           SCHMIDT_O2 = a + SSTT* (-b + SSTT * (c + SSTT * (-d)))
         ELSEWHERE
            SCHMIDT_O2 = 1e30
          END WHERE
          return
          END subroutine CSCHMIDT_O2

  !*****************************************************************************

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE O2SATU   ****************************
!
!-------------------------------------------------------------------------

          subroutine O2SATU(SSTT, SSSS,landmask, O2SAT,istr,
     &         iend,jstr,jend)
        implicit none
! needed for c0
#include "param_bec2.h"
       integer::istr,iend,jstr,jend
       LOGICAL,DIMENSION(istr:iend,jstr:jend) :: landmask

     !---------------------------------------------------------------------------
     !
     !   Computes oxygen saturation concentration at 1 atm total pressure
     !   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
     !   in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
     !
     !   FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
     !   THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
     !
     !   *** NOTE: THE "A_3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
     !   *** IT SHOULDN'T BE THERE.                                ***
     !
     !   O2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
     !   0 permil <= S <= 42 permil
     !   CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil,
     !   O2SAT = 282.015 mmol/m^3
     !
     !---------------------------------------------------------------------------
 
     !---------------------------------------------------------------------------
     !   arguments
     !---------------------------------------------------------------------------

       REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: 
     &    SSTT,  ! sea surface temperature (C)
     &    SSSS    ! sea surface salinity (psu)

    !---------------------------------------------------------------------------
    !   result declaration
    !---------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend) :: O2SAT

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend,jstr:jend) :: TS

    !---------------------------------------------------------------------------
    !   coefficients in expansion
    !---------------------------------------------------------------------------

        REAL, PARAMETER :: 
     &    a_0 = 2.00907, 
     &    a_1 = 3.22014, 
     &    a_2 = 4.05010, 
     &    a_3 = 4.94457, 
     &    a_4 = -2.56847E-1, 
     &    a_5 = 3.88767, 
     &    b_0 = -6.24523E-3, 
     &    b_1 = -7.37614E-3, 
     &    b_2 = -1.03410E-2, 
     &    b_3 = -8.17083E-3, 
     &    c_0 = -4.88682E-7

        WHERE (LANDMASK)
           TS = LOG( ((T0_Kelvin+25.0) - SSTT) / (T0_Kelvin + SSTT) )

         O2SAT = EXP(a_0+TS*(a_1+TS*(a_2+TS*(a_3+TS*(a_4+TS*a_5)))) + 
     &         SSSS*( (b_0+TS*(b_1+TS*(b_2+TS*b_3))) + SSSS*c_0 ))
        ELSEWHERE
          O2SAT = c0
       END WHERE
    !---------------------------------------------------------------------------
    !   Convert from ml/l to mmol/m^3
    !---------------------------------------------------------------------------

        O2SAT = O2SAT * 44.6596

        return
        END subroutine O2SATU

  !*****************************************************************************

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE CSCHMIDT_CO2   ****************************
!
!-------------------------------------------------------------------------

          subroutine CSCHMIDT_CO2(SSTT,landmask,SCHMIDT_CO2,istr,
     &                     iend,jstr,jend )
        implicit none
! needed for c0
#include "param_bec2.h"
       integer::istr,iend,jstr,jend
         LOGICAL,DIMENSION(istr:iend,jstr:jend) :: landmask

    !---------------------------------------------------------------------------
    !   Compute Schmidt number of CO2 in seawater as function of SST
    !   where LAND_MASK is true. Give zero where LAND_MASK is false.
    !
    !   ref : Wanninkhof, J. Geophys. Res, Vol. 97, No. C5,
    !   pp. 7373-7382, May 15, 1992
    !---------------------------------------------------------------------------
    !DL: over land, do NOT assign zero to the result, but use very large value instead
    ! in order to prevent division by zero.

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend), INTENT(IN) :: SSTT

    !---------------------------------------------------------------------------
    !   result declaration
    !---------------------------------------------------------------------------

          REAL, DIMENSION(istr:iend,jstr:jend), INTENT(OUT) :: SCHMIDT_CO2

    !---------------------------------------------------------------------------
    !   coefficients in expansion
    !---------------------------------------------------------------------------

         REAL, PARAMETER :: 
     &    a = 2073.1, 
     &    b = 125.62, 
     &    c = 3.6276,
     &    d = 0.043219

         WHERE (LANDMASK)
           SCHMIDT_CO2 = a + SSTT * (-b + SSTT * (c + SSTT * (-d)))
         ELSEWHERE
           SCHMIDT_CO2 = 1e30
         END WHERE
         return
        END subroutine  CSCHMIDT_CO2

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE CO2CALC_ROW   ***************************
!-------------------------------------------------------------------------

       SUBROUTINE co2calc_row(mask, t, s, dic_in, ta_in, pt_in,
     &    sit_in, phlo, phhi, ph, xco2_in, atmpres, co2star,
     &    dco2star, pco2oc, dpco2, istr,iend)


    !---------------------------------------------------------------------------
    !   SUBROUTINE CO2CALC
    !
    !   PURPOSE : Calculate delta co2*, etc. from total alkalinity, total CO2,
    !             temp, salinity (s), etc.
    !---------------------------------------------------------------------------
        implicit none
       real,parameter::c0=0.0, c1=1.0, c10=10.0, 
     &   c1000=1000.0, T0_Kelvin=273.16,rho_sw=4.1/3.996
       integer::istr,iend
    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------
       REAL, PARAMETER :: xacc = 1e-10
       INTEGER, PARAMETER :: maxit = 100

  !-----------------------------------------------------------------------------
  !   declarations for function coefficients & species concentrations
  !-----------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend) :: 
     &  k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff, 
     &  bt, st, ft, dic, ta, pt, sit

         LOGICAL, DIMENSION(istr:iend), INTENT(IN) :: mask
         REAL, DIMENSION(istr:iend), INTENT(IN) :: 
     &    t,         ! temperature (degrees C)
     &    s,         ! salinity (PSU)
     &    dic_in,    ! total inorganic carbon (mmol/m^3)
     &    ta_in,     ! total alkalinity (meq/m^3)
     &    pt_in,     ! inorganic phosphate (mmol/m^3)
     &    sit_in,    ! inorganic silicate (mmol/m^3)
     &    phlo,      ! lower limit of pH range
     &    phhi,      ! upper limit of pH range
     &    xco2_in,   ! atmospheric mole fraction CO2 in dry air (ppmv)
     &    atmpres     ! atmospheric pressure (atmosphere)


    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

       REAL, DIMENSION(istr:iend), INTENT(OUT) :: 
     &    ph,        ! computed ph values, for initial guess on next time step
     &    co2star,   ! CO2*water (mmol/m^3)
     &    dco2star,  ! delta CO2 (mmol/m^3)
     &    pco2oc,    ! oceanic pCO2 (ppmv)
     &    dpco2      ! Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)


    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

       INTEGER :: i

        REAL :: 
     &    mass_to_vol,   ! (mol/kg) -> (mmol/m^3)
     &    vol_to_mass,   ! (mmol/m^3) -> (mol/kg)
     &    tk,            ! temperature (K)
     &    is,            ! ionic strength
     &    scl,           ! chlorinity
     &    co2starair,    ! co2star saturation
     &    tk100, tk1002, invtk, dlogtk, is2, sqrtis, 
     &    s2, sqrts, s15, htotal2



        REAL, DIMENSION(istr:iend) :: 
     &    xco2,          ! atmospheric CO2 (atm)
     &    htotal,        ! free concentration of H ion
     &    x1, x2          ! bounds on htotal for solver

    !---------------------------------------------------------------------------
    !   check for existence of ocean points
    !---------------------------------------------------------------------------


       IF (COUNT(mask) == 0) THEN
         ph          = c0
         co2star     = c0
         dco2star    = c0
         pCO2oc      = c0
         dpCO2       = c0
       RETURN
      END IF


    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

       mass_to_vol = 1e6 * rho_sw
       vol_to_mass = c1 / mass_to_vol

    !---------------------------------------------------------------------------
    !   convert tracer units to per mass  xco2 from uatm to atm
    !---------------------------------------------------------------------------

        DO i = istr,iend
          IF (mask(i)) THEN
            dic(i)  = dic_in(i)  * vol_to_mass
            ta(i)   = ta_in(i)   * vol_to_mass
            pt(i)   = pt_in(i)   * vol_to_mass
            sit(i)  = sit_in(i)  * vol_to_mass
           xco2(i) = xco2_in(i) * 1e-6

          !---------------------------------------------------------------------
          !   Calculate all constants needed to convert between various
          !   measured carbon species. References for each equation are
          !   noted in the code.  Once calculated, the constants are stored
          !   and passed in the common block "const". The original version
          !   of this code was based on the code by Dickson in Version 2 of
          !   "Handbook of Methods for the Analysis of the Various Parameters
          !   of the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3,
          !   p25-26).
          !   Derive simple terms used more than once
          !---------------------------------------------------------------------

            tk       = T0_Kelvin + t(i)
            tk100    = tk * 1e-2
            tk1002   = tk100 * tk100
            invtk    = c1 / tk
            dlogtk   = LOG(tk)

            is       = 19.924 * s(i) / (c1000 - 1.005 * s(i))
            is2      = is * is
            sqrtis   = SQRT(is)
            sqrts    = SQRT(s(i))
            s15      = s(i) ** 1.5
            s2       = s(i) ** 2
            scl      = s(i) / 1.80655

          !---------------------------------------------------------------------
          !   f = k0(1-pH2O)*correction term for non-ideality
          !   Weiss  Price (1980, Mar. Chem., 8, 347-359;
          !                 Eq 13 with table 6 values)
          !---------------------------------------------------------------------

            ff(i) = EXP(-162.8301 + 218.2968/tk100 + 
     &          90.9241*LOG(tk100) - 
     &          1.47696*tk1002 + s(i)*(.025695 - .025225*tk100 + 
     &          0.0049867*tk1002))

          !---------------------------------------------------------------------
          !   K0 from Weiss 1974
          !---------------------------------------------------------------------

            k0(i) = EXP(93.4517/tk100 - 60.2409 + 23.3585*LOG(tk100) + 
     &          s(i)*(0.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))

          !---------------------------------------------------------------------
          !   k1 = [H][HCO3]/[H2CO3]
          !   k2 = [H][CO3]/[HCO3]
          !   Millero p.664 (1995) using Mehrbach et al. data on seawater scale
          !---------------------------------------------------------------------

          k1(i) = 10**(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk - 
     &         0.0118*s(i) + 0.000116*s2))

          k2(i) = 10**(-1*(1394.7*invtk + 4.777 - 
     &           0.0184*s(i) + 0.000118*s2))
 
          !---------------------------------------------------------------------
          !   kb = [H][BO2]/[HBO2]
          !   Millero p.669 (1995) using data from Dickson (1990)
          !---------------------------------------------------------------------

          kb(i) = EXP((-8966.90 - 2890.53*sqrts - 77.942*s(i) + 
     &          1.728*s15 - 0.0996*s2)*invtk + 
     &          (148.0248 + 137.1942*sqrts + 1.62142*s(i)) + 
     &          (-24.4344 - 25.085*sqrts - 0.2474*s(i)) * 
     &          dlogtk + 0.053105*sqrts*tk)

          !---------------------------------------------------------------------
          !   k1p = [H][H2PO4]/[H3PO4]
          !   DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
          !---------------------------------------------------------------------

          k1p(i) = EXP(-4576.752*invtk + 115.525 - 18.453 * dlogtk + 
     &         (-106.736*invtk + 0.69171) * sqrts + 
     &          (-0.65643*invtk - 0.01844) * s(i))

          !---------------------------------------------------------------------
          !   k2p = [H][HPO4]/[H2PO4]
          !   DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
          !---------------------------------------------------------------------

          k2p(i) = EXP(-8814.715*invtk + 172.0883 - 27.927 * dlogtk + 
     &          (-160.340*invtk + 1.3566) * sqrts + 
     &          (0.37335*invtk - 0.05778) * s(i))

          !---------------------------------------------------------------------
          !   k3p = [H][PO4]/[HPO4]
          !   DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
          !---------------------------------------------------------------------

          k3p(i) = EXP(-3070.75*invtk - 18.141 +  
     &          (17.27039*invtk + 2.81197) * sqrts + 
     &          (-44.99486*invtk - 0.09984) * s(i))

          !---------------------------------------------------------------------
          !   ksi = [H][SiO(OH)3]/[Si(OH)4]
          !   Millero p.671 (1995) using data from Yao and Millero (1995)
          !---------------------------------------------------------------------

          ksi(i) = EXP(-8904.2*invtk + 117.385 - 19.334 * dlogtk + 
     &         (-458.79*invtk + 3.5913) * sqrtis + 
     &         (188.74*invtk - 1.5998) * is + 
     &          (-12.1652*invtk + 0.07871) * is2 + 
     &          LOG(1.0-0.001005*s(i)))

          !---------------------------------------------------------------------
          !   kw = [H][OH]
          !   Millero p.670 (1995) using composite data
          !---------------------------------------------------------------------

          kw(i) = EXP(-13847.26*invtk + 148.9652 - 23.6521 * dlogtk + 
     &        (118.67*invtk - 5.977 + 1.0495 * dlogtk) * 
     &         sqrts - 0.01615 * s(i))

          !---------------------------------------------------------------------
          !   ks = [H][SO4]/[HSO4]
          !   Dickson (1990, J. chem. Thermodynamics 22, 113)
          !---------------------------------------------------------------------

          ks(i) = EXP(-4276.1*invtk + 141.328 - 23.093*dlogtk + 
     &          (-13856*invtk + 324.57 - 47.986*dlogtk) * 
     &          sqrtis + 
     &          (35474*invtk - 771.54 + 114.723*dlogtk) * is - 
     &          2698*invtk*is**1.5 + 1776*invtk*is2 + 
     &          LOG(1.0 - 0.001005*s(i)))

          !---------------------------------------------------------------------
          !   kf = [H][F]/[HF]
          !   Dickson and Riley (1979) -- change pH scale to total
          !---------------------------------------------------------------------

          kf(i) = EXP(1590.2*invtk - 12.641 + 1.525*sqrtis + 
     &         LOG(1.0 - 0.001005*s(i)) +  
     &          LOG(1.0 + (0.1400/96.062)*(scl)/ks(i)))

          !---------------------------------------------------------------------
          !   Calculate concentrations for borate, sulfate, and fluoride
          !   bt : Uppstrom (1974)
          !   st : Morris  Riley (1966)
          !   ft : Riley (1965)
          !---------------------------------------------------------------------

          bt(i) = 0.000232 * scl/10.811
          st(i) = 0.14 * scl/96.062
          ft(i) = 0.000067 * scl/18.9984

          x1(i) = c10 ** (-phhi(i))
          x2(i) = c10 ** (-phlo(i))

        END IF ! if mask

       END DO ! i loop


    !---------------------------------------------------------------------------
    !   If DIC and TA are known then either a root finding or iterative
    !   method must be used to calculate htotal. In this case we use
    !   the Newton-Raphson "safe" method taken from "Numerical Recipes"
    !   (function "rtsafe.f" with error trapping removed).
    !
    !   As currently set, this procedure iterates about 12 times. The
    !   x1 and x2 values set below will accomodate ANY oceanographic
    !   values. If an initial guess of the pH is known, then the
    !   number of iterations can be reduced to about 5 by narrowing
    !   the gap between x1 and x2. It is recommended that the first
    !   few time steps be run with x1 and x2 set as below. After that,
    !   set x1 and x2 to the previous value of the pH +/- ~0.5.
    !---------------------------------------------------------------------------

        CALL drtsafe_row(mask, x1, x2, xacc, htotal,istr,
     &    iend, k0, k1, k2, 
     &    kw, kb, ks, kf, k1p, k2p, k3p, ksi, 
     &    ff, bt, st, ft, dic, ta, pt, sit)
    !---------------------------------------------------------------------------
    !   Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
    !   ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
    !---------------------------------------------------------------------------

        DO i = istr,iend
          IF (mask(i)) THEN

          htotal2 = htotal(i) ** 2
          co2star(i) = dic(i) * htotal2 / 
     &          (htotal2 + k1(i)*htotal(i) + k1(i)*k2(i))
          co2starair = xco2(i) * ff(i) * atmpres(i)
          dco2star(i) = co2starair - co2star(i)
          ph(i) = -LOG10(htotal(i))

          !---------------------------------------------------------------------
          !   Add two output arguments for storing pco2oc
          !   Should we be using K0 or ff for the solubility here?
          !---------------------------------------------------------------------

          pco2oc(i) = co2star(i) / ff(i)
          dpCO2(i)    = pco2oc(i) - xco2(i) * atmpres(i)

          !---------------------------------------------------------------------
          !   Convert units of output arguments
          !   Note: pco2oc and dpCO2 are calculated in atm above.
          !---------------------------------------------------------------------

          co2star(i)  = co2star(i) * mass_to_vol
          dco2star(i) = dco2star(i) * mass_to_vol

          pco2oc(i) = pco2oc(i) * 1e6
          dpCO2(i)    = dpCO2(i) * 1e6

        ELSE ! if mask

          ph(i)       = c0
          co2star(i)  = c0
          dco2star(i) = c0
          pco2oc(i) = c0
          dpCO2(i)    = c0


        END IF ! if mask
       END DO ! i loop

      END SUBROUTINE co2calc_row

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE TALK_ROW   ****************************
!
!-------------------------------------------------------------------------

       SUBROUTINE talk_row(mask, x, fn, df,istr,iend, k0, k1, k2,
     &  kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff, bt, st, ft, dic, ta, pt, sit)

        implicit none
       REAL, PARAMETER :: xacc = 1e-10
       INTEGER, PARAMETER :: maxit = 100

  !-----------------------------------------------------------------------------
  !   declarations for function coefficients & species concentrations
  !-----------------------------------------------------------------------------
         integer :: istr,iend

        REAL, DIMENSION(istr:iend) :: 
     &  k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff, 
     &  bt, st, ft, dic, ta, pt, sit

    !---------------------------------------------------------------------------
    !   This routine computes TA as a function of DIC, htotal and constants.
    !   It also calculates the derivative of this function with respect to
    !   htotal. It is used in the iterative solution for htotal. In the call
    !   "x" is the input value for htotal, "fn" is the calculated value for
    !   TA and "df" is the value for dTA/dhtotal.
    !---------------------------------------------------------------------------
         real,parameter::c1=1.0, c2=2.0, c3=3.0 

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------
        LOGICAL, DIMENSION(istr:iend), INTENT(IN) :: mask
        REAL, DIMENSION(istr:iend), INTENT(IN) :: x

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

         REAL, DIMENSION(istr:iend), INTENT(OUT) :: fn, df

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

         INTEGER :: i

          REAL :: 
     &    x1, x2, x3, k12, k12p, k123p, a, a2, da, b, b2, db, c

    !---------------------------------------------------------------------------

      DO i = istr,iend
        IF (mask(i)) THEN
          x1 = x(i)
          x2 = x1 * x1
          x3 = x2 * x1
          k12 = k1(i) * k2(i)
          k12p = k1p(i) * k2p(i)
          k123p = k12p * k3p(i)
          a = x3 + k1p(i) * x2 + k12p * x1 + k123p
          a2 = a * a
          da = c3 * x2 + c2 * k1p(i) * x1 + k12p
          b = x2 + k1(i) * x1 + k12
          b2 = b * b
          db = c2 * x1 + k1(i)
          c = c1 + st(i)/ks(i)
          !---------------------------------------------------------------------
          !   fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
          !---------------------------------------------------------------------

          fn(i) = k1(i) * x1 * dic(i)/b + 
     &          c2 * dic(i) * k12/b + 
     &          bt(i)/(c1 + x1/kb(i)) + 
     &          kw(i)/x1 + 
     &          pt(i) * k12p * x1/a + 
     &          c2 * pt(i) * k123p/a + 
     &          sit(i)/(c1 + x1/ksi(i)) - 
     &          x1/c - 
     &          st(i)/(c1 + ks(i)/x1/c) - 
     &          ft(i)/(c1 + kf(i)/x1) - 
     &          pt(i) * x3/a - 
     &          ta(i)

          !---------------------------------------------------------------------
          !   df = d(fn)/dx
          !---------------------------------------------------------------------

          df(i) = ((k1(i)*dic(i)*b) - k1(i)*x1*dic(i)*db)/b2 - 
     &          c2 * dic(i) * k12 * db/b2 - 
     &          bt(i)/kb(i)/(c1+x1/kb(i)) ** 2 - 
     &          kw(i)/x2 + 
     &          (pt(i) * k12p * (a - x1 * da))/a2 - 
     &          c2 * pt(i) * k123p * da/a2 - 
     &          sit(i)/ksi(i)/(c1+x1/ksi(i)) ** 2 - 
     &          c1/c + 
     &          st(i) * (c1 + ks(i)/x1/c)**(-2) * (ks(i)/c/x2) + 
     &          ft(i) * (c1 + kf(i)/x1)**(-2) * kf(i)/x2 - 
     &          pt(i) * x2 * (c3 * a - x1 * da)/a2

        END IF ! if mask
        END DO ! i loop

      END SUBROUTINE talk_row

!-------------------------------------------------------------------------
!
!     ************   SUBROUTINE DRTSAFE_ROW   ****************************
!
!-------------------------------------------------------------------------

        SUBROUTINE drtsafe_row(mask_in, x1, x2, xacc, soln,istr,
     &   iend, k0, k1, k2, 
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi, 
     &    ff, bt, st, ft, dic, ta, pt, sit)
        implicit none
       INTEGER, PARAMETER :: maxit = 100
       integer::istr,iend
  !-----------------------------------------------------------------------------
  !   declarations for function coefficients & species concentrations
  !-----------------------------------------------------------------------------

        REAL, DIMENSION(istr:iend) :: 
     &  k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff, 
     &  bt, st, ft, dic, ta, pt, sit

    !---------------------------------------------------------------------------
    !   Vectorized version of drtsafe, which was a modified version of
    !      Numerical Recipes algorithm.
    !   Keith Lindsay, Oct 1999
    !
    !   Algorithm comment :
    !      Iteration from Newton's method is used unless it leaves
    !      bracketing interval or the dx is > 0.5 the previous dx.
    !      In that case, bisection method is used.
    !---------------------------------------------------------------------------
        real,parameter::c0=0.0, c2=2.0
  
    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------
  
        LOGICAL, DIMENSION(istr:iend), INTENT(IN) :: mask_in
        REAL, DIMENSION(istr:iend), INTENT(IN) :: x1, x2
        REAL, INTENT(IN) :: xacc

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

       REAL, DIMENSION(istr:iend), INTENT(OUT) :: soln

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

       LOGICAL :: leave_bracket, dx_decrease
       LOGICAL, DIMENSION(istr:iend) :: mask
       INTEGER ::  i, it
       REAL :: temp
       REAL, DIMENSION(istr:iend) :: xlo, xhi, flo, fhi, f, 
     &     df, dxold, dx

    !---------------------------------------------------------------------------
    !   bracket root at each location and set up first iteration
    !---------------------------------------------------------------------------

       mask = mask_in

       CALL talk_row(mask, x1, flo, df,istr,iend, k0, k1, k2, 
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi, 
     &    ff, bt, st, ft, dic, ta, pt, sit)
       CALL talk_row(mask, x2, fhi, df,istr,iend, k0, k1, k2, 
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi, 
     &    ff, bt, st, ft, dic, ta, pt, sit)

       DO i = istr,iend
         IF (mask(i)) THEN
           IF (flo(i) .LT. c0) THEN
              xlo(i) = x1(i)
              xhi(i) = x2(i)
           ELSE
             xlo(i) = x2(i)
             xhi(i) = x1(i)
             temp = flo(i)
             flo(i) = fhi(i)
             fhi(i) = temp
          END IF
          soln(i) = 0.5 * (xlo(i) + xhi(i))
          dxold(i) = ABS(xlo(i) - xhi(i))
          dx(i) = dxold(i)
         END IF
        END DO

        CALL talk_row(mask, soln, f, df,istr,iend, k0, k1, k2, 
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi, 
     &    ff, bt, st, ft, dic, ta, pt, sit)

    !---------------------------------------------------------------------------
    !   perform iterations, zeroing mask when a location has converged
    !---------------------------------------------------------------------------

       DO it = 1,maxit
        DO i = istr,iend
           IF (mask(i)) THEN
             leave_bracket = ((soln(i)-xhi(i))*df(i)-f(i)) * 
     &             ((soln(i)-xlo(i))*df(i)-f(i)) .GE. 0
             dx_decrease = ABS(c2 * f(i)) .LE. ABS(dxold(i) * df(i))
             IF (leave_bracket .OR. .NOT. dx_decrease) THEN
                dxold(i) = dx(i)
                dx(i) = 0.5 * (xhi(i) - xlo(i))
                soln(i) = xlo(i) + dx(i)
                IF (xlo(i) .EQ. soln(i)) mask(i) = .FALSE.
             ELSE
                dxold(i) = dx(i)
                dx(i) = -f(i) / df(i)
                temp = soln(i)
                soln(i) = soln(i) + dx(i)
                IF (temp .EQ. soln(i)) mask(i) = .FALSE.
             END IF
             IF (ABS(dx(i)) .LT. xacc) mask(i) = .FALSE.
          END IF
        END DO

        IF (.NOT. ANY(mask)) RETURN

        CALL talk_row(mask, soln, f, df,istr,iend, k0, k1, k2, 
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi, 
     &    ff, bt, st, ft, dic, ta, pt, sit)

        DO i = istr,iend
          IF (mask(i)) THEN
             IF (f(i) .LT. c0) THEN
                xlo(i) = soln(i)
                flo(i) = f(i)
             ELSE
                xhi(i) = soln(i)
                fhi(i) = f(i)
             END IF
          END IF
        END DO

       END DO ! iteration loop

       END SUBROUTINE drtsafe_row


#else  /* BIOLOGY_BEC2 */
       subroutine ecosys_bec2_tile_empty ()
       
       end
#endif  /* BIOLOGY_BEC2 */
