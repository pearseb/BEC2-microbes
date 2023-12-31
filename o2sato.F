!  ---------------------------------------------------------------------
!  $Log: o2sato.F,v $
!  Revision 1.1  2001/10/19 21:22:14  frenzel
!  Routines from the OCMIP project for the ROMS carbonized ecosystem model.
!  Some of these routines were adapted for the use in a parallel code.
!
! 
!  Revision 1.1  1999/03/22 12:57:48  orr
!  Initial revision
! 
!  ---------------------------------------------------------------------
!  
#include "cppdefs.h"
#if defined BIOLOGY_NPZDOC && defined OXYGEN && defined OCMIP_OXYGENSAT
      function o2sato(T,S)
!
! ********************************************************************
!                                     
! Computes the oxygen saturation concentration at 1 atm total pressure
! in mol/m^3 given the temperature (t, in deg C) and the salinity (s,
! in permil). 
!
! FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
! THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
!
! *** NOTE: THE "A3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
! *** IT SHOULDN'T BE THERE.                                ***
!
! o2sato IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
! 0 permil <= S <= 42 permil
! C
! CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil, 
! o2sato = 0.282015 mol/m^3
!
! The unit for the return value is now [mmol/m^3], as needed by ROMS
!
! ********************************************************************
!
      implicit none
      real :: o2sato

! input variables
      real, intent(in) :: T   ! Temperature [C]
      real, intent(in) :: S   ! Salinity [PSU]

      real, parameter :: A0 = 2.00907
      real, parameter :: A1 = 3.22014
      real, parameter :: A2 = 4.05010
      real, parameter :: A3 = 4.94457
      real, parameter :: A4 = -2.56847E-1
      real, parameter :: A5 = 3.88767
      real, parameter :: B0 = -6.24523E-3
      real, parameter :: B1 = -7.37614E-3
      real, parameter :: B2 = -1.03410E-2
      real, parameter :: B3 = -8.17083E-3
      real, parameter :: C0 = -4.88682E-7

! other variables
      real :: TT      
      real :: TK        ! Temperature [K]
      real :: TS, TS2, TS3, TS4, TS5
      real :: CO
!
      TT  = 298.15-T
      TK  = 273.15+T
      TS  = LOG(TT/TK)
      TS2 = TS**2
      TS3 = TS**3
      TS4 = TS**4
      TS5 = TS**5
      CO  = A0 + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5
     $     + S*(B0 + B1*TS + B2*TS2 + B3*TS3)
     $     + C0*(S*S)
      o2sato = EXP(CO)
!
!  Convert from ml/l to mmol/m^-3 for ROMS
!
      o2sato = o2sato/22.3916*1000.0
      return
      end
#else
      subroutine o2sato_empty()
      return
      end
#endif /* BIOLOGY_NPZDOC && OXYGEN && OCMIP_OXYGENSAT */
