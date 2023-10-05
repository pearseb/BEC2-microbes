#include "cppdefs.h"


      subroutine check_switches1(ierr)

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if new CPP-switches  !!!!!!
!!!!!! were introduced into "cppdefs.h".  NO ACTION IS NEEDED  !!!!!!
!!!!!! if changes in "cppdefs.h" are limited to activation or  !!!!!!
!!!!!! deactivation of previously known switches.              !!!!!!
!!!!!! To refresh this file compile and execute "cppcheck.F"   !!!!!!
!!!!!! as an independent program, or use commands              !!!!!!
!!!!!! "make checkdefs" or "make depend".                      !!!!!!
!!!!!! Total number of CPP-switches:                     61 !!!!!!


      implicit none
      integer ierr, is,ie, iexample
#include "param.h"
#include "strings.h"
#ifdef MPI
# include "scalars.h"
#endif
      do is=1,max_opt_size
        cpps(is:is)=' '
      enddo
      iexample=0
!
! cppdefs_UP.h                                                    
!
      is=1 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='<cppdefs_UP.h>'
      cpps(ie+1:ie+1)=' '

#ifdef SOLVE3D
      iexample=iexample+1
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SOLVE3D'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef UV_ADV
      is=ie+2 ; ie=is+5
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='UV_ADV'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef UV_COR
      is=ie+2 ; ie=is+5
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='UV_COR'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef NONLIN_EOS
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='NONLIN_EOS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPLIT_EOS
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SPLIT_EOS'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SALINITY
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SALINITY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SFLX_CORR
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SFLX_CORR'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DIURNAL_SRFLUX
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='DIURNAL_SRFLUX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef QCORRECTION
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='QCORRECTION'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef T_FRC_BRY
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='T_FRC_BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef Z_FRC_BRY
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='Z_FRC_BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef M3_FRC_BRY
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='M3_FRC_BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef M2_FRC_BRY
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='M2_FRC_BRY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPONGE
      is=ie+2 ; ie=is+5
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SPONGE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPONGE_WIDTH
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SPONGE_WIDTH'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef UV_VIS2
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='UV_VIS2'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef TS_DIF2
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='TS_DIF2'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_MIXING
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_MIXING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_KPP
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_KPP'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_NONLOCAL
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_NONLOCAL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_RIMIX
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_RIMIX'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_CONVEC
      is=ie+2 ; ie=is+9
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_CONVEC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef LMD_BKPP
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='LMD_BKPP'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef CURVGRID
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='CURVGRID'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPHERICAL
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SPHERICAL'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef MASKING
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='MASKING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef MASK_LAND_DATA
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='MASK_LAND_DATA'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_M2FLATHER
      is=ie+2 ; ie=is+12
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_M2FLATHER'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_M3ORLANSKI
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_M3ORLANSKI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_TORLANSKI
      is=ie+2 ; ie=is+12
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_TORLANSKI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BIOLOGY_BEC2
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BIOLOGY_BEC2'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BIOLOGY
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BIOLOGY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DAILYPAR_BEC
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='DAILYPAR_BEC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BIOLOGY_NPZDOC
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BIOLOGY_NPZDOC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DAILYPAR_PHOTOINHIBITION
      is=ie+2 ; ie=is+23
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='DAILYPAR_PHOTOINHIBITION'
      cpps(ie+1:ie+1)=' '
#endif
!
! cppdefs.h                                                       
!
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='<cppdefs.h>'
      cpps(ie+1:ie+1)=' '

#ifdef GRID_SIZE
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='GRID_SIZE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DOMAIN_TILING
      is=ie+2 ; ie=is+12
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='DOMAIN_TILING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_WEST
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_WEST'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_NORTH
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_NORTH'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_SOUTH
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_SOUTH'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_M3ORLANSKI
      is=ie+2 ; ie=is+13
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_M3ORLANSKI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef OBC_TORLANSKI
      is=ie+2 ; ie=is+12
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='OBC_TORLANSKI'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SPONGE
      is=ie+2 ; ie=is+5
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SPONGE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SFLX_CORR
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SFLX_CORR'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef QCORRECTION
      is=ie+2 ; ie=is+10
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='QCORRECTION'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef VFLX_CORR
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='VFLX_CORR'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef AVERAGES
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='AVERAGES'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef SLICE_AVG
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='SLICE_AVG'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef CALENDAR
      is=ie+2 ; ie=is+7
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='CALENDAR'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef STARTDATE
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='STARTDATE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BIOLOGY_BEC2
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BIOLOGY_BEC2'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BIOLOGY
      is=ie+2 ; ie=is+6
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BIOLOGY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DEFAULT_BRY_VALUES
      is=ie+2 ; ie=is+17
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='DEFAULT_BRY_VALUES'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef BEC2_DIAG
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='BEC2_DIAG'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef Ncycle_SY
      is=ie+2 ; ie=is+8
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='Ncycle_SY'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef EXPLICIT_MICROBES
      is=ie+2 ; ie=is+16
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='EXPLICIT_MICROBES'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef FACULTATIVE_MICROBES
      is=ie+2 ; ie=is+19
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='FACULTATIVE_MICROBES'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef IODINE
      is=ie+2 ; ie=is+5
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='IODINE'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef PCO2AIR_FORCING
      is=ie+2 ; ie=is+14
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='PCO2AIR_FORCING'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef DAILYPAR_BEC
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='DAILYPAR_BEC'
      cpps(ie+1:ie+1)=' '
#endif
#ifdef WRITE_DEPTHS
      is=ie+2 ; ie=is+11
      if (ie >= max_opt_size) goto 99
      cpps(is:ie)='WRITE_DEPTHS'
      cpps(ie+1:ie+1)=' '
#endif
      return
  99  mpi_master_only write(*,'(/1x,2A/14x,A)')
     &  '### ERROR: checkdefs :: Insufficient size of string "cpps" ',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
