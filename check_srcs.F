#include "cppdefs.h"

      subroutine check_srcs

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if the new files     !!!!!!
!!!!!! were introduced into or deleted from the list of source !!!!!!
!!!!!! codes SRCS in the Makefile.                             !!!!!!

      implicit none
      integer i
#include "strings.h"
      do i=1,max_opt_size
        srcs(i:i)=' '
      enddo
      srcs(   1:   5)='IOSRC'
      srcs(   6:   7)=' :'
      srcs(   8:  19)=' checkdims.F'
      srcs(  20:  33)=' find_record.F'
      srcs(  34:  47)=' insert_node.F'
      srcs(  48:  65)=' ncdf_read_write.F'
      srcs(  66:  76)=' get_grid.F'
      srcs(  77:  87)=' get_init.F'
      srcs(  88:  98)=' wrt_grid.F'
      srcs(  99: 116)=' put_global_atts.F'
      srcs( 117: 126)=' def_rst.F'
      srcs( 127: 136)=' wrt_rst.F'
      srcs( 137: 146)=' def_his.F'
      srcs( 147: 156)=' wrt_his.F'
      srcs( 157: 166)=' set_avg.F'
      srcs( 167: 176)=' wrt_avg.F'
      srcs( 177: 188)=' def_statn.F'
      srcs( 189: 200)=' wrt_statn.F'
      srcs( 201: 213)=' get_forces.F'
      srcs( 214: 226)=' get_stflux.F'
      srcs( 227: 240)=' get_bry_all.F'
      srcs( 241: 253)=' set_forces.F'
      srcs( 254: 266)=' set_smflux.F'
      srcs( 267: 279)=' set_stflux.F'
      srcs( 280: 293)=' nf_read_bry.F'
      srcs( 294: 311)=' bulk_flux_COAMS.F'
      srcs( 312: 324)=' set_seaice.F'
      srcs( 325: 336)=' get_tides.F'
      srcs( 337: 348)=' set_tides.F'
      srcs( 349: 361)=' set_srflux.F'
      srcs( 362: 374)=' get_tclima.F'
      srcs( 375: 387)=' get_uclima.F'
      srcs( 388: 398)=' closecdf.F'
      srcs( 399: 412)=' init_agedye.F'
      srcs( 413: 425)=' get_agedye.F'
      srcs( 426: 437)=' set_cycle.F'
      srcs( 438: 451)=' get_tsource.F'
      srcs( 452: 458)=' KRNSRC'
      srcs( 459: 460)=' :'
      srcs( 461: 470)=' $(DRVSRC)'
      srcs( 471: 481)=' $(STEPSRC)'
      srcs( 482: 491)=' $(BRYSRC)'
      srcs( 492: 501)=' $(EOSSRC)'
      srcs( 502: 511)=' $(VMXSRC)'
      srcs( 512: 521)=' $(CHKSRC)'
      srcs( 522: 531)=' $(MPISRC)'
      srcs( 532: 536)=' SRCS'
      srcs( 537: 538)=' :'
      srcs( 539: 548)=' $(KRNSRC)'
      srcs( 549: 557)=' $(IOSRC)'
      srcs( 558: 567)=' $(BIOSRC)'
      return
      end
