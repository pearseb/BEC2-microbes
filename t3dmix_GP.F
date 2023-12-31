#include "cppdefs.h"
#if defined SOLVE3D && (defined TS_DIF2 || defined TS_DIF4 )

      subroutine t3dmix (tile)
      implicit none
      integer tile, itrc
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      do itrc=1,NT
        call t3dmix_tile (istr,iend,jstr,jend, itrc, A3d(1,1),A2d(1,1),
     &                                     A2d(1,2), A2d(1,3),A2d(1,5),
     &                       A2d(1,7), A2d(1,9), A2d(1,11), A2d(1,13))
      enddo
      end

      subroutine t3dmix_tile (istr,iend,jstr,jend, itrc, LapT,FX,FE,
     &                               FSC,dTdz, dTdx,dTde,dZdx,dZde)
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, itrc, i,j,k,k1,k2
# ifdef TS_DIF4
     &                           , imin,imax,jmin,jmax
# endif
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: LapT
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX,FE
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,2) ::  FSC,dTdz,
     &                                    dTdx,dTde, dZdx,dZde
      real cff, cff1, cff2, cff3, cff4

# include "grid.h"
# include "ocean3d.h"
# include "mixing.h"
# include "climat.h"
# include "scalars.h"
# ifdef PHYS_FLUX_ANALYSIS
#  include "physflux.h"
# endif
# include "compute_auxiliary_bounds.h"

# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
# ifdef TS_DIF2
                                             !!! WARNING: RECURSIVE
! Compute horizontal (along geopotential     !!!   BLOCKING SEQUENCE
! surfaces) Laplacian diffusion terms.       !!!

      k2=1            ! vertical   dTdx,dTde(:,:,k2) k+1   rho-points
      do k=0,N,+1     ! placement    dTdz,FSC(:,:,k2) k+1/2   W-points
        k1=k2         ! of scratch dTdx,dTde(:,:,k1) k     rho-points
        k2=3-k1       ! variables    dTdz,FSC(:,:,k1) k-1/2   W-points
        if (k<N) then
          do j=jstr,jend
            do i=istr,iend+1
              cff=0.5*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
              dZdx(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i-1,j,k+1))
              dTdx(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)
     &                               -t(i-1,j,k+1,nrhs,itrc)
#  ifdef CLIMAT_TS_MIXH
     &               -tclm(i,j,k+1,itrc)+tclm(i-1,j,k+1,itrc)
#  endif
     &                                                      )
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
              dZde(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i,j-1,k+1))
              dTde(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)
     &                               -t(i,j-1,k+1,nrhs,itrc)
#  ifdef CLIMAT_TS_MIXH
     &               -tclm(i,j,k+1,itrc)+tclm(i,j-1,k+1,itrc)
#  endif
     &                                                      )
            enddo
          enddo
        endif

        if (k==0 .or. k==N) then
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              dTdz(i,j,k2)=0.0     ! Set no-flux bottom or
              FSC(i,j,k2)=0.0       ! top boundary conditions.
            enddo
          enddo
        else
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              dTdz(i,j,k2)=(t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc)
#  ifdef CLIMAT_TS_MIXH
     &                          -tclm(i,j,k+1,itrc)+tclm(i,j,k,itrc)
#  endif
     &                                   )/(z_r(i,j,k+1)-z_r(i,j,k))
            enddo
          enddo
        endif

        if (k>0) then
          do j=jstr,jend
            do i=istr,iend+1
              FX(i,j)=0.25*(diff2(i,j,itrc)+diff2(i-1,j,itrc))
     &                     *(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &         dTdx(i  ,j,k1)-0.5*(
     &             min(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k1)+dTdz(i,j,k2))
     &            +max(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              FE(i,j)=0.25*(diff2(i,j,itrc)+diff2(i,j-1,itrc))
     &                     *(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &         dTde(i,j,k1)-0.5*(
     &             min(dZde(i,j,k1),0.)*(dTdz(i,j-1,k1)+dTdz(i,j,k2))
     &            +max(dZde(i,j,k1),0.)*(dTdz(i,j-1,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          if (k<N) then
            do j=jstr,jend
              do i=istr,iend
c**             cff1=min(dZdx(i,j,k1)+dZdx(i+1,j,k2),0.)
c**             cff2=max(dZdx(i,j,k2)+dZdx(i+1,j,k1),0.)
c**
c**             cff3=min(dZde(i,j,k1)+dZde(i,j+1,k2),0.)
c**             cff4=max(dZde(i,j,k2)+dZde(i,j+1,k1),0.)
c**
c**             FSC(i,j,k2)=0.25*diff2(i,j,itrc)*(
c**  &                             dTdz(i,j,k2)*( (cff1+cff2)**2
c**  &                                           +(cff3+cff4)**2)
c**
c**  &                           -cff1*(dTdx(i,j,k1)+dTdx(i+1,j,k2))
c**  &                           -cff2*(dTdx(i,j,k2)+dTdx(i+1,j,k1))
c**
c**  &                           -cff3*(dTde(i,j,k1)+dTde(i,j+1,k2))
c**  &                           -cff4*(dTde(i,j,k2)+dTde(i,j+1,k1))
c**  &                                                             )

                cff1=min(dZdx(i  ,j,k1),0.)
                cff2=min(dZdx(i+1,j,k2),0.)

                cff3=max(dZdx(i  ,j,k2),0.)
                cff4=max(dZdx(i+1,j,k1),0.)

                FSC(i,j,k2)=0.5*diff2(i,j,itrc)*(
     &                       cff1*(cff1*dTdz(i,j,k2)-dTdx(i  ,j,k1))
     &                        +cff2*(cff2*dTdz(i,j,k2)-dTdx(i+1,j,k2))

     &                        +cff3*(cff3*dTdz(i,j,k2)-dTdx(i  ,j,k2))
     &                        +cff4*(cff4*dTdz(i,j,k2)-dTdx(i+1,j,k1))
     &                                                               )
                cff1=min(dZde(i,j  ,k1),0.)
                cff2=min(dZde(i,j+1,k2),0.)

                cff3=max(dZde(i,j  ,k2),0.)
                cff4=max(dZde(i,j+1,k1),0.)

                FSC(i,j,k2)=FSC(i,j,k2)+0.5*diff2(i,j,itrc)*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTde(i,j  ,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTde(i,j+1,k2))

     &                       +cff3*(cff3*dTdz(i,j,k2)-dTde(i,j  ,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTde(i,j+1,k1))
     &                                                              )
              enddo
            enddo
          endif
          do j=jstr,jend
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*( pm(i,j)*
     &               pn(i,j)*( FX(i+1,j)-FX(i,j) +FE(i,j+1)-FE(i,j) )
     &                                       +FSC(i,j,k2)-FSC(i,j,k1)
     &                                                   )/Hz(i,j,k)
            enddo
          enddo
#  ifdef FULL_PHYS_FLUX_ANALYSIS
#if !defined SELECTED_FLUX_ANALYSIS /* AH */
        do j=Jstr,Jend
           do i=Istr,Iend+1
              ! unit: mmol/s
              HorXMixFlux(i,j,k,itrc) = -FX(i,j)
#   ifdef MASKING
     &              * rmask(i,j)
#   endif
           end do
        end do
        do j=Jstr,Jend+1
           do i=Istr,Iend+1
              ! unit: mmol/s
              HorYMixFlux(i,j,k,itrc) = -FE(i,j)
#   ifdef MASKING
     &              * rmask(i,j)
#   endif
           end do
        end do
        do j=Jstr,Jend
           do i=Istr,Iend
              ! unit: mmol/(m2 s)
              ! reverse sign so that fluxes are positive upwards
              VertMixFlux(i,j,k,itrc) = -FSC(i,j,k2)
#   ifdef MASKING
     &             * rmask(i,j)
#   endif
            enddo
          enddo
       else ! <-- k.gt.0
          do j=Jstr,Jend+1
             do i=Istr,Iend+1
                VertMixFlux(i,j,k,itrc) = 0.0
             end do
          end do
#endif /* SELECTED_FLUX_ANALYSIS */
#  endif /* FULL_PHYS_FLUX_ANALYSIS */
        endif  ! <-- k>0
      enddo   ! --> k
# endif /* T_DIF2 */
# ifdef TS_DIF4
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        imin=istr
      else
        imin=istr-1
      endif
      if (EASTERN_EDGE) then
        imax=iend
      else
        imax=iend+1
      endif
#  else
      imin=istr-1
      imax=iend+1
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if (NORTHERN_EDGE) then
        jmax=jend
      else
        jmax=jend+1
      endif
#  else
      jmin=jstr-1
      jmax=jend+1
#  endif

! Compute horizontal biharmonic diffusion      !!! WARNING: RECURSIVE
! along geopotential surfaces. The biharmonic  !!!  BLOCKING SEQUENCE
! operator is computed by applying the rotated !!!
! Laplacian operator twice.

! THE FIRST LAPLACIAN: ! vertical   dTdx,dTde(:,:,k2)  k+1   rho-points
!                      ! placement    dTdz,FSC(:,:,k2) k+1/2   W-points
      k2=1             ! of scratch dTdx,dTde(:,:,k1)  k     rho-points
      do k=0,N,+1      ! variables    dTdz,FSC(:,:,k1) k-1/2   W-points
        k1=k2
        k2=3-k1
        if (k<N) then
          do j=jmin,jmax
            do i=imin,imax+1
              cff=0.5*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
              dZdx(i,j,k2)=cff*(z_r(i,j,k+1)-z_r(i-1,j,k+1))
              dTdx(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-
     &                               t(i-1,j,k+1,nrhs,itrc))
            enddo
          enddo
          do j=jmin,jmax+1
            do i=imin,imax
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
              dZde(i,j,k2)=cff*(z_r(i,j,k+1)-z_r(i,j-1,k+1))
              dTde(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-
     &                               t(i,j-1,k+1,nrhs,itrc))
            enddo
          enddo
        endif

        if (k==0 .or. k==N) then
          do j=jmin-1,jmax+1
            do i=imin-1,imax+1
              dTdz(i,j,k2)=0.0
              FSC(i,j,k2)=0.0
            enddo
          enddo
        else
          do j=jmin-1,jmax+1
            do i=imin-1,imax+1
              dTdz(i,j,k2)=(t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc))
     &                                    /(z_r(i,j,k+1)-z_r(i,j,k))
            enddo
          enddo
        endif

        if (k>0) then
          do j=jmin,jmax
            do i=imin,imax+1
              FX(i,j)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &                                                dTdx(i  ,j,k1)
     &      -0.5*( min(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k1)+dTdz(i,j,k2))
     &            +max(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          do j=jmin,jmax+1
            do i=imin,imax
              FE(i,j)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &                                                  dTde(i,j,k1)
     &      -0.5*( min(dZde(i,j,k1),0.)*(dTdz(i,j-1,k1)+dTdz(i,j,k2))
     &            +max(dZde(i,j,k1),0.)*(dTdz(i,j-1,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          if (k<N) then
            do j=jmin,jmax
              do i=imin,imax
c**             cff1=min(dZdx(i,j,k1)+dZdx(i+1,j,k2),0.)
c**             cff2=max(dZdx(i,j,k2)+dZdx(i+1,j,k1),0.)
c**
c**             cff3=min(dZde(i,j,k1)+dZde(i,j+1,k2),0.)
c**             cff4=max(dZde(i,j,k2)+dZde(i,j+1,k1),0.)
c**
c**             FSC(i,j,k2)=0.25*( dTdz(i,j,k2)*( (cff1+cff2)**2
c**  &                                           +(cff3+cff4)**2)
c**
c**  &                           -cff1*(dTdx(i,j,k1)+dTdx(i+1,j,k2))
c**  &                           -cff2*(dTdx(i,j,k2)+dTdx(i+1,j,k1))
c**
c**  &                           -cff3*(dTde(i,j,k1)+dTde(i,j+1,k2))
c**  &                           -cff4*(dTde(i,j,k2)+dTde(i,j+1,k1))
c**  &                                                             )

                cff1=min(dZdx(i  ,j,k1),0.)
                cff2=min(dZdx(i+1,j,k2),0.)

                cff3=max(dZdx(i  ,j,k2),0.)
                cff4=max(dZdx(i+1,j,k1),0.)

                FSC(i,j,k2)=0.5*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTdx(i  ,j,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTdx(i+1,j,k2))

     &                       +cff3*(cff3*dTdz(i,j,k2)-dTdx(i  ,j,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTdx(i+1,j,k1))
     &                                                              )
                cff1=min(dZde(i,j  ,k1),0.)
                cff2=min(dZde(i,j+1,k2),0.)

                cff3=max(dZde(i,j  ,k2),0.)
                cff4=max(dZde(i,j+1,k1),0.)

                FSC(i,j,k2)=FSC(i,j,k2)+0.5*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTde(i,j  ,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTde(i,j+1,k2))

     &                       +cff3*(cff3*dTdz(i,j,k2)-dTde(i,j  ,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTde(i,j+1,k1))
     &                                                              )
              enddo
            enddo
          endif
          do j=jmin,jmax
            do i=imin,imax
              LapT(i,j,k)=( pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                       +FE(i,j+1)-FE(i,j) )
     &                                   +FSC(i,j,k2)-FSC(i,j,k1)
     &                                               )/Hz(i,j,k)
            enddo
          enddo
        endif   ! <-- k>0
      enddo   ! --> k

#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then            ! Apply lateral boundary
          do k=1,N                        ! conditions to the first
            do j=jmin,jmax                ! Laplacian (in cases other
#   ifdef WESTERN_WALL
              LapT(istr-1,j,k)=0.         ! than periodic: closed
#   else
              LapT(istr-1,j,k)=LapT(istr,j,k)      !  or gradient)
#   endif
            enddo
          enddo
        endif
        if (EASTERN_EDGE) then
          do k=1,N
            do j=jmin,jmax
#   ifdef EASTERN_WALL
              LapT(iend+1,j,k)=0.
#   else
              LapT(iend+1,j,k)=LapT(iend,j,k)
#   endif
            enddo
          enddo
        endif
#  endif /* !EW_PERIODIC */
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do k=1,N
            do i=imin,imax
#   ifdef SOUTHERN_WALL
              LapT(i,jstr-1,k)=0.
#   else
              LapT(i,jstr-1,k)=LapT(i,jstr,k)
#   endif
            enddo
          enddo
        endif
        if (NORTHERN_EDGE) then
          do k=1,N
            do i=imin,imax
#   ifdef NORTHERN_WALL
              LapT(i,jend+1,k)=0.
#   else
              LapT(i,jend+1,k)=LapT(i,jend,k)
#   endif
            enddo
          enddo
        endif
#  endif /* !NS_PERIODIC */
!
! THE SECOND LAPLACIAN
!
      k2=1
      do k=0,N,+1
       k1=k2
       k2=3-k1

        if (k<N) then
          do j=jstr,jend
            do i=istr,iend+1
              cff=0.5*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
              dZdx(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i-1,j,k+1))
              dTdx(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i-1,j,k+1))
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
              dZde(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i,j-1,k+1))
              dTde(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i,j-1,k+1))
            enddo
          enddo
        endif

        if (k==0 .or. k==N) then
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              dTdz(i,j,k2)=0.0
              FSC(i,j,k2)=0.0
            enddo
          enddo
        else
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              dTdz(i,j,k2)=(LapT(i,j,k+1)-LapT(i,j,k))
     &                       /(z_r(i,j,k+1)-z_r(i,j,k))
            enddo
          enddo
        endif

        if (k>0) then
          do j=jstr,jend
            do i=istr,iend+1
              FX(i,j)=-0.25*(diff4(i,j,itrc)+diff4(i-1,j,itrc))
     &                     *(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &         dTdx(i  ,j,k1)-0.5*(
     &             min(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k1)+dTdz(i,j,k2))
     &            +max(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              FE(i,j)=-0.25*(diff4(i,j,itrc)+diff4(i,j-1,itrc))
     &                     *(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &         dTde(i,j,k1)-0.5*(
     &             min(dZde(i,j,k1),0.)*(dTdz(i,j-1,k1)+dTdz(i,j,k2))
     &            +max(dZde(i,j,k1),0.)*(dTdz(i,j-1,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          if (k<N) then
            do j=jstr,jend
              do i=istr,iend
c**             cff1=min(dZdx(i,j,k1)+dZdx(i+1,j,k2),0.)
c**             cff2=max(dZdx(i,j,k2)+dZdx(i+1,j,k1),0.)
c**
c**             cff3=min(dZde(i,j,k1)+dZde(i,j+1,k2),0.)
c**             cff4=max(dZde(i,j,k2)+dZde(i,j+1,k1),0.)
c**
c**             FSC(i,j,k2)=-0.25*diff4(i,j,itrc)*(
c**  &                             dTdz(i,j,k2)*( (cff1+cff2)**2
c**  &                                           +(cff3+cff4)**2)
c**
c**  &                           -cff1*(dTdx(i,j,k1)+dTdx(i+1,j,k2))
c**  &                           -cff2*(dTdx(i,j,k2)+dTdx(i+1,j,k1))
c**
c**  &                           -cff3*(dTde(i,j,k1)+dTde(i,j+1,k2))
c**  &                           -cff4*(dTde(i,j,k2)+dTde(i,j+1,k1))
c**  &                                                             )

                cff1=min(dZdx(i  ,j,k1),0.)
                cff2=min(dZdx(i+1,j,k2),0.)

                cff3=max(dZdx(i  ,j,k2),0.)
                cff4=max(dZdx(i+1,j,k1),0.)

                FSC(i,j,k2)=-0.5*diff4(i,j,itrc)*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTdx(i  ,j,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTdx(i+1,j,k2))

     &                       +cff3*(cff3*dTdz(i,j,k2)-dTdx(i  ,j,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTdx(i+1,j,k1))
     &                                                              )
                cff1=min(dZde(i,j  ,k1),0.)
                cff2=min(dZde(i,j+1,k2),0.)

                cff3=max(dZde(i,j  ,k2),0.)
                cff4=max(dZde(i,j+1,k1),0.)

                FSC(i,j,k2)=FSC(i,j,k2)-0.5*diff4(i,j,itrc)*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTde(i,j  ,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTde(i,j+1,k2))

     &                       +cff3*(cff3*dTdz(i,j,k2)-dTde(i,j  ,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTde(i,j+1,k1))
     &                                                              )
              enddo
            enddo
          endif
          do j=jstr,jend
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) + dt*( pm(i,j)
     &              *pn(i,j)*( FX(i+1,j)-FX(i,j) +FE(i,j+1)-FE(i,j) )
     &                        +FSC(i,j,k2)-FSC(i,j,k1)    )/Hz(i,j,k)
            enddo
          enddo
#  ifdef FULL_PHYS_FLUX_ANALYSIS
#if !defined SELECTED_FLUX_ANALYSIS /* AH */
        do j=Jstr,Jend
           do i=Istr,Iend+1
              ! unit: mmol/s
              HorXMixFlux(i,j,k,itrc) = -FX(i,j)
#   ifdef MASKING
     &              * rmask(i,j)
#   endif
           end do
        end do
        do j=Jstr,Jend+1
           do i=Istr,Iend+1
              ! unit: mmol/s
              HorYMixFlux(i,j,k,itrc) = -FE(i,j)
#   ifdef MASKING
     &              * rmask(i,j)
#   endif
           end do
        end do
        do j=Jstr,Jend
           do i=Istr,Iend
              ! unit: mmol/(m2 s)
              ! reverse sign so that fluxes are positive upwards
              VertMixFlux(i,j,k,itrc) = -FSC(i,j,k2)
#   ifdef MASKING
     &             * rmask(i,j)
#   endif
            enddo
          enddo
       else ! <-- k.gt.0
          do j=Jstr,Jend+1
             do i=Istr,Iend+1
                VertMixFlux(i,j,k,itrc) = 0.0
             end do
          end do
#endif /* SELECTED_FLUX_ANALYSIS */
#  endif /* FULL_PHYS_FLUX_ANALYSIS */
        endif  ! <-- k>0
      enddo   ! --> k
# endif /* TS_DIF4 */

# ifdef EXCHANGE
        call exchange_tile (istr,iend,jstr,jend,
     &                      t(START_2D_ARRAY,1,nnew,itrc),N)
# endif
      end
#else
      subroutine t3dmix_empty
      end
#endif /* SOLVE3D */
