!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!     This subroutine calculates the analytical solutions if available and
!     writes the results in the output files.
!     MODIFIED 16-JUL-1993 to include the analytical solution of test 'RCSR'
!     MODIFIED 26-JAN-1994 to include the analytical solution of test 'RBWI'
!     MODIFIED 30-OCT-1995 to include analytical solutions for the tracer
!===============================================================================
      subroutine pltout( basenm, nsdim, nx, ny, nz, mnx1, mny5, mnz1, &
      igeomx, nstep, gamma_ad, dt, timej,                &
      xznl, xzn, xznr, yznl, yzn, yznr, zzn, physic,     &
      rstfil, outfil, tdlfil, datfil, suffix, suffixd,   &
      nout1, tout1 )
      USE constantes
      USE tipos
      USE workarrays
      implicit none
# include "types.h"
!-------------------------------------------------------Input Variables
      
      character(4):: basenm

      integer(is):: nsdim, nx, ny, nz, mnx1, mny5, mnz1, igeomx, nstep

      real(rs)   :: gamma_ad, dt, timej

      real(rs)   :: xznl(*), xzn(*), xznr(*), &
      yznl(*), yzn(*), yznr(*), &
      zzn(*)

!type(Tphysic):: physic(0:mnx1,-4:mny5,0:mnz1)
      type(Tphysic):: physic(-4:mnx1,-4:mny5,-4:mnz1)

!------------------------------------------------------Output Variables
      character(7):: rstfil, outfil, tdlfil
      character(8):: datfil
      character(2):: suffix
      character(3):: suffixd

      integer(is):: nout1

      real(rs)   :: tout1

!-------------------------------------------------------Local Variables
      character(256)::fich

      character(4),external :: int2ch4

      integer(is):: i, j, k, ierr

      real(rs)   :: g, dvlme, dmas, dmasn, den, tmas, tmasn, ten, &
      ekinn, eintn, tenn, aux

      real(rs), pointer :: pa3(:,:), rhoa3(:,:), vela3(:,:)

!--------------------------------------------------------------------------BEGIN

! mnx,mny vectors
      pa3   => wa400
      rhoa3 => wa401
      vela3 => wa402


      g     = gamma_ad
      tmas  = 0.0
      tmasn = 0.0
      ten   = 0.0
      ekinn = 0.0
      eintn = 0.0
      tenn  = 0.0

      IF (nsdim == 1) THEN
         DO i= 1, nx
            if (igeomx == 2) then
               dvlme = 4.0*pi*(xznr(i)**3-xznl(i)**3)/3.0
            else if (igeomx==1) then
               dvlme = pi*(xznr(i)**2-xznl(i)**2)
            else
               dvlme = xznr(i)-xznl(i)
            endif

            dmas  = physic(i,1,1)%denstw*dvlme
            dmasn = physic(i,1,1)%densty*dvlme
            den   = physic(i,1,1)%energy*dvlme
            tmas  = tmas  + dmas
            tmasn = tmasn + dmasn
            ten   = ten   + den
            ekinn = ekinn + 0.5*dmasn*physic(i,1,1)%velx*physic(i,1,1)%velx
            eintn = eintn + dmasn*physic(i,1,1)%eps
         ENDDO
         tenn = ekinn + eintn

      ELSE
         DO j= 1, ny
            do i= 1, nx
               if (igeomx==2) then
                  dvlme = 2.0*pi*(xznr(i)**3-xznl(i)**3)*     &
                  (-cos(yznr(j))+cos(yznl(j)))/3.0
               else if (igeomx==1) then
                  dvlme = pi*(xznr(i)**2-xznl(i)**2)*(yznr(j)-yznl(j))
               else
                  dvlme = (xznr(i)-xznl(i))*(yznr(j)-yznl(j))
               endif

               dmas  = physic(i,j,1)%denstw*dvlme
               dmasn = physic(i,j,1)%densty*dvlme
               den   = physic(i,j,1)%energy*dvlme
               tmas  = tmas  + dmas
               tmasn = tmasn + dmasn
               ten   = ten   + den
               ekinn = ekinn + 0.5*dmasn*(physic(i,j,1)%velx*physic(i,j,1)%velx &
               + physic(i,j,1)%vely*physic(i,j,1)%vely)
               eintn = eintn + dmasn*physic(i,j,1)%eps
            enddo
         ENDDO
         tenn = ekinn + eintn
      ENDIF


      IF (basenm == 'GSCH') THEN
!Output for SCHN test in 3D along the diagonal
         fich = TRIM(outRoot)//outfil//'_'//int2ch4(cartrank)
         open( 12, file= fich, form= 'formatted', iostat= ierr )
         if (ierr /= 0) then
            call ERROR( 'pltout: Error openning ', fich, __FILE__, __LINE__  )
         endif

         do k= 1, nz
            do j= 1, ny
               do i= 1, nx
                  write(12,'(6(f25.17,1x))')                        &
                  physic(i,j,k)%denstw, physic(i,j,k)%momenx, &
                  physic(i,j,k)%momeny, physic(i,j,k)%momenz, &
                  physic(i,j,k)%energy, physic(i,j,k)%tracer
               enddo
            enddo
         enddo
         close(12)
      ENDIF

      IF (basenm == 'ESSR') THEN
!Output for ESSR test in 3D.
         fich = TRIM(outRoot)//outfil//'_'//int2ch4(cartrank)
         open( 16, file= fich, form= 'formatted', iostat= ierr )
         if (ierr /= 0) then
            call ERROR( 'pltout: Error openning ', fich, __FILE__, __LINE__  )
         endif
         write(16,'(A,I6,A,E12.3,A,E12.3,A,E12.3,A,E12.3,A,E12.3,A,E12.3)')      &
         'n =',nstep,'   time = ', timej, '   dt = ', dt,            &
         '   cmass = ',  tmas, '   nmass = ',tmasn, '   cen = ',ten, &
         '   nen = ', tenn

         do k= 1, nz
            do j= 1, ny
               do i= 1, nx
                  aux = SQRT( physic(i,j,k)%velx*physic(i,j,k)%velx &
                  + physic(i,j,k)%vely*physic(i,j,k)%vely &
                  + physic(i,j,k)%velz*physic(i,j,k)%velz )

                  write(16,'(3(f8.4,1x),2(f30.16,1x),6(f30.16,1x))') &
                  xzn(i), yzn(j), zzn(k), physic(i,j,k)%pres,  &
                  pa3(i,j), physic(i,j,k)%densty, rhoa3(i,j),  &      
                  aux, vela3(i,j)
               enddo
            enddo
         enddo
         close(16)
      ENDIF

      IF (basenm == 'RSSR') THEN
!Output for RSSR test in 3D.
         fich = TRIM(outRoot)//outfil//'_'//int2ch4(cartrank)
         open( 16, file= fich, form= 'formatted', iostat= ierr )
         if (ierr /= 0) then
            call ERROR( 'pltout: Error openning ', fich, __FILE__, __LINE__  )
         endif
         write(16,'(A,I6,A,E12.3,A,E12.3,A,E12.3,A,E12.3,A,E12.3,A,E12.3)')      &
         'n =',nstep,'   time = ', timej, '   dt = ', dt,            &
         '   cmass = ',  tmas, '   nmass = ',tmasn, '   cen = ',ten, &
         '   nen = ', tenn

         do j= 1, ny
            do i= 1, nx
               aux = SQRT( physic(i,j,1)%velx*physic(i,j,1)%velx  &
               + physic(i,j,1)%vely*physic(i,j,1)%vely  &
               + physic(i,j,1)%velz*physic(i,j,1)%velz )

               write(16,'(3(f8.4,1x),2(f30.16,1x),6(f30.16,1x))') &
               xzn(i), yzn(j), physic(i,j,1)%pres,          &
               pa3(i,j),physic(i,j,1)%densty, rhoa3(i,j),   &
               aux, vela3(i,j)
            enddo
         enddo
         close(16)
      ENDIF
      nout1 = 0
      tout1 = 0.0

      call filnam( 'a', basenm, datfil, rstfil, outfil, tdlfil, suffixd, suffix )
!----------------------------------------------------------------------------END
      END subroutine pltout
