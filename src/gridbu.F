!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!     This subroutine computes the grid coordinates.
!     Adapted from Prometheus.
!     MODIFIED 16-JUL-1993 to include boundary conditions for test 'RCSR'
!     MODIFIED 08-JUL-1998 to include the generation of a geometrically
!     spanning grid in y-direction
!===============================================================================
      subroutine grid( basenm, igeomx, nx, nx2, ny, nyh, nys, nz, nz2,    &
      nsdim, ixsym, izsym, iextrax, iextraz,             &
      delx, dely, delz, small,                           &
      gridlx, gridlx2, gridly, gridlz, gridlz2,          &
      bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &
      xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, &
      xCPU, yCPU, zCPU )
      USE tipos
      USE parallel
      USE memoria
      USE workarrays

      implicit none
# include "types.h"
!-------------------------------------------------------Input Variables
      character, intent(in):: basenm*(*)


!     ny is the original size of the homogeneous grid
!     nys is the total (adding the extended grid in the last part) size of the grid
!     nyh is the global size
      integer(is):: igeomx, nx, nx2, ny, nyh, nys, nz, nz2,  &
      nsdim, ixsym, izsym, iextrax, iextraz

      real(rs):: delx, dely, delz, small,                 &
      gridlx, gridlx2, gridly, gridlz, gridlz2
      integer :: xCPU, yCPU, zCPU

!------------------------------------------------------Output Variables
      integer(is):: bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz

      real(rs):: xznl(*), xzn(*), xznr(*), &
      yznl(*), yzn(*), yznr(*), &
      zznl(*), zzn(*), zznr(*)

!-------------------------------------------------------Local Variables
      integer(is):: ii, nnx, nnz, cont, limite, proc

      real(rs)   :: lx, ly, lz, fac, ofac, ffac, dfdfac, errfac, facm

!--------------------------------------------------------------------------BEGIN
!     Type of coordinate
!     igeom = 0 ===> planar geometry
!     igeom = 1 ===> cylindrical geometry (radial)
!     igeom = 2 ===> spherical geometry (radial)
!     igeom = 4 ===> spherical geometry (angular-theta)
!     
!     Boundary conditions
!     bndm.. = 1 ===> reflecting boundary (min+max)
!     bndm.. = 2 ===> flow out boundary   (min+max)
!     bndm.. = 3 ===> flow in boundary    (min+max)
!     bndm.. = 4 ===> periodic boundary   (min+max)
!     bndm.. = 5 ===> jet boundary (min); analytical boundary (max)
!     bndm.. = 6 ===> fixed decreasing atmosphere boundary    (max)
!     bndm.. = 7 ===> exhaustive boundary for jets            (max)
!     bndm.. = 8 ===> periodic bndry for half wavelength of antisym. perturbations (min+max)
!     bndm.. = 6 ===> perturbation at the inlet (min)


!Boundary conditions
      if (basenm=='RST_'.OR.basenm=='RBW_'.OR.basenm=='RBWI'.OR. &
      basenm=='RRST'.OR.basenm=='RRBW'.OR.basenm=='RSCH'.OR. &
      basenm=='PRST') then
      bndmnx = 2
      bndmxx = 2
      bndmny = 2
      bndmxy = 2
      if (nsdim==3) then
         bndmnz = 2
         bndmxz = 2
      endif

      else if (basenm=='RSR_') then
         bndmnx = 1
         bndmxx = 2

      else if (basenm=='RCSR') then
         bndmnx = 1
         bndmxx = 5

      else if (basenm=='CRST') then
         bndmnx = 1
         bndmxx = 2
         bndmny = 1
         bndmxy = 2

      else if (basenm=='RSSR') then
         bndmnx = 7
         bndmxx = 8
         bndmny = 7
         bndmxy = 8

      else if (basenm=='GSCH') then
         bndmnx = 2
         bndmxx = 2
         bndmny = 2
         bndmxy = 2
         bndmnz = 2
         bndmxz = 2

      else if (basenm=='ESSR') then
         bndmnx = 7
         bndmxx = 8
         bndmny = 7
         bndmxy = 8
         bndmnz = 7
         bndmxz = 8

      else if (basenm(1:1)=='H') then
         if (nsdim==2) then
            if (ixsym==1) then
               bndmnx = 1
               bndmxx = 2
            else
               bndmnx = 2
               bndmxx = 2
            endif
         else if (nsdim==3) then
            bndmnx = 2
            bndmxx = 2
            bndmnz = 2
            bndmxz = 2
         endif
         bndmny = 4
         bndmxy = 4

      else if (basenm(1:1)=='J'.and. basenm(3:3) == 'O') then
         bndmnx = 2
         bndmxx = 2
         bndmnz = 2
         bndmxz = 2
         bndmny = 9
         bndmxy = 7

      else if (basenm(1:1)=='J'.and. basenm(2:2) /= 'M'.and. basenm(3:3) /= 'O'  &
               .and. basenm(3:3) /= 'C') then
         if (nsdim==2) then
            bndmnx = 1
            bndmxx = 2
         else if (nsdim==3) then
            bndmnx = 2
            bndmxx = 2
            bndmnz = 2
            bndmxz = 2
         endif
         bndmny = 5
         bndmxy = 7

      else if (basenm(1:1)=='J' .and. basenm(2:2) == 'M') then
         bndmnx = 2
         bndmxx = 2
         bndmnz = 2
         bndmxz = 2
         bndmny = 6
         bndmxy = 7

      else if (basenm(1:1)=='J' .and. basenm(3:3) == 'C') then
         bndmnx = 2
         bndmxx = 2
         bndmnz = 2
         bndmxz = 2
         bndmny = 5
         bndmxy = 7

      else if (basenm(1:1)=='M'.and. basenm(2:2) /= 'L') then
         bndmnx = 2
         bndmxx = 9
         bndmnz = 2
         bndmxz = 2
         bndmny = 8
         bndmxy = 2

      else if (basenm(1:1)=='M' .and. basenm(2:2) == 'L') then
         bndmnx = 2
         bndmxx = 2
         bndmnz = 2
         bndmxz = 2
         bndmny = 2
         bndmxy = 2

      else if (basenm == 'SCH3') then
         bndmnx = 2
         bndmxx = 2
         bndmnz = 2
         bndmxz = 2
         bndmny = 2
         bndmxy = 2

      endif

      if ( (bndmnx==4.AND.bndmxx/=4).OR.(bndmny==4.AND.bndmxy/=4).OR.  &
      (bndmnz==4.AND.bndmxz/=4) ) then
      call ERROR( 'grid', 'Periodic boundary not correct', __FILE__, __LINE__ )
      endif

!      if ((basenm(1:1)=='J' .AND. basenm(2:2)/='M') .AND. (bndmny/=5 .and. bndmny /= 9)) then
!         call ERROR( 'grid', 'Jet boundary not correct', __FILE__, __LINE__ )
!      endif

!Set up x-grid
!gridlx is the length of the x-grid
!nx is the number of cells in the x-direction


!     -GHANGES BY 3D PARALLELIZATION
      
      delx    = gridlx/(REAL(nx)*xCPU)
      !print*, 'GRIDBU gridlx, delx, nx: ', gridlx, delx, nx
    
      call global2local( nx, limite, proc, 1 )
      !print*, 'GRIDBU, nx, limite, proc: ', nx, limite, proc


            
      if (mpi_coords(1) < proc) then
         limite = nx
      else if (proc > mpi_coords(1)) then
         limite = -1
      endif
      
      do ii= 1, limite+1
         xznl(ii) = REAL(iniTramo(1,mpi_coords(1))+ii-2)*delx - gridlx/2.
      enddo  

      do ii= 1, limite
         xznr(ii) = xznl(ii+1)
      enddo

      do ii= 1, limite
         xzn(ii) = 0.5*(xznr(ii)+xznl(ii))
      enddo

      ! print*, 'GRIDBU, nx, xzn(1), xzn(nx): ', nx, xzn(1), xzn(nx)
!     - END CHANGES 


!     Add an extra x-grid of length gridlx2 and nx2 cells geometrically
!     spaced.
      IF (iextrax==1) THEN
         nnx    = nx+nx2
         lx     = gridlx2/delx
         ofac   = 1.2
         cont  = 0
         errfac = 1.0
         do while (errfac >= 1.0E-8)
            cont  = cont + 1
            ffac   = ofac**(nx2 + 1)-(lx + 1.0)*ofac + lx
            dfdfac = REAL(nx2 + 1)*ofac**(nx2)-(lx + 1.0)
            fac    = ofac - ffac/dfdfac
            errfac = ABS(1.0-fac/ofac)
            ofac   = MAX(fac,1.0+small)

            if (cont>=10000) then
               call ERROR( 'grid', &
               'Non-convergence in x-grid generation', &
               __FILE__, __LINE__ )
            endif
         enddo

         do ii= nx+2, nnx+1
            xznl(ii)   = xznl(ii-1) + fac*(xznl(ii-1)-xznl(ii-2))
            xznr(ii-1) = xznl(ii)
            xzn(ii-1)  = 0.5*(xznl(ii-1)+xznr(ii-1))
         enddo

         nx     = nnx
         gridlx = gridlx + gridlx2

      ENDIF
!+++  
      
!     Duplicate the x-grid around the x=0 axis in the case of
!     non-symetric problems
!       IF (ixsym/=1) THEN
!          nnx = 2*nx
          
!          do ii= 1, nx+1
!             xznl(nx+ii) = xznl(ii)
!          enddo

!          do ii= 1, nx
!             xznr(nx+ii) = xznl(ii+1)
!             xzn(nx+ii)  = 0.5*(xznl(ii)+xznr(ii))
!          enddo
!
!          do ii= 1, nx
!             xznl(ii) = - xznl(2*nx+2-ii)
!          enddo
!          
!          do ii= 1, nx
!             xznr(ii) = xznl(ii+1)
!             xzn(ii)  = 0.5*(xznl(ii)+xznr(ii))       
!          enddo
!          nx = 2*nx
!
!          
!       ENDIF
      

! gridly is the length of the y-grid
! ny is the number of cells in the y-direction
      IF (nsdim>=2) THEN
!Split y-grid through processors
!Set up y-grid
         dely = gridly/(REAL(nyh)*yCPU)
         
         call global2local( nyh, limite, proc, 2 )
         
         if (mpi_coords(2) < proc) then
            limite = ny
         else if (proc > mpi_coords(2)) then
            limite = -1
         endif
         
         do ii= 1, limite+1
            yznl(ii) = REAL(iniTramo(2,mpi_coords(2))+ii-2)*dely
         enddo  

         do ii= 1, limite
            yznr(ii) = yznl(ii+1)
         enddo

         do ii= 1, limite
            yzn(ii) = 0.5*(yznr(ii)+yznl(ii))
         enddo

      ELSE
         yznl(1) = 0.0
         yznr(1) = 0.0
      ENDIF

! Add an extra y-grid of length gridly/2 and ny/5 cells geometrically
! spaced.
     ! IF (basenm(1:1)=='J'.OR.basenm(1:1)=='3') THEN
      IF (basenm(1:1)=='3') THEN

         call global2local( nyh*yCPU, limite, proc, 2)
         
         
         if (mpi_coords(2) == proc) then
            ly    = (gridly/dely)*0.5
            facm  = ((5*ly+5)/(nyh+5))**(5/nyh)
            ofac  = 1.2
            cont  = 0

            errfac = 1.0
            do while (errfac>=1.0E-8)
               cont  = cont + 1
               ffac   = ofac**(nyh/5 + 1)-(ly + 1.0)*ofac + ly
               dfdfac = (nyh/5 + 1)*ofac**(nyh/5)-(ly + 1.0)
               fac    = ofac - ffac/dfdfac
               errfac = ABS(1.0-fac/ofac)
               ofac   = MAX(fac,facm+small)
               if (cont>=10000) then
                  call ERROR( 'grid', &
                  'Non-convergence in y-grid generation', &
                  __FILE__, __LINE__ )
               endif
            enddo

            do ii= limite+2, ny+1
               yznl(ii)   = yznl(ii-1) + fac*(yznl(ii-1)-yznl(ii-2))
               yznr(ii-1) = yznl(ii)
               yzn(ii-1)  = 0.5*(yznl(ii-1)+yznr(ii-1))

            enddo

            if (mpi_coords(2).ne.mpi_dims(2) - 1) then
               call sendYgrid( proc, ny, fac, yznl )
            endif
            
         elseif (mpi_coords(2) > proc) then
            call recieveYgrid( proc, ny, yznl, yzn, yznr )
         endif

      ENDIF

! Add an extra y-grid of length gridly/3.5 and ny/4 cells geometrically
! spaced.
      IF (basenm(1:1)=='M'.and.basenm(2:2)/='L') THEN
        if (mpi_coords(2).eq.mpi_dims(2) - 1) then

         ly    = (gridly/nbproc/dely)*0.286
         facm  = ((5*ly+5)/(ny+5))**(4/ny)
         ofac  = 1.2
         cont  = 0

         errfac = 1.0
         do while (errfac>=1.0E-8)
            cont  = cont + 1
            ffac   = ofac**(ny/4 + 1)-(ly + 1.0)*ofac + ly
            dfdfac = (ny/4 + 1)*ofac**(ny/4)-(ly + 1.0)
            fac    = ofac - ffac/dfdfac
            errfac = ABS(1.0-fac/ofac)
            ofac   = MAX(fac,facm+small)
            if (cont>=10000) then
               call ERROR( 'grid', 'Non-convergence in y-grid generation', &
               __FILE__, __LINE__ )
            endif
         enddo

         do ii= ny+2, nys+1
            yznl(ii)   = yznl(ii-1) + fac*(yznl(ii-1)-yznl(ii-2))
            yznr(ii-1) = yznl(ii)
            yzn(ii-1)  = 0.5*(yznl(ii-1)+yznr(ii-1))
         enddo


      endif
      ENDIF


!Set up z-grid
      zznl(1) = 0.0
      zznr(1) = 0.0

!gridlz is the length of the z-grid
!nz is the number of cells in the z-direction
      IF (nsdim==3) THEN
         
!     -GHANGES BY 3D PARALLELIZATION
         
         delz    = gridlz/(REAL(nz)*zCPU)
        
         call global2local( nz, limite, proc, 3 )
        
         if (mpi_coords(3) < proc) then
            limite = nz
         else if (proc > mpi_coords(3)) then
            limite = -1
         endif
         
         do ii= 1, limite+1
            zznl(ii) = REAL(iniTramo(3,mpi_coords(3))+ii-2)*delz -gridlz/2.
         enddo  

         do ii= 1, limite
            zznr(ii) = zznl(ii+1)
         enddo
         
         do ii= 1, limite
            zzn(ii) = 0.5*(zznr(ii)+zznl(ii))
         enddo
!print*,'GRID, delz, limite, nz, zznl(1), zznl(limite+1), cartrank', delz, limite, nz, zznl(1), zznl(limite+1), cartrank,  mpi_coords(1), mpi_coords(2), mpi_coords(3)
!call flush(6)
         
!     - END CHANGES 



!     Add an extra z-grid of length gridlz2 and nz2 cells
!     geometrically space d.
         if (iextraz==1) then
            nnz    = nz+nz2
            lz     = gridlz2/delz
            ofac   = 1.2
            cont   = 0
            errfac = 1.0

            do while (errfac>=1.0E-8)
               cont  = cont + 1
               ffac   = ofac**(nz2 + 1)-(lz + 1.0)*ofac + lz
               dfdfac = REAL(nz2 + 1)*ofac**(nz2)-(lz + 1.0)
               fac    = ofac - ffac/dfdfac
               errfac = ABS(1.0-fac/ofac)
               ofac   = MAX(fac,1.0+small)
               if (cont>=10000) then
                  call ERROR( 'grid', 'Non-convergence in z-grid generation', &
                  __FILE__, __LINE__ )
               endif
            enddo

            do ii= nz+2, nnz+1
               zznl(ii)   = zznl(ii-1) + fac*(zznl(ii-1)-zznl(ii-2))
               zznr(ii-1) = zznl(ii)
               zzn(ii-1)  = 0.5*(zznl(ii-1)+zznr(ii-1))
               
            enddo

            nz     = nnz
            gridlz = gridlz + gridlz2
         endif

!Duplicate the z-grid around the z=0 axis
!          if (izsym/=1) then
!             nnz = 2*nz
!             do ii= 1, nz+1
!                zznl(nz+ii) = zznl(ii)
!             enddo
!             do ii= 1, nz
!                zznr(nz+ii) = zznl(ii+1)
!                zzn(nz+ii)  = 0.5*(zznl(ii)+zznr(ii))
!             enddo
!
!             do ii= 1, nz
!                zznl(ii) = - zznl(2*nz+2-ii)
!             enddo
!            
!             do ii= 1, nz
!                zznr(ii) = zznl(ii+1)
!                zzn(ii)  = 0.5*(zznl(ii)+zznr(ii))      
!             enddo
!             nz = 2*nz
!          endif
      ENDIF
!----------------------------------------------------------------------------END
      END subroutine grid
