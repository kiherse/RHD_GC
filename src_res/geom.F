!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!     This subroutine calculates the sections of the cells normal to
!     the direction of the sweep and the cell volumes.
!     Adapted from Prometheus.
!===============================================================================

      subroutine geom( fznp, nznp, igeom, xyswp, xl, xr, dx, xzn, &
      area, areal, arear, dvol )
      implicit none
# include "types.h"

!-------------------------------------------------------Input Variables

      integer(is):: j, fznp, nznp, igeom, xyswp

! size -4:mn5
      real(rs)   :: xl(-4:*), xr(-4:*), dx(-4:*)

! size mnx1
      real(rs)   :: xzn(*)

!------------------------------------------------------Output Variables
! size -4:mn5
      real(rs)   :: area(-4:*), areal(-4:*), arear(-4:*), dvol(-4:*)
      
!-------------------------------------------------------Local Variables

      integer(is):: i, ig

      real(rs)   :: tt

!--------------------------------------------------------------------------BEGIN

      if ((igeom==1 .OR. igeom==2) .AND. xyswp/=1) then
         call ERROR( 'GEOM', 'Incorrect geometry', __FILE__, __LINE__ )
      endif

      ig = igeom+1

      IF (ig == 1) THEN
         do i= fznp, nznp+1
            areal(i) = 1.0
            arear(i) = 1.0
            area(i)  = 1.0
            dvol(i)  = dx(i)
         enddo

      ELSE IF (ig == 2) THEN
         do i= fznp, nznp+1
            areal(i) = ABS(xl(i))
            arear(i) = ABS(xr(i))
            area(i)  = 0.5*(areal(i) + arear(i))
            dvol(i)  = area(i)*dx(i)
         enddo

      ELSE IF (ig == 3) THEN
         tt = 1.0 / 3.0
         do i= fznp, nznp+1
            areal(i) = xl(i)*xl(i)
            arear(i) = xr(i)*xr(i)
            dvol(i)  = (xr(i)*arear(i)-xl(i)*areal(i))*tt
            area(i)  = dvol(i)/dx(i)
         enddo

      ELSE IF (ig == 4) THEN
         call ERROR( 'geom', 'angular cylindrical coordinate not implemented', &
         __FILE__, __LINE__ )

!     ELSE IF (ig == 5) THEN
!     do i= fznp, nznp+1
!     areal(i) = SIN(xl(i))
!     area(i)  = SIN(0.5*(xl(i)+xr(i)))
!     arear(i) = SIN(xr(i))
!     dvol(i)  = area(i)*dx(i)*xzn(j)
!     enddo

      ELSE IF (ig == 6) THEN
         call ERROR( 'geom', 'angular spherical coordinate not implemented', &
         __FILE__, __LINE__ )
      ENDIF

!----------------------------------------------------------------------------END
      END subroutine geom
