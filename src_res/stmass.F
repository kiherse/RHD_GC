!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!===============================================================================
      subroutine stmass( igeomx, igeomy, nx, ny, nz, mnx1, &
      mnx5, mny1, mny5, mnz1, mnz5, alpn,  &
      betn, gamn, rbn, m_0, R_b, c, pc, shift, rho_0a, q0, &
      xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,   &
      physic )
      USE tipos
      USE workarrays
      implicit none
# include "types.h"
!-------------------------------------------------------Input Variables
      integer(is):: igeomx, igeomy, nx, ny, nz, mnx1, mnx5, mny1, mny5, mnz1, mnz5

      real(rs)   :: alpn, betn, gamn, rbn, m_0, &
      R_b, c, pc, rho_0a, q0, shift

! sizes mnx1, mny1, mnz1
      real(rs)   :: xznl(mnx1), xzn(mnx1), xznr(mnx1), &
      yznl(mny1), yzn(mny1), yznr(mny1), &
      zznl(mnz1), zzn(mnz1), zznr(mnz1)

!------------------------------------------------------Output Variables
      type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)

!-------------------------------------------------------Local Variables
      integer(is):: i, j, k
      
      real(rs)   :: radi, rhost, mloss0, mloss  

!-------------------------------------------------------------------------BEGIN
!     The stellar density in the galaxy follows a Nuker profile

!     q0 is the mass-loss of stars in the core (in g/yr/pc^3). It is changed to 
!     (g/s/cm^3)
      q0 = q0/3600.0/24.0/365.25/pc**3

!     mloss0 is the mass-loss per unit time and volume of the code in the core
      mloss0 = q0/rho_0a*R_b/c


!$OMP PARALLEL DEFAULT(none) &
!$OMP PRIVATE(j, i, k, radi) &
!$OMP SHARED(ny, physic, alpn, betn, gamn, rbn) &
!$OMP SHARED(nx, nz, mloss0, shift)&
!$OMP SHARED(xzn, yzn, zzn)
      do k=1, nz     
         do j=1, ny
            do i=1, nx  
!     shift is the shift of the jet injection to the nucleus
!     the stellar density follows a Nuker profile (Lauer et al. 2007)
!     alpn=2.4, betn=1.43, gamn=0.0, rbn=260 pc (numbers from NGC4261)
               radi  = sqrt(xzn(i)**2+(yzn(j)+shift)**2+zzn(k)**2)
               physic(i,j,k)%stems = mloss0*(rbn/radi)**gamn* &
               (1.+(radi/rbn)**alpn)**((gamn-betn)/alpn) 
            enddo
         enddo
      enddo
!$OMP END PARALLEL
!----------------------------------------------------------------------------END
      END subroutine stmass
