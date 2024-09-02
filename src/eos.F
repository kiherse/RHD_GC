!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL Variables
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine eos( m, n, g, u, rho, p, dpdrh, dpdu, h, cs )
  implicit none
# include "types.h"
!------------------------------------------------------- Input Variables
  integer(is):: m, n

  real(rs)   :: g

  ! size -4:mn6
  real(rs)   :: u(-4:*), rho(-4:*)

!------------------------------------------------------ Output Variables
  ! size -4:mn6
  real(rs)   :: p(-4:*), dpdrh(-4:*), dpdu(-4:*), h(-4:*), cs(-4:*)

!------------------------------------------------------- Local Variables
  integer(is):: i

  real(rs)   :: gam1

!------------------------------------------------------------------------- BEGIN
  gam1 = g - 1.0

  do i= m, n
    p(i)     = gam1*rho(i)*u(i)
    dpdrh(i) = gam1*u(i)
    dpdu(i)  = gam1*rho(i)
    h(i)     = 1.0 + u(i) + p(i)/rho(i)
    cs(i)    = SQRT((dpdrh(i) + p(i)*dpdu(i)/rho(i)/rho(i))/h(i))
      if (p(i)<=0.0) then 
         write(*,*) i, p(i), gam1, rho(i), u(i)
         call ERROR( 'eos', 'p(i)<=0.0', __FILE__, __LINE__ )
      endif
  enddo

!----------------------------------------------------------------------------END
END subroutine eos

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL Variables
! DESCRIPTION
! SEE ALSO
!===============================================================================
subroutine eos_synge( input, tid, m, n, me, mp, rho, xe, tem, eps,   & 
                  p, gam, h, cs, dpdrho, dpdeps )

  USE workarrays 
  implicit none
# include "types.h"
!------------------------------------------------------- Input Variables
  integer(is):: input, m, n, tid

  real(rs)   ::  me, mp

  ! size -4:mn6
  real(rs)   :: tem(-4:*), rho(-4:*), xe(-4:*)

!------------------------------------------------------ Output Variables
  ! size -4:mn6
  real(rs)   :: dpdrho(-4:*), dpdeps(-4:*), h(-4:*), cs(-4:*), gam(-4:*), &
                p(-4:*), eps(-4:*)

!------------------------------------------------------- Local Variables
  integer(is):: i, it, itmax

  real(rs)   ::  errlim , errmax

  real(rs), pointer   ::  ze(:), zp(:), meff1(:),              &
                          ge(:), dge(:), gp(:), dgp(:),        &
                          f(:), df(:), oldtem(:),              &
                          relerr(:)

!------------------------------------------------------------------------- BEGIN

!==================== e q u a t i o n   o f   s t a t e =======================

!      for a ideal, nondegenerate electron - positron - proton - gas

!      The state variables are calculated at indices m to n

!      rho    :    density
!      xe     :    electron + positron mass fraction (rho_e- + rho_e+) / rho
!      tem    :    temperature
!      eps    :    internal energy / mass
!      p      :    pressure
!      gam    :    d(ln p) / d(ln rho)  at constant entropy
!      h      :    enthalpy / mass
!      cs     :    speed of sound
!      dpdrho :    dp / drho at constant eps
!      dpdeps :    dp / deps at constant rho

!      units: c = 1, R_b = 1, rho_0a = 1

!      input :  (1)  rho, xe, tem    ==>    eps, p, ...
!               (2)  rho, xe, eps    ==>    tem, p, ...
 
  errlim = 1.e-7
  itmax = 20

  ze     => wa(tid)%wa277
  zp     => wa(tid)%wa278
  meff1  => wa(tid)%wa279               
  ge     => wa(tid)%wa280
  dge    => wa(tid)%wa281
  gp     => wa(tid)%wa282 
  dgp    => wa(tid)%wa283       
  f      => wa(tid)%wa284
  df     => wa(tid)%wa285
  oldtem => wa(tid)%wa286 
  relerr => wa(tid)%wa287



!-----[ rho, xe and tem given ]-------------------------------------------------
      if (input .eq. 1)  then
      do i = m, n
        ze(i) = me / tem(i)
	zp(i) = mp / tem(i)
	meff1(i) = (1.0-xe(i))/mp + xe(i)/me
	dpdrho(i) = tem(i)*meff1(i)
	p(i) = rho(i)*dpdrho(i)
        ge(i)=5.0/2.0/ze(i)+sqrt(9.0/4.0/ze(i)**2+1.0)
        gp(i)=5.0/2.0/zp(i)+sqrt(9.0/4.0/zp(i)**2+1.0)
        dge(i)=-5.0/2.0/ze(i)**2-9.0/4.0/ze(i)**3/ &
               sqrt(9.0/4.0/ze(i)**2+1.0)
        dgp(i)=-5.0/2.0/zp(i)**2-9.0/4.0/zp(i)**3/ &
               sqrt(9.0/4.0/zp(i)**2+1.0)
      enddo
      
      do i = m, n
        eps(i) = (1.0-xe(i))*gp(i) + xe(i)*ge(i) - dpdrho(i) - 1.0
	h(i) = 1.0 + eps(i) + dpdrho(i)
	gam(i) = 1.0 - 1.0/(1.0+( (1.0-xe(i))*dgp(i)*zp(i)**2/mp + & 
                            xe(i)*dge(i)*ze(i)**2/me )/meff1(i))
        dpdeps(i) = rho(i)*(gam(i)-1.0)
	cs(i) = sqrt( gam(i)*p(i) / (h(i)*rho(i)) )
      enddo

      endif
!-------------------------------------------------------------------------------
    
!-----[ rho, xe and eps given ]-------------------------------------------------
!     (tem is used as start value for Newton-Raphson-iteration)
      if (input .eq. 2)  then

      do i = m, n
	meff1(i) = (1-xe(i))/mp + xe(i)/me
      enddo

!     --- start of iteration ---
      do it = 1, itmax      
        do i = m, n
          ze(i) = me / tem(i)
	  zp(i) = mp / tem(i)
	  dpdrho(i) = tem(i)*meff1(i)
          ge(i)=5.0/2.0/ze(i)+sqrt(9.0/4.0/ze(i)**2+1.0)
          gp(i)=5.0/2.0/zp(i)+sqrt(9.0/4.0/zp(i)**2+1.0)
          dge(i)=-5.0/2.0/ze(i)**2-9.0/4.0/ze(i)**3/  &
                  sqrt(9.0/4.0/ze(i)**2+1.0)
          dgp(i)=-5.0/2.0/zp(i)**2-9.0/4.0/zp(i)**3/  &
                  sqrt(9.0/4.0/zp(i)**2+1.0)
       enddo
  
        do i = m, n
	  f(i) = (1.0-xe(i))*gp(i) + xe(i)*ge(i) - dpdrho(i) - 1.0 - eps(i)
	  df(i) = -(1.0-xe(i))*dgp(i)*zp(i)**2/mp -   &
                  xe(i)*dge(i)*ze(i)**2/me - meff1(i)
          oldtem(i) = tem(i)
	  tem(i) = tem(i) - f(i)/df(i)
	  relerr(i) = abs((tem(i)-oldtem(i))/oldtem(i))
	enddo
	do i = m, n
	   if(tem(i).le.0.0) then	      
	      tem(i) = oldtem(i)/2.0
	      relerr(i) = abs((tem(i)-oldtem(i))/oldtem(i))
!	      print *, '::: T<0 at i=', i
	   endif
	enddo
	errmax = 0.0  
        do i = m, n
	  if( relerr(i) .gt. errmax ) then
	    errmax = relerr(i)
	  endif
	enddo
	if( errmax .le. errlim ) then
	  goto 10
	endif
      enddo
 10   continue
!     TODO: gather & scatter after first iteration?
!      print *, "#ITERATIONS:", it, "ERROR:", errmax, "Xe=", xe(1)
!     --- end of iteration ---

      do i = m, n
	p(i) = rho(i)*dpdrho(i)
	h(i) = 1 + eps(i) + dpdrho(i)
	gam(i) = 1 - 1/(1+( (1-xe(i))*dgp(i)*zp(i)**2/mp +    &
                            xe(i)*dge(i)*ze(i)**2/me )/meff1(i))
        dpdeps(i) = rho(i)*(gam(i)-1)
	cs(i) = sqrt( gam(i)*p(i) / (h(i)*rho(i)) )
      enddo

      endif

!----------------------------------------------------------------------------END
END subroutine eos_synge
