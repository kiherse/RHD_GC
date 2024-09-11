!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!     This subroutine provides the boundary conditions.
!     Adapted from Prometheus.
!     MODIFIED 16-JUL-1993 to include special boundary condition for
!     'RCSR' test.
!     MODIFIED 10-FEB-1995 to include sinusoidal variations of the beam
!     flow velocities at the beam boundary in jet simulations.
!     MODIFIED 30-NOV-1995 to eliminate de sinusoidal jet conditions and
!     to include boundary conditions for the tracer
!     MODIFIED 28-FEB-1996 to include exhaustive boundary conditions
!     (Wilson & Falle,  1985,  MNRAS,  216,  971) for jets
!     MODIFIED 20-DEC-1996 to change the perturbation data input (now
!     through common /PERT/)
!===============================================================================
      subroutine bndry( fznp, nznp, j, l, nsdim, bndmin, bndmax, xyswp, ipr,    &
      gamma_ad, timej, mny5, dt, t0, t1, vptbd, rptbd, pptbd, &
      ms, pb, rhob, rhoeb, velb, mm1, mm2, nn, dk1,  dk2,     &
      psh1, psh2, kb, m_0, cl, me, mp, delx, uwi, vwi,        &
      shift, fg, rm, tem0, tem1, rorb, nmod, A0, gridly,      &
      fra, frs, pha, pa, rhoa, xzn, yzn, zzn, p, rho, rhoe,   &
      rhopr, vel, velt, veltt, f, x, xl, xr, dx )

      USE parallel
      USE workarrays
      implicit none
# include "types.h"
!-------------------------------------------------------Input Variables
      integer(is):: fznp, nznp, j, l, mny5

      integer(is):: nsdim, bndmin, bndmax, xyswp, ms, nmod, ipr

      real(rs)   :: gamma_ad, timej, dt, delx,  t0, t1, vptbd, rptbd, pptbd,  & 
      pb, rhob, rhoeb, velb, mm1, mm2, nn, cl, dk1, dk2,        &
      psh1, psh2, me, mp, shift, fg, rm, tem0, tem1, uwi, vwi,  &
      rorb, A0, gridly, kb, m_0, fra(25), frs(25), pha(25)  

      real(rs)   :: pa(-4:mny5), rhoa(-4:mny5)

! sizes mnx1, mny1, mnz1
      real(rs)   :: xzn(*), yzn(*), zzn(*) 

!------------------------------------------------------Output Variables
! size -4:mn5
      real(rs)   :: p(-4:*), rho(-4:*), rhoe(-4:*), vel(-4:*), velt(-4:*), &
      veltt(-4:*),  rhopr(-4:*), f(-4:*)

! size -4:mn5
      real(rs)   :: x(-4:*), xl(-4:*), xr(-4:*), dx(-4:*)

!-------------------------------------------------------Local Variables
      integer(is):: i, ln

      real(rs)   :: vboun, rboun, reboun, pboun, vi, epstmp, radius, gamm1, irho, & 
      tt1, tt2

      real(rs)   :: scr4, scr5, scr6, scr8, scr9, xt, xtt, theta1, theta2, &
      p0, pi, z1, z2, alpha, theta

      real(rs)   :: radi, rhoa3

      real(rs)   :: rho_sy(1),xe_s(1),tem_s(1),   &
      eps_sy(1),p_sy(1),gam_s(1),   &
      h_s(1), cs_s(1),              &
      dpdrho_s(1), dpdeps_s(1),     &
      ze(1),zp(1),meff1(1),dpdrho(1)

!---------------------------------------------------------------------------BEGIN
      IF ((xyswp == 1 .and. mpi_coords(1) == 0) .or. &
      (xyswp == 2 .and. mpi_coords(2) == 0) .or. &
      (xyswp == 3 .and. mpi_coords(3) == 0)) THEN 

! LEFT BOUNDARY
      IF (bndmin == 1) THEN
!Reflecting boundary
         do i= -4, 0
            rho(fznp-1+i)  = rho(fznp-i)
            rhoe(fznp-1+i)  = rhoe(fznp-i)
            rhopr(fznp-1+i) = rho(fznp-1+i) - rhoe(fznp-1+i)
            vel(fznp-1+i)  = -vel(fznp-i)
            velt(fznp-1+i) = velt(fznp-i)
            veltt(fznp-1+i)= veltt(fznp-i)
            p(fznp-1+i)    = p(fznp-i)
            f(fznp-1+i)    = f(fznp-i)
            dx(fznp-1+i)   = dx(fznp-i)
         enddo

         do i=fznp-1, fznp-5, -1
            xl(i) = xl(i+1) - dx(i)
            xr(i) = xr(i+1) - dx(i+1)
            x(i)  = 0.5*(xl(i) + xr(i))
         enddo

      ELSE IF (bndmin == 2) THEN
!Flow-out boundary
         do i= fznp-5, fznp-1
            rho(i)  = rho(fznp)
            rhoe(i)  = rhoe(fznp)
            rhopr(i) = rho(i) - rhoe(i)
            vel(i)  = vel(fznp)
            velt(i) = velt(fznp)
            veltt(i)= veltt(fznp)
            p(i)    = p(fznp)
            f(i)    = f(fznp)
            dx(i)   = dx(fznp)
         enddo

         do i=fznp-1, fznp-5, -1
            xl(i) = xl(i+1) - dx(i)
            xr(i) = xr(i+1) - dx(i+1)
            x(i)  = 0.5*(xl(i) + xr(i))
         enddo

      ELSE IF (bndmin == 4) THEN
!Periodic boundary
         do i= -4, 0
            rho  (fznp-1+i) = rho  (nznp+i)
            rhoe  (fznp-1+i) = rhoe  (nznp+i)
            rhopr(fznp-1+i) = rho  (fznp-1+i) - rhoe (fznp-1+i)
            vel  (fznp-1+i) = vel  (nznp+i)
            velt (fznp-1+i) = velt (nznp+i)
            veltt(fznp-1+i) = veltt(nznp+i)
            p    (fznp-1+i) = p    (nznp+i)
            f    (fznp-1+i) = f    (nznp+i)
            dx   (fznp-1+i) = dx   (nznp+i)
         enddo


      ELSE IF (bndmin == 5) THEN
!Inflow jet boundary + Fixed boundary
         do i= -4, 0

            if ((xzn(l)*xzn(l)+zzn(j)*zzn(j)) <= 1.0) then
               if (timej <= t0 .OR. timej >= t1) then
                  vboun = velb
                  rboun = rhob
                  reboun= rhoeb
                  pboun = pb
               else
                  vboun = vptbd
                  rboun = rptbd
                  reboun= rptbd/rhob * rhoeb
                  pboun = pptbd
               endif
               rho(i)  = rboun
               irho    = 1.0 /rho(i) 
               rhoe(i) = reboun
               rhopr(i)= rho(i) - rhoe(i)
               vel(i)  = vboun

!              Helicoidal PERTURBATIONS of the velocity 
!              (to break simmetry and reproduce observations)
               velt(i) = vel(i)*0.00025*(cos(0.0001*timej)+cos(0.05*timej)+cos(0.1*timej)+cos(0.5*timej))
               veltt(i)= vel(i)*0.00025*(sin(0.0001*timej)+sin(0.05*timej)+sin(0.1*timej)+sin(0.5*timej))
               p(i)    = pboun
               f(i)    = 1.0 !tracer
               dx(i)   = dx(fznp)

            else
                  rho(i)  = rho(fznp)
                  rhoe(i)  = rhoe(fznp)
                  rhopr(i)= rho(i) - rhoe(i)
                  vel(i)  = 0.0
                  velt(i) = 0.0
                  veltt(i)= 0.0
                  p(i)    = p(fznp)
                  f(i)    = 0.0
                  dx(i)   = dx(fznp)
!               endif
            endif
         enddo

      ELSE IF (bndmin == 6) THEN
!     Inflow perturbed jet boundary + Fixed boundary

         pi = acos(-1.)  
         tt2 = sqrt(xzn(l)*xzn(l) + zzn(j)*zzn(j))
         tt1 = tt2**ms
         tt1 = 1.0 / cosh( tt1 )    

         do i=fznp-5,fznp-1
            if (abs(tt1) <= 100.0) then
               gamm1 = 1.0 / (gamma_ad - 1.0)

               rho(i)  = 1.0 - (1.0-rhob)*tt1
               rhoe(i) = 1.0/1836.11 - (1.0/1836.11 -rhoeb)*tt1
               rhopr(i)= rho(i) - rhoe(i)
               irho    = 1.0 / rho(i)          
               vel(i)  = velb*tt1
               velt(i) = 0.0
               veltt(i)= 0.0
               p(i)    = pb
               f(i)    = f(fznp)
  
               if ( abs(tt2) <= 1.0 ) then 
                  p0 = A0*pb 
                  ln = 1

                  if (xzn(l) > 1.e-6) then
                     p(i)= p(i)+p0*cos(frs(ln)*timej                 &
                     +2.0*atan(zzn(j)/xzn(l))+pha(ln)*2.0*pi)  &
                     *(sin(pi*sqrt(xzn(l)*xzn(l)+zzn(j)*zzn(j))))**2
                  else if (abs(xzn(l)) < 1.e-6) then
                     p(i)= p(i)+p0*cos(frs(ln)*timej    &
                     +pi+pha(ln)*2.0*pi)          &
                     *(sin(pi*sqrt(xzn(l)*xzn(l)+zzn(j)*zzn(j))))**2
                  else 
                     p(i)=p(i)+p0*(cos(frs(ln)*timej              &
                     +2.0*(atan(zzn(j)/xzn(l))-pi)           &
                     +pha(ln)*2.0*pi))                       &
                     *(sin(pi*sqrt(xzn(l)*xzn(l)+zzn(j)*zzn(j))))**2
                  endif              
               endif
            else
               rho(i)  = rho(fznp)
               rhoe(i) = rhoe(fznp)
               rhopr(i)= rho(i) - rhoe(i)
               velt(i) = velt(fznp)
               veltt(i)= veltt(fznp)
               p(i)    = p(fznp)
               f(i)    = f(fznp)
               if (vel(fznp) <= 0.0) then
! outflow
                  vel(i) = vel(fznp)
               else
! reflecting
                  vel(i) = -vel(fznp)
               endif
            endif
         enddo
         

      ELSE IF (bndmin == 7) THEN
!     Cartesian relativistic shock reflexion (2D,3D):
!     Analytic conditions

         if (xyswp == 1) then
            xt  = yzn(j)
            xtt = zzn(l)
         else if (xyswp == 2) then
            xt  = zzn(j)
            xtt = xzn(l)
         else if (xyswp == 3) then
            xt  = xzn(j)
            xtt = yzn(l)
         endif

         vi     = 0.9999
         epstmp = 1.0E-5/SQRT(1.0-vi*vi)

         gamm1 = gamma_ad-1.0

         scr4   = vi*(timej-dt)
         scr5   = -vi*xt
         scr6   = -vi*xtt
         scr8   = gamm1*epstmp
         scr9   = 1.0 + epstmp + scr8

         do i= fznp-1, fznp-5, -1
            xl(i) = xl(i+1) - delx
            xr(i) = xr(i+1) -  delx
            x(i)  = 0.5*(xl(i) + xr(i))
         enddo

         gamm1 = (gamma_ad-1.0)*epstmp

         do i= fznp-5, fznp-1
            radius   = 1.0 / SQRT( x(i)*x(i) + xt*xt + xtt*xtt )
            tt1   = vi*radius

!that's the correct density for 2D or 3D
            rho(i)   = (1.0 + tt1*(timej-dt))**(nsdim-1)
            rhoe(i) = rhoeb
            rhopr(i) = rho(i) - rhoe(i)
            irho      = 1.0 / rho(i)
            vel(i)   = -tt1*x(i)
            velt(i)  = -tt1*xt
            veltt(i) = -tt1*xtt
            p(i)     = gamm1*rho(i)
            f(i)     = f(fznp)
            dx(i)    = delx
         enddo

      ELSE IF (bndmin == 8) THEN
!     Inflow perturbed jet boundary + Fixed boundary
         do i= -4, 0

            tt1 = sqrt(xzn(l)*xzn(l) + zzn(j)*zzn(j))

            if (tt1 <= 1.0) then
               gamm1 = 1.0 / (gamma_ad - 1.0)

               rho(i)  = rhob
               rhoe(i) = rhoeb
               rhopr(i) = rho(i) - rhoe(i)
               irho    = 1.0 / rho(i)
               vel(i)  = velb
               velt(i) = 0.0
               veltt(i)= 0.0
               p(i)    = pb
               f(i)    = 1.0

            else
               rho(i)  = 1.0
               rhoe(i) = rhoeb/rhob
               rhopr(i) = rho(i) - rhoe(i)
               theta1 = atan(zzn(j)/       &
               sqrt(shift**2+(rorb-xzn(l))**2))
               theta2 = atan(10.0/(500.0-xzn(l)))
               vel(i)  =  vwi*cos(theta1)*sin(theta2)   
               velt(i) =  vwi*sin(theta1)
               veltt(i)= -vwi*cos(theta1)*cos(theta2)
               p(i)    = (gamma_ad - 1.0)*rho(i)*uwi
               f(i)    = 0.0
            endif

         enddo  

      ELSE IF (bndmin == 9) THEN

!Inflow jet boundary + Fixed boundary                                                                       
         do i= -4, 0

            if ((xzn(l)*xzn(l)+zzn(j)*zzn(j)) <= 1.0) then
               if (timej <= t0 .OR. timej >= t1) then
                  vboun = velb
                  rboun = rhob
                  reboun= rhoeb
                  pboun = pb
               else
                  vboun = vptbd
                  rboun = rptbd
                  reboun= rptbd/rhob * rhoeb
                  pboun = pptbd
               endif

!               gamm1 = 1.0 / (gamma_ad - 1.0)                                                              
               rho(i)  = rboun
               irho    = 1.0 /rho(i)
               rhoe(i) = reboun
               rhopr(i)= rho(i) - rhoe(i)
               vel(i)  = vboun!/(1.+timej-5784.5)**2                                                        
!               velt(i) = 0.0                                                                               
!               veltt(i)= 0.0                                                                               
!               vel(i)  = 0.5/(timej-12000.)**2                                                             

               if (xzn(l) > 0. .and. zzn(j) > 0.) then
                  velt(i) = 0.1/sqrt(2.)*velb*sin(atan(zzn(j)/xzn(l))) !+ &
!                        vel(i)*0.000025*(cos(0.0001*timej)+cos(0.05*timej)+cos(0.1*timej)+cos(0.5*timej))   
                  veltt(i)= 0.1/sqrt(2.)*velb*cos(atan(zzn(j)/xzn(l))) !+ &
!                        vel(i)*0.000025*(sin(0.0001*timej)+sin(0.05*timej)+sin(0.1*timej)+sin(0.5*timej))     
               else if (xzn(l) < 0. .and. zzn(j) > 0.) then
                  velt(i) = 0.1/sqrt(2.)*velb*sin(atan(abs(zzn(j)/xzn(l)))) !+ &                                                                             
!                        vel(i)*0.000025*(cos(0.0001*timej)+cos(0.05*timej)+cos(0.1*timej)+cos(0.5*timej))                                                   
                  veltt(i)= -0.1/sqrt(2.)*velb*cos(atan(abs(zzn(j)/xzn(l)))) !+ &                                                                             
!                        vel(i)*0.000025*(sin(0.0001*timej)+sin(0.05*timej)+sin(0.1*timej)+sin(0.5*timej))      
               else if (xzn(l) < 0. .and. zzn(j) < 0.) then
                  velt(i) = -0.1/sqrt(2.)*velb*sin(atan(abs(zzn(j)/xzn(l)))) !+ &                                                 
!                        vel(i)*0.000025*(cos(0.0001*timej)+cos(0.05*timej)+cos(0.1*timej)+cos(0.5*timej))                                                   
                  veltt(i)= -0.1/sqrt(2.)*velb*cos(atan(abs(zzn(j)/xzn(l)))) !+ &                                                                           
!                        vel(i)*0.000025*(sin(0.0001*timej)+sin(0.05*timej)+sin(0.1*timej)+sin(0.5*timej))    
               else if (xzn(l) > 0. .and. zzn(j) < 0.) then
                  velt(i) = -0.1/sqrt(2.)*velb*sin(atan(abs(zzn(j)/xzn(l)))) !+ &   
!                        vel(i)*0.000025*(cos(0.0001*timej)+cos(0.05*timej)+cos(0.1*timej)+cos(0.5*timej))                                                   
                  veltt(i)= 0.1/sqrt(2.)*velb*cos(atan(abs(zzn(j)/xzn(l)))) !+ &                            
!                        vel(i)*0.000025*(sin(0.0001*timej)+sin(0.05*timej)+sin(0.1*timej)+sin(0.5*timej))                                                    
               endif

               p(i)    = p(fznp)
               f(i)    = f(fznp)
               dx(i)   = dx(fznp)
!               endif                                                                                       
 
            else

               rho(i)  = rho(fznp)
               rhoe(i)  = rhoe(fznp)
               rhopr(i)= rho(i) - rhoe(i)
               vel(i)  = -vel(fznp)
               velt(i) = velt(fznp)
               veltt(i)= veltt(fznp)
               p(i)    = p(fznp)
               f(i)    = f(fznp)
               dx(i)   = dx(fznp)
!               endif                                                                                                                      
            endif
         enddo

!! Wind inflow boundary conditions

!         do i=nznp+1, nznp+5
!            xl(i) = xl(i-1) + dx(nznp)
!            xr(i) = xr(i-1) + dx(nznp)
!            x(i)  = 0.5*(xl(i) + xr(i))
     
 !            radi  = sqrt(xzn(i)**2+(yzn(j)+shift)**2+zzn(k)**2)
 !            rhoa3 = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
 !                       + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)

 !            rho_sy(1)=rhoa3
                   
 !            xe_s(1)  = 1.0/1836.110

 !            if (radi < rm) then 
 !              tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
 !              if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
 !            else
 !                tem_s(1)=tem1*kb/m_0/c**2
 !            endif

          
 !            call eos_synge( 1, tid, 1, 1, me, mp, rho_sy, xe_s, tem_s, eps_sy,   &
 !                          p_sy, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)          

  !          vel(i)  =  vwi
  !          velt(i) =  0.0
  !          veltt(i)=  0.0

  !           idum=-3

  !           rhoco = ran1(idum)*0.2-0.1


  !           physic(i,j,k)%tracer = 0.0

  !          physic(i,j,k)%densty  = rho_sy(1)+rho_sy(1)*rhoco
  !          physic(i,j,k)%denstye = xe_s(1)*physic(i,j,k)%densty
  !          physic(i,j,k)%pres    = p_sy(1)
  !          physic(i,j,k)%opres   = physic(i,j,k)%pres 
            
  !          epsco = p_sy(1)/(gam_s(1)-1.0)/physic(i,j,k)%densty - eps_sy(1)

  !          physic(i,j,k)%eps     = eps_sy(1)+epsco
     
   !         f(i)    =  0.0
   !         dx(i)   = dx(nznp)
   ! enddo



      ELSE IF (bndmin == 10) THEN
!Flow-out boundary with approx... 
         do i= fznp-1, fznp-5, -1
            rho(i)  = (rho(fznp+1)-rho(fznp))/(x(fznp+1)-x(fznp))*(x(i)-x(fznp))+rho(fznp)
!            rhoe(i)  = (rhoe(fznp+1)-rhoe(fznp))/(x(fznp+1)-x(fznp))*(x(i)-x(fznp))+rhoe(fznp) 
            rhoe(i)  = rho(i)/1836.11
            rhopr(i) = rho(i) - rhoe(i)
!            if (vel(fznp) < 0.) then
!               vel(i)  = vel(fznp)
!            else
!              vel(i)  = -vel(fznp)
!            endif
            vel(i) = 0.
            velt(i) = velt(fznp)
            veltt(i)= veltt(fznp)
            p(i)    = (p(fznp+1)-p(fznp))/(x(fznp+1)-x(fznp))*(x(i)-x(fznp))+p(fznp)
            f(i)    = f(fznp)
            dx(i)   = dx(fznp)
         enddo

         do i=fznp-1, fznp-5, -1
            xl(i) = xl(i+1) - dx(i)
            xr(i) = xr(i+1) - dx(i+1)
            x(i)  = 0.5*(xl(i) + xr(i))
         enddo


      ELSE
# ifdef DEBUG
         write(utrace,*) 'Unexpected value for bndmin: ', bndmin
# endif
         if (xyswp == 1) then
            call ERROR( 'bndry','Unexpected x-min boundary', __FILE__, __LINE__ )
         else if (xyswp == 2) then
            call ERROR( 'bndry','Unexpected y-min boundary', __FILE__, __LINE__ )
         else
            call ERROR( 'bndry','Unexpected z-min boundary', __FILE__, __LINE__ )
         endif

      ENDIF


      do i= fznp-1, fznp-5, -1
         xl(i) = xl(i+1) - dx(i)
         xr(i) = xr(i+1) - dx(i+1)
         x(i)  = 0.5*(xl(i) + xr(i))
      enddo




      ENDIF
      

      IF ((xyswp == 1 .and. mpi_coords(1) == mpi_dims(1)-1) .or. &
      (xyswp == 2 .and. mpi_coords(2) == mpi_dims(2)-1) .or. &
      (xyswp == 3 .and. mpi_coords(3) == mpi_dims(3)-1)) THEN 
!RIGHT BOUNDARY
      IF (bndmax == 1) THEN
!Reflecting boundary
         do i= 1, 5
            rho(nznp+i)  =  rho(nznp+1-i)
            rhoe(nznp+i) =  rhoe(nznp+1-i)
            rhopr(nznp+i)=  rho(nznp+i) - rhoe(nznp+i)
            vel(nznp+i)  = -vel(nznp+1-i)
            velt(nznp+i) =  velt(nznp+1-i)
            veltt(nznp+i)=  veltt(nznp+1-i)
            p(nznp+i)    =  p(nznp+1-i)
            f(nznp+i)    =  f(nznp+1-i)
            dx(nznp+i)   =  dx(nznp+1-i)
         enddo

      ELSE IF (bndmax == 2) THEN
!Flow-out boundary
         do i= nznp+1, nznp+5
            rho(i)  = rho(nznp)
            rhoe(i) = rhoe(nznp)
            rhopr(i)= rho(i) - rhoe(i) 
            vel(i)  = vel(nznp)
            velt(i) = velt(nznp)
            veltt(i)= veltt(nznp)
            p(i)    = p(nznp)
            f(i)    = f(nznp)
            dx(i)   = dx(nznp)
         enddo


      ELSE IF (bndmax == 4) THEN
!Periodic boundary
         do i= 1, 5
            rho (nznp+i) = rho (i)
            rhoe (nznp+i) = rhoe (i)
            rhopr(nznp+i) = rho (nznp+i) - rhoe (nznp+i)
            vel (nznp+i) = vel (i)
            velt(nznp+i) = velt(i)
            veltt(nznp+i)= veltt(i)
            p   (nznp+i) = p   (i)
            f   (nznp+i) = f   (i)
            dx  (nznp+i) = dx  (i)
         enddo


      ELSE IF (bndmax == 5) THEN
!     Special boundary Inflow condition for the cylindrical radial
!     shock reflection in 1d
         do i= nznp+1, nznp+5
            dx(i)   = dx(nznp)
            xl(i)   = xl(i-1) + dx(i-1)
            xr(i)   = xr(i-1) + dx(i)
            x(i)    = 0.5*(xl(i)+xr(i))
            rho(i)  = 1.0+ABS(vel(nznp))*timej/x(i)
            rhoe(i) = rho(i) 
            rhopr(i) = rho(i) - rhoe(i)
            vel(i)  = vel(nznp)
            velt(i) = velt(nznp)
            veltt(i)= veltt(nznp)
            p(i)    = p(nznp)
            f(i)    = f(nznp)
         enddo


      ELSE IF (bndmax == 6) THEN
! Special boundary
! Fixed decreasing atmosphere
         gamm1 = 1.0 / (gamma_ad - 1.0)

         do i= 1, 5
            p(nznp+i)    = pb/dk1/(1.0+(yzn(nznp+i)/psh1)**nn)**(mm1/nn)
            rho(nznp+i)  = p(nznp+i)/(pb/dk1)*1.0
            rhoe(nznp+i) = rho(nznp+i)/1836.11
            rhopr(nznp+i) = rho(nznp+i) - rhoe(nznp+i)
            vel(nznp+i)  = 0.0
            velt(nznp+i) = 0.0
            veltt(nznp+i)= 0.0
            f(nznp+i)    = f(nznp)
            dx(nznp+i)    = dx(nznp)
         enddo

      ELSE IF (bndmax == 7) THEN
!Outflow + Exhaustive boundary for jets

!         gamm1 = 1.0 / (gamma_ad - 1.0)

         do i= 1, 5
!            if ((f(nznp) >= 0.01) .OR. (vel(nznp) >= 0.0)) then
               rho(nznp+i)  = rho(nznp)
               rhoe(nznp+i) = rhoe(nznp)
               rhopr(nznp+i)= rho(nznp) - rhoe(nznp)
!               vel(nznp+i)  = vel(nznp)
               vel(nznp+i)  = 0.
               velt(nznp+i) = velt(nznp)
               veltt(nznp+i)= veltt(nznp)
               p(nznp+i)    = p(nznp)
               f(nznp+i)    = f(nznp)
               dx(nznp+i)    = dx(nznp)
!            else   
!               if (ipr == 1) then
!                  p(nznp+i)    = pb/dk1/(1.+yzn(nznp+i))**mm1
!                  rho(nznp+i)  = p(nznp+i)/(pb/dk1)*1.0
!                  rhoe(nznp+i) = rho(nznp+i)/1836.11 
!                  rhopr(nznp+i)= rho(nznp+i) - rhoe(nznp+i)
!                  vel(nznp+i)  = vel(nznp)
!                  velt(nznp+i) = velt(nznp)
!                  veltt(nznp+i)= veltt(nznp)
!                  f(nznp+i)    = f(nznp)
!                  dx(nznp+i)    = dx(nznp)
!               else
!                  p(nznp+i)    = pa(nznp+i)
!                  rho(nznp+i)  = rhoa(nznp+i)
!                  rhoe(nznp+i) = rho(nznp+i)/1836.11 
!                  rhopr(nznp+i)= rho(nznp+i) - rhoe(nznp+i)
!                  vel(nznp+i)  = vel(nznp)
!                  velt(nznp+i) = velt(nznp)
!                  veltt(nznp+i)= veltt(nznp)
!                  f(nznp+i)    = f(nznp)   
!                  dx(nznp+i)    = dx(nznp)
!               endif
!            endif
         enddo

      ELSE IF (bndmax == 8) THEN
!     spherical shock reflection test in spherical or cylindrical
!     coordinates copy left boundary values and reflect velocity
!     this is a trick to get perfect symmetry!

         do i=nznp+1, nznp+5
            rho(i)   =  rho(nznp+1-i)
            rhoe(i)  =  rhoe(nznp+1-i)
            rhopr(i) =  rho(i) - rhoe(i)
            vel(i)   = -vel(nznp+1-i)
            velt(i)  =  velt(nznp+1-i)
            veltt(i) =  veltt(nznp+1-i)
            p(i)     =  p(nznp+1-i)
            f(i)     =  f(nznp+1-i)
            dx(i)    =  delx
            xl(i)    = -xr(nznp+1-i)
            x(i)     = -x(nznp+1-i)
            xr(i)    = -xl(nznp+1-i)            
         enddo

      ELSE IF (bndmax == 9) THEN
!     Wind inflow boundary conditions

         do i=nznp+1, nznp+5
            rho(i)  = 1.0
            rhoe(i) = rho(i)/1836.11
            rhopr(i) = rho(i) - rhoe(i)
            theta1 = atan(zzn(l)/       &
            sqrt((yzn(j)+shift)**2+(rorb-xzn(nznp))**2))
            theta2 = atan((yzn(j)+shift)/(rorb-xzn(nznp)))
            vel(i)  = -vwi*cos(theta1)*cos(theta2)   
            velt(i) =  vwi*cos(theta1)*sin(theta2) 
            veltt(i)=  vwi*sin(theta1)
            p(i)    = (gamma_ad - 1.0)*rho(i)*uwi
            f(i)    =  0.0
            dx(i)   = dx(nznp)
         enddo

      ELSE IF (bndmin == 10) THEN
!Flow-out boundary with approx... 
         do i= nznp+1, nznp+5
            rho(i)  = (rho(nznp)-rho(nznp-1))/(x(nznp)-x(nznp-1))*(x(i)-x(nznp))+rho(nznp)
!            rhoe(i)  = (rhoe(nznp)-rhoe(nznp-1))/(x(nznp)-x(nznp-1))*(x(i)-x(nznp))+rhoe(nznp) 
            rhoe(i)  = rho(i)/1836.11
            rhopr(i) = rho(i) - rhoe(i)
!            if (vel(nznp) > 0.) then
!               vel(i)  = vel(nznp)
!            else
!               vel(i)  = -vel(nznp)
!            endif
            vel(i) = 0.
            velt(i) = velt(nznp)
            veltt(i)= veltt(nznp)
            p(i)    = (p(nznp)-p(nznp-1))/(x(nznp)-x(nznp-1))*(x(i)-x(nznp))+p(nznp)
            f(i)    = f(nznp)
            dx(i)   = dx(nznp)
         enddo


      ENDIF

      do i=nznp+1, nznp+5
         xl(i) = xl(i-1) + dx(i)
         xr(i) = xr(i-1) + dx(i)
         x(i)  = 0.5*(xl(i) + xr(i))
      enddo

      ENDIF

!-----------------------------------------------------------------------------END
      END subroutine bndry

