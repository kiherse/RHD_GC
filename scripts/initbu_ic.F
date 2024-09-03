!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!     This subroutine defines the initial model
!     MODIFIED 16-JUL-1993 to include the initial model for test 'RCSR'
!     MODIFIED 01-DEC-1995 to include initial conditions for the tracer
!     MODIFIED 21-FEB-1996 to include the input parameters defining the
!     gradient of the ambient pressure in jet
!     simulations
!     MODIFIED 20-DEC-1996 to change initial data input for jet
!     simulations
!     (now defined in routine INPUT_J)
!     MODIFIED 08-JUL-1998 to include new initial models for jet
!     simulations
!===============================================================================

      subroutine init( basenm, nmod, nstep,                   &
      nx, ny, nyh, nz, mnx1, mnx5, mny1, mny5, mnz1, mnz5,    &
      mn1, mn5, mn6, nyold, nxtot, nxo, nztot, nzo,           & 
      dtini, gamma_ad, pb, velb, rhob,                        &
      rhoeb, temb, machb, gamb, dk1, dk2, psh1, psh2,         &
      mm1, mm2, nn, fg, rm, brt, tem0, tem1,                  &
      ipr, oang, ms, gridlx, gridly, timej, dt,               &
      gamma_e, gamma_c, gamma_ce, eps_e, eps_c,               &
      rho_e, rho_c, p_e, p_c, w_c, nu_ce, xd, alpha,          &
      delx, dely, delz, xznl, xzn, xznr, yznl, yzn,           & 
      yznr, zznl, zzn, zznr, pha, frs, fra, rhoa, pa,         &
      dpady, A0, uwi, vwi, vwp, rhowp, uwp, vk, shift,        &
      rorb, rho_0a, R_b, pc, c, kb, m_0, me, mp, q0,          &
      mloss0, ncl, nz2, ic, jc, kc, physic )

      USE constantes
      USE parallel
      USE tipos
      USE workarrays

#if defined(PARALELO) && defined(MN) 
      USE mpi
#endif
      IMPLICIT NONE
#if defined(PARALELO) &&        !defined(MN)
      include "mpif.h"
#endif

# include "types.h"

!-------------------------------------------------------Input Variables

      character(4), intent(in):: basenm

      integer(is):: nmod, nstep, nx, ny, nyh, nz, nyt, &
      mnx1, mnx5, mny1, mny5, mnz1, mnz5, ipr, mn1, mn5, mn6, nz2, brt   

      integer(is):: nyold, nxtot, nxo, nztot, nzo

      real(rs)   :: dtini, pb, velb, rhob, rhoeb, temb, machb, &
                gamb, dk1, dk2, psh1, psh2, mm1, mm2, nn, oang,      &
                fg, rm, tem0, tem1, pc, c, kb, rho_0a, m_0, me, mp,  &
                q0, mloss0, R_b, gamma_ad

      real(rs)   :: gridlx, gridly, timej, dt

      real(rs)   :: gamma_e, gamma_c, gamma_ce, &
                eps_e, eps_c, rho_e, rho_c, &
                p_e, p_c, w_c, nu_ce, xd,   &
                alpha, delx, dely, delz

      real(rs)   :: xznl(mnx1), xzn(mnx1), xznr(mnx1), &
                yznl(mny1), yzn(mny1), yznr(mny1), &
                zznl(mnz1), zzn(mnz1), zznr(mnz1)
      
      real(rs)   :: A0, pha(25), frs(25), fra(25)

!------------------------------------------------------Output Variables

      real(rs)   :: pa(-4:mny5), dpady(-4:mny5), rhoa(-4:mny5)
      
      real(rs)   :: uwi, vwi, shift, rorb, vk, vwp, rhowp, uwp, csb2    

      type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)

!-------------------------------------------------------Local Variables

      character(256)::fich, fich2

      character(4),external :: int2ch4

      integer(is):: i, j, l, k, w, ierr, nhar,    &
      index1, index2, point, process, ms, ll,  &
      idum, ii, jj, iic, jjc, kkc, &
      nrad, ci, cj, ck, ncl, tid, kk

      integer(is):: ic(2000), jc(2000), kc(2000)

      integer(is):: ini, ifi, jni, jfi, kni, kfi

      real(rs):: scr1, scrc1, scr2, scrc2, scr3, scrc3, scr4, scrc4,      &
             scr5, scr6, qscr6, qscr7, qscr8, qr, qvx,                &
             vi, eps_tmp, tmp1, velz_e, velz_c, p0, eps_i, gamm1, tt, &
             s, gamma_L, nu, Gamma_i, Gamma_ex, rho_ir0, rho_er0,     &
             c1, c2, c3,  k0, velx_e, vely_e,                         &
             velx_c, vely_c, tmp, ar, fr, dfdr, erro, mass_c

     real(rs), allocatable :: scrch1(:), scrch2(:), scrch3(:), scrch4(:)
     real(rs), allocatable :: xid(:), rad(:), rho2(:), vel2(:)             


     real(rs):: theta1, theta2, yz, orho, WL, xc, yc, zc, ac, bc, cc, dist, &
              distc, rhocl, delyh, rhot, rangl, sangl, radj, xci, theta1p, &
              theta2p, lfwp, rhowp2, z1, z2, rstar, rr, rhoco, epsco, ran1, &
              ran0

     real(rs):: tem_min, tem_max, hb, epsb, etars, radi,          & 
               rc, tema, teme, he, game, pe, epse, rhoe, rhoa3 

     real(rs), pointer :: rho_sy(:),xe_s(:),tem_s(:),   &
                         eps_sy(:),p_sy(:),gam_s(:),   &
                         h_s(:), cs_s(:),              &
                         dpdrho_s(:), dpdeps_s(:)

     real(rs):: prp, erp, vrp, csp, rcrit

     real(rs):: yzncl(2000), x1(70), x2(70), xsh(70), rhosh(70), &
               rhosp(70), t1sh(70), t1sp(70), tem1sh(70), tem1sp(70)

     integer(is):: npsh(70), npsp(70)

     character(72) :: text

     external ran1, ran0

     real(4), pointer :: scr(:,:,:)

!########################################### TRYIN THINGS!!!!!!!!!

      real(8), dimension(:), allocatable :: ic_radius, ic_density, ic_pressure, ic_spec_int_energy
      real(8), allocatable :: ic_density_interp, ic_spec_int_energy_interp, ic_pressure_interp

      integer :: unit_number, n

      integer :: num_lines, iostat
      character(len=100) :: line

!########################################### TRYIN THINGS!!!!!!!!!

     !--------- CODI FELIX
      integer                     :: total_grups, restes
      real(rs)                    :: temporal
      type(TplanePhysic), pointer :: xplane(:,:), zplane(:,:)
      integer(is)                 :: lonx(mpi_dims(2)), lonz(mpi_dims(2)), &
      offsetx(mpi_dims(2)), offsetz(mpi_dims(2))
      type(TplanePhysic)          :: xplane0(nz,nyh), zplane0(nx,nyh)
     !--------- FI CODI FELIX

!------------------------------------------------------------------BEGIN

      !-4:mn6 vectors
      !scrch1 => wa(1)%wa100
      !scrch2 => wa(1)%wa101
      !scrch3 => wa(1)%wa102
      !scrch4 => wa(1)%wa103
      allocate (scrch1(-4:mn6), stat= ierr)
      allocate (scrch2(-4:mn6), stat= ierr)
      allocate (scrch3(-4:mn6), stat= ierr)
      allocate (scrch4(-4:mn6), stat= ierr)
      !print*,'INITBU, after allocate scrch'

      !mn1 vectors
      !xid  => wa(1)%wa300
      allocate (xid(mn1), stat= ierr)
      !print*,'INITBU, after allocate xid'

      !-4:mn5 vectors
      !rad  => wa(1)%wa201
      !rho2 => wa(1)%wa202
      !vel2 => wa(1)%wa203
      allocate (rad(-4:mn5), stat= ierr)
      allocate (rho2(-4:mn5), stat= ierr)
      allocate (vel2(-4:mn5), stat= ierr)
      !print*,'INITBU, after allocate rad, rho2, vel2'
      
      ! -4:mn5 vectors
      tid = 1

      rho_sy => wa(tid)%wa204
      xe_s   => wa(tid)%wa205
      tem_s  => wa(tid)%wa206
      eps_sy => wa(tid)%wa207
      p_sy   => wa(tid)%wa208
      gam_s  => wa(tid)%wa209 
      h_s    => wa(tid)%wa210
      cs_s   => wa(tid)%wa211
      dpdrho_s => wa(tid)%wa212
      dpdeps_s => wa(tid)%wa213

      ! INITIAL MODEL PARAMETERS

      !Model number: nstep = 0
      !Initial time: timej = 0.0
      !Initial timestep: dt = dtini

      c = 2.99792458e10 !cgs

!---------------------------------------------- 3D SCHNEIDER TEST 

      !print*,'INITBU, before if basenm==SCH3'                                                                                                                
      if (basenm == 'SCH3') then

         gamma_ad = 5.0/3.0 !adiabatic index

         do k= 1, nz
            do i= 1, ny
               do j= 1, nx
                  if (yzn(i) <= gridly/2.0) then
                     physic(j,i,k)%velx   =  0.0
                     physic(j,i,k)%vely   =  0.0
                     physic(j,i,k)%velz   =  0.0
                     physic(j,i,k)%densty = 10.0
                     physic(j,i,k)%denstye = 10.0                     
                     !physic(j,i,k)%eps    =  2.0
                     tem_s(1) = 4.e12*kb/(m_0*c**2)
                     physic(j,i,k)%tracer =  1.0
                     scrch1(j)            =  1.0
                  else
                     physic(j,i,k)%velx   = 0.0
                     physic(j,i,k)%vely   = 0.0
                     physic(j,i,k)%velz   = 0.0
                     physic(j,i,k)%densty = 1.0
                     physic(j,i,k)%denstye= 1.0
                     !physic(j,i,k)%eps    = 1.0E-6
                     tem_s(1) = 2.0e6*kb/(m_0*c**2)                     
                     physic(j,i,k)%tracer = 0.0
                     scrch1(j)            = 1.0
                  endif

                  rho_sy(1) = physic(j,i,k)%densty
                  xe_s(1)  = physic(j,i,k)%denstye/physic(j,i,k)%densty
        
                  call eos_synge( 1, tid, 1, 1, me, mp, rho_sy, xe_s, tem_s, eps_sy,  &
                       p_sy, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)
                                             
                  physic(j,i,k)%pres = p_sy(1)
                
                  physic(j,i,k)%denstye = xe_s(1)*rho_sy(1)
                  physic(j,i,k)%pres    = p_sy(1)
                  physic(j,i,k)%opres   = physic(j,i,k)%pres 
                  physic(j,i,k)%eps     = eps_sy(1)
                  physic(j,i,k)%denstw  = physic(j,i,k)%densty * scrch1(j)
                  physic(j,i,k)%denstew = physic(j,i,k)%denstye * scrch1(j)
                  physic(j,i,k)%momenx  = physic(j,i,k)%denstw * h_s(1) &
                  * scrch1(j) * physic(j,i,k)%velx
                  physic(j,i,k)%momeny  = physic(j,i,k)%denstw * h_s(1) &
                  * scrch1(j) * physic(j,i,k)%vely
                  physic(j,i,k)%momenz  = physic(j,i,k)%denstw * h_s(1) &
                  * scrch1(j) * physic(j,i,k)%velz
                  physic(j,i,k)%energy  = physic(j,i,k)%denstw * h_s(1) &
                  * scrch1(j) - physic(j,i,k)%pres   &
                  - physic(j,i,k)%denstw                                                                        

               enddo
            enddo
         enddo
      endif

!----------------------------------------------JET SIMULATIONS

  !    The units of the code are chosen so that next ones are equal to 1:
  !       - RB     (beam radius)
  !       - C      (light speed)
  !       - RHOEXT (density of the external medium)
  !
  !    Length Unit is RB,
  !    Time   Unit is RB/C,
  !    Mass   Unit is RHOEXT*RB**3.
  !
  !    The beam/jet conditions are determined by
  !       - VELB:  beam velocity           (standard values: 0     - 1   )
  !       - RHOB:  beam density            (                 0.001 - 10  )
  !       - MACHB: beam proper Mach number (                 1     - 10  )
  !       - GAMMA: adiabatic index         (                 1.33  - 1.67)
  !
  !    The external medium is determined by
  !       - DK: beam-to-external pressure ratio (K=1 ==> Pressure matched jets)
  !       - The velocity of the external medium which is fixed to 0.
  !       - The parameters defining the pressure gradient:
  !                 p(y) = pb/dk/[1+(y/PSH)^NN]^(MM/NN)

       
! J*** basenm
      IF (basenm(1:1) == 'J' .AND. basenm(2:2) /= '3') THEN

            tem_s(1) = temb
            xe_s(1)  = rhoeb/rhob
            rho_sy(1)= rhob      

         call eos_synge( 1, tid, 1, 1, me, mp, rho_sy, xe_s, tem_s, eps_sy, & 
                        p_sy, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)	
               
         temb = tem_s(1)
         epsb = eps_sy(1)
         pb   = p_sy(1)
         hb   = h_s(1)
         gamb = gam_s(1)
         ! csb2  = (velb/machb)**2
         csb2 = cs_s(1)*cs_s(1)
         machb = velb/ cs_s(1) 
               
         ! MEDIUM PARAMETERS
         rho_sy(1) = 1.0
         xe_s(1)  = 1.0/1836.110
         tem_s(1) = tem0
               
         call eos_synge( 1, tid, 1, 1, me, mp, rho_sy, xe_s, tem_s, eps_sy,   &
                           p_sy, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)          

         do l= 1, nz
            do j= 1, ny
            do i= 1, nx

               physic(i,j,l)%velx = 0.0
               physic(i,j,l)%vely = 0.0
               physic(i,j,l)%velz = 0.0

               scrch1(i) = 1.0

               physic(i,j,l)%tracer = 0.0

               physic(i,j,l)%densty  = rho_sy(1)
               physic(i,j,l)%denstye = xe_s(1)*rho_sy(1)
               physic(i,j,l)%pres    = p_sy(1)
               physic(i,j,l)%opres   = physic(i,j,l)%pres 
               physic(i,j,l)%eps     = eps_sy(1)
               physic(i,j,l)%denstw  = physic(i,j,l)%densty * scrch1(i)
               physic(i,j,l)%denstew = physic(i,j,l)%denstye * scrch1(i)
               physic(i,j,l)%momenx  = physic(i,j,l)%denstw * h_s(1) &
                                       * scrch1(i) * physic(i,j,k)%velx
               physic(i,j,l)%momeny  = physic(i,j,l)%denstw * h_s(1) &
                                       * scrch1(i) * physic(i,j,k)%vely
               physic(i,j,l)%momenz  = physic(i,j,l)%denstw * h_s(1) &
                                       * scrch1(i) * physic(i,j,k)%velz
               physic(i,j,l)%energy  = physic(i,j,l)%denstw * h_s(1) &
                                       * scrch1(i) - physic(i,j,l)%pres   &
                                       - physic(i,j,l)%denstw        

               !Isothermal external medium
               !physic(i,j,l)%densty = physic(i,j,l)%pres/(pb/dk) * 1.0
               !physic(i,j,l)%denstye= xe_s(1)*physic(i,j,l)%densty

            enddo
            enddo
         enddo
         
         print*,'INITBU: physic(1,1,1)%pres: ', physic(1,1,1)%pres
         print*,'INITBU: physic(nx, ny,nz)%pres: ', physic(nx, ny, nz)%pres
         print*,'INITBU: physic(1,1,1)%densty: ', physic(1,1,1)%densty
         print*,'INITBU: physic(nx, ny,nz)%densty: ', physic(nx, ny, nz)%densty
         print*,'INITBU: physic(1,1,1)%momenx: ', physic(1,1,1)%momenx
         print*,'INITBU: physic(nx, ny,nz)%momenx: ', physic(nx, ny, nz)%momenx
         print*,'INITBU: physic(1,1,1)%momeny: ', physic(1,1,1)%momeny
         print*,'INITBU: physic(nx, ny,nz)%momeny: ', physic(nx, ny, nz)%momeny
         print*,'INITBU: physic(1,1,1)%momenz: ', physic(1,1,1)%momenz
         print*,'INITBU: physic(nx, ny,nz)%momenz: ', physic(nx, ny, nz)%momenz
         print*,'INITBU: physic(1,1,1)%energy: ', physic(1,1,1)%energy
         print*,'INITBU: physic(nx, ny,nz)%energy: ', physic(nx, ny, nz)%energy

         tema= physic(nx/2,1,nz/2)%pres/ (physic(nx/2,1,nz/2)%denstye/me &
         + (physic(nx/2,1,nz/2)%densty-physic(nx/2,1,nz/2)%denstye)/mp )
         etars = (rhob*hb/h_s(1))/(1-velb**2)
         
         if (cartrank == 0) then
            print *, '---------- beam and medium parameters ----------'
            print *, '  vb [c]        = ', velb
            print *, '  Mb            = ', velb/sqrt(csb2)
            print *, '  Tb [K]        = ', temb*m_0*c**2/kb
            print *, '  Tm [K]        = ', tema*m_0*c**2/kb
            print *, '  rhob [g*cm-3] = ', rhob*rho_0a
            print *, '  rhom [g*cm-3] = ', rho_0a
            print *, '  eta           = ', rhob
            print *, '  Xeb           = ', rhoeb/rhob
            print *, '  epsb [c2]     = ', epsb
            print *, '  epsm [c2]     = ', physic(nx/2,1,nz/2)%eps 
            print *, '  pb [c2]       = ', pb
            print *, '  pm [c2]       = ', physic(nx/2,1,nz/2)%pres
            print *, '  pb/pm         = ', pb/physic(nx/2,1,nz/2)%pres
            print *, '  gamb          = ', gamb
            print *, '  gamm          = ', gam_s(1)
            print *, '  vjet1d [c]    = ', velb*sqrt(etars)/    &
                     (1.0+sqrt(etars))
            print *, '------------------------------------------------'
         endif      
         
      ENDIF
      !print*,'INITBU, BEFORE SETUP'

! J3** basenm     
      IF (basenm(1:1) == 'J' .AND. basenm(2:2) == '3') THEN

         tem_s(1) = temb
         xe_s(1)  = rhoeb/rhob
         rho_sy(1)= rhob      

         !print*,'INITBU, BEFORE SYNGE'
            call eos_synge( 1, tid, 1, 1, me, mp, rho_sy, xe_s, tem_s, eps_sy, & 
                        p_sy, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)	
         !print*,'INITBU, AFTER SYNGE'
                  
         temb = tem_s(1)
         epsb = eps_sy(1)
         pb   = p_sy(1)
         hb   = h_s(1)
         gamb = gam_s(1)
         !csb2  = (velb/machb)**2
         csb2 = cs_s(1)*cs_s(1)
         machb = velb/ cs_s(1) 
               
         ! MEDIUM PARAMETERS
                  rho_sy(1) = 1.0
                  !idum = -3
         ! RANDOM SEED
                  !idum = int(secnds(0.0))*cartrank-1

! J300 basenm 
            if (basenm(3:3) /= 'C') then 

            ! Initial conditions profile
            open(unit=unit_number, file='perseus_ics.txt', status='unknown', iostat=ierr)
            n = count_lines(unit_number) - 1

            allocate(ic_radius(n), ic_density(n), ic_pressure(n), ic_spec_int_energy(n))

            read(unit_number,*) ! Skip header (first line)
            do w = 1, n
                  read(unit_number, *) ic_radius(w), ic_density(w), ic_pressure(w), ic_spec_int_energy(w)        
            end do
            close(unit_number)

            do k = 1, nz
               do j = 1, ny
                  do i = 1, nx

                     ! Radi distance of the point in CGS
                     radi  = sqrt(xzn(i)**2+(yzn(j)+shift)**2+zzn(k)**2) * R_b !R_B in cm
                                                                              !########################################### TRYIN THINGS!!!!!!!!!
                     ! INTERPOLATION (also in CGS units)
                     allocate(ic_density_interp, ic_spec_int_energy_interp, ic_pressure_interp)
                     ic_density_interp = interpolate(radi, ic_radius, ic_density)
                     ic_spec_int_energy_interp = interpolate(radi, ic_radius, ic_spec_int_energy)
                     ic_pressure_interp = interpolate(radi, ic_radius, ic_pressure)

                     ! Print test of the interpolation
                     if (i==500 .and. j==500 .and. k==500) then
                        print*,'Original radius', ic_radius(:120)
                        print*,'Original density', ic_density(:120)
                        print*,'Interpolation radius', radi
                        print*,'Interpolated density', ic_density_interp
                        print*,''
                        print*,'x,y,shift,z',xzn(i),yzn(i),shift,zzn(i)
                        print*,'gridlx,gridly,nx,ny',gridlx,gridly,nx,ny
                        print*,'meff',(xe_s(1)/me+(1-xe_s(1))/mp)**(-1)
                     endif

                     ! Interpolated profiles (including change to code units)
                     rhoa3 = ic_density_interp / rho_0a
                     rho_sy(1) = rhoa3     
                     xe_s(1)  = 1.0/1836.110

                     gamma_ad = 5.0/3.0
                     tem_s(1) = (ic_pressure_interp/ic_density_interp*(xe_s(1)/me+(1-xe_s(1))/mp)**(-1)*m_0/kb)*kb/(m_0*c**2)

                     !print*,'INITBU, BEFORE SYNGE',i,j,k
                     call eos_synge( 1, tid, 1, 1, me, mp, rho_sy, xe_s, tem_s, eps_sy,   &
                                    p_sy, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)          
                     !print*,'INITBU, AFTER SYNGE', i,j,k

                     physic(i,j,k)%velx = 0.0
                     physic(i,j,k)%vely = 0.0
                     physic(i,j,k)%velz = 0.0

                     scrch1(i)   = 1.0/SQRT(1.0 - physic(i,j,k)%velx   &
                                                * physic(i,j,k)%velx &
                                                - physic(i,j,k)%vely   &
                                                * physic(i,j,k)%vely &
                                                - physic(i,j,k)%velz   &
                                                * physic(i,j,k)%velz)

                     if (basenm(3:3) /= '0' .and. basenm(3:3) /= '8') then                            
                        rhoco = ran0(idum)*0.2 - 0.1
                     else
                        rhoco = 0.0
                     endif

                     physic(i,j,k)%tracer = 0.0
                     physic(i,j,k)%densty  = rho_sy(1)+rho_sy(1)*rhoco
                     physic(i,j,k)%denstye = xe_s(1)*physic(i,j,k)%densty
                     physic(i,j,k)%pres    = p_sy(1)
                     physic(i,j,k)%opres   = physic(i,j,k)%pres 
                     
                     epsco = p_sy(1)/(gam_s(1)-1.0)/physic(i,j,k)%densty !- eps_sy(1)

                     physic(i,j,k)%eps     = epsco
                     physic(i,j,k)%denstw  = physic(i,j,k)%densty * scrch1(i)
                     physic(i,j,k)%denstew = physic(i,j,k)%denstye * scrch1(i)
                     physic(i,j,k)%momenx  = physic(i,j,k)%denstw * h_s(1) &
                                             * scrch1(i) * physic(i,j,k)%velx
                     physic(i,j,k)%momeny  = physic(i,j,k)%denstw * h_s(1) &
                                             * scrch1(i) * physic(i,j,k)%vely
                     physic(i,j,k)%momenz  = physic(i,j,k)%denstw * h_s(1) &
                                             * scrch1(i) * physic(i,j,k)%velz
                     physic(i,j,k)%energy  = physic(i,j,k)%denstw * h_s(1) &
                                             * scrch1(i) - physic(i,j,k)%pres   &
                                             - physic(i,j,k)%denstw  
                     !print*,'INITBU, after ran2', rhoco,rho_sy(1),cartrank                      
                     !stop

                     deallocate(ic_density_interp, ic_spec_int_energy_interp, ic_pressure_interp)

                  enddo
               enddo
            enddo
      
! J3C* basenm            
            else if (basenm(3:3) == 'C') then
               allocate (scr(512,512,512), stat= ierr)

               ini=mpi_coords(1)*nx
               ifi=ini+nx
               jni=mpi_coords(2)*ny
               jfi=jni+ny
               kni=mpi_coords(3)*nz
               kfi=kni+nz

            if (basenm(4:4) == '0') then   
               fich= '/storage/scratch/lv72/lv72805/D02/dens_D02_091'
            elseif (basenm(4:4) == '1') then
               fich= '/storage/scratch/lv72/lv72805/D10/dens_D10_091'
            elseif (basenm(4:4) == '3') then
               fich= '/storage/scratch/lv72/lv72805/D03/dens_D03_091'
            elseif (basenm(4:4) == '6') then
               fich= '/storage/scratch/lv72/lv72805/D06/dens_D06_091'
            endif

            print*,'hello 1',mpi_coords(1),mpi_coords(2),mpi_coords(3)
            print*,ini,ifi,jni,jfi,kni,kfi

               open( 1, file=fich, form='unformatted', status='old', iostat=ierr )
               if (ierr/=0) then
               call ERROR( 'initbu: Error openning ', fich, __FILE__, __LINE__ )
               endif

               read(1) (((scr(i,j,k), i=1,512), j=1,512), k=1,512)             
               close(1)

               print*,'hello 2'

            do kk=1,nz
               do jj=1,ny
                  do ii=1,nx
                     physic(ii,jj,kk)%densty = scr(ii+ini,jj+jni,kk+kni)*2.*0.72**2 /(3.18)**3/rho_0a*1.e-20
                  enddo
               enddo
            enddo

            print*, rho_0a,physic(1,1,1)%densty, physic(nx,ny,nz)%densty

! J3C0 basenm            
            if (basenm(4:4) == '0') then
               fich= '/storage/scratch/lv72/lv72805/D02/temp_D02_091'
! J3C1 basenm 
            elseif (basenm(4:4) == '1') then
               fich= '/storage/scratch/lv72/lv72805/D10/temp_D10_091'
! J3C3 basenm 
            elseif (basenm(4:4) == '3') then
               fich= '/storage/scratch/lv72/lv72805/D03/temp_D03_091'
! J3C6 basenm 
            elseif (basenm(4:4) == '6') then
               fich= '/storage/scratch/lv72/lv72805/D06/temp_D06_091'
            endif

               open( 2, file=fich, form='unformatted', status='old', iostat=ierr )                                                       
               if (ierr/=0) then
               call ERROR( 'initbu: Error openning ', fich, __FILE__, __LINE__ )                                                      
               endif

               read(2) (((scr(i,j,k), i=1,512), j=1,512), k=1,512)

            close(2)

            print*,'till here works'

            do k=1,nz
               do j=1,ny
                  do i=1,nx
                     tem_s(1) = scr(i+ini,j+jni,k+kni)*kb/m_0/c**2
                  
                     rho_sy(1) = physic(i,j,k)%densty

      !      if (basenm(4:4) == '0') then
      !         fich= '/storage/scratch/lv72/lv72805/D02/pres_D02_091'
      !      elseif (basenm(4:4) == '1') then
      !         fich= '/storage/scratch/lv72/lv72805/D10/pres_D10_091'
      !      endif

      !         open( 2, file=fich, form='unformatted', status='old', iostat=ierr )                                                       
      !         if (ierr/=0) then
      !           call ERROR( 'initbu: Error openning ', fich, __FILE__, __LINE__ )                                        !                
      !         endif

      !         read(2) (((scr(i,j,k), i=1,1024), j=1,1024), k=1,1024)

      !        close(2)

      !        do kk=1,nz
      !           do jj=1,ny
      !              do ii=1,nx
      !                 physic(ii,jj,kk)%pres = scr(ii+ini,jj+jni,kk+kni)
      !              enddo
      !           enddo
      !        enddo


            xe_s(1)  = 1.0/1836.110

                  call eos_synge( 1, tid, 1, 1, me, mp, rho_sy, xe_s, tem_s, eps_sy,   &
                                 p_sy, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)          
                  !print*,'INITBU, AFTER SYNGE', i,j,k

                  physic(i,j,k)%velx = 0.0
                  physic(i,j,k)%vely = 0.0
                  physic(i,j,k)%velz = 0.0

                  scrch1(i)   = 1.0/SQRT(1.0 - physic(i,j,k)%velx   &
                                             * physic(i,j,k)%velx &
                                             - physic(i,j,k)%vely   &
                                             * physic(i,j,k)%vely &
                                             - physic(i,j,k)%velz   &
                                             * physic(i,j,k)%velz)

                  physic(i,j,k)%tracer = 0.0
                  physic(i,j,k)%densty  = rho_sy(1)+rho_sy(1)*rhoco
                  physic(i,j,k)%denstye = xe_s(1)*physic(i,j,k)%densty
                  physic(i,j,k)%pres    = p_sy(1)
                  physic(i,j,k)%opres   = physic(i,j,k)%pres 
                  
                  epsco = p_sy(1)/(gam_s(1)-1.0)/physic(i,j,k)%densty !- eps_sy(1)

                  physic(i,j,k)%eps     = epsco
                  physic(i,j,k)%denstw  = physic(i,j,k)%densty * scrch1(i)
                  physic(i,j,k)%denstew = physic(i,j,k)%denstye * scrch1(i)
                  physic(i,j,k)%momenx  = physic(i,j,k)%denstw * h_s(1) &
                                          * scrch1(i) * physic(i,j,k)%velx
                  physic(i,j,k)%momeny  = physic(i,j,k)%denstw * h_s(1) &
                                          * scrch1(i) * physic(i,j,k)%vely
                  physic(i,j,k)%momenz  = physic(i,j,k)%denstw * h_s(1) &
                                          * scrch1(i) * physic(i,j,k)%velz
                  physic(i,j,k)%energy  = physic(i,j,k)%denstw * h_s(1) &
                                          * scrch1(i) - physic(i,j,k)%pres   &
                                          - physic(i,j,k)%denstw  
                 !print*,'INITBU, after ran2', rhoco,rho_sy(1),cartrank                      
                 !stop

               enddo
               enddo
            enddo
            
            deallocate ( scr)
            !print*,'INITBU',cartrank,'init done'

            endif

            if (cartrank == 0) then
            print *, '---------- beam and medium parameters ----------'
            print *, '  vb [c]        = ', velb
            print *, '  Mb            = ', velb/sqrt(csb2)
            print *, '  Tb [K]        = ', temb*m_0*c**2/kb
            tema= physic(nx/2,1,nz/2)%pres/ (physic(nx/2,1,nz/2)%denstye/me &
                  + (physic(nx/2,1,nz/2)%densty-physic(nx/2,1,nz/2)%denstye)/mp )
            print *, '  Tm [K]        = ', tema*m_0*c**2/kb
            print *, '  rhob [g*cm-3] = ', rhob*rho_0a
            print *, '  rhom [g*cm-3] = ', rho_0a
            print *, '  eta           = ', rhob
            print *, '  Xeb           = ', rhoeb/rhob
            print *, '  epsb [c2]     = ', epsb
            print *, '  epsm [c2]     = ', physic(nx/2,1,nz/2)%eps 
            print *, '  pb [c2]       = ', pb
            print *, '  pm [c2]       = ', physic(nx/2,1,nz/2)%pres
            print *, '  pb/pm         = ', pb/physic(nx/2,1,nz/2)%pres
            print *, '  gamb          = ', gamb
            print *, '  gamm          = ', gam_s(1)
            etars = (rhob*hb/h_s(1))/(1-velb**2)
            print *, '  vjet1d [c]    = ', velb*sqrt(etars)/    &
                        (1.0+sqrt(etars))
            print *, '------------------------------------------------'
            endif

         ENDIF

!************************************************************************************
      
      contains

      FUNCTION interpolate(x, xp, fp) result(fp_interp)

      real(8), intent(in) :: x
      real(8), dimension(:), intent(in) :: xp, fp
      real(8) :: fp_interp
      integer :: j
   
      j = 1
      do while (xp(j) < x .and. j < size(xp))
         j = j + 1
      end do
      if (j == 1) then
         fp_interp = fp(1)
      else if (j == size(xp)) then
         fp_interp = fp(size(xp))
      else
         fp_interp = fp(j-1) + (x - xp(j-1)) * (fp(j) - fp(j-1)) / (xp(j) - xp(j-1))
      end if
      END

      !---------------------------------------------------------------------------------

      FUNCTION count_lines(unit_number) result(num_lines)

      integer, intent(in) :: unit_number
      integer :: num_lines, iostat
      character(len=100) :: line
   
      num_lines = 0
   
      do
         read(unit_number, '(a)', iostat=iostat) line
         if (iostat /= 0) exit
         num_lines = num_lines + 1
      end do
      rewind(unit_number)
      END

!************************************************************************************ 

!----------------------------------------------------------------------------END

      END subroutine init

!***********************************************************************
      FUNCTION ran1(idum)
!***********************************************************************  
         
#ifdef PARALELO
#   include "mpif.h"
#endif

# include "types.h"

      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL(rs) ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
            NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      !      �\200\234Minimal�\200\235 random number generator of Park and Miller with Bays-Durham shuffle and
      !      added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
      !      the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
      !      alter idum between successive deviates in a sequence. RNMX should approximate the largest
      !      floating value that is less than 1.
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then

            idum=max(-idum,1)

            do j=NTAB+8,1,-1
                  k=idum/IQ
                  idum=IA*(idum-k*IQ)-IR*k
                  if (idum.lt.0) idum=idum+IM
                  if (j.le.NTAB) iv(j)=idum
            enddo
            iy=iv(1)
      endif

      k=idum/IQ

      idum=IA*(idum-k*IQ)-IR*k

      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END

!***********************************************************************
      FUNCTION RAN0(idum)
!***********************************************************************
      INTEGER(is) idum,IA,IM,IQ,IR,MASK
      REAL(rs) ran0,AM
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
                 MASK=123459876)
      INTEGER k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      return
      END 
