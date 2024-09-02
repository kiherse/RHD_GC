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

      subroutine init_ext( basenm, nmod, nstep,               &
      nx, ny, nyh, nz, mnx1, mnx5, mny1, mny5, mnz1, mnz5,    &
      mn1, mn5, mn6,                                          &
      nyold, nxtot, nxo, nztot, nzo,                          & 
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

#if defined(PARALELO) 
      USE mpi
#endif

      IMPLICIT NONE

#if defined(PARALELO) && defined(THIN)
      include "mpif.h"
#endif

# include "types.h"

!-------------------------------------------------------Input Variables
      character(4), intent(in):: basenm

      integer(is):: nmod, nstep, nx, ny, nyh, nz, nyt, &
      mnx1, mnx5, mny1, mny5, mnz1, mnz5, ipr, mn1, mn5, mn6, nz2, brt   

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

      integer(is):: i, j, l, k, ierr, nhar,    &
      index1, index2, point, process, ms, ll,  &
      idum, ii, jj, iic, jjc, kkc, &
      nrad, ci, cj, ck, ncl, tid, kk, & 
      inix, fix, iniy, fiy, iniz, fiz, &
      dxoff, dyoff, dzoff,             &
      nyold, nxtot, nxo, nztot, nzo


      integer(is):: ic(2000), jc(2000), kc(2000)

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


!real(rs), pointer :: scrch1(:), scrch2(:), scrch3(:), scrch4(:)
!real(rs), pointer :: xid(:), rad(:), rho2(:), vel2(:)

!     --------- CODI FELIX
      integer                     :: total_grups, restes
      real(rs)                    :: temporal
      type(TplanePhysic), pointer :: xplane(:,:), zplane(:,:)
      integer(is)                 :: lonx(mpi_dims(2)), lonz(mpi_dims(2)), &
      offsetx(mpi_dims(2)), offsetz(mpi_dims(2))
      type(TplanePhysic)          :: xplane0(nz,nyh), zplane0(nx,nyh)
!     --------- FI CODI FELIX

!------------------------------------------------------------------BEGIN

!-4:mn6 vectors
!scrch1 => wa(1)%wa100
!scrch2 => wa(1)%wa101
!scrch3 => wa(1)%wa102
!scrch4 => wa(1)%wa103

! print*,'INITBU, basenm(1:1): ',basenm(1:1)

      allocate (scrch1(-4:mn6), stat= ierr)
      allocate (scrch2(-4:mn6), stat= ierr)
      allocate (scrch3(-4:mn6), stat= ierr)
      allocate (scrch4(-4:mn6), stat= ierr)


!print*,'INITBU, after allocate scrch'

!mn1 vectors
!xid  => wa(1)%wa300
      allocate (xid(mn1), stat= ierr)

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

                               


  ! JET SIMULATIONS
  !    The units of the code are chosen so that
  !       - RB     (beam radius)
  !       - C      (light speed)
  !       - RHOEXT (density of the external medium)
  !    are equal to 1.
  !
  !    This means that the Length Unit is RB,
  !                        Time   Unit is RB/C,
  !                        Mass   Unit is RHOEXT*RB**3.
  !
  !    The beam conditions are determined by
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

          
  IF (basenm(1:1) == 'J' .AND. basenm(2:2) /= '3') THEN
!       print*,'INITBU: dentro del IF, mm1:', mm1


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
!     csb2  = (velb/machb)**2
      csb2 = cs_s(1)*cs_s(1)
      machb = velb/ cs_s(1) 
          
!---  MEDIUM PARAMETERS	---
         rho_sy(1) = 1.0

!    if (mm == 0.0) then
!      do j= 1, ny
!        pa(j)    = pb/dk
!        dpady(j) = 0.0
!      enddo
!    else
!      do j= 1, ny
!         pa(j)    = pb/dk/(1.0 + (yzn(j)/psh)**nn)**(mm/nn)
!         dpady(j) = -(mm*pb/dk/psh**nn)*yzn(j)**(nn-1.0)/  &
!                      (1.0+(yzn(j)/psh)**nn)**(mm/nn+1.0)
!      enddo
!    endif
!    print*,'INITBU: pa(1)',pa(1), 'pa(ny)', pa(ny)
     

      xe_s(1)  = 1.0/1836.110

!      if (radi < rm) then 
!         tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
!         if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
!      else
!         tem_s(1)=tem1*kb/m_0/c**2
!      endif

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


          !isothermal external medium
!          physic(i,j,l)%densty = physic(i,j,l)%pres/(pb/dk) * 1.0
!          physic(i,j,l)%denstye= xe_s(1)*physic(i,j,l)%densty 
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


! print*,'INITBU, BEFORE SETUP'

  IF (basenm(1:1) == 'J' .AND. basenm(2:2) == '3') THEN

    tem_s(1) = temb
    xe_s(1)  = rhoeb/rhob
    rho_sy(1)= rhob      

!   print*,'init_ext',temb,rhob,tem0,tem1


! print*,'INITBU, BEFORE SYNGE'

      call eos_synge( 1, tid, 1, 1, me, mp, rho_sy, xe_s, tem_s, eps_sy, & 
                    p_sy, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)	

! print*,'INITBU, AFTER SYNGE'

            
    temb = tem_s(1)
    epsb = eps_sy(1)
    pb   = p_sy(1)
    hb   = h_s(1)
    gamb = gam_s(1)
!    csb2  = (velb/machb)**2
    csb2 = cs_s(1)*cs_s(1)
    machb = velb/ cs_s(1) 
          
!---  MEDIUM PARAMETERS	---
         rho_sy(1) = 1.0

!         idum = -3
! RANDOM SEED
              dxoff  = mpi_coords(1) * nx      
              dyoff  = mpi_coords(2) * ny       
              dzoff  = mpi_coords(3) * nz     

!              print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),dxoff,dyoff,dzoff

!                  if (dyoff < nyold .and. dyoff+ny > nyold) then
!                    iniy = nyold - dyoff + 1
!                    fiy = ny
!                 else
                    iniy = 1
                    fiy = ny
!                 endif
!                 if (dxoff < (nxtot-nxo)/2 .and. nx > (nxtot-nxo)/2) then
!                    inix = 1
!                    fix = (nxtot-nxo)/2
!                 else if (dxoff < (nxtot+nxo)/2 .and. nx > (nxtot+nxo)/2) then
!                    inix = (nxtot + nxo) /2 +1
!                    fix = nx
!                 else
                    inix = 1
                    fix = nx
!                 endif
!
!                 if (dzoff < (nztot-nzo)/2 .and. nz > (nztot-nzo)/2) then
!                    iniz = 1
!                    fiz = (nztot-nzo)/2
!                 else if (dzoff < (nztot+nzo)/2 .and. nz > (nztot+nzo)/2) then
!                    iniz = (nztot + nzo) /2 + 1
!                    fiz = nz      
!                 else
                    iniz = 1
                    fiz = nz
!                 endif

!!         idum = int(secnds(0.0))*cartrank-1


!                 print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),dxoff,dyoff,dzoff
!                 print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),inix,fix,iniy,fiy,iniz,fiz


        
      if (basenm(3:3) /= 'C') then 

        do k = iniz, fiz
          do j = iniy, fiy
!            idum = int(secnds(0.0)*j)*2-1
            do i = inix, fix

               if (k <=(nztot-nzo)/2 .or. k+dzoff >=(nztot + nzo) /2 + 1 .or.    &
                   j + dyoff > nyold .or.                                        &
                   i <= (nxtot-nxo)/2 .or. i+dxoff >= (nxtot + nxo) /2 + 1) then  

             radi  = sqrt(xzn(i)**2+(yzn(j)+shift)**2+zzn(k)**2)
             rhoa3 = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                        + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)

             rho_sy(1)=rhoa3
                   
             xe_s(1)  = 1.0/1836.110

             if (radi < rm) then 
               tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
               if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
             else
                 tem_s(1)=tem1*kb/m_0/c**2
             endif

!             print*,'INITBU, BEFORE SYNGE'
          

             call eos_synge( 1, tid, 1, 1, me, mp, rho_sy, xe_s, tem_s, eps_sy,   &
                           p_sy, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)          


!             print*,'INITBU, AFTER SYNGE', i,j,k


            physic(i,j,k)%velx = 0.0
            physic(i,j,k)%vely = 0.0
            physic(i,j,k)%velz = 0.0


            scrch1(i)   = 1.0/SQRT(1.0 - physic(i,j,k)%velx   &
                                         * physic(i,j,k)%velx &
                                       - physic(i,j,k)%vely   &
                                         * physic(i,j,k)%vely &
                                       - physic(i,j,k)%velz   &
                                         * physic(i,j,k)%velz)


!             print*,'INITBU, before ran', idum,ran1(idum)
             

             idum = int(secnds(0.0))*2-1


             if (basenm(3:3) /= '0' .and. basenm(3:3) /= '8' ) then                            
                rhoco = ran0(idum)*0.2 - 0.1
                physic(i,j,k)%densty  = rho_sy(1)+rho_sy(1)*rhoco
                physic(i,j,k)%denstye = xe_s(1)*physic(i,j,k)%densty
                physic(i,j,k)%pres    = p_sy(1)
                physic(i,j,k)%opres   = physic(i,j,k)%pres

                epsco = p_sy(1)/(gam_s(1)-1.0)/physic(i,j,k)%densty !- eps_sy(1) 

             else
                rhoco = 0.0
                physic(i,j,k)%densty  = rho_sy(1)
                physic(i,j,k)%denstye = xe_s(1)*physic(i,j,k)%densty
                physic(i,j,k)%pres    = p_sy(1)
                physic(i,j,k)%opres   = physic(i,j,k)%pres

                epsco = eps_sy(1)   

             endif

!             if (i==1) print*,'INITBU, after ran', j,k,rhoco,rho_sy(1),cartrank

!             if 

             physic(i,j,k)%tracer = 0.0

 
!            physic(i,j,k)%eps     = eps_sy(1)+epsco

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

!           print*,'INITBU, after ran2', rhoco,rho_sy(1),cartrank                      
            endif
!           stop
           enddo
         enddo
        enddo
 
      else if (basenm(3:3) == 'C') then

      fich= TRIM(path_wkd)//'/'//'prof3D_D2_den_091_0.d'
      
      open( 1, file=fich, form='formatted', status='old', iostat=ierr )
      if (ierr/=0) then
         call ERROR( 'initbu: Error openning ', fich, __FILE__, __LINE__ )
      endif

      read (1,'(a72)') text
      print*,text
      read (1,'(a72)') text
      print*,text
      read (1,'(a72)') text
      print*,text
      read (1,'(a72)') text
      print*,text
      read (1,'(a72)') text

      do kk = 1, 70
         read(1,'(1x,3(f6.1,3x),1x,I5,2x,I6,2x,2(f9.3,2x))')     &
                 x1(kk), x2(kk), xsh(kk),                      &
                 npsh(kk), npsp(kk), rhosh(kk), rhosp(kk)
         print*,xsh(kk),rhosh(kk)
      enddo   
      close(1)


      fich= TRIM(path_wkd)//'/'//'prof3D_D2_temp_091_0.d'
      
      open( 2, file=fich, form='formatted', status='old', iostat=ierr )
      if (ierr/=0) then
         call ERROR( 'initbu: Error openning ', fich, __FILE__, __LINE__ )
      endif

      read (2,'(a72)') text
      print*,text
      read (2,'(a72)') text
      print*,text
      read (2,'(a72)') text
      print*,text

      do kk = 1, 70
         read(2,'(1x,2(f6.1,3x),f6.1,2x,f6.3,2x,f5.3,1x,f6.3,2x,f5.3)')     &
                 x1(kk), x2(kk), xsh(kk),                      &
                 t1sh(kk), t1sp(kk), tem1sh(kk), tem1sp(kk)
                 print*,xsh(kk),tem1sh(kk)
      enddo   
      close(2)
      do kk=1,70
         print*,xsh(kk),rhosh(kk),tem1sh(kk)
      enddo

      stop
      endif


!      print*,'INITBU',cartrank,'init done'

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




!----------------------------------------------------------------------------END
         END subroutine init_ext


