!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!     CREATED 21-FEB-96 to compute the pseudogravity term balancing the
!     ambient pressure gradient given by
!     p(y) = pb/dk/[1+(y/psh)^nn]^(mm/nn)
!     UPDATED TO 3D 23-FEB-06
!===============================================================================

      subroutine psgrav( basenm, igodu, nx, ny, nz, mnx5, mny5, mnz5, ipr,         &
      nn, psh1, psh2, mm1, mm2, gamma_ad, dk1, dk2, pb,         &
      small, epsiln, omg1, omg2, me, mp,                        &
      fg, shift, rm, kb, m_0, c, tem0, tem1,                    &
      xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,        &
      boundPosX, boundPosY, boundPosZ, rhoa, pa, grav, physic )

      USE tipos
      USE workarrays
      implicit none
# include "types.h"

!-------------------------------------------------------Input Variables

      character, intent(in):: basenm*(*)

      integer(is):: igodu, nx, ny, nz, mnx5, mny5, mnz5, ipr

      real(rs)   :: nn, psh1, psh2, mm1, mm2, gamma_ad, dk1, dk2, pb, &
      small, epsiln, omg1, omg2, me, mp, shift, fg, rm, &
      tem0, tem1, m_0, kb, c, radi, rhoa3

      ! sizes mnx1, mny1, mnz1
      real(rs)   :: xznl(*), xzn(*), xznr(*), &
      yznl(*), yzn(*), yznr(*), &
      zznl(*), zzn(*), zznr(*)

      real(rs) :: rhoa(-4:mny5), pa(-4:mny5)

      type(TboundPosX) :: boundPosX(-4:5)
      type(TboundPosY) :: boundPosY(-4:5)
      type(TboundPosZ) :: boundPosZ(-4:5)

!------------------------------------------------------Output Variables

      real(rs)     :: grav(-4:mny5)
      type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)

!-------------------------------------------------------Local Variables

      integer(is):: tid, j, i, k, fzn, nzn, ini, ifi, jni, jfi, kni, kfi
      integer(is):: iniz, finz, deltax, deltay
      integer(is):: iniy, finy
      integer(is):: ierr
      character(256)::fich

!-------------------------------------------------------OMP Variables

      integer(is) :: nthreads, curth, tmp1, tmp2

      real(rs)   :: tt1, tt2, tt3, tt4, tt5, tt6, tt7, gamm1

      real(rs), pointer :: coeff1(:), coeff2(:), coeff3(:), &
      coeff4(:), coeff5(:)
      real(rs), pointer :: p(:), rho(:), vel(:), velt(:), &
      veltt(:), rhoe(:), f(:)
      real(rs), pointer :: rflx(:), mflx(:), mtflx(:), &
      mttflx(:), eflx(:)
      real(rs), pointer :: pm(:), pp(:), rhom(:), rhop(:), &
      rhoem(:), rhoep(:), velm(:), velp(:), veltm(:),     &
      veltp(:), velttm(:), velttp(:), &
      fm(:), fp(:)
      real(rs), pointer :: delp(:), delrho(:), delvel(:), &
      delvlt(:), delvltt(:), delrhoe(:)
      real(rs), pointer :: dp(:), p6(:), drho(:), rho6(:), &
      drhoe(:), rhoe6(:), dvel(:), vel6(:), dvelt(:),     &
      velt6(:), dveltt(:), veltt6(:)
      real(rs), pointer :: flatn(:), flatn1(:)
      real(rs), pointer :: scrch1(:), scrch2(:)

      real(rs), pointer :: pl(:), pr(:), rhol(:), rhor(:), &
      rhoel(:), rhoer(:),             &
      xel(:), xer(:), teml(:),        &
      temr(:),                        &
      vell(:), velr(:), veltl(:),     &
      veltr(:), velttl(:), velttr(:), &
      ul(:), ur(:), csl(:),           &
      csr(:), hl(:), hr(:),           &
      dpdrhl(:), dpdrhr(:),           &
      dpdul(:), dpdur(:),             &
      gaml(:), gamr(:),               &
      rl(:), rr(:), ml(:),            &
      mr(:), mtl(:), mtr(:),          &
      mttl(:), mttr(:), el(:),        &
      er(:), fl(:), fr(:)

      real(rs), pointer :: x(:), xl(:), xr(:), dx(:)
      real(rs), pointer :: wm(:), wp(:), wl(:), wr(:)
      real(rs), pointer :: w(:), wvel(:), wvelt(:), wveltt(:)

      real(rs)          :: rho_sy(1),xe_s(1),tem_s(1),   &
      eps_sy(1),p_sy(1),gam_s(1),   &
      h_s(1), cs_s(1),              &
      dpdrho_s(1), dpdeps_s(1),     &
      ze(1),zp(1),meff1(1),dpdrho(1)

      real(4), pointer :: scr(:,:,:)

!-------------------------------------------------------------------------BEGIN
      
      nthreads = omp_th

      ! SIMULATIONS ON COSMOLOGICAL BOXES
      if (basenm(3:3) == 'C') then
         ! allocate (scr(512,512,512), stat = ierr)
         ! allocate (scr(1024,1024,1024), stat = ierr)
         allocate (scr(2048,2048,2048), stat = ierr)

         if (basenm(4:4) == '0') then
            ! fich= '/storage/scratch/lv72/lv72805/D02/512_256kpc/temp_D02_091'
            ! fich= '/storage/scratch/lv72/lv72805/D02/1024_256kpc/temp_D02_091'
            fich= '/storage/scratch/lv72/lv72805/D02/2048_256kpc/temp_D02_091'
         elseif (basenm(4:4) == '3') then
            fich= '/storage/scratch/lv72/lv72805/D03/2048_256kpc/temp_D03_091'
         elseif (basenm(4:4) == '6') then
            fich= '/storage/scratch/lv72/lv72805/D06/2048_256kpc/temp_D06_091'
         elseif (basenm(4:4) == '1') then
            fich= '/storage/scratch/lv72/lv72805/D10/2048_256kpc/temp_D10_091'
         endif

         open( 2, file=fich, form='unformatted', status='old', iostat=ierr )                                                                           
         if (ierr/=0) then
            call ERROR( 'initbu: Error openning ', fich, __FILE__, __LINE__ )                                                                           
         endif
         ! read(2) (((scr(i,j,k), i=1,512), j=1,512), k=1,512)
         ! read(2) (((scr(i,j,k), i=1,1024), j=1,1024), k=1,1024)
         read(2) (((scr(i,j,k), i=1,2048), j=1,2048), k=1,2048)
         close(2)

         deltax = 640 ! Cell displacement in the x and z direction
         deltay = 16 ! Cell displacement in the y direction 
         ! equivalent to 2kpc when 4 cells/Rb resolution 

         ini=mpi_coords(1)*nx+deltax
         ifi=ini+nx
         jni=mpi_coords(2)*ny+deltay
         jfi=jni+ny
         kni=mpi_coords(3)*nz+deltax
         kfi=kni+nz
      endif

!$OMP PARALLEL DEFAULT(none) &
!$OMP PRIVATE(tid, j, i, k, fzn, nzn) &
!$OMP PRIVATE(tt1, tt2, tt3, tt4, tt5, tt6, tt7, gamm1) &
!$OMP PRIVATE(coeff1, coeff2, coeff3, coeff4, coeff5) &
!$OMP PRIVATE(p, rho, rhoe, vel, velt, veltt, f) &
!$OMP PRIVATE(rflx, mflx, mtflx, mttflx, eflx) &
!$OMP PRIVATE(pm, pp, rhom, rhop, rhoem, rhoep, velm, velp, veltm, veltp, velttm, velttp) &
!$OMP PRIVATE(fm, fp, delp, delrho, delrhoe, delvel, delvlt, delvltt) &
!$OMP PRIVATE(dp, p6, drho, rho6, drhoe, rhoe6, dvel, vel6, dvelt, velt6, dveltt, veltt6)&
!$OMP PRIVATE(flatn, flatn1, scrch1, scrch2) &
!$OMP PRIVATE(pl, pr, rhol, rhor, rhoel, rhoer, vell, velr, veltl, veltr, velttl, velttr) &
!$OMP PRIVATE(ul, ur, csl, csr,  xel, xer, teml, temr) &
!$OMP PRIVATE(hl, hr, dpdul, dpdur, dpdrhl, dpdrhr, gaml, gamr) &
!$OMP PRIVATE(rl, rr, ml, mr, mtl, mtr, mttl, mttr) &
!$OMP PRIVATE(er, el, fl, fr) &
!$OMP PRIVATE(tmp1, tmp2, nthreads, curth) &
!$OMP PRIVATE(iniz, finz, iniy, finy, ipr, rhoa3) &
!$OMP PRIVATE( x, xl, xr, dx, wm, wp, wl, wr, w, wvel, wvelt, wveltt) &
!$OMP PRIVATE(rho_sy,xe_s,tem_s,eps_sy,p_sy,gam_s,h_s,cs_s,dpdrho_s,dpdeps_s)&
!$OMP PRIVATE(radi,ze,zp,dpdrho,meff1) &
!$OMP SHARED(epsiln, igodu, omg2, omg1, wa, basenm ) &
!$OMP SHARED(nx, ny, nz, small, physic, nn, psh1, grav, mm1, mm2, psh2) &
!$OMP SHARED(boundPosX, boundPosY,boundPosZ ) &
!$OMP SHARED(gamma_ad, dk1, dk2, pb, me, mp, rhoa, pa, cartrank) & 
!$OMP SHARED(xzn, xznl, xznr, yzn, yznl, yznr, zzn, zznl,zznr ) &
!$OMP SHARED(shift,fg,rm,kb,m_0,c,tem0,tem1,ini,jni,kni,scr)&
!$OMP SHARED(mpi_coords, mpi_dims)

# ifdef OPENMP
      call my_omp_get_thread_num(curth)
      tid      = curth+1

      iniz = 1
      finz = nz
      iniy = 1
      finy = ny

# else
      tid  = 1
      iniz = 1
      finz = nz
      iniy = 1
      finy = ny
# endif

! -4:mn6 vectors
      rflx   => wa(tid)%wa100
      mflx   => wa(tid)%wa101
      mtflx  => wa(tid)%wa102
      mttflx => wa(tid)%wa103
      eflx   => wa(tid)%wa104

      scrch1 => wa(tid)%wa105
      scrch2 => wa(tid)%wa106

      pl     => wa(tid)%wa107
      pr     => wa(tid)%wa108
      rhol   => wa(tid)%wa109
      rhor   => wa(tid)%wa110
      rhoel  => wa(tid)%wa188
      rhoer  => wa(tid)%wa189
      vell   => wa(tid)%wa111
      velr   => wa(tid)%wa112
      veltl  => wa(tid)%wa113
      veltr  => wa(tid)%wa114
      velttl => wa(tid)%wa115
      velttr => wa(tid)%wa116
      ul     => wa(tid)%wa117
      ur     => wa(tid)%wa118
      csl    => wa(tid)%wa119
      csr    => wa(tid)%wa120
      hl     => wa(tid)%wa190
      hr     => wa(tid)%wa191  
      dpdrhl => wa(tid)%wa192
      dpdrhr => wa(tid)%wa193        
      dpdul  => wa(tid)%wa194
      dpdur  => wa(tid)%wa195          
      gaml   => wa(tid)%wa196
      gamr   => wa(tid)%wa197
      rl     => wa(tid)%wa121
      rr     => wa(tid)%wa122
      ml     => wa(tid)%wa123
      mr     => wa(tid)%wa124
      mtl    => wa(tid)%wa125
      mtr    => wa(tid)%wa126
      mttl   => wa(tid)%wa127
      mttr   => wa(tid)%wa128
      el     => wa(tid)%wa129
      er     => wa(tid)%wa130
      fl     => wa(tid)%wa131
      fr     => wa(tid)%wa132

      wm     => wa(tid)%wa133
      wp     => wa(tid)%wa134
      wl     => wa(tid)%wa135
      wr     => wa(tid)%wa136
      xel    => wa(tid)%wa198
      xer    => wa(tid)%wa199
      teml   => wa(tid)%wa1101
      temr   => wa(tid)%wa1102

! -4:mn5 vectors
      coeff1 => wa(tid)%wa200
      coeff2 => wa(tid)%wa201
      coeff3 => wa(tid)%wa202
      coeff4 => wa(tid)%wa203
      coeff5 => wa(tid)%wa204

      p     => wa(tid)%wa205
      rho   => wa(tid)%wa206
      rhoe  => wa(tid)%wa268
      vel   => wa(tid)%wa207
      velt  => wa(tid)%wa208
      veltt => wa(tid)%wa209
      f     => wa(tid)%wa210

      pm     => wa(tid)%wa211
      pp     => wa(tid)%wa212
      rhom   => wa(tid)%wa213
      rhop   => wa(tid)%wa214
      rhoem  => wa(tid)%wa269
      rhoep  => wa(tid)%wa270
      velm   => wa(tid)%wa215
      velp   => wa(tid)%wa216
      veltm  => wa(tid)%wa217
      veltp  => wa(tid)%wa218
      velttm => wa(tid)%wa219
      velttp => wa(tid)%wa220
      fm     => wa(tid)%wa221
      fp     => wa(tid)%wa222

      delp    => wa(tid)%wa223
      delrho  => wa(tid)%wa224
      delrhoe => wa(tid)%wa271
      delvel  => wa(tid)%wa225
      delvlt  => wa(tid)%wa226
      delvltt => wa(tid)%wa227

      dp     => wa(tid)%wa228
      p6     => wa(tid)%wa229
      drho   => wa(tid)%wa230
      rho6   => wa(tid)%wa231
      drhoe  => wa(tid)%wa272
      rhoe6  => wa(tid)%wa273
      dvel   => wa(tid)%wa232
      vel6   => wa(tid)%wa233
      dvelt  => wa(tid)%wa234
      velt6  => wa(tid)%wa235
      dveltt => wa(tid)%wa236
      veltt6 => wa(tid)%wa237

      flatn  => wa(tid)%wa238
      flatn1 => wa(tid)%wa239

      x  => wa(tid)%wa240
      xl => wa(tid)%wa241
      xr => wa(tid)%wa242
      dx => wa(tid)%wa243

      w      => wa(tid)%wa244
      wvel   => wa(tid)%wa245
      wvelt  => wa(tid)%wa246
      wveltt => wa(tid)%wa247

!-------------------------------------------------------------------------

      IF (basenm(1:1) == 'J'.and.basenm(2:2) == '3') THEN

      !---------------------------------------------------- FIRST SWEEP: GRAVX. 

         fzn = 1
         nzn = nx
         
!$OMP DO SCHEDULE(static)
         DO k = iniz, finz
            do i = 1, ny

               do j= fzn, nzn
                  x(j)    = xzn(j)
                  xl(j)   = xznl(j)
                  xr(j)   = xznr(j)
                  dx(j)   = xr(j) - xl(j)

                  if (basenm(3:3) == 'C') then
                     rho(j) = physic(j,i,k)%densty
                     tem_s(1) = scr(j+ini,i+jni,k+kni)*kb/m_0/c**2
                     rhoe(j) =  rho(j)/1836.110
                                    
                  else 
                     radi  = sqrt(xzn(j)**2+(yzn(i)+shift)**2+zzn(k)**2)
                     rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                        + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)

                     rhoe(j) =  rho(j)/1836.110

                     if (radi < rm) then 
                        tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                        if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                     else
                        tem_s(1)=tem1*kb/m_0/c**2
                     endif
                  endif

                  ze = me / tem_s(1)
                  zp = mp / tem_s(1)
                  meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                  dpdrho(1) = tem_s(1)*meff1(1)
                  p(j) = rho(j)*dpdrho(1)
                 
                  vel(j)  = 0.0
                  velt(j) = 0.0
                  veltt(j)= 0.0
               enddo

               if (mpi_coords(1) /= 0) then
                  do j= fzn-1, fzn-5, -1
                     xl(j)   = boundPosX(j)%xznl
                     x(j)    = boundPosX(j)%xzn
                     xr(j)   = boundPosX(j)%xznr
                     dx(j)   = boundPosX(j)%dx

                     if (basenm(3:3) == 'C') then
                        rho(j) = physic(j,i,k)%densty
                        tem_s(1) = scr(j+ini,i+jni,k+kni)*kb/m_0/c**2
                        rhoe(j) =  rho(j)/1836.110
                                 
                     else 
                        radi  = sqrt(xzn(j)**2+(yzn(i)+shift)**2+zzn(k)**2)
                        rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                           + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)

                        rhoe(j) =  rho(j)/1836.110

                        if (radi < rm) then 
                           tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                           if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                        else
                           tem_s(1)=tem1*kb/m_0/c**2
                        endif
                     endif

                     ze = me / tem_s(1)
                     zp = mp / tem_s(1)
                     meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                     dpdrho(1) = tem_s(1)*meff1(1)
                     p(j) = rho(j)*dpdrho(1)
                  
                     vel(j)  = 0.0
                     velt(j) = 0.0
                     veltt(j)= 0.0
                  enddo       
               endif

               if  (mpi_coords(1) /= mpi_dims(1)-1) then
                  do j= nzn+1, nzn+5
                     xl(j)   = boundPosX(j-nzn)%xznl
                     x(j)    = boundPosX(j-nzn)%xzn
                     xr(j)   = boundPosX(j-nzn)%xznr
                     dx(j)   = boundPosX(j-nzn)%dx

                     if (basenm(3:3) == 'C') then
                        rho(j) = physic(j,i,k)%densty
                        tem_s(1) = scr(j+ini,i+jni,k+kni)*kb/m_0/c**2
                        rhoe(j) =  rho(j)/1836.110
                                       
                     else
                        radi  = sqrt(xzn(j)**2+(yzn(i)+shift)**2+zzn(k)**2)
                        rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                           + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)

                        rhoe(j) =  rho(j)/1836.110

                        if (radi < rm) then 
                           tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                           if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                        else
                           tem_s(1)=tem1*kb/m_0/c**2
                        endif
                     endif

                     ze = me / tem_s(1)
                     zp = mp / tem_s(1)
                     meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                     dpdrho(1) = tem_s(1)*meff1(1)
                     p(j) = rho(j)*dpdrho(1)
                  
                     vel(j)  = 0.0
                     velt(j) = 0.0
                     veltt(j)= 0.0
                  enddo
               endif

               if (mpi_coords(1) == 0) then
                  do j= fzn-1, fzn-5, -1
                     ! GHANGES BY 3D PARALLELIZATION: copy X boundaries 
                     dx(j)   = dx(fzn)
                     xl(j)   = xl(j+1)-dx(j)
                     xr(j)   = xr(j+1)-dx(j+1)
                     x(j)    = 0.5*(xl(j)+xr(j))

                     if (basenm(3:3) == 'C') then                                               
                        rho(j)  = rho(fzn) ! This is what we should change if we want to be more correct with the bndrybu
                        rhoe(j) =  rho(j)/1836.110
                        tem_s(1) = scr(fzn+ini,i+jni,k+kni)
                        tem_s(1) =  tem_s(1) *kb/m_0/c**2

                     else
                        radi  = sqrt(x(j)**2+(yzn(i)+shift)**2+zzn(k)**2)
                        rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                             + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)
                        rhoe(j) =  rho(j)/1836.110

                        if (radi < rm) then 
                           tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                           if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                        else
                           tem_s(1)=tem1*kb/m_0/c**2
                        endif
                     endif

                     ze = me / tem_s(1)
                     zp = mp / tem_s(1)
                     meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                     dpdrho(1) = tem_s(1)*meff1(1)
                     p(j) = rho(j)*dpdrho(1)

                     vel(j)  = 0.0
                     velt(j) = 0.0
                     veltt(j)= 0.0
                  enddo
               endif
               
               if  (mpi_coords(1) == mpi_dims(1)-1) then
                  do j= nzn+1, nzn+5
                     ! GHANGES BY 3D PARALLELIZATION: copy X boundaries
                     dx(j)   = dx(nzn)
                     xl(j)   = xl(j-1) + dx(j-1)
                     xr(j)   = xr(j-1) + dx(j)
                     x(j)    = 0.5*(xl(j)+xr(j))
                 
                     if (basenm(3:3) == 'C') then                                         
                        rho(j)  =  rho(nzn)
                        rhoe(j) =  rho(j)/1836.110
                        tem_s(1)= scr(nzn+ini,i+jni,k+kni)
                        tem_s(1) =  tem_s(1) *kb/m_0/c**2

                     else
                        radi  = sqrt(x(j)**2+(yzn(i)+shift)**2+zzn(k)**2)
                        rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                             + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)
                        rhoe(j) =  rho(j)/1836.110

                        if (radi < rm) then 
                           tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                           if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                        else
                           tem_s(1)=tem1*kb/m_0/c**2
                        endif
                     endif
                   
                     ze = me / tem_s(1)
                     zp = mp / tem_s(1)
                     meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                     dpdrho(1) = tem_s(1)*meff1(1)
                     p(j) = rho(j)*dpdrho(1)
                     vel(j)  = 0.0
                     velt(j) = 0.0
                     veltt(j)= 0.0
                  enddo
               endif

               call coef( tid, fzn, nzn, dx, coeff1, coeff2, coeff3, coeff4, coeff5 )

               call interp( tid, fzn, nzn, small,                    &
               coeff1, coeff2, coeff3, coeff4, coeff5,  &
               p, pm, pp, delp )

               call interp( tid, fzn, nzn, small,                    &
               coeff1, coeff2, coeff3, coeff4, coeff5,  &
               rho, rhom, rhop, delrho )

               call interp( tid, fzn, nzn, small,                    &
               coeff1, coeff2, coeff3, coeff4, coeff5,  &
               rhoe, rhoem, rhoep, delrhoe )

               call flaten( tid, fzn, nzn, igodu, epsiln, omg1, omg2,&
               p, vel, flatn, flatn1 ) 

               do j= fzn-1, nzn+1
                  tt1 = flatn(j)*rho(j)
                  tt7 = flatn(j)*rhoe(j)
                  tt5 = flatn(j)*p(j)
               
                  rhom(j)  = tt1 + flatn1(j)*rhom(j)
                  rhop(j)  = tt1 + flatn1(j)*rhop(j)
                  rhoem(j) = tt7 + flatn1(j)*rhoem(j)
                  rhoep(j) = tt7 + flatn1(j)*rhoep(j)  
                  pm(j)    = tt5 + flatn1(j)*pm(j)
                  pp(j)    = tt5 + flatn1(j)*pp(j)
               enddo

               call monot( tid, fzn, nzn, small, p,    pm,     pp,     dp,     p6     )
               call monot( tid, fzn, nzn, small, rho,  rhom,   rhop,   drho,   rho6   )
               call monot( tid, fzn, nzn, small, rhoe, rhoem,  rhoep,  drhoe,  rhoe6  )

               do j= fzn, nzn+1

                  rhol(j)  = rhop(j-1)
                  rhoel(j) = rhoep(j-1)
                  pl(j)    = pp(j-1)
                  vell(j)  = 0.
                  veltl(j) = 0.
                  velttl(j)= 0.
                  wl(j)    = 1.0
                  xel(j)   = rhoel(j)/rhol(j)
                  teml(j)  = pl(j) / (rhol(j)*( xel(j)/me + (1.0-xel(j))/mp )) 


                  rhor(j)  = rhom(j)
                  rhoer(j) = rhoem(j)
                  pr(j)    = pm(j)
                  velr(j)  = 0.
                  veltr(j) = 0.
                  velttr(j)= 0.
                  wr(j)    = 1.
                  xer(j)   = rhoer(j)/rhor(j)
                  temr(j)  = pr(j) / (rhor(j)*( xer(j)/me + (1.0-xer(j))/mp )) 

                  ! if (mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1 .and. &
                  ! j == 324 .and. i == 217 .and. k == 217) then
                  !    print*,'------------------------------ Psgrav x: before Synge -----------------------------'
                  !    print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell ',j,i,k
                  !    print*,'me, mp, rhol, xel, teml, ul: ',me, mp, rhol(j), xel(j), teml(j), ul(j)
                  !    print*,'pl, gaml, hl, csl, dpdrhl, dpdul: ',pl(j), gaml(j), hl(j), csl(j), dpdrhl(j), dpdul(j)                                   
                  !    print*,'-----------------------------------------------------------------------------------------'
                  ! endif

               enddo

               call eos_synge( 1, tid, fzn, nzn+1, me, mp, rhol, xel, teml, ul,   &
                                 pl, gaml, hl, csl, dpdrhl, dpdul)	
               call eos_synge( 1, tid, fzn, nzn+1, me, mp, rhor, xer, temr, ur,   &
                                 pr, gamr, hr, csr, dpdrhr, dpdur)

               do j= fzn, nzn+1  
                  rl(j)    = rhol(j)*wl(j)
                  tt2      = rl(j)*hl(j)*wl(j)
                  ml(j)    = tt2*vell(j)
                  mtl(j)   = tt2*veltl(j)
                  el(j)    = tt2 - rl(j) - pl(j)
                  fl(j)    = fp(j-1)

                  rr(j)    = rhor(j)*wr(j)
                  tt4      = rr(j)*hr(j)*wr(j)
                  mr(j)    = tt4*velr(j)
                  mtr(j)   = tt4*veltr(j)
                  er(j)    = tt4 - rr(j) - pr(j)
                  fr(j)    = fm(j)

                  ! if (mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1 .and. &
                  ! j == 324 .and. i == 217 .and. k == 217) then
                  !    print*,'------------------------------ Psgrav x: Riemann parameters -----------------------------'
                  !    print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell ',j,i,k
                  !    print*,'vell, velr, veltl, veltr, velttl, velttr, wl, wr: ', vell(j), velr(j), veltl(j), veltr(j), velttl(j), velttr(j), wl(j), wr(j) 
                  !    print*,'csl, csr, pl, pr, hl, hr, dpdrhl, dpdrhr, dpdul: ',csl(j), csr(j), pl(j), pr(j), hl(j), hr(j), dpdrhl(j), dpdrhr(j), dpdul(j)
                  !    print*,'dpdur, rl, rr, ml, mr, mtl, mtr: ',dpdur(j), rl(j), rr(j), ml(j), mr(j), mtl(j), mtr(j)
                  !    print*,'mttl, mttr, el, er: ',mttl(j), mttr(j), el(j), er(j)                                       
                  !    print*,'-----------------------------------------------------------------------------------------'
                  ! endif

               enddo

               ! Approximate Riemann solver
               call riem3dp( tid, fzn, nzn+1, rhol, rhor, ul, ur,              &
                              vell, velr, veltl, veltr, velttl, velttr, wl, wr, &
                              csl, csr, pl, pr, hl, hr, dpdrhl, dpdrhr, dpdul,  &
                              dpdur, rl, rr, ml, mr, mtl, mtr,                  &
                              mttl, mttr, el, er,                               &
                              rflx, mflx, mtflx, mttflx, eflx )

               do j= fzn, nzn
               physic(j,i,k)%gravx = (mflx(j+1) - mflx(j))/abs(dx(j))

               ! if (mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1 .and. &
               ! j == 324 .and. i == 217 .and. k == 217) then
               !    print*,'---------- Psgrav x: after Riemann ------------'
               !    print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell ',j,i,k
               !    print*,'Fluxes j+1 and j: ',mflx(j+1),mflx(j) 
               !    print*,'dx(j): ',dx(j)   
               !    print*,'Gravx: ',physic(j,i,k)%gravx                                              
               !    print*,'----------------------------------------------------------'
               ! endif

               enddo

            enddo
         ENDDO

!$OMP END DO
      ENDIF
!$OMP BARRIER

! Second SWEEP. GRAVY
      fzn = 1
      nzn = ny
!$OMP DO SCHEDULE(static)
      DO k = iniz, finz
         do i = 1, nx

            do j= fzn, nzn
               x(j)    = yzn(j)
               xl(j)   = yznl(j)
               xr(j)   = yznr(j)
               dx(j)   = xr(j) - xl(j)
               vel(j)  = 0.0
               velt(j) = 0.0
               veltt(j)= 0.0
            enddo

            if (mpi_coords(2) /= 0) then
               do j= fzn-1, fzn-5, -1 
                  dx(j)   = boundPosY(j)%dx
                  xl(j)   = boundPosY(j)%yznl
                  xr(j)   = boundPosY(j)%yznr
                  x(j)    = boundPosY(j)%yzn
                  vel(j)  = 0.0
                  velt(j) = 0.0
                  veltt(j)= 0.0
               enddo
            else
               do j= fzn-1, fzn-5, -1 
                  dx(j)   = dx(fzn)
                  xl(j)   = xl(j+1) - dx(j)
                  xr(j)   = xr(j+1) - dx(j)
                  x(j)    = x(j+1) - dx(j)
                  vel(j)  = 0.0
                  velt(j) = 0.0
                  veltt(j)= 0.0
               enddo               
            endif
            if (mpi_coords(2) /= mpi_dims(2)-1) then 
               do j= nzn+1, nzn+5
                  dx(j)   = boundPosY(j-nzn)%dx
                  xl(j)   = boundPosY(j-nzn)%yznl
                  xr(j)   = boundPosY(j-nzn)%yznr
                  x(j)    = boundPosY(j-nzn)%yzn
                  vel(j)  = 0.0
                  velt(j) = 0.0
                  veltt(j)= 0.0
               enddo
            else
               do j= nzn+1, nzn+5 
                  dx(j)   = dx(nzn)
                  xl(j)   = xl(j-1) + dx(j)
                  xr(j)   = xr(j-1) + dx(j)
                  x(j)    = x(j-1) + dx(j)
                  vel(j)  = 0.0
                  velt(j) = 0.0
                  veltt(j)= 0.0
               enddo
            endif


            IF (basenm(1:1) == 'J' .and. basenm(2:2) == '3') THEN

               if (basenm(3:3) == 'C') then
                  do j=fzn,nzn
                     rho(j) = physic(i,j,k)%densty
                     tem_s(1) = scr(i+ini,j+jni,k+kni)*kb/m_0/c**2
                     rhoe(j) =  rho(j)/1836.110
                     ze = me / tem_s(1)
                     zp = mp / tem_s(1)
                     meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                     dpdrho(1) = tem_s(1)*meff1(1)
                     p(j) = rho(j)*dpdrho(1)
                  enddo
                        
                  if (mpi_coords(2) /= 0) then  
                     do j = fzn-5,fzn-1
                        rho(j) = physic(i,j,k)%densty
                        tem_s(1) = scr(i+ini,j+jni,k+kni)*kb/m_0/c**2
                        rhoe(j) =  rho(j)/1836.110
                        ze = me / tem_s(1)
                        zp = mp / tem_s(1)
                        meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                        dpdrho(1) = tem_s(1)*meff1(1)
                        p(j) = rho(j)*dpdrho(1)
                     enddo
                  else 
                     do j = fzn-5,fzn-1
                        rho(j) = rho(fzn)
                        rhoe(j) =  rho(j)/1836.110
                        tem_s(1) = scr(i+ini,fzn+jni,k+kni)
                        tem_s(1) =  tem_s(1) *kb/m_0/c**2
                        ze = me / tem_s(1)
                        zp = mp / tem_s(1)
                        meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                        dpdrho(1) = tem_s(1)*meff1(1)
                        p(j) = rho(j)*dpdrho(1)
                     enddo
                  endif

                  if (mpi_coords(2) == mpi_dims(2)-1) then
                     do j = nzn+1,nzn+5 
                        rho(j) = rho(nzn)
                        rhoe(j) =  rho(j)/1836.110
                        tem_s(1) = scr(i+ini,nzn+jni,k+kni) 
                        tem_s(1) =  tem_s(1) *kb/m_0/c**2
                        ze = me / tem_s(1)
                        zp = mp / tem_s(1)
                        meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                        dpdrho(1) = tem_s(1)*meff1(1)
                        p(j) = rho(j)*dpdrho(1)
                     enddo
                  else
                     do j = nzn+1,nzn+5
                           rho(j) = physic(i,j,k)%densty
                           tem_s(1) = scr(i+ini,j+jni,k+kni)*kb/m_0/c**2
                           rhoe(j) =  rho(j)/1836.110
                           ze = me / tem_s(1)
                           zp = mp / tem_s(1)
                           meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                           dpdrho(1) = tem_s(1)*meff1(1)
                           p(j) = rho(j)*dpdrho(1)
                        enddo
                  endif

               else

                  do j = fzn-5,nzn+5
         
                     radi  = sqrt(xzn(i)**2+(x(j)+shift)**2+zzn(k)**2)
                     rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                        + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)
                     rhoe(j) =  rho(j)/1836.110

                     if (radi < rm) then 
                        tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                        if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                     else
                        tem_s(1)=tem1*kb/m_0/c**2
                     endif
                     ze = me / tem_s(1)
                     zp = mp / tem_s(1)
                     meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                     dpdrho(1) = tem_s(1)*meff1(1)
                     p(j) = rho(j)*dpdrho(1)
                  enddo   
               endif

            ELSE
               if (ipr == 1) then
                  do j= fzn, nzn 
                     ! p(j)    = pb/dk/(1.0+(x(j)/psh)**nn)**(mm/nn)
                     p(j)     = pb/dk1/(1.+x(j))**mm1
                     rho(j)  = p(j)/(pb/dk1)*1.0
                     rhoe(j) = rho(j)/1836.11
                  enddo

                  do j= fzn-5, fzn-1
                     p(j)    = pb/dk1/(1.+x(j))**mm1
                     rho(j)  = p(j)/(pb/dk1)*1.0
                     rhoe(j) = rho(j)/1836.11
                  enddo

                  do j= nzn+1, nzn+5
                     ! p(j)    = pb/dk/(1.0+(x(j)/psh)**nn)**(mm/nn)
                     p(j)     = pb/dk1/(1.+x(j))**mm1
                     rho(j)  = p(j)/(pb/dk1)*1.0
                     rhoe(j) = rho(j)/1836.11
                  enddo
               else
                  do j= fzn, nzn 
                     p(j)    = pa(j)
                     rho(j)  = rhoa(j)
                     rhoe(j) = rho(j)/1836.11
                  enddo           
                  do j= fzn-5, fzn-1
                     p(j)    = p(fzn)
                     rho(j)  = rhoa(j)
                     rhoe(j) = rho(j)/1836.11
                  enddo
                  do j= nzn+1, nzn+5
                     p(j)    = pa(j)
                     rho(j)  = rhoa(j)
                     rhoe(j) = rho(j)/1836.11
                  enddo
               endif  
            ENDIF
            
            call coef( tid, fzn, nzn, dx, coeff1, coeff2, coeff3, coeff4, coeff5 )

            call interp( tid, fzn, nzn, small,                    &
               coeff1, coeff2, coeff3, coeff4, coeff5,  &
               p, pm, pp, delp )

            call interp( tid, fzn, nzn, small,                    &
               coeff1, coeff2, coeff3, coeff4, coeff5,  &
               rho, rhom, rhop, delrho )

            call interp( tid, fzn, nzn, small,                    &
               coeff1, coeff2, coeff3, coeff4, coeff5,  &
               rhoe, rhoem, rhoep, delrhoe )

            call flaten( tid, fzn, nzn, igodu, epsiln, omg1, omg2,&
               p, vel, flatn, flatn1 ) 

            do j= fzn-1, nzn+1
               tt1 = flatn(j)*rho(j)
               tt7 = flatn(j)*rhoe(j)
               tt5 = flatn(j)*p(j)
               
               rhom(j)  = tt1 + flatn1(j)*rhom(j)
               rhop(j)  = tt1 + flatn1(j)*rhop(j)
               rhoem(j) = tt7 + flatn1(j)*rhoem(j)
               rhoep(j) = tt7 + flatn1(j)*rhoep(j)  
               pm(j)    = tt5 + flatn1(j)*pm(j)
               pp(j)    = tt5 + flatn1(j)*pp(j)
            enddo

            call monot( tid, fzn, nzn, small, p,    pm,     pp,     dp,     p6     )
            call monot( tid, fzn, nzn, small, rho,  rhom,   rhop,   drho,   rho6   )
            call monot( tid, fzn, nzn, small, rhoe, rhoem,  rhoep,  drhoe,  rhoe6  )

            do j= fzn, nzn+1
               rhol(j)  = rhop(j-1)
               rhoel(j) = rhoep(j-1)
               pl(j)    = pp(j-1)
               vell(j)  = 0.
               veltl(j) = 0.
               velttl(j)= 0.
               wl(j)    = 1.0
               xel(j)   = rhoel(j)/rhol(j)
               teml(j)  = pl(j) / (rhol(j)*( xel(j)/me + (1.0-xel(j))/mp )) 


               rhor(j)  = rhom(j)
               rhoer(j) = rhoem(j)
               pr(j)    = pm(j)
               velr(j)  = 0.
               veltr(j) = 0.
               velttr(j)= 0.
               wr(j)    = 1.
               xer(j)   = rhoer(j)/rhor(j)
               temr(j)  = pr(j) / (rhor(j)*( xer(j)/me + (1.0-xer(j))/mp )) 
            enddo


            call eos_synge( 1, tid, fzn, nzn+1, me, mp, rhol, xel, teml, ul,   &
                              pl, gaml, hl, csl, dpdrhl, dpdul)	
            call eos_synge( 1, tid, fzn, nzn+1, me, mp, rhor, xer, temr, ur,   &
                              pr, gamr, hr, csr, dpdrhr, dpdur)


            do j= fzn, nzn+1  
               rl(j)    = rhol(j)*wl(j)
               tt2      = rl(j)*hl(j)*wl(j)
               ml(j)    = tt2*vell(j)
               mtl(j)   = tt2*veltl(j)
               el(j)    = tt2 - rl(j) - pl(j)
               fl(j)    = fp(j-1)

               rr(j)    = rhor(j)*wr(j)
               tt4      = rr(j)*hr(j)*wr(j)
               mr(j)    = tt4*velr(j)
               mtr(j)   = tt4*veltr(j)
               er(j)    = tt4 - rr(j) - pr(j)
               fr(j)    = fm(j)
            enddo

            ! Approximate Riemann solver
            call riem3dp( tid, fzn, nzn+1, rhol, rhor, ul, ur,              &
                           vell, velr, veltl, veltr, velttl, velttr, wl, wr, &
                           csl, csr, pl, pr, hl, hr, dpdrhl, dpdrhr, dpdul,  &
                           dpdur, rl, rr, ml, mr, mtl, mtr,                  &
                           mttl, mttr, el, er,                               &
                           rflx, mflx, mtflx, mttflx, eflx )


            if (basenm(1:1) == 'J' .and. basenm(2:2) == '3') THEN
               do j= fzn, nzn
                  physic(i,j,k)%gravy = (mflx(j+1) - mflx(j))/dx(j)
               enddo
            else
                do j= fzn, nzn
                  grav(j) = (mflx(j+1) - mflx(j))/dx(j)
               enddo
            endif
         enddo
      ENDDO
!$OMP END DO
!$OMP BARRIER

      IF (basenm(1:1) == 'J' .and. basenm(2:2) == '3') THEN
! THIRD SWEEP. GRAVZ.
         fzn = 1
         nzn = nz
!$OMP DO SCHEDULE(static)
         DO k = iniy, finy
            do i = 1, nx

               do j= fzn, nzn

                  x(j)    = zzn(j)
                  xl(j)   = zznl(j)
                  xr(j)   = zznr(j)
                  dx(j)   = xr(j) - xl(j)

                 if (basenm(3:3) == 'C') then
                     rho(j) = physic(i,k,j)%densty
                     tem_s(1) = scr(i+ini,k+jni,j+kni)
                     tem_s(1) =  tem_s(1) *kb/m_0/c**2
                     rhoe(j) =  rho(j)/1836.110
                     
                  else

                     radi  = sqrt(xzn(i)**2+(yzn(k)+shift)**2+zzn(j)**2)
                     rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                          + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)
                     rhoe(j) =  rho(j)/1836.110

                     if (radi < rm) then 
                        tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                        if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                     else
                        tem_s(1)=tem1*kb/m_0/c**2
                     endif
                  endif

                  ze = me / tem_s(1)
                  zp = mp / tem_s(1)
                  meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                  dpdrho(1) = tem_s(1)*meff1(1)
                  p(j) = rho(j)*dpdrho(1)
                 
                  vel(j)  = 0.0
                  velt(j) = 0.0
                  veltt(j)= 0.0
               enddo

               if (mpi_coords(3) /= 0) then
                  do j= fzn-1, fzn-5, -1
                  ! GHANGES BY 3D PARALLELIZATION: copy Z boundaries  

                     dx(j)   = boundPosZ(j)%dx
                     xl(j)   = boundPosZ(j)%zznl
                     xr(j)   = boundPosZ(j)%zznr
                     x(j)    = boundPosZ(j)%zzn

                     if (basenm(3:3) == 'C') then
                        rho(j) = physic(i,k,j)%densty
                        tem_s(1) = scr(i+ini,k+jni,j+kni)*kb/m_0/c**2
                        rhoe(j) =  rho(j)/1836.110
                     else
                        radi  = sqrt(xzn(i)**2+(yzn(k)+shift)**2+x(j)**2)
                        rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                              + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)
                        rhoe(j) =  rho(j)/1836.110

                        if (radi < rm) then 
                           tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                           if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                        else
                           tem_s(1)=tem1*kb/m_0/c**2
                        endif
                     endif
                     ze = me / tem_s(1)
                     zp = mp / tem_s(1)
                     meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                     dpdrho(1) = tem_s(1)*meff1(1)
                     p(j) = rho(j)*dpdrho(1)
                     vel(j)  = 0.0
                     velt(j) = 0.0
                     veltt(j)= 0.0
                  enddo

               else
                  do j= fzn-1, fzn-5, -1
                     dx(j)   = dx(fzn)
                     xl(j)   = xl(j+1)-dx(j)
                     xr(j)   = xr(j+1)-dx(j+1)
                     x(j)    = 0.5*(xl(j)+xr(j))
 
                     if (basenm(3:3) == 'C') then         
                        rho(j)  = rho(fzn)
                        rhoe(j) =  rho(j)/1836.110         
                        tem_s(1) = scr(i+ini,k+jni,fzn+kni)
                        tem_s(1) =  tem_s(1) *kb/m_0/c**2
                     else
                        radi  = sqrt(xzn(i)**2+(yzn(k)+shift)**2+x(j)**2)
                        rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                              + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)
                        rhoe(j) =  rho(j)/1836.110

                        if (radi < rm) then 
                           tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                           if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                        else
                           tem_s(1)=tem1*kb/m_0/c**2
                        endif
                     endif

                     ze = me / tem_s(1)
                     zp = mp / tem_s(1)
                     meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                     dpdrho(1) = tem_s(1)*meff1(1)
                     p(j) = rho(j)*dpdrho(1)
                  
                     vel(j)  = 0.0
                     velt(j) = 0.0
                     veltt(j)= 0.0
                  enddo
               endif

               if (mpi_coords(3) /= mpi_dims(3)-1) then
                  do j= nzn+1, nzn+5
                     dx(j)   = boundPosZ(j-nzn)%dx
                     xl(j)   = boundPosZ(j-nzn)%zznl
                     xr(j)   = boundPosZ(j-nzn)%zznr
                     x(j)    = boundPosZ(j-nzn)%zzn

                     if (basenm(3:3) == 'C') then
                        rho(j) = physic(i,k,j)%densty
                        tem_s(1) = scr(i+ini,k+jni,j+kni)*kb/m_0/c**2
                        rhoe(j) =  rho(j)/1836.110
                     else
                        radi  = sqrt(xzn(i)**2+(yzn(k)+shift)**2+x(j)**2)
                        rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                              + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)
                        rhoe(j) =  rho(j)/1836.110

                        if (radi < rm) then 
                           tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                           if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                        else
                           tem_s(1)=tem1*kb/m_0/c**2
                        endif
                     endif

                     ze = me / tem_s(1)
                     zp = mp / tem_s(1)
                     meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                     dpdrho(1) = tem_s(1)*meff1(1)
                     p(j) = rho(j)*dpdrho(1)
                  
                     vel(j)  = 0.0
                     velt(j) = 0.0
                     veltt(j)= 0.0
                  enddo

               else     
                  do j= nzn+1, nzn+5
                     ! GHANGES BY 3D PARALLELIZATION: copy  boundaries  
                     dx(j)   = dx(nzn)
                     xl(j)   = xl(j-1) + dx(j-1)
                     xr(j)   = xr(j-1) + dx(j)
                     x(j)    = 0.5*(xl(j)+xr(j))

                     if (basenm(3:3) == 'C') then                                                 
                        rho(j)  = rho(nzn)
                        rhoe(j) =  rho(j)/1836.110              
                        tem_s(1) = scr(i+ini,k+jni,nzn+kni) 
                        tem_s(1) = tem_s(1) *kb/m_0/c**2
                            
                     else
                        radi  = sqrt(xzn(i)**2+(yzn(k)+shift)**2+x(j)**2)
                        rho(j) = 1./(1.+(radi/psh1)**nn)**(mm1/nn)      & 
                             + fg*1./(1.+(radi/psh2)**nn)**(mm2/nn)
                        rhoe(j) =  rho(j)/1836.110

                        if (radi < rm) then 
                           tem_s(1)=(tem0+(tem1-tem0)*radi/rm)*kb/m_0/c**2
                           if (tem_s(1) > tem1*kb/m_0/c**2) tem_s(1)=tem1*kb/m_0/c**2
                        else
                           tem_s(1)=tem1*kb/m_0/c**2
                        endif
                     endif

                     ze = me / tem_s(1)
                     zp = mp / tem_s(1)
                     meff1(1) = (1.0-1./1836.110)/mp + 1./1836.110/me
                     dpdrho(1) = tem_s(1)*meff1(1)
                     p(j) = rho(j)*dpdrho(1)
                     
                     vel(j)  = 0.0
                     velt(j) = 0.0
                     veltt(j)= 0.0
                  enddo
               endif
                              
               call coef( tid, fzn, nzn, dx, coeff1, coeff2, coeff3, coeff4, coeff5 )

               call interp( tid, fzn, nzn, small,                    &
               coeff1, coeff2, coeff3, coeff4, coeff5,  &
               p, pm, pp, delp )

               call interp( tid, fzn, nzn, small,                    &
               coeff1, coeff2, coeff3, coeff4, coeff5,  &
               rho, rhom, rhop, delrho )

               call interp( tid, fzn, nzn, small,                    &
               coeff1, coeff2, coeff3, coeff4, coeff5,  &
               rhoe, rhoem, rhoep, delrhoe )

               call flaten( tid, fzn, nzn, igodu, epsiln, omg1, omg2,&
               p, vel, flatn, flatn1 ) 

               do j= fzn-1, nzn+1
                  tt1 = flatn(j)*rho(j)
                  tt7 = flatn(j)*rhoe(j)
                  tt5 = flatn(j)*p(j)
               
                  rhom(j)  = tt1 + flatn1(j)*rhom(j)
                  rhop(j)  = tt1 + flatn1(j)*rhop(j)
                  rhoem(j) = tt7 + flatn1(j)*rhoem(j)
                  rhoep(j) = tt7 + flatn1(j)*rhoep(j)  
                  pm(j)    = tt5 + flatn1(j)*pm(j)
                  pp(j)    = tt5 + flatn1(j)*pp(j)
               enddo

               call monot( tid, fzn, nzn, small, p,    pm,     pp,     dp,     p6     )
               call monot( tid, fzn, nzn, small, rho,  rhom,   rhop,   drho,   rho6   )
               call monot( tid, fzn, nzn, small, rhoe, rhoem,  rhoep,  drhoe,  rhoe6  )

               do j= fzn, nzn+1
                  rhol(j)  = rhop(j-1)
                  rhoel(j) = rhoep(j-1)
                  pl(j)    = pp(j-1)
                  vell(j)  = 0.
                  veltl(j) = 0.
                  velttl(j)= 0.
                  wl(j)    = 1.0
                  xel(j)   = rhoel(j)/rhol(j)
                  teml(j)  = pl(j) / (rhol(j)*( xel(j)/me + (1.0-xel(j))/mp )) 


                  rhor(j)  = rhom(j)
                  rhoer(j) = rhoem(j)
                  pr(j)    = pm(j)
                  velr(j)  = 0.
                  veltr(j) = 0.
                  velttr(j)= 0.
                  wr(j)    = 1.
                  xer(j)   = rhoer(j)/rhor(j)
                  temr(j)  = pr(j) / (rhor(j)*( xer(j)/me + (1.0-xer(j))/mp )) 
               enddo


               call eos_synge( 1, tid, fzn, nzn+1, me, mp, rhol, xel, teml, ul,   &
                                 pl, gaml, hl, csl, dpdrhl, dpdul)	
               call eos_synge( 1, tid, fzn, nzn+1, me, mp, rhor, xer, temr, ur,   &
                                 pr, gamr, hr, csr, dpdrhr, dpdur)


               do j= fzn, nzn+1  
                  rl(j)    = rhol(j)*wl(j)
                  tt2      = rl(j)*hl(j)*wl(j)
                  ml(j)    = tt2*vell(j)
                  mtl(j)   = tt2*veltl(j)
                  el(j)    = tt2 - rl(j) - pl(j)
                  fl(j)    = fp(j-1)

                  rr(j)    = rhor(j)*wr(j)
                  tt4      = rr(j)*hr(j)*wr(j)
                  mr(j)    = tt4*velr(j)
                  mtr(j)   = tt4*veltr(j)
                  er(j)    = tt4 - rr(j) - pr(j)
                  fr(j)    = fm(j)
               enddo

               ! Approximate Riemann solver
               call riem3dp( tid, fzn, nzn+1, rhol, rhor, ul, ur,              &
                              vell, velr, veltl, veltr, velttl, velttr, wl, wr, &
                              csl, csr, pl, pr, hl, hr, dpdrhl, dpdrhr, dpdul,  &
                              dpdur, rl, rr, ml, mr, mtl, mtr,                  &
                              mttl, mttr, el, er,                               &
                              rflx, mflx, mtflx, mttflx, eflx )

               do j= fzn, nzn
                  physic(i,k,j)%gravz = (mflx(j+1) - mflx(j))/abs(dx(j))
               enddo

            enddo
         ENDDO
!$OMP END DO
      ENDIF

!$OMP END PARALLEL 

      deallocate (scr)

      if (cartrank == 0) then
         print*,'-------------------------- GRAV -------------------------------'
         print*,mpi_coords(1), mpi_coords(2),mpi_coords(3),physic(nx/2,1,nz/2)%gravx,physic(nx/2,1,nz/2)%gravy,physic(nx/2,1,nz/2)%gravz
         print*,mpi_coords(1), mpi_coords(2),mpi_coords(3),physic(nx/2,ny,nz/2)%gravx,physic(nx/2,ny,nz/2)%gravy,physic(nx/2,ny,nz/2)%gravz
         print*,'---------------------------------------------------------------'
      endif
      
!----------------------------------------------------------------------------END
      END subroutine psgrav
