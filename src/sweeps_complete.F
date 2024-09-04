!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!===============================================================================
      subroutine sweeps(basenm, rkin, igodu, igeomx, igeomy, igeomz,            &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,       &
               nx, nz, nsdim, ny, nys, xyswp, mnx5, mny5, mnz5,      &
               alpha_R, beta_R, gamma_ad, dt,                        &
               timej, t0, t1, mm1, mm2, nn, brt,  vptbd, rptbd,      &
               pptbd, ms, delx, dk1, dk2, pb, psh1, psh2, rhob,      &
               rhoeb, velb, small, smlrho,                           &
               epsiln, omg1, omg2, uwi, vwi, shift, rorb,            &
               R_b, meff, rho_0a, me, mp, KB1,                       &
               KB2, cl, pc, m_0, kb, mloss0, tem0, tem1, rm, fg,     &
               nmod, ipr, A0, gridly, fra, frs, pha, pa, rhoa,       &
               xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,    &
               boundPosX, boundPosY, boundPosZ, grav, physic, curth, &
               ic_radius, ic_density, ic_pressure, nstep, nbegin, &
               fxbounds, nxbounds, fzbounds, nzbounds)

      USE tipos
      USE omp_lib
      USE workarrays
      implicit none
# include "types.h"

!-------------------------------------------------------Input Variables

      character, intent(in):: basenm*(*)

      integer(is):: rkin, igodu, igeomx, igeomy, igeomz, ms,         &
      bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,  &
      nx, nz, nsdim, ny, nys, xyswp, mnx5, mny5, mnz5, &
      nmod, ipr, brt, ierr

      real(rs):: alpha_R, beta_R, gamma_ad, dt,                &
      timej, t0, t1, mm1, mm2, nn, vptbd, rptbd, pptbd,    &
      delx, dk1, dk2, pb, psh1, psh2, rhob, rhoeb, velb,   &
      small, smlrho, epsiln, omg1, omg2, me, mp, KB1, KB2, &
      cl, rho_0a, m_0, kb, R_b, pc, mloss0, tem0, tem1,    &
      rm, fg, uwi, vwi, shift, rorb, meff, gridly  
      
      real(rs)   :: A0, pha(25), frs(25), fra(25)
      
      real(rs) :: rhoa(-4:mny5), pa(-4:mny5)

! sizes mnx1, mny1, mnz1
      real(rs):: xznl(*), xzn(*), xznr(*), &
      yznl(*), yzn(*), yznr(*), &
      zznl(*), zzn(*), zznr(*)

      type(TboundPosX) :: boundPosX(-4:5)
      type(TboundPosY) :: boundPosY(-4:5)
      type(TboundPosZ) :: boundPosZ(-4:5)
      
      real(8), dimension(10000) :: ic_radius, ic_density, ic_pressure

!------------------------------------------------------Output Variables
      
      real(rs)   :: grav(-4:mny5)

      type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)

!-------------------------------------------------------Local Variables

      integer(is):: i, j, k, fznp, nznp, nxp, nyp, nzp, &
      igeom, bndmin, bndmax, fzn, nzn

      integer :: tid

!-------------------------------------------------------OMP Variables

      integer(is) :: nthreads, curth


      real(rs)   :: a, b, c, aux, gam1, dedt, tem, tem_s, xe_s

      real(rs), pointer :: vel(:), velt(:), veltt(:), u(:),      &
      rho(:), rhoe(:), rhopr(:), p(:), h(:), &
      x(:), xl(:), xr(:), dx(:),             &
      area(:), areal(:), arear(:), dvol(:)

      real(rs), pointer :: r(:), re(:), rp(:), m(:), mt(:), mtt(:), &
      e(:), f(:), fict(:)

      real(rs), pointer :: rflx(:), reflx(:), mflx(:), mtflx(:), &
      mttflx(:), eflx(:), fflx(:)

      real(rs), pointer :: trac(:,:,:)

      type(Tbounds) :: fxbounds(-4:0), nxbounds(1:5)
      type(Tbounds) :: fybounds(-4:0), nybounds(1:5)
      type(Tbounds) :: fzbounds(-4:0), nzbounds(1:5)
      integer(is) :: nstep, nbegin

!---------------------------------------------variables de timing

      real(4) :: temps
      integer(4) :: estat

!--------------------------------------------------------------------------BEGIN

      nthreads = omp_th

!      trac => wa409                                                                                                                                           
      allocate (trac(mnx5,mny5,mnz5), stat= ierr)

!$OMP PARALLEL DEFAULT(shared)& 
!$OMP PRIVATE(tid,curth,i,j,k,fznp,nznp,nxp,nyp,nzp,igeom,bndmin,bndmax,fzn,nzn)&
!$OMP PRIVATE(a,b,c,aux,gam1)&
!$OMP PRIVATE(vel,velt,veltt,u,rho,rhoe,rhopr,p,x,xl,xr,dx)&
!$OMP PRIVATE(area,areal,arear,dvol,r,re,rp,m,mt,mtt,e,f,fict)&
!$OMP PRIVATE(rflx,reflx,mflx,mtflx,mttflx,eflx,fflx)&
!$OMP PRIVATE(xe_s,meff,tem,dedt)&
!$OMP PRIVATE(temps, estat)

# ifdef OPENMP
      call my_omp_get_thread_num(curth)
      tid = curth + 1
# else
      tid = 1
# endif
      
! -4:mn6 vectors
      rflx   => wa(tid)%wa100
      reflx  => wa(tid)%wa188
      mflx   => wa(tid)%wa101
      mtflx  => wa(tid)%wa102
      mttflx => wa(tid)%wa103
      eflx   => wa(tid)%wa104
      fflx   => wa(tid)%wa105

! -4:mn5 vectors
      vel   => wa(tid)%wa200
      velt  => wa(tid)%wa201
      veltt => wa(tid)%wa202
      u     => wa(tid)%wa203
      rho   => wa(tid)%wa204
      rhoe  => wa(tid)%wa205
      rhopr => wa(tid)%wa206
      p     => wa(tid)%wa207
      x     => wa(tid)%wa208
      xl    => wa(tid)%wa209
      xr    => wa(tid)%wa210
      dx    => wa(tid)%wa211
      area  => wa(tid)%wa212
      areal => wa(tid)%wa213
      arear => wa(tid)%wa214
      dvol  => wa(tid)%wa215

      r    => wa(tid)%wa216
      re   => wa(tid)%wa217
      m    => wa(tid)%wa218
      mt   => wa(tid)%wa219
      mtt  => wa(tid)%wa220
      e    => wa(tid)%wa221
      f    => wa(tid)%wa222
      fict => wa(tid)%wa223

!      trac => wa409


! print*,'SWEEPS X-sweep '
!print*,cartrank,'x-sweep'
!X-sweep
      if (rkin==1) then
         a = 1.0
         b = 0.0
         c = 1.0
      else if (rkin==2) then
         a = beta_R/alpha_R
         b = 1.0/alpha_R
         c = b
      else
         a = 1.0/3.0
         b = 2.0/3.0
         c = b
      endif

!      gam1   = gamma_ad
      xyswp  = 1
      bndmin = bndmnx
      bndmax = bndmxx
      igeom  = igeomx


      fzn    = 1
      nzn    = nx
      fznp   = fzn
      nznp   = nzn
      nyp    = ny
      nzp    = nz

! Check working OMP threads      
!     if (cartrank==0 .and. rkin==1) then
!        print*,'SWEEPS WITH',omp_get_num_threads() 
!     endif
 
!$OMP DO SCHEDULE(static)
       do k= 1, nzp
         do j= 1, nyp
            do i= fznp, nznp

               !print*,'SWEEP LOOP (i,j,k) =', i, j, k
               !print*,'THREAD NUMBER',omp_get_thread_num()
               
               vel(i)  = physic(i,j,k)%velx
               velt(i) = physic(i,j,k)%vely
               veltt(i)= physic(i,j,k)%velz

               rho(i)  = physic(i,j,k)%densty
               rhoe(i) = physic(i,j,k)%denstye
               rhopr(i)= physic(i,j,k)%densty-physic(i,j,k)%denstye
               p(i)    = physic(i,j,k)%pres

               f(i)    = physic(i,j,k)%tracer
               
               trac(i,j,k) = physic(i,j,k)%tracer

               xl(i)   = xznl(i)
               x(i)    = xzn(i)
               xr(i)   = xznr(i)
               dx(i)   = xr(i) - xl(i)
               fict(i) = 0.0
               
            enddo
            
            if ( mpi_coords(1) == 0 .or.   &
            mpi_coords(1) == (mpi_dims(1)-1) ) then

               ! if (j==1 .and. k==1) then
               !    print*, 'Sweeps beggining with nstep = ', nstep, 'and nbegin = ', nbegin
               ! endif

               call bndry(basenm, fznp, nznp, j, k, nsdim, bndmin, bndmax, xyswp, ipr,    &
               gamma_ad, timej, mny5, dt, t0, t1, vptbd, rptbd, pptbd, &
               ms, pb, rhob, rhoeb, velb, mm1, mm2, nn, dk1,  dk2,     &
               psh1, psh2, kb, m_0, cl, me, mp, delx, uwi, vwi,        &
               shift, fg, rm, tem0, tem1, rorb, nmod, A0, gridly,      &
               fra, frs, pha, pa, rhoa, xzn, yzn, zzn, p, rho, rhoe,   &
               rhopr, vel, velt, veltt, f, x, xl, xr, dx,              &
               ic_radius, ic_density, ic_pressure, fxbounds, nxbounds, nstep, nbegin)

            endif


!     -GHANGES BY 3D PARALLELIZATION: added this part 
         
         if (mpi_coords(1) /= 0) then
! Left
! print*,'SWEEPS X-sweep - left'
            
            do i= fznp-5, fznp-1
               vel(i)  = physic(i,j,k)%velx
               velt(i) = physic(i,j,k)%vely
               veltt(i)= physic(i,j,k)%velz

               rho(i)  = physic(i,j,k)%densty
               rhoe(i) = physic(i,j,k)%denstye
               rhopr(i) = rho(i) - rhoe(i)
               p(i)    = physic(i,j,k)%pres
 
               f(i)    = physic(i,j,k)%tracer
               
               xl(i)   = boundPosX(i)%xznl 
               x(i)    = boundPosX(i)%xzn  
               xr(i)   = boundPosX(i)%xznr 
               dx(i)   = boundPosX(i)%dx   
               fict(i) = 0.0
              
            enddo
         endif   

         if (mpi_coords(1) /= mpi_dims(1)-1) then        
            
! Right
! print*,'SWEEPS X-sweep -  right'
            do i= nznp+1, nznp+5
               vel(i)  = physic(i,j,k)%velx
               velt(i) = physic(i,j,k)%vely
               veltt(i)= physic(i,j,k)%velz

               rho(i)  = physic(i,j,k)%densty
               rhoe(i) = physic(i,j,k)%denstye
               rhopr(i) = rho(i) - rhoe(i)           
               p(i)    = physic(i,j,k)%pres
               f(i)    = physic(i,j,k)%tracer

               xl(i)   = boundPosX(i-nznp)%xznl 
               x(i)    = boundPosX(i-nznp)%xzn  
               xr(i)   = boundPosX(i-nznp)%xznr 
               dx(i)   = boundPosX(i-nznp)%dx   
               fict(i) = 0.0
              
            enddo
         endif         
!     - END GHANGES BY 3D PARALLELIZATION

         call geom( fznp, nznp, igeom, xyswp, xl, xr, dx, xzn, &
         area, areal, arear, dvol )


         call hydrow3( tid, fznp, nznp, igodu,                           &
         small, smlrho, epsiln, omg1, omg2, me, mp,   &
         x, dx, areal, p, rho, rhoe, rhopr, vel,      &
         velt, veltt, f, rflx, reflx, mflx, mtflx, mttflx, eflx, fflx )
        

         do i= fznp, nznp
           
            aux = dt / dvol(i)
             
            physic(i,j,k)%tracer = a * physic(i,j,k)%odenstw      &
            * physic(i,j,k)%otracer      &
            + b * physic(i,j,k)%denstw       &
            * physic(i,j,k)%tracer       &
            - c * (fflx(i+1)-fflx(i))*aux

            physic(i,j,k)%denstw = a*physic(i,j,k)%odenstw  &
            + b*physic(i,j,k)%denstw   &
            - c*(rflx(i+1)-rflx(i))*aux

            physic(i,j,k)%denstew = a*physic(i,j,k)%odenstew      &
            + b*physic(i,j,k)%denstew   &
            - c*(reflx(i+1)-reflx(i))*aux

            physic(i,j,k)%momenx = a*physic(i,j,k)%omomenx  &
            + b*physic(i,j,k)%momenx   &
            - c*(mflx(i+1)-mflx(i))*aux

            physic(i,j,k)%momeny = a*physic(i,j,k)%omomeny    &
            + b*physic(i,j,k)%momeny     &
            - c*(mtflx(i+1)-mtflx(i))*aux


            physic(i,j,k)%momenz = a*physic(i,j,k)%omomenz      &
            + b*physic(i,j,k)%momenz       &
            - c*(mttflx(i+1)-mttflx(i))*aux

            physic(i,j,k)%energy = a*physic(i,j,k)%oenergy  &
            + b*physic(i,j,k)%energy   &
            - c*(eflx(i+1)-eflx(i))*aux
          
         enddo
        
         do i= fznp, nznp
            physic(i,j,k)%momenx = physic(i,j,k)%momenx + c*dt*physic(i,j,k)%gravx
            physic(i,j,k)%energy = physic(i,j,k)%energy +   &
                               c*dt*physic(i,j,k)%gravx*physic(i,j,k)%velx
         enddo


      enddo
      enddo
!$OMP END DO


!$OMP BARRIER
 !print*,'SWEEPS Y-sweep '
! print*,cartrank,'y-sweep'
!Y-sweep
      if (rkin==1) then
         b = 1.0
         c = 1.0
      else if (rkin==2) then
         b = 1.0
         c = 1.0/alpha_R
      else
         b = 1.0
         c = 2.0/3.0
      endif

      IF (nsdim > 1) THEN
        
         xyswp  = 2
         bndmin = bndmny
         bndmax = bndmxy
         igeom  = igeomy
         fzn    = 1
         nzn    = ny
         fznp   = fzn
         nznp   = nzn
         nxp    = nx
         nzp    = nz
         
!$OMP DO SCHEDULE(static)
         do k= 1, nzp
            do i= 1, nxp
! Middle
               do j= fznp, nznp
                  vel(j)  = physic(i,j,k)%vely
                  velt(j) = physic(i,j,k)%velz
                  veltt(j)= physic(i,j,k)%velx

                  rho(j)   = physic(i,j,k)%densty
                  rhoe(j)  = physic(i,j,k)%denstye
                  rhopr(j) = physic(i,j,k)%densty - physic(i,j,k)%denstye

                  p(j)    = physic(i,j,k)%pres
                  f(j)    = trac(i,j,k)

                  xl(j)   = yznl(j)
                  x(j)    = yzn(j)
                  xr(j)   = yznr(j)
                  dx(j)   = xr(j)-xl(j)
                  fict(j) = 0.0
               
               enddo
!     -GHANGES BY 3D PARALLELIZATION: nuproc --> mpi_coords(2), nbproc --> mpi_dims(2)

               if (bndmny /= 4 .and. ( mpi_coords(2) == 0 .or.   &
               mpi_coords(2) == (mpi_dims(2)-1) )) then

                  call bndry(basenm, fznp, nznp, k, i, nsdim, bndmin, bndmax, xyswp, ipr,    &
                  gamma_ad, timej, mny5, dt, t0, t1, vptbd, rptbd, pptbd, &
                  ms, pb, rhob, rhoeb, velb, mm1, mm2, nn, dk1,  dk2,     &
                  psh1, psh2, kb, m_0, cl, me, mp, delx, uwi, vwi,        &
                  shift, fg, rm, tem0, tem1, rorb, nmod, A0, gridly,      &
                  fra, frs, pha, pa, rhoa, xzn, yzn, zzn, p, rho, rhoe,   &
                  rhopr, vel, velt, veltt, f, x, xl, xr, dx,              &
                  ic_radius, ic_density, ic_pressure, fybounds, nybounds, nstep, nbegin)
   
               endif
            


!     -GHANGES BY 3D PARALLELIZATION: nuproc --> mpi_coords(2), nbproc --> mpi_dims(2)

            if (mpi_coords(2) /= 0) then
! Left
               do j= fznp-5, fznp-1
                  vel(j)  = physic(i,j,k)%vely
                  velt(j) = physic(i,j,k)%velz
                  veltt(j)= physic(i,j,k)%velx

                  rho(j)  = physic(i,j,k)%densty
                  rhoe(j) = physic(i,j,k)%denstye
                  rhopr(j) = rho(j)-rhoe(j)
                  p(j)    = physic(i,j,k)%pres
              
                  f(j)    = physic(i,j,k)%tracer
                  
                  xl(j)   = boundPosY(j)%yznl !yznl(j)
                  x(j)    = boundPosY(j)%yzn !yzn(j)
                  xr(j)   = boundPosY(j)%yznr !yznr(j)
                  dx(j)   = boundPosY(j)%dx !xr(j)-xl(j)
                  fict(j) = 0.0
               enddo
            endif   

            if (mpi_coords(2) /= mpi_dims(2)-1) then        
               
! Right
               do j= nznp+1, nznp+5
                  vel(j)  = physic(i,j,k)%vely
                  velt(j) = physic(i,j,k)%velz
                  veltt(j)= physic(i,j,k)%velx

                  rho(j)  = physic(i,j,k)%densty
                  rhoe(j) = physic(i,j,k)%denstye
                  rhopr(j) = rho(j)-rhoe(j)
                  p(j)    = physic(i,j,k)%pres

                  f(j)    = physic(i,j,k)%tracer

                  xl(j)   = boundPosY(j-nznp)%yznl !yznl(j)
                  x(j)    = boundPosY(j-nznp)%yzn !yzn(j)
                  xr(j)   = boundPosY(j-nznp)%yznr !yznr(j)
                  dx(j)   = boundPosY(j-nznp)%dx !xr(j)-xl(j)
                  fict(j) = 0.0
               enddo
            endif         
            
            call geom( fznp, nznp, igeom, xyswp, xl, xr, dx, xzn, &
            area, areal, arear, dvol )


!            if (i == 1 .and. rkin == 2) then
!               print*,cartrank,tid,'hydrow',p(1),i,k
!               print*,i,k, vel(0), velt(0), veltt(0),rho(0),rhoe(0)
!               print*,i,k, vel(ny+1), velt(ny+1), veltt(ny+1),rho(ny+1),rhoe(ny+1)
!            endif
               

            call hydrow3( tid, fznp, nznp, igodu,                           &
            small, smlrho, epsiln, omg1, omg2, me, mp,   &
            x, dx, areal, p, rho, rhoe, rhopr, vel,      &
            velt, veltt, f, rflx, reflx, mflx, mtflx, mttflx, eflx, fflx )

!            if (i == 1 .and. rkin == 2) then
!               print*,cartrank,tid,'flux',p(1),i,k
!               print*,i,k, mflx(1), mtflx(1), mttflx(1),rflx(1),eflx(1)
!               print*,i,k, mflx(ny+1), mtflx(ny+1), mttflx(ny+1),rflx(ny+1),eflx(ny+1)                            
!            endif


           
            do j= fznp, nznp
               aux = dt / dvol(j)

               physic(i,j,k)%tracer = b*physic(i,j,k)%tracer       &
               - c*(fflx(j+1)-fflx(j))*aux
               physic(i,j,k)%denstw = b*physic(i,j,k)%denstw       &
               - c*(rflx(j+1)-rflx(j))*aux
               physic(i,j,k)%denstew = b*physic(i,j,k)%denstew       &
               - c*(reflx(j+1)-reflx(j))*aux
               physic(i,j,k)%momenx = b*physic(i,j,k)%momenx       &
               - c*(mttflx(j+1)-mttflx(j))*aux
               physic(i,j,k)%momeny = b*physic(i,j,k)%momeny       &
               - c*(mflx(j+1)-mflx(j))*aux
               physic(i,j,k)%momenz = b*physic(i,j,k)%momenz       &
               - c*(mtflx(j+1)-mtflx(j))*aux
               physic(i,j,k)%energy = b*physic(i,j,k)%energy       &
               - c*(eflx(j+1)-eflx(j))*aux

            enddo


!            if (i == 1 .and. rkin == 2) then
!               print*,cartrank,tid,'grav',p(1),i,k
!               print*,i,k, physic(1,1,k)%denstw,physic(1,1,k)%energy,physic(1,1,k)%gravy
!               print*,i,k, physic(1,ny,k)%denstw,physic(1,ny,k)%energy,physic(1,ny,k)%gravy
!            endif


!            print*,cartrank,tid,'grav',physic(i,1,k)%energy,i,k


            if (basenm(1:1) == 'J'.and.basenm(2:2)=='3') then
                do j= fznp, nznp
                  physic(i,j,k)%momeny = physic(i,j,k)%momeny + c*dt*physic(i,j,k)%gravy
                  physic(i,j,k)%energy = physic(i,j,k)%energy +   &
                                         c*dt*physic(i,j,k)%gravy*physic(i,j,k)%vely
               enddo

            else if (basenm(1:1) == 'J'.and.basenm(2:2)/='3') then
               do j= fznp, nznp
                  physic(i,j,k)%momeny = physic(i,j,k)%momeny + c*dt*grav(j)
                  physic(i,j,k)%energy = physic(i,j,k)%energy + c*dt*grav(j)*physic(i,j,k)%vely
               enddo
            endif
         enddo
      enddo
!$OMP END DO
      ENDIF

!$OMP BARRIER
! print*,cartrank,'z-sweep '
!Z-sweep
      IF (nsdim > 2) THEN
         if (rkin==1) then
            b = 1.0
            c = 1.0
         else if (rkin==2) then
            b = 1.0
            c = 1.0/alpha_R
         else
            b = 1.0
            c = 2.0/3.0
         endif

         xyswp  = 3
         bndmin = bndmnz
         bndmax = bndmxz
         igeom  = igeomz
         fzn    = 1
         nzn    = nz
         fznp   = fzn
         nznp   = nzn
         nxp    = nx
         nyp    = ny
         

!$OMP DO SCHEDULE(static)
         do j= 1, nyp
            do i= 1, nxp
               do k= fznp, nznp    
                  vel(k)  = physic(i,j,k)%velz
                  velt(k) = physic(i,j,k)%velx
                  veltt(k)= physic(i,j,k)%vely

                  rho(k)  = physic(i,j,k)%densty
                  rhoe(k) = physic(i,j,k)%denstye
                  rhopr(k) = physic(i,j,k)%densty - physic(i,j,k)%denstye 
                  p(k)    = physic(i,j,k)%pres
                  f(k)    = trac(i,j,k)

                  xl(k)   = zznl(k)
                  x(k)    = zzn(k)
                  xr(k)   = zznr(k)
                  dx(k)   = xr(k)-xl(k)
                  fict(k) = 0.0
               enddo
             

               if ( mpi_coords(3) == 0 .or.   &
               mpi_coords(3) == (mpi_dims(3)-1) ) then
         
                  call bndry(basenm, fznp, nznp, i, j, nsdim, bndmin, bndmax, xyswp, ipr,    &
                  gamma_ad, timej, mny5, dt, t0, t1, vptbd, rptbd, pptbd, &
                  ms, pb, rhob, rhoeb, velb, mm1, mm2, nn, dk1,  dk2,     &
                  psh1, psh2, kb, m_0, cl, me, mp, delx, uwi, vwi,        &
                  shift, fg, rm, tem0, tem1, rorb, nmod, A0, gridly,      &
                  fra, frs, pha, pa, rhoa, xzn, yzn, zzn, p, rho, rhoe,   &
                  rhopr, vel, velt, veltt, f, x, xl, xr, dx,              &
                  ic_radius, ic_density, ic_pressure, fzbounds, nzbounds, nstep, nbegin)
      
               endif

!     -GHANGES BY 3D PARALLELIZATION: added this part 
               if (mpi_coords(3) /= 0) then
! Left
                  do k= fznp-5, fznp-1
                     vel(k)  = physic(i,j,k)%velz
                     velt(k) = physic(i,j,k)%velx
                     veltt(k)= physic(i,j,k)%vely
                     
                     rho(k)  = physic(i,j,k)%densty
                     rhoe(k) = physic(i,j,k)%denstye
                     rhopr(k) = rho(k) -rhoe(k) 
                     p(k)    = physic(i,j,k)%pres
                     
                     f(k)    = physic(i,j,k)%tracer
                     
                     xl(k)   = boundPosZ(k)%zznl 
                     x(k)    = boundPosZ(k)%zzn  
                     xr(k)   = boundPosZ(k)%zznr 
                     dx(k)   = boundPosz(k)%dx   
                     fict(k) = 0.0
                  enddo
               endif   
               
            if (mpi_coords(3) /= mpi_dims(3)-1) then        
!     Right
                  do k= nznp+1, nznp+5
                     vel(k)  = physic(i,j,k)%velz
                     velt(k) = physic(i,j,k)%velx
                     veltt(k)= physic(i,j,k)%vely
                     
                     rho(k)  = physic(i,j,k)%densty
                     rhoe(k) = physic(i,j,k)%denstye
                     rhopr(k) = rho(k) -rhoe(k) 
                     p(k)    = physic(i,j,k)%pres
                     
                     f(k)    = physic(i,j,k)%tracer
                     
                     xl(k)   = boundPosZ(k-nznp)%zznl 
                     x(k)    = boundPosZ(k-nznp)%zzn  
                     xr(k)   = boundPosZ(k-nznp)%zznr 
                     dx(k)   = boundPosZ(k-nznp)%dx   
                     fict(k) = 0.0
                  enddo
               endif         
!     - END GHANGES BY 3D PARALLELIZATION


               call geom( fznp, nznp, igeom, xyswp, xl, xr, dx, xzn, &            
                          area, areal, arear, dvol )

               
               call hydrow3( tid, fznp, nznp, igodu,                           &
               small, smlrho, epsiln, omg1, omg2, me, mp,   &
               x, dx, areal, p, rho, rhoe, rhopr, vel,      &
               velt, veltt, f, rflx, reflx, mflx, mtflx, mttflx, eflx, fflx )
               

               do k= fznp, nznp
                  aux = dt/dvol(k)              

                  physic(i,j,k)%denstw = b*physic(i,j,k)%denstw         &
                  - c*(rflx(k+1)-rflx(k))*aux
                  physic(i,j,k)%denstew = b*physic(i,j,k)%denstew         &
                  - c*(reflx(k+1)-reflx(k))*aux
                  physic(i,j,k)%momenx = b*physic(i,j,k)%momenx         &
                  - c*(mtflx(k+1)-mtflx(k))*aux
                  physic(i,j,k)%momeny = b*physic(i,j,k)%momeny         &
                  - c*(mttflx(k+1)-mttflx(k))*aux
                  physic(i,j,k)%momenz = b*physic(i,j,k)%momenz         &
                  - c*(mflx(k+1)-mflx(k))*aux
                  physic(i,j,k)%energy = b*physic(i,j,k)%energy         &
                  - c*(eflx(k+1)-eflx(k))*aux
                  physic(i,j,k)%tracer = (b*physic(i,j,k)%tracer        &
                  - c*(fflx(k+1)-fflx(k))*aux) / physic(i,j,k)%denstw
                 


                  if ((basenm(1:1) == 'M' .and. basenm(2:2) == 'C').or. &
                  (basenm(1:1) == 'P'.and. basenm(2:2) == 'W').or.      &
                  (basenm(1:1) == 'J'.and. basenm(2:2) == '3')) then ! .and. basenm(4:4) /= '0')) then
                     xe_s  =  physic(i,j,k)%denstye/physic(i,j,k)%densty
                     meff= 1.0/9.1095e-28 * (xe_s + (1.0-xe_s)/1836.11)
                     tem = physic(i,j,k)%pres/physic(i,j,k)%densty*9.e20/meff/kb
                  
                     if (tem .lt. 1e4) then
                        dedt = 0.
                     elseif (tem .ge. 1.e4 .and. tem .lt. 1.e5) then
                        dedt = (physic(i,j,k)%densty*rho_0a*meff)**2 * 7.e-27*tem 
                     elseif (tem .ge. 1.e5 .and. tem .lt. 4.e7) then
                        dedt = (physic(i,j,k)%densty*rho_0a*meff)**2 * 7.e-19*tem**(-0.6)
                     elseif (tem .ge. 4.e7) then
                        dedt = (physic(i,j,k)%densty*rho_0a*meff)**2 * 3.e-27*tem**(0.5)
                     endif
                  
                     dedt = dedt/rho_0a/9.e20/3.e10*R_b
                     physic(i,j,k)%energy = physic(i,j,k)%energy - c*dedt*dt
                     
!     if (nuproc == 0 .and. physic(i,j,k)%energy <= 0.) then
!     print*,nuproc,i,j,k,physic(i,j,k)%energy+c*dedt*dt,physic(i,j,k)%energy,dedt,dt
!     endif
                  
                  endif

                  if (basenm(3:3) == 'M') then                     
                     if (physic(i,j,k)%tracer > 0.0) then
                        physic(i,j,k)%denstw  = physic(i,j,k)%denstw   &
                        + c*dt*physic(i,j,k)%stems                         
                        physic(i,j,k)%denstew = physic(i,j,k)%denstew   &
                        + c*dt*physic(i,j,k)%stems/1836.11
                     endif
                  endif 

              enddo

              if ((basenm(1:1) == 'J'.and. basenm(2:2) == '3')) then
                 do k = fznp, nznp
                     physic(i,j,k)%momenz = physic(i,j,k)%momenz + c*dt*physic(i,j,k)%gravz
                     physic(i,j,k)%energy = physic(i,j,k)%energy +   &
                                c*dt*physic(i,j,k)%gravz*physic(i,j,k)%velz
                 enddo
              endif 

           enddo
         enddo
!$OMP END DO
      
      ENDIF
!$OMP END PARALLEL


 deallocate ( trac )

! print*,'END SWEEPS'
!----------------------------------------------------------------------------END
      END subroutine sweeps
