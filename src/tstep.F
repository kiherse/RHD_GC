!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!     This subroutine calculates the new timestep value.
!===============================================================================

      subroutine tstep( basenm, nx, ny, nz, mnx5, mny5, mnz5, nsdim,   &
      igeomy, itstp, nstep, timej, cfl,                          &
      meff, R_b, rho_0a, kb, me, mp, cl, m_0,                    &
      xznl, xznr, xzn, yznl, yznr, zznl, zznr, physic, dt)
      
      USE tipos
      USE parallel
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

      integer(is):: nx, ny, nz, mnx5, mny5, mnz5, &
      nsdim, igeomy, itstp, nstep, xCPUs, yCPUs, zCPUs 

      real(rs):: timej, cfl, meff, R_b, rho_0a, kb, me, mp, cl, m_0

      real(rs):: xznl(*), xznr(*), xzn(*), &
      yznl(*), yznr(*),         &
      zznl(*), zznr(*)

      type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------------------Output Variables

      real(rs) :: dt

!-------------------------------------------------------Local Variables

      integer(is):: tid, i, j, k, ic, jc, kc, fx, pmin, reduceMinDT, on

      integer(is):: ict(64), jct(64), kct(64)
      
      real(rs):: dtt(64), soundt(64), dtcc_omp(64), dtc_omp(64)

# ifdef OPENMP
      integer :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
# endif

      real(rs):: lambd1, lambd4, lambda, lamb1x, lamb4x, lamb1y, lamb4y, &
      lamb1z, lamb5x, lamb5y, lamb5z,                         &
      dtcc, dtestmax, sound, v, olddt, tt1, tt2, tt3, tt4

      real(rs):: tem, dedt, dtc, dtco, temc, dtcmin, temcmin, xe

      real(rs), pointer :: lambdx(:), lambdy(:), lambdz(:), &
      scrch1(:), scrch2(:), scrch3(:)

      real(rs), pointer :: rho_s(:),xe_s(:),tem_s(:),         & 
      eps_s(:),p_s(:),gam_s(:),          &
      h_s(:),cs_s(:),                    &
      dpdrho_s(:),dpdeps_s(:)

      real(rs), allocatable :: dtcc_thread(:, :, :), dtc_thread(:, :, :)

      allocate( dtcc_thread(nx,ny,nz), dtc_thread(nx,ny,nz) )

!--------------------------------------------------------------------------BEGIN

      dtcc = 0.d0
      dtc  = 1.d0
      ic    = 0
      jc    = 0
      kc    = 0
      olddt = dt

!$OMP PARALLEL DEFAULT(shared)& 
!$OMP PRIVATE(tid,i,j,k,tt1,tt2,tt3,tt4)&
!$OMP PRIVATE(xe_s,tem_s,rho_s,eps_s,p_s)&
!$OMP PRIVATE(gam_s,h_s,cs_s,dpdrho_s,dpdeps_s)&
!$OMP PRIVATE(scrch1,scrch2,scrch3,lamb1x,lamb4x,lamb1y,lamb4y)&
!$OMP PRIVATE(lamb1z,lamb5x,lamb5y,lamb5z)&
!$OMP PRIVATE(lambdx,lambdy,lambdz)&
!$OMP PRIVATE(xe,tem,dedt,dtco,dtcmin,temcmin)&
!$OMP PRIVATE(dtestmax)
      !If OMP reduction works: !$OMP REDUCTION(max:dtcc)&
      !If OMP reduction works:  !$OMP REDUCTION(min:dtc)

# ifdef OPENMP
     tid = omp_get_thread_num()+1
# else
     tid = 1
# endif

  ! -4:mn5 vectors
      rho_s => wa(tid)%wa201
      xe_s  => wa(tid)%wa202
      tem_s => wa(tid)%wa203
      eps_s => wa(tid)%wa204
      p_s   => wa(tid)%wa205
      gam_s => wa(tid)%wa206 
      h_s   => wa(tid)%wa207
      cs_s  => wa(tid)%wa208
      dpdrho_s => wa(tid)%wa209
      dpdeps_s => wa(tid)%wa210

! -4:mn6
      scrch1 => wa(tid)%wa100
      scrch2 => wa(tid)%wa101
      scrch3 => wa(tid)%wa102

! mn vectors
      lambdx => wa(tid)%wa300
      lambdy => wa(tid)%wa301
      lambdz => wa(tid)%wa302

      fx    = 1.d0
      dtcmin= 10.d0
      temcmin = 10.d0

      dtcc_omp(tid) = 0.d0
      dtc_omp(tid) = 1.d0


      IF (nsdim == 3) THEN
!$OMP DO SCHEDULE(static)
         do k=1, nz  
            do j= 1, ny
               do i= 1, nx
                  xe_s(i)   = physic(i,j,k)%denstye/physic(i,j,k)%densty
                  tem_s(i)  = physic(i,j,k)%pres / (physic(i,j,k)%densty*( xe_s(i)/me + (1.0-xe_s(i))/mp ))
                  rho_s(i)  = physic(i,j,k)%densty
                  eps_s(i)  = physic(i,j,k)%eps
               enddo
               
               call eos_synge( 1, tid, 1, nx, me, mp, rho_s, xe_s, tem_s, eps_s,   &
               p_s, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)     
               
               do i=1,nx
                  tt1 = physic(i,j,k)%velx*physic(i,j,k)%velx
                  tt2 = physic(i,j,k)%vely*physic(i,j,k)%vely
                  tt3 = physic(i,j,k)%velz*physic(i,j,k)%velz
                  tt4 = cs_s(i)*cs_s(i)
                  
                  scrch1(i)  = tt1 + tt2 + tt3
                  scrch2(i)  = 1.0 - tt1 - (scrch1(i)-tt1)*tt4
                  scrch2(i)  = SQRT( (1.0-scrch1(i))*scrch2(i) )
                  scrch3(i)  = 1.0-scrch1(i)*tt4
                  
                  lamb1x     = physic(i,j,k)%velx*(1.0-tt4) &
                  - cs_s(i)*scrch2(i)
                  lamb5x     = physic(i,j,k)%velx*(1.0-tt4) &
                  + cs_s(i)*scrch2(i)
                  
                  scrch2(i)  = 1.0 - tt2 - (scrch1(i)-tt2)*tt4
                  scrch2(i)  = SQRT( (1.0-scrch1(i))*scrch2(i) )
                  
                  lamb1y     = physic(i,j,k)%vely*(1.0-tt4) &
                  - cs_s(i)*scrch2(i)
                  lamb5y     = physic(i,j,k)%vely*(1.0-tt4) &
                  + cs_s(i)*scrch2(i)
                  
                  scrch2(i)  = 1.0 - tt3 - (scrch1(i)-tt3)*tt4
                  scrch2(i)  = SQRT( (1.0-scrch1(i))*scrch2(i) )
                  
                  lamb1z     = physic(i,j,k)%velz*(1.0-tt4) &
                  - cs_s(i)*scrch2(i)
                  lamb5z     = physic(i,j,k)%velz*(1.0-tt4) &
                  + cs_s(i)*scrch2(i)

                  lambdx(i) = MAX( ABS(lamb1x), ABS(lamb5x) )
                  lambdy(i) = MAX( ABS(lamb1y), ABS(lamb5y) )
                  lambdz(i) = MAX( ABS(lamb1z), ABS(lamb5z) )
               enddo
               
               do i= 1, nx
                  dtestmax = MAX( lambdx(i)/(xznr(i)-xznl(i)),               &
                  lambdy(i)/(yznr(j)-yznl(j)),               &
                  lambdz(i)/(zznr(k)-zznl(k)) ) / scrch3(i)
                  
                  if (dtestmax > dtcc_omp(tid)) then
                     ict(tid)   = i
                     jct(tid)   = j
                     kct(tid)   = k
                     dtcc_omp(tid) = dtestmax
                     dtt(tid) = dtestmax
                     soundt(tid) = cs_s(i)
                  endif  
                  
                  dtcc_thread(i,j,k) = max(dtcc, dtcc_omp(tid))
                  ! if (dtcc == dtcc_omp(tid)) then
                  !    print*,'dtcc:',cfl/dtcc/SQRT(3.0),' in ',i,j,k
                  ! endif

                  if ((basenm(1:1) == 'M' .and. basenm(2:2) == 'C').or.basenm(1:1) == 'J') then
                     xe  =  physic(i,j,k)%denstye/physic(i,j,k)%densty
                     tem =  physic(i,j,k)%pres*m_0*cl**2/ (physic(i,j,k)%densty*    &
                            (xe/me + (1.0-xe)/mp))/kb

                    if (tem < 1.e4) then
                       dedt = 0.
                       dtco  = 1.

                    elseif (tem >= 1.e4 .and. tem < 1.e5) then
                       dedt = (physic(i,j,k)%densty*rho_0a*                     &
                               (xe/me + (1.0-xe)/mp))**2/m_0**2 * 7.e-27*tem
                       dedt = dedt/rho_0a/9.e20/3.e10*R_b
                       dtco = 0.1*physic(i,j,k)%energy/dedt
                     !   if (dtco < 1.d0) then
                     !    print*, '---------------- DTCO -------------------'
                     !    print*,'dtco:',dtco
                     !    print*,'i,j,k:',i,j,k
                     !   endif
                    elseif (tem >= 1.e5 .and. tem < 4.e7) then
                       dedt = (physic(i,j,k)%densty*rho_0a*                   &
                               (xe/me + (1.0-xe)/mp))**2/m_0**2 *                   &
                               7.e-19*tem**(-0.6)
                       dedt = dedt/rho_0a/9.e20/3.e10*R_b
                       dtco = 0.1*physic(i,j,k)%energy/dedt
                     !   if (dtco < 1.d0) then
                     !    print*, '---------------- DTCO -------------------'
                     !    print*,'dtco:',dtco
                     !    print*,'i,j,k:',i,j,k
                     !   endif
                    elseif (tem .ge. 4.e7) then
                       dedt = (physic(i,j,k)%densty*rho_0a*                     &
                               (xe/me + (1.0-xe)/mp))**2/m_0**2 *                 & 
                               3.e-27*tem**(0.5)
                       dedt = dedt/rho_0a/9.e20/3.e10*R_b
                       dtco = 0.1*physic(i,j,k)%energy/dedt
                    endif

                    if (dtco < dtcmin) then
                       dtcmin = dtco
                       temcmin = tem
                    endif

                    if (dtco < dtc_omp(tid)) then
                       dtc_omp(tid)  = dtco
                    endif

                    dtc_thread(i,j,k) = min(dtc, dtc_omp(tid))

                 endif
              enddo
           enddo
        ENDDO

!$OMP END DO
!$OMP MASTER

         ! dt MIN FOR THE WORKING NODE
         call par_barrier()

         dtcc = MAXVAL(RESHAPE(dtcc_thread, (/SIZE(dtcc_thread)/)))
         dtc = MINVAL(RESHAPE(dtc_thread, (/SIZE(dtcc_thread)/)))

         dt = MIN( cfl/dtcc/SQRT(3.0), 1.05*olddt, dtc)

         if (mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1 .and. mod(nstep,50) == 0) then
            write(6,*) '------------- TSTEP ----------------'
            write(6,*) '1.05*olddt,dtcc,dtc: ',1.05*olddt,cfl/dtcc/SQRT(3.0),dtc
            write(6,*) 'dt: ',dt
            write(6,*)  '-------------------------------------'
         endif

         ! Find the cell that gives the dtcc
         on = omp_get_num_threads()
         do k=1,on
            if (dtt(k) == dtcc) then
               ic = ict(k)
               jc = jct(k)
               kc = kct(k)
               v  = SQRT( physic(ic,jc,kc)%velx*physic(ic,jc,kc)%velx   &
               + physic(ic,jc,kc)%vely*physic(ic,jc,kc)%vely   &
               + physic(ic,jc,kc)%velz*physic(ic,jc,kc)%velz )
               sound = soundt(k)   
            endif
         enddo
!$OMP END MASTER         
      ENDIF
!$OMP END PARALLEL

      ! dt MIN AND BROADCAST FOR ALL NODES

      if (mpi_dims(1)*mpi_dims(2)*mpi_dims(3) > 1) then
         pmin  = reduceMinDT( ic, jc, kc, dt, sound, v )
      endif

      IF (mod(nstep, itstp) == 0) THEN
         if (mpi_dims(1)*mpi_dims(2)*mpi_dims(3) > 1) then
            if (cartrank == 0) then
               write(6,'(i5,2x,i4,2x,e9.2,3(2x,i5),e11.2,2x,e11.2,2x,e11.2)') &
               nstep, pmin, dt, ic, jc, kc, sound, v, timej
            endif
         else
            write(6,'(i5,2x,e9.2,3(2x,i5),e11.2,2x,e11.2,2x,e11.2)')      &
            nstep, dt, ic, jc, kc, sound, v, timej
         endif
      ENDIF

      deallocate(dtcc_thread,dtc_thread)

!----------------------------------------------------------------------------END
      END subroutine tstep
