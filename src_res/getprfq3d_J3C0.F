!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL Variables
! DESCRIPTION
! SEE ALSO
!===============================================================================

subroutine getprfq3d( basenm, nx, ny, nz, mnx1, mnx5, mny5, mnz1, mnz5,     &
                      small, smallp, vwi, vk, timej,  &
		                vwp, rhowp, uwp, yc, me, mp,    &
                      rkin, xzn, yzn, zzn, physic )

  USE tipos
  USE parallel
  USE workarrays

   implicit none
# include "types.h"

!------------------------------------------------------- Input Variables

  character, intent(in):: basenm*(*)

  integer(is):: nx, ny, nz, mnx1, mnx5, mny5, mnz1, mnz5, rkin

  real(rs)   :: small, smallp, me, mp, vwi, vk, timej, vwp, rhowp, uwp 

  real(rs)   :: xzn(*), yzn(*), zzn(*)

!------------------------------------------------------ Output Variables

  type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------------------- Local Variables

  integer(is):: tid, i, j, k, ll, it, itmax, nl, nb_nan
      
# ifdef OPENMP
    integer :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads,curth
# endif

  integer(is):: lista(3,maxele)

!     input/output for eos_synge
  real(rs)   :: op_s(mnx1),w(mnx1),v2(mnx1),                  &
                denstw_s(mnx1),energy_s(mnx1),momenx_s(mnx1), &
                momeny_s(mnx1),momenz_s(mnx1),pmin(mnx1),     &
                f(mnx1), df(mnx1), relerr(mnx1)
                
  real(rs)   :: errmax, errlim, xe_so, rho_so, tem_so, eps_so, p_so, &
                gam_so, h_so, cs_so, dpdrho_so, dpdeps_so

  real(rs)   :: nve, small2, HT, WT, temi1, temi2, temj1, temj2, temk1, temk2
  
  real(rs)   :: ze,zp,meff1,ge,gp,dge,dgp

  real(rs)   :: yc, theta1, theta2, distc, distc2

  real(rs), pointer :: P(:), rho_s(:),xe_s(:),tem_s(:),   & 
                       eps_s(:),p_s(:),gam_s(:),          &
                       h_s(:),cs_s(:),                    &
                       dpdrho_s(:),dpdeps_s(:)

   integer, allocatable :: nele_thread(:)
   integer :: num_threads

!------------------------------------------------------------------------- BEGIN

  small2 = SMALL*SMALL

  errlim = 1.e-6
  itmax = 100

  nele = 0

!$OMP PARALLEL
!$OMP SINGLE
   num_threads = omp_get_num_threads()
!$OMP END SINGLE
!$OMP END PARALLEL
   allocate(nele_thread(num_threads))
   nele_thread = 0

!$OMP PARALLEL DEFAULT(shared)& 
!$OMP PRIVATE(tid,curth,i,j,k,P,rho_s,xe_s,tem_s,eps_s)&
!$OMP PRIVATE(p_s,gam_s,h_s,cs_s,dpdrho_s,dpdeps_s)&
!$OMP PRIVATE(op_s,w,v2,denstw_s,energy_s,momenx_s)&
!$OMP PRIVATE(momeny_s,momenz_s)&
!$OMP PRIVATE(pmin,f,df,relerr,errmax,it,lista,nl)
      !If OMP reduction works: !$OMP REDUCTION(+:nele)

# ifdef OPENMP
    call my_omp_get_thread_num(curth)
    tid = curth+1
# else
    tid = 1
    num_threads = 1
# endif

  ! -4:mn5 vectors
  P => wa(tid)%wa200
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

!-------------------------------------------------------------------------

  nl=1
  lista(1,1) = 0
  lista(2,1) = 0
  lista(3,1) = 0

!$OMP DO

   do k= 1, nz

      do j= 1, ny

         do i= 1, nx

            p_s(i)      = physic(i,j,k)%pres
            denstw_s(i) = physic(i,j,k)%denstw
            energy_s(i) = physic(i,j,k)%energy
            momenx_s(i) = physic(i,j,k)%momenx
            momeny_s(i) = physic(i,j,k)%momeny
            momenz_s(i) = physic(i,j,k)%momenz
            xe_s(i)     = physic(i,j,k)%denstew/denstw_s(i)
            pmin(i)     = max( sqrt(momenx_s(i)**2+momeny_s(i)**2+momenz_s(i)**2)   &  
                           - energy_s(i) - denstw_s(i), smallp )
            tem_s(i)    =  physic(i,j,k)%pres / (physic(i,j,k)%densty*              &
                           (xe_s(i)/me + (1.0-xe_s(i))/mp))             

            if (pmin(i).gt.(1.6666667d0-1.0)*energy_s(i)) then
               if (i /= lista(1,nl) .or. j /= lista(2,nl) .or. k /= lista(3,nl)) then 
                  nele_thread(tid) = nele_thread(tid)+1
                  nl=nele_thread(tid)

                  lista(1,nl) = i
                  lista(2,nl) = j
                  lista(3,nl) = k
               endif

               print*,'-------------------------- Getprfq: PMIN PROBLEM --------------------------------'
               print*,'Momentum = ',momenx_s(i),momeny_s(i),momenz_s(i)
               print*,'Energy = ', energy_s(i)
               print*,'Conserved density  = ', denstw_s(i)
               print*,'Pmin_expression = ', sqrt(momenx_s(i)**2+momeny_s(i)**2+momenz_s(i)**2) - energy_s(i) - denstw_s(i)
               print*,'Smallp = ', smallp
               print*,'Pmin = ', pmin(i)
               print*,'(1.67-1)Energy = ', (1.6666667d0-1.0)*energy_s(i) 
               print*,'Pressure [code]:', p_s(i)
               print*,'Temperature [code, K]:', tem_s(i)
               print*,'mpi_coords:',mpi_coords(1), mpi_coords(2), mpi_coords(3)
               print*,'i,j,k:',i,j,k
               print*,'nele, tid:', nele_thread(tid), tid
               print*,'------------------------------------------------------------------------' 
            endif
         enddo
      
         it = 0
         errmax = 1. 

         DO WHILE ((errmax>errlim) .AND. (it < itmax))  

            it=it+1

            do i = 1, nx

               v2(i) = energy_s(i) + denstw_s(i) + p_s(i)
               v2(i) = ( momenx_s(i)*momenx_s(i) + momeny_s(i)*momeny_s(i) + momenz_s(i)*momenz_s(i)) & 
                        / (v2(i)*v2(i))
               w(i) = 1.0/sqrt(1.0-v2(i))
               rho_s(i) = denstw_s(i)/w(i)
               eps_s(i) = ( energy_s(i) + denstw_s(i)*(1.0-w(i)) +  &
                           p_s(i)*(1.0-w(i)*w(i)) ) / (denstw_s(i)*w(i))
               op_s(i) = p_s(i)

               if (eps_s(i).lt.0.0) then
                  if (i /= lista(1,nl) .or. j /= lista(2,nl) .or. k /= lista(3,nl)) then 
                     nele_thread(tid) = nele_thread(tid)+1
                     nl=nele_thread(tid)

                     lista(1,nl) = i
                     lista(2,nl) = j
                     lista(3,nl) = k
                  endif

                  print*,'----------- Getprfq:  Negative EPS problem ----------------------------'
                  print*, 'Cells ', i, j, 'in mpi_coords ',mpi_coords(1), mpi_coords(2), mpi_coords(3)
                  print*,'nele, it, eps_s: ', nele_thread(tid),it,eps_s(i)
                  print*,'------------------------------------------------------------------------' 
               endif
            enddo

            call eos_synge( 2, tid, 1, nx, me, mp, rho_s, xe_s, tem_s, eps_s,   &
                              p_s, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)

            errmax = 0.0

            ! NEWTON-RAPHSON for the pressure
            do i=1, nx 

               f(i) = p_s(i) - op_s(i)
               df(i) = v2(i)*cs_s(i)*cs_s(i) - 1.0
               p_s(i) = op_s(i) - f(i)/df(i)
               relerr(i) = abs((p_s(i)-op_s(i))/op_s(i))   

               if (p_s(i).le.pmin(i)) then
                  if (i /= lista(1,nl) .or. j /= lista(2,nl) .or. k /= lista(3,nl)) then 
                     nele_thread(tid) = nele_thread(tid)+1
                     nl=nele_thread(tid)

                     lista(1,nl) = i
                     lista(2,nl) = j
                     lista(3,nl) = k
                  endif

                  if (it < itmax) then
                     p_s(i) = op_s(i)
                  else
                     p_s(i) = smallp   
                  endif

                  print*,'---------------- Getprfq:  p < pmin problem ----------------------------'
                  print*, 'Cells ', i, j, 'in mpi_coords ',mpi_coords(1), mpi_coords(2), mpi_coords(3)
                  print*, 'nele,it,p_s(i): ', nele_thread(tid),it,p_s(i)
                  print*,'------------------------------------------------------------------------'

               endif

               if (relerr(i) > errmax) then
                  if (i /= lista(1,nl) .or. j /= lista(2,nl) .or. k /= lista(3,nl)) then
                     errmax = relerr(i)
                  endif
               endif

               if( errmax > errlim .and. it == itmax ) then
                  if (i /= lista(1,nl) .or. j /= lista(2,nl) .or. k /= lista(3,nl)) then
                     nele_thread(tid) = nele_thread(tid)+1
                     nl=nele_thread(tid)

                     lista(1,nl) = i
                     lista(2,nl) = j
                     lista(3,nl) = k
                  endif

                  print*,'------------- Getprfq3d (non convergence): errmax and itmax ---------------------'
                  print*, 'nele,it,errmax,errlim,itmax: ', nele_thread(tid),it,errmax,errlim,itmax
                  print*,'Cell',i,j,k, 'at', mpi_coords(1), mpi_coords(2), mpi_coords(3)
                  print*,'---------------------------------------------------------------------------------'
                  
                  p_s(i) = smallp

               endif

            enddo

         ENDDO

         !---------------------- end of iteration -------------------------
   
         do i=1,nx
            physic(i,j,k)%pres = p_s(i)
            physic(i,j,k)%opres = p_s(i)
            physic(i,j,k)%densty = rho_s(i)
            physic(i,j,k)%denstye = xe_s(i)*rho_s(i)
            physic(i,j,k)%eps = eps_s(i)
            physic(i,j,k)%velx = momenx_s(i)/(energy_s(i)+denstw_s(i)+p_s(i))
            physic(i,j,k)%vely = momeny_s(i)/(energy_s(i)+denstw_s(i)+p_s(i))
            physic(i,j,k)%velz = momenz_s(i)/(energy_s(i)+denstw_s(i)+p_s(i))
         enddo

    enddo
   enddo

!$OMP END DO
!$OMP END PARALLEL

   nele = sum(nele_thread)

!-------------------------------------------------------------------------

   ! if (cartrank == 0) then

      if (nele >= maxele) then
         print*,'GETPRFQ: max elements error, list full',cartrank        
         call ERROR( 'getprfq3d', 'list full', __FILE__, __LINE__  )
      endif

      if (nele > 0) then

         call intercambiaBoundPhysicPlane( nx, ny, nz, mnx5, mny5, mnz5, physic )
         call par_barrier()

!$OMP PARALLEL DEFAULT(shared)& 
!$OMP PRIVATE(tid,i,j,nve,temi1,temi2,temj1,temj2)&
!$OMP PRIVATE(rho_so,xe_so,tem_so,eps_so)&
!$OMP PRIVATE(p_so,ze,zp,meff1,dpdrho_so,ge,gp,dge,dgp)&
!$OMP PRIVATE(WT,HT)

!$OMP DO SCHEDULE(static)
         
         do ll= 1, nele

            i = lista(1,ll)
            j = lista(2,ll)
            k = lista(3,ll)

            if (ll == 1) then
               print*,'------------------------- Getprfq: recovery ----------------------------'
               print*,'mpi_coords:',mpi_coords(1), mpi_coords(2), mpi_coords(3)
               print*,'nele: ', nele
               print*,'------------------------------------------------------------------------' 
            endif

            physic(i,j,k)%velx    = 0.0
            physic(i,j,k)%vely    = 0.0
            physic(i,j,k)%densty  = 0.0
            physic(i,j,k)%denstye = 0.0
            physic(i,j,k)%eps     = 0.0
            physic(i,j,k)%pres    = 0.0
      
            nve = 0.0
            if (physic(i-1,j,k)%eps /= 0.0) then
               nve= nve+1.0
            endif
            if (physic(i+1,j,k)%eps /= 0.0) then
               nve= nve+1.0
            endif
            if (physic(i,j-1,k)%eps /= 0.0) then
               nve= nve+1.0
            endif
            if (physic(i,j+1,k)%eps /= 0.0) then
               nve= nve+1.0
            endif
            if (physic(i,j,k-1)%eps /= 0.0) then
               nve= nve+1.0
            endif
            if (physic(i,j,k+1)%eps /= 0.0) then
               nve= nve+1.0
            endif
      
            if (nve == 0) then
               nve = 1.0
            else
               nve = 1.0 / nve
            endif
      
            physic(i,j,k)%velx = (physic(i-1,j,k)%velx+physic(i+1,j,k)%velx  &
            + physic(i,j-1,k)%velx+physic(i,j+1,k)%velx  &
            + physic(i,j,k-1)%velx+physic(i,j,k+1)%velx)*nve
      
            physic(i,j,k)%vely = (physic(i-1,j,k)%vely+physic(i+1,j,k)%vely &  
            + physic(i,j-1,k)%vely+physic(i,j+1,k)%vely &
            + physic(i,j,k-1)%vely+physic(i,j,k+1)%vely)*nve
            
            physic(i,j,k)%velz = (physic(i-1,j,k)%velz+physic(i+1,j,k)%velz &
            + physic(i,j-1,k)%velz+physic(i,j+1,k)%velz &
            + physic(i,j,k-1)%velz+physic(i,j,k+1)%velz)*nve
            
            physic(i,j,k)%densty = (physic(i-1,j,k)%densty+physic(i+1,j,k)%densty &
            + physic(i,j-1,k)%densty+physic(i,j+1,k)%densty &
            + physic(i,j,k-1)%densty+physic(i,j,k+1)%densty)*nve
            
            physic(i,j,k)%denstye = (physic(i-1,j,k)%denstye+physic(i+1,j,k)%denstye &
            + physic(i,j-1,k)%denstye+physic(i,j+1,k)%denstye  &
            + physic(i,j,k-1)%denstye+physic(i,j,k+1)%denstye)*nve
            
            physic(i,j,k)%eps = (physic(i-1,j,k)%eps+physic(i+1,j,k)%eps &
            + physic(i,j-1,k)%eps+physic(i,j+1,k)%eps &
            + physic(i,j,k-1)%eps+physic(i,j,k+1)%eps)*nve
      
            physic(i,j,k)%tracer = (physic(i-1,j,k)%tracer+physic(i+1,j,k)%tracer &
            + physic(i,j-1,k)%tracer+physic(i,j+1,k)%tracer  &
            + physic(i,j,k-1)%tracer+physic(i,j,k+1)%tracer)*nve

            if (physic(i-1,j,k)%eps /= 0.0) then
               temi1 =  physic(i-1,j,k)%pres / (physic(i-1,j,k)%densty*              &
               (physic(i-1,j,k)%denstye/physic(i-1,j,k)%densty/me +                    &
               (1.0-physic(i-1,j,k)%denstye/physic(i-1,j,k)%densty)/mp))
            else
               temi1 = 0.0
            endif
            if (physic(i+1,j,k)%eps /= 0.0) then
               temi2 =  physic(i+1,j,k)%pres / (physic(i+1,j,k)%densty*              &
               (physic(i+1,j,k)%denstye/physic(i+1,j,k)%densty/me +                    &
               (1.0-physic(i+1,j,k)%denstye/physic(i+1,j,k)%densty)/mp))
            else
               temi2 = 0.0
            endif
            if (physic(i,j-1,k)%eps /= 0.0) then
               temj1 =  physic(i,j-1,k)%pres / (physic(i,j-1,k)%densty*              &
               (physic(i,j-1,k)%denstye/physic(i,j-1,k)%densty/me +                    &
               (1.0-physic(i,j-1,k)%denstye/physic(i,j-1,k)%densty)/mp))
            else
               temj1 = 0.0
            endif
            if (physic(i,j+1,k)%eps /= 0.0) then
               temj2 =  physic(i,j+1,k)%pres / (physic(i,j+1,k)%densty*              &
               (physic(i,j+1,k)%denstye/physic(i,j+1,k)%densty/me +                    &
               (1.0-physic(i,j+1,k)%denstye/physic(i,j+1,k)%densty)/mp))
            else
               temj2 = 0.0
            endif
            if (physic(i,j,k-1)%eps /= 0.0) then
               temk1 =  physic(i,j,k-1)%pres / (physic(i,j,k-1)%densty*              &
               (physic(i,j,k-1)%denstye/physic(i,j,k-1)%densty/me +                    &
               (1.0-physic(i,j,k-1)%denstye/physic(i,j,k-1)%densty)/mp))
            else
               temk1 = 0.0
            endif
            if (physic(i,j,k+1)%eps /= 0.0) then
               temk2 =  physic(i,j,k+1)%pres / (physic(i,j,k+1)%densty*              &
               (physic(i,j,k+1)%denstye/physic(i,j,k+1)%densty/me +                    &
               (1.0-physic(i,j,k+1)%denstye/physic(i,j,k+1)%densty)/mp))
            else
               temk2 = 0.0
            endif
      
            rho_so =  physic(i,j,k)%densty 
            xe_so  =  physic(i,j,k)%denstye/physic(i,j,k)%densty
            eps_so =  physic(i,j,k)%eps 
      
            tem_so  = (temi1+temi2+temj1+temj2+temk1+temk2)*nve

            ze = me / tem_so
            zp = mp / tem_so
            meff1 = (1.0-xe_so)/mp + xe_so/me
            dpdrho_so = tem_so*meff1
            p_so = rho_so*dpdrho_so
            ge=5.0/2.0/ze+sqrt(9.0/4.0/ze**2+1.0)
            gp=5.0/2.0/zp+sqrt(9.0/4.0/zp**2+1.0)
            dge=-5.0/2.0/ze**2-9.0/4.0/ze**3/ &
               sqrt(9.0/4.0/ze**2+1.0)
            dgp=-5.0/2.0/zp**2-9.0/4.0/zp**3/ &
               sqrt(9.0/4.0/zp**2+1.0)
            
            eps_so = (1.0-xe_so)*gp + xe_so*ge - dpdrho_so - 1.0
      
            physic(i,j,k)%pres  = p_so
      
            physic(i,j,k)%eps   = eps_so
      
            physic(i,j,k)%opres = physic(i,j,k)%pres
            
            WT  = 1.0/ SQRT(1.0 - (physic(i,j,k)%velx*physic(i,j,k)%velx+ &
            physic(i,j,k)%vely*physic(i,j,k)%vely  + &
            physic(i,j,k)%velz*physic(i,j,k)%velz))
            
            HT  = 1.0 + physic(i,j,k)%eps + physic(i,j,k)%pres/physic(i,j,k)%densty
            
            physic(i,j,k)%denstw  = physic(i,j,k)%densty * WT
            physic(i,j,k)%denstew = physic(i,j,k)%denstye * WT
            physic(i,j,k)%momenx  = physic(i,j,k)%denstw*HT*WT*physic(i,j,k)%velx
            physic(i,j,k)%momeny  = physic(i,j,k)%denstw*HT*WT*physic(i,j,k)%vely
            physic(i,j,k)%momenz  = physic(i,j,k)%denstw*HT*WT*physic(i,j,k)%velz
            physic(i,j,k)%energy  = physic(i,j,k)%denstw*HT*WT                  &
            - physic(i,j,k)%pres - physic(i,j,k)%denstw

         enddo
!$OMP END DO
!$OMP END PARALLEL
      endif
   !sendif

   deallocate(nele_thread)

!----------------------------------------------------------------------------END
   
      END subroutine getprfq3d
