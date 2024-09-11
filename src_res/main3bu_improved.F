!===============================================================================
!     NAME F90 SPECIFICATION ARGUMENTS GLOBAL Variables DESCRIPTION SEE
!     ALSO This is the main program for two-dimensional simulations of
!     relativistic flows. The numerical method used belongs to the class
!     of high-resolution shock-capturing techniques and is based on an
!     approximate Riemann solver. A parabolic reconstruction technique
!     is used to get third order accuracy in space.  Advance in time is
!     performed simultaneously in the two spatial directions Adapted
!     from Prometheus.  MODIFIED 30-OCT-1995 to include the integration
!     of a passive scalar field tracing the jet particles.  MODIFIED
!     22-NOV-1995 to change GETPRFQ to GETPRFQ2D MODIFIED 21-FEB-1996 to
!     introduce commons /beam/, /ambnt/ and /gravty/ and add the call to
!     PSGRAV MODIFIED 20-DEC-1996 to include the call to INPUT_J
!     containing initial data for jet simulations MODIFIED 20-MAY-1998
!     to include a call to routine CLEANATM after writing restarting
!     files MODIFIED 3-MAR-2003. Started change to 3D. Manolo.  MODIFIED
!     22-FEB-2006. Moving and cleaning to Fortran90. BSC MODIFIED
!     22-MAR-2006. Paralelizing 3D. BSC MODIFIED JUL-2010. Paralelizing
!     3D. Carmen Aloy
!===============================================================================

      PROGRAM RATPENAT

      ! Including external modules and libraries
      
      USE tipos
      USE memoria
      USE parallel
      USE workarrays

      ! Preprocessing

#ifdef TIEMPO
      USE tiempo
#endif

#if defined(PARALELO) && defined(MN) 
      USE mpi
#endif
      IMPLICIT NONE

# include "types.h"

!-------------------------------------------------------Local Variables

      logical(ls):: NoSalir, NoEixir

      integer :: part,ndirs,freqpl,final,H_openr,H_opend,rst_file,rst_rec, rst_filed

      character(20) :: ndirsch
      character(8)  :: datfil
      character(7)  :: rstfil, outfil, tdlfil
      character(4)  :: basenm
      character(3)  :: suffixd
      character(2)  :: suffix

      !  Local sizes for nx, ny, nz
      !  nyh -> Global homogeneous size
      !  nyt -> Global size
      !  nyold -> used for restart and increase of grid
      integer(is):: nx, nx2, nxo, &
      nz, nz2, nzo, &
      ny, nny, ny0, nyold, nyh, nys, nyt, &
      imin, imax, jmin, jmax, kmin, kmax, &
      nxtot, nytot, nztot

      integer(is):: dxoff, dyoff, dzoff

      !  Local & maximum sizes for workarrays (+1, +5, and +6)
      integer(is):: mnx1, mnx5, mnx6, &
      mny1, mny5, mny6, &
      mnz1, mnz5, mnz6, &
      mn1,  mn5,  mn6

      integer(is):: ierr, i, j, k,                    &
      xyswp, nbegin, iwrtd, rkin,                     &
      igeomx, igeomy, igeomz,                         &
      bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz, &
      nend, irstrt, nrstrt, nout, iapprs, irk, igodu, &
      itstp, niter, nsdim, nmod,                      &
      ixsym, izsym, iextrax, iextraz,                 &
      nstep, nout1, nrst, itdelay, ipr, ms, mclouds,  &
      brt
      
      real(rs)   :: alpha_R, beta_R, timed, dtimed,    &
      gamma_ad, alpha, tmax, trstrt, tout, cfl, dtini, &
      small, smlrho, smallp, smallu,                   &
      gridlx, gridlx2, gridly, gridlz,                 &
      gridlz2, epsiln, omg1, omg2,                     &
      delx, dely, delz, temb,                          &
      pb, rhob, rhoeb, velb, machb, gamb, csb2,        &
      dk1, dk2, psh1, psh2, mm1, mm2,                  &
      nn, oang, t0, t1, vptbd, rptbd, pptbd,           &
      gamma_e, gamma_c, gamma_ce,                      &
      eps_e, eps_c, rho_e, rho_c, p_e, p_c, w_c,       &
      nu_ce, xd, timej, dt, tout1, trst,               &
      alpn, betn, gamn, rbn, radius, cdenst,           & 
      q0, mloss0, fg, rm, tem0, tem1, epsb, hb, lofb,  &
      lk, me, mp
      
      real(rs)   :: pc, me_cgs, c, mp_cgs, R_b, rho_0a, m_0, kb, meff, &
      KB1, KB2
      
      real(rs)   :: A0, pha(25), frs(25), fra(25)
      
      ! Definitions for microquasar and wind simulations 
      real(rs)   :: uwi, vwi, vk, shift, rorb, rad, sangl, rangl, delc, &
      dist, ran1, vwp, rhowp, uwp, yc, rstar, csp

      integer(is) :: idum, ncl, ic(2000), jc(2000), kc(2000)
      integer(is) :: ii, jj
      
      real(rs), allocatable:: xznl(:), xzn(:), xznr(:), &
      yznl(:), yzn(:), yznr(:), &
      zznl(:), zzn(:), zznr(:), &
      grav(:), rhoa(:),         &
      pa(:),   dpady(:)

      type(Tphysic), allocatable:: physic(:,:,:)
      
      type(TboundPosX), allocatable:: boundPosX(:)
      type(TboundPosY), allocatable:: boundPosY(:)
      type(TboundPosZ), allocatable:: boundPosZ(:)

      integer :: ticks, ticks_second, ticks_max
      integer :: seconds_actual, seconds_base, seconds_limit, salva_rst

      integer :: omp_get_num_threads, omp_get_max_threads, curth

      integer :: ofznx, onznx, ofzny, onzny, ofznz, onznz

      integer :: numArgs, xCPUs, yCPUs, zCPUs 

      character(len=256) :: program_name, tmpStr

      type(Tbounds) :: fxcortex(:,:,:), nxcortex(:,:,:), fycortex(:,:,:), nycortex(:,:,:), fzcortex(:,:,:), nzcortex(:,:,:)
      
!--------------------------------------------------------------------------BEGIN
      
      ! Process the command line
      outRoot =''
      path_wkd =''
      numArgs = command_argument_count()

      call get_command_argument(0, program_name)
      if (numArgs.ne.5) then
         print*, 'Usage: ', TRIM(program_name), &
         'OutRoot, Working-Path xCPUS yCPUS zCPUS'
         STOP ' '
      endif
      call get_command_argument(1, outRoot)
      call get_command_argument(2, path_wkd)
      call get_command_argument(3, tmpstr)
      read(tmpstr, *) xCPUs
      call get_command_argument(4, tmpstr)
      read(tmpstr, *) yCPUs
      call get_command_argument(5, tmpstr)
      read(tmpstr, *) zCPUs

      ! Constants
      pc = 1 / 3.24254215e-19   
      me_cgs = 9.1095e-28
      c = 2.99792458e10
      kb = 1.38066e-16
      mp_cgs = 1836.11 * me_cgs
      KB1 = 7.5e20
      KB2 = 4.4e-10
      kb= 1.38066e-16
      meff= 1.0/9.1095e-28 * 2.0/1836.11

      !  Read initial parameters
      call input(  suffix, suffixd,                & 
      basenm, rst_file, rst_rec,                   &
      nend, irstrt, nrstrt, nout, itstp, iapprs,   &
      irk, igodu, niter, nsdim, ixsym, izsym,      &
      iextrax, iextraz, nx, nx2, nxo, ny, nyold, nyh,   &
      nz, nz2, nzo, mnx1, mnx5, mnx6, mnz1, mnz5, mnz6, &
      igeomx, igeomy, igeomz,                      &
      xd, alpha, tmax, trstrt, tout, cfl, dtini,   &
      small, smlrho, smallp, smallu,               &
      gridlx, gridlx2, gridlz, gridlz2, gridly,    &
      epsiln, omg1, omg2, freqpl )
      itdelay = 0

      nxtot = nx
      nytot = ny 
      nztot = nz  
      
      ! Read initial parameters for the environment and the jet
      if (basenm(1:1) == 'J' .OR. basenm(1:1) == '3'.OR. basenm(1:1) == 'M') then

         call input_j(  basenm,                                           &
                  suffixd, tdlfil, itdelay, gamma_ad, velb, rhob, rhoeb,  &
                  temb, machb, radius, cdenst, shift, q0, fg, ms,         &
                  dk1, dk2, rm, psh1, psh2, nn, mm1, mm2, ipr,            &
                  alpn, betn, gamn, rbn, brt, tem0, tem1, oang, csb2,     &
                  pb, mclouds, t0, t1, vptbd, rptbd, pptbd)
 
         ! Code units
         R_b    = radius * pc
         rho_0a = cdenst*mp_cgs
         m_0    = rho_0a * R_b**3
         me     = me_cgs / m_0
         mp     = mp_cgs / m_0
         psh1   = psh1*pc/R_b
         psh2   = psh2*pc/R_b 
         rm     = rm*pc/R_b
         rbn    = rbn*pc/R_b
         shift  = shift*pc/R_b    
         temb   = temb*kb/(m_0*c**2)
      else
         R_b    = pc
         rho_0a = cdenst*mp_cgs
         m_0    = rho_0a * R_b**3
         me     = me_cgs / m_0
         mp     = mp_cgs / m_0
      endif

      ! MPI & OMP initialitation
      call mpi_init3D (xCPUs, yCPUs, zCPUs, nsdim)
      !   print*, ' MAIN:  after MPI_INIT'
      !    call flush(6)
      call omp_init
      !   print*, ' MAIN:  after OMP_INIT'
      !    call flush(6)
      call distribute_ompth (nsdim) 
      !   print*, ' MAIN:  after DISTRIBUTE OMPTH'
      !    !call flush(6)
      call decompose_mpi_domain (xCPUs, yCPUs, zCPUs, nx, ny, nz, nx2, nz2, nyold, nyh)
      !   print*, ' MAIN:  after DECOMPOSE_MPI_DOMAIN'
      !    !call flush(6)

      ! ny is going to be the total size of a local grid
      ! ny0 is the homogeneous size of the local grid (used with growing jets)
      
      ! Split of the grid between the parallelized nodes
      call reparto_trabajo( basenm, nx*xCPUs, ny*yCPUs, nz*zCPUs, mn1, mn5, mn6, mnx1, mnx5, mnx6, &
      mny1, mny5, mny6, mnz1, mnz5, mnz6, ny0 )

      if (cartrank==0) then
         print*,'--------------------------- MAIN ----------------------------------'
         print*,'Whole grid nx, ny, nz:',nxtot, nytot, nztot
         print*,'After reparto_trabajo nx, ny, nz, mnx5, mny5, mnz5, cartrank:'
         print*,nx,ny,nz, mnx5, mny5, mnz5,cartrank
         print*,'-------------------------------------------------------------------'
      endif

      ! Allocation of grid, boundary and physical variables
      allocate( xznl(mnx1), xzn(mnx1), xznr(mnx1),   &
      yznl(mny1), yzn(mny1), yznr(mny1),             &
      zznl(mnz1), zzn(mnz1), zznr(mnz1),             &
      grav(-4:mny5), rhoa(-4:mny5),                  &
      pa(-4:mny5), dpady(-4:mny5), stat=ierr )

      ! Feasible errors in allocation
      if (ierr /= 0) then
         print*, 'Main ERROR. Allocation of main variables xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, grav, rhoa, pa, dpady'
      endif

      allocate( physic(-4:mnx5,-4:mny5,-4:mnz5), stat=ierr )

      if (ierr /= 0) then
         print*, 'Main ERROR. Allocation of main variable physic', ierr, cartrank
      endif

      allocate( boundPosX(-4:5), boundPosY(-4:5), boundPosZ(-4:5), stat=ierr )
      
      if (ierr /= 0) then
         print*, 'Main ERROR. Allocation of main variables boundPosX, boundPosY, boundPosZ'
      endif

#ifdef HDF
      call rrhdf5_init(H_openr, H_opend)
#endif  

!-------------------------------------------- STARTING/RESTARTING THE SIMULATION
      
      imin  = mpi_coords(1) * nx       
      jmin  = mpi_coords(2) * ny       
      kmin  = mpi_coords(3) * nz      

      IF (irstrt == 0 .or. &
         (irstrt == 1 .and. (ny*yCPUs > nyold .or. nx*xCPUs > nxo .or. nz*zCPUs > nzo)) ) THEN

         ! CONSTRUCT A NEW MODEL
         ! NEW simulation
         if (irstrt == 0) then

            timej  = 0.0
            nstep  = 0
            dt     = dtini
            nbegin = 1
            
            rst_file  = 0
            rst_rec   = 0
            rst_filed = 0

            if (suffix(2:2) /= 'A') then
               call ERROR( 'MAIN', 'Check input file suffix', __FILE__, __LINE__ )
            endif

            if (cartrank == 0) then
               write(6,*) ' '
               write(6,*) 'MAIN: Constructing new initial model'
               write(6,*) ' '
            endif  

            call grid( basenm, igeomx, nx, nx2, ny0, nyh, ny, nz, nz2,    &                                     
               nsdim, ixsym, izsym, iextrax, iextraz,             &                            
               delx, dely, delz, small,                           &                                
               gridlx, gridlx2, gridly, gridlz, gridlz2,          &                                   
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &                                      
               xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, &                                            
               xCPUs, yCPUs, zCPUs )                                                                                            
                                                                                                                                  
            call getWAMemorySpace( mn6, mn5, mn1, mnx1, mny1, mnz1, nx, ny, nyt, nz, mpi_dims(1), mpi_dims(2), mpi_dims(3) )     

            if (cartrank==0) then
               print*,'MAIN: grid done'
            endif

         ! RESTART simulation
         else if (irstrt == 1) then

            ! Read restart file (R)
            rstfil  = basenm//'R'//suffix
            datfil  = basenm//'D'//suffixd

            if (cartrank == 0) then
               write(6,*)
               write(6,*) 'MAIN: Restarting from ', rstfil
               write(6,*)
            endif
         
            call grid(basenm, igeomx, nx, nx2, ny0, nyh, ny, nz, nz2,    &
                        nsdim, ixsym, izsym, iextrax, iextraz,             &
                        delx, dely, delz, small,                           &
                        gridlx, gridlx2, gridly, gridlz, gridlz2,          &
                        bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &
                        xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, &
                        xCPUs, yCPUs, zCPUs )

            call getWAMemorySpace( mn6, mn5, mn1, mnx1, mny1, mnz1, nx, ny, nyt, nz, mpi_dims(1), mpi_dims(2), mpi_dims(3) )

#ifdef HDF         
            call restrt( basenm, 1, mn6, mn5, mn1, mnx1, mnx5, &
            mny1, mny5, mnz1, mnz5, freqpl, dtini,             &
            H_openr, H_opend, rst_file, rst_rec, rst_filed,    &
            datfil, rstfil, outfil, tdlfil, suffixd, suffix,   &
            nrst, nout1,                                       & 
            nx, nxo, ny0, ny, nyold, nyh, nz, nzo, nyt,        &
            nxtot, nytot, nztot,                               &
            nsdim, igodu, igeomx, igeomy, igeomz, nstep,       &
            bndmxz, bndmny, bndmxy, bndmnz, bndmxx, bndmnx,    &
            small, epsiln, omg1, omg2, boundPosY,              &
            dt, dk1, dk2, psh1, psh2, temb, machb, csb2,       &
            nn, mm1, mm2, tem0, tem1,                          &          
            oang, ipr, ms, nmod, A0, fra, frs, pha,            &
            dely, delz, tout1, trst, gamma_ad,                 &
            gridlx, gridly, gridlz, pb, rhob, rhoeb, velb,     &
            timej, grav, rhoa, pa, dpady, uwi, vwi, vk,        &
            radius, cdenst, shift, fg, rm, brt, rorb,          &
            t0, t1, vptbd, rptbd, pptbd, alpn, betn, gamn,     &
            rbn, m_0, R_b, c, pc, rho_0a, q0, me, mp,          &
            xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, &
            delx, iextrax, iextraz, gridlx2, gridlz2,          &
            nx2, nz2, ixsym, izsym, ncl, ic, jc, kc, mclouds,  &
            physic, xCPUs, yCPUs, zCPUs )
#endif 
         endif
      
         dxoff  = mpi_coords(1) * nx      
         dyoff  = mpi_coords(2) * ny       
         dzoff  = mpi_coords(3) * nz     

         ! NEW simulation
         if (irstrt == 0) then

            call init( basenm, nmod, nstep,                  &
                  nx, ny, nyt, nz, mnx1, mnx5, mny1, mny5, mnz1, mnz5, &
                  mn1, mn5, mn6,                                       &
                  nyold, nxtot, nxo, nztot, nzo,                       &
                  dtini, gamma_ad, pb, velb, rhob, rhoeb, temb,        &
                  machb, gamb, dk1, dk2, psh1, psh2, mm1, mm2,         &
                  nn, fg, rm, brt, tem0, tem1,                         &
                  ipr, oang, ms, gridlx, gridly, timej, dt,            &
                  gamma_e, gamma_c, gamma_ce, eps_e, eps_c,            &
                  rho_e, rho_c, p_e, p_c, w_c, nu_ce, xd, alpha,       &
                  delx, dely, delz,xznl, xzn, xznr, yznl, yzn, yznr,   &
                  zznl, zzn, zznr, pha, frs, fra, rhoa, pa, dpady,     &                       
                  A0, uwi, vwi, vwp, rhowp, uwp, vk, shift, rorb,      &                                                            
                  rho_0a, R_b, pc, c, kb, m_0, me, mp, q0,             &
                  mloss0, ncl, nz2, ic, jc, kc, physic,                &
                  fxcortex, nxcortex, fycortex, nycortex, fzcortex, nzcortex )

            if (cartrank==0) then
               print*, '------------------ MAIN: after init ------------------'
               print*, 'mpi_coords:', mpi_coords(1),mpi_coords(2),mpi_coords(3)
               print*, 'cdenst,fg,physic(1,1,1)%pres,physic(1,1,1)%densty:'
               print*,cdenst,fg,physic(1,1,1)%pres,physic(1,1,1)%densty
               print*, '------------------------------------------------------'
            endif

         ! RESTART simulation with GRID INCREASE
         else if (dyoff >= nyold .or. dxoff < (nxtot-nxo)/2 .or. dxoff >= (nxtot-nxo)/2 + nxo  &
                     .or. dzoff < (nztot-nzo)/2 .or. dzoff >= (nztot-nzo)/2 + nzo                   & 
                     .or. (dyoff + ny) > nyold .or. (dxoff + nx) > (nxtot-nxo)/2 + nxo                 &
                     .or. (dzoff + nz) > (nztot-nzo)/2 + nzo) then

            call grid( basenm, igeomx, nx, nx2, ny0, nyh, ny, nz, nz2,    &
                        nsdim, ixsym, izsym, iextrax, iextraz,             &
                        delx, dely, delz, small,                           &
                        gridlx, gridlx2, gridly, gridlz, gridlz2,          &
                        bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &
                        xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, &
                        xCPUs, yCPUs, zCPUs )
            
            call init_ext( basenm, nmod, nstep,                  &
            nx, ny, nyt, nz, mnx1, mnx5, mny1, mny5, mnz1, mnz5, &
            mn1, mn5, mn6,                                       &
            nyold, nxtot, nxo, nztot, nzo,                       & 
            dtini, gamma_ad, pb, velb, rhob, rhoeb, temb,        &
            machb, gamb, dk1, dk2, psh1, psh2, mm1, mm2,         & 
            nn, fg, rm, brt, tem0, tem1,                         &
            ipr, oang, ms, gridlx, gridly, timej, dt,            &
            gamma_e, gamma_c, gamma_ce, eps_e, eps_c,            &
            rho_e, rho_c, p_e, p_c, w_c, nu_ce, xd, alpha,       &
            delx, dely, delz,xznl, xzn, xznr, yznl, yzn, yznr,   &
            zznl, zzn, zznr, pha, frs, fra, rhoa, pa, dpady,     &              
            A0, uwi, vwi, vwp, rhowp, uwp, vk, shift, rorb,      &                      
            rho_0a, R_b, pc, c, kb, m_0, me, mp, q0,             &
            mloss0, ncl, nz2, ic, jc, kc, physic )

         endif

         ! Comunicate positions bounds between processors
         call intercambiaBoundPos( nx, xzn, xznl, xznr, ny, yzn, yznl, yznr, nz, zzn, zznl, zznr, boundPosX, boundPosY, boundPosZ )
         call par_barrier( ) 

         ! Restarted simulation time step change
         if (irstrt == 1) then
            nbegin = nstep + 1
         endif      
         
         ! Pseudo-gravity for steady jets in pressure decreasing atmospheres
         if (basenm(1:1) == 'J' .OR. basenm(1:1) == 'P') then
         
            ! Comunicate physics bounds between processors
            call intercambiaBoundPhysic( physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy  )                                                            

            ! Pseudo-gravity
            call psgrav( basenm, igodu, nx, ny, nz, mnx5, mny5, mnz5, ipr, &
            nn, psh1, psh2, mm1, mm2, gamma_ad, dk1, dk2, pb,         &
            small, epsiln, omg1, omg2, me, mp,                        &
            fg, shift, rm, kb, m_0, c, tem0, tem1,                    &
            xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,        &
            boundPosX, boundPosY, boundPosZ, rhoa, pa, grav, physic )
            
         endif

         ! Mass injection from stars
         if (basenm(1:1) == 'J' .and. basenm(3:3) == 'M') then
            call stmass( igeomx, igeomy, nx, ny, nz, mnx1, &
            mnx5, mny1, mny5, mnz1, mnz5, alpn,  &
            betn, gamn, rbn, m_0, R_b, c, pc, shift, rho_0a, q0, &
            xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,   &
            physic )
         endif
         
         ! Initial time step
         if (nstep /= 0) then
            nbegin = nstep + 1
            nend   = nstep + nend
         endif

         ! File names
         suffixd = 'A'//suffix
         outfil  = basenm//'O'//suffix
         rstfil  = basenm//'R'//suffix
         datfil  = basenm//'D'//suffixd

      ! COMPLETE RESTART
      ELSE 
         ! Read restart file (R)
         rstfil  = basenm//'R'//suffix
         datfil  = basenm//'D'//suffixd

         if (cartrank == 0) then
            write(6,*)
            write(6,*) 'MAIN: Restarting from ', rstfil
            write(6,*)
         endif
         
         call grid( basenm, igeomx, nx, nx2, ny0, nyh, ny, nz, nz2,    &
                     nsdim, ixsym, izsym, iextrax, iextraz,             &
                     delx, dely, delz, small,                           &
                     gridlx, gridlx2, gridly, gridlz, gridlz2,          &
                     bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &
                     xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, &
                     xCPUs, yCPUs, zCPUs )

         call getWAMemorySpace( mn6, mn5, mn1, mnx1, mny1, mnz1, nx, ny, nyt, nz, mpi_dims(1), mpi_dims(2), mpi_dims(3) )


#ifdef HDF         
         call restrt( basenm, 1, mn6, mn5, mn1, mnx1, mnx5, &
         mny1, mny5, mnz1, mnz5, freqpl, dtini,             &
         H_openr, H_opend, rst_file, rst_rec, rst_filed,    &
         datfil, rstfil, outfil, tdlfil, suffixd, suffix,   &
         nrst, nout1,                                       & 
         nx, nxo, ny0, ny, nyold, nyh, nz, nzo, nyt,        &
         nxtot, nytot, nztot,                               &
         nsdim, igodu, igeomx, igeomy, igeomz, nstep,       &
         bndmxz, bndmny, bndmxy, bndmnz, bndmxx, bndmnx,    &
         small, epsiln, omg1, omg2, boundPosY,              &
         dt, dk1, dk2, psh1, psh2, temb, machb, csb2,       &
         nn, mm1, mm2, tem0, tem1,                          &          
         oang, ipr, ms, nmod, A0, fra, frs, pha,            &
         dely, delz, tout1, trst, gamma_ad,                 &
         gridlx, gridly, gridlz, pb, rhob, rhoeb, velb,     &
         timej, grav, rhoa, pa, dpady, uwi, vwi, vk,        &
         radius, cdenst, shift, fg, rm, brt, rorb,          &
         t0, t1, vptbd, rptbd, pptbd, alpn, betn, gamn,     &
         rbn, m_0, R_b, c, pc, rho_0a, q0, me, mp,          &
         xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, &
         delx, iextrax, iextraz, gridlx2, gridlz2,          &
         nx2, nz2, ixsym, izsym, ncl, ic, jc, kc, mclouds,  &
         physic, xCPUs, yCPUs, zCPUs ) 
#endif 
         call par_barrier( )


         ! Comunicate positions bounds between processors
         call intercambiaBoundPos( nx, xzn, xznl, xznr, ny, yzn, yznl, yznr, nz, zzn, zznl, zznr, boundPosX, boundPosY, boundPosZ )
         call par_barrier( ) 
         nbegin = nstep + 1

         ! ! Print initial conditions but when restarting (except gamma)
         ! if (cartrank == 0) then
         !    print*,' '
         !    print*,'---------------------------------------------------------------------'
         !    print*,'RESTART: variables as in input_j.dat'
         !    print*,gamma_ad,velb,rhob,rhoeb,temb,machb,radius,cdenst,shift,ms,q0,fg,rm
         !    print*,oang,dk1,dk2,psh1,psh2,nn,mm1,mm2,ipr,tem0,tem1,mclouds,alpn,betn,gamn
         !    print*,rbn,brt,t0,t1,vptbd,rptbd,pptbd,itdelay,suffixd
         !    print*,' '

         !    epsb = pb/(gamb - 1.)/rhob
         !    ! the gamma is computed in the eosynge and hence cannot be used here,
         !    ! we would need to redefine it with the eosynge of smth for getting the input gamma
         !    hb   = 1.+ gamb*epsb
         !    lofb = 1./sqrt(1.-velb*velb)
         !    ! lk = (hb*lofb-1.)*rhob*rho_0a*lofb*acos(-1.)*(R_b)**2*velb*c**3 
         !    lk = hb*lofb**2*rhob*rho_0a*acos(-1.)*(R_b)**2*velb*c**3 ! including rest mass
         !    print*,' '
         !    print*,'RESTART: jet physical input'
         !    print*,'gamb,epsb,hb,lofb,rhob,rho_0a,R_b,velb:'
         !    print*,gamb,epsb,hb,lofb,rhob,rho_0a,R_b,velb
         !    print*,' '
         !    print*,'Jet kinetic power (cgs)', lk
         !    print*,'---------------------------------------------------------------------'
         !    print*,' '
         ! endif

         ! Mass loading from stars  
         if (basenm(1:1) == 'J' .and. basenm(3:3) == 'M') then
            call stmass( igeomx, igeomy, nx, ny, nz, mnx1, &
            mnx5, mny1, mny5, mnz1, mnz5, alpn,  &
            betn, gamn, rbn, m_0, R_b, c, pc, shift, rho_0a, q0, &
            xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,   &
            physic )
         endif
         
      ENDIF

      IF (irstrt == 0 .and. cartrank == 0) THEN
         if (basenm(1:1)=='J'.or.basenm(1:1) =='M') then
            call showParameters( suffixd, rstfil, itdelay, &
            gamma_ad, velb, rhob, machb, pb, dk1, dk2, psh1, psh2, &
            nn, mm1, mm2, t0, t1, vptbd, rptbd, pptbd )

            epsb = pb/(gamb - 1.)/rhob
            hb   = 1.+ gamb*epsb
            lofb = 1./sqrt(1.-velb*velb)
            lk = (hb*lofb-1.)*rhob*rho_0a*lofb*acos(-1.)*(R_b)**2*velb*c**3
            print*,'----------------- MAIN: jet properties --------------------'
            print*,'gamb,epsb,hb,lofb,rhob,rho_0a,R_b,velb:'
            print*,gamb,epsb,hb,lofb,rhob,rho_0a,R_b,velb
            print*,'Jet kinetic power (cgs)', lk
            print*,'------------------------------------------------------------'
         endif
      ENDIF

      if (mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1) then
         print*,'------------ MAIN: before main loop ---------------'
         print*,'Time step: ',nstep
         print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell 324 217 217'                                                        
         print*,'Density and conserved density: ',physic(324,217,217)%densty,physic(324,217,217)%denstw
         print*,'Pressure: ',physic(324,217,217)%pres
         print*,'Internal energy: ',physic(324,217,217)%eps
         print*,'Velocity: ',physic(324,217,217)%velx,physic(324,217,217)%vely,physic(324,217,217)%velz
         print*,'Energy: ',physic(324,217,217)%energy
         print*,'Moment: ',physic(324,217,217)%momenx,physic(324,217,217)%momeny,physic(324,217,217)%momenz
         print*,'----------------------------------------------------'
      endif

!-------------------------------------------------------------- TIME EVOLUTION

      ! Time interval between time delay data files
      iwrtd = 0
      if (itdelay == 1) then
         iwrtd = 1
         dtimed = gridly/REAL(nyh)
         if (cartrank == 0) then
            write(6,*) 'Time interval between time delay files: ', dtimed
         endif
      endif

      timed = timej + dt
      nout1 = 0
      tout1 = 0.0
      nrst  = 0
      trst  = 0.0

      if (irk==2) then
         alpha_R = 2.0
         beta_R  = 1.0
      else
         alpha_R = 4.0
         beta_R  = 3.0
      endif

      nstep   = nbegin
      NoSalir = .TRUE.
      final   = 0
      nys = ny

      call par_barrier( )

      !------------------------------------------ START MAIN LOOP

      DO WHILE (nstep <= nend .AND. NoSalir )

         timej = timej + dt
         nout1 = nout1 + 1
         tout1 = tout1 + dt
         nrst  = nrst  + 1
         trst  = trst  + dt

         IF (timej > tmax) THEN
            NoSalir = .FALSE.
            final = 1
         ELSE
            NoSalir = .TRUE.

            do k= 1, nz
               do j= 1, ny
                  do i= 1, nx
                     physic(i,j,k)%omomenx = physic(i,j,k)%momenx
                     physic(i,j,k)%omomeny = physic(i,j,k)%momeny
                     physic(i,j,k)%omomenz = physic(i,j,k)%momenz
                     physic(i,j,k)%oenergy = physic(i,j,k)%energy
                     physic(i,j,k)%odenstw = physic(i,j,k)%denstw
                     physic(i,j,k)%odenstew = physic(i,j,k)%denstew
                     physic(i,j,k)%opres   = physic(i,j,k)%pres
                     physic(i,j,k)%otracer = physic(i,j,k)%tracer
                     
                     if (physic(i,j,k)%denstw <=0. .or. physic(i,j,k)%energy <= 0. .or.  &
                          physic(i,j,k)%pres <= 0.) then
                        print*,mpi_coords(1),mpi_coords(2),mpi_coords(3),i,j,k
                        print*,physic(i,j,k)%densty,physic(i,j,k)%denstye,physic(i,j,k)%pres
                        stop
                     endif
                  enddo
               enddo
            enddo

            ! Runge-Kutta loop
            DO rkin= 1, irk

               call intercambiaBoundPhysic( physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy  )

               if (mod(nstep,50) == 0  .and. rkin==1 .and. mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1) then
                  print*,'------------- MAIN: before sweeps ----------------'
                  print*,'Time step:', nstep
                  print*,'Rkin:', rkin
                  print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell 324 217 217'                                                        
                  print*,'Density and conserved density: ',physic(324,217,217)%densty,physic(324,217,217)%denstw
                  print*,'Pressure: ',physic(324,217,217)%pres
                  print*,'Internal energy: ',physic(324,217,217)%eps
                  print*,'Velocity: ',physic(324,217,217)%velx,physic(324,217,217)%vely,physic(324,217,217)%velz
                  print*,'Energy and conserved energy: ',physic(324,217,217)%energy
                  print*,'Moment: ',physic(324,217,217)%momenx,physic(324,217,217)%momeny,physic(324,217,217)%momenz
                  print*,'----------------------------------------------------'
               endif

               call sweeps( basenm, rkin, igodu, igeomx, igeomy, igeomz,            &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,       &
               nx, nz, nsdim, ny, nys, xyswp, mnx5, mny5, mnz5,      &
               alpha_R, beta_R, gamma_ad, dt,                        &
               timej, t0, t1, mm1, mm2, nn, brt,  vptbd, rptbd,      &
               pptbd, ms, delx, dk1, dk2, pb, psh1, psh2, rhob,      &
               rhoeb, velb, small, smlrho,                           &
               epsiln, omg1, omg2, uwi, vwi, shift, rorb,            &
               R_b, meff, rho_0a, me, mp, KB1,                       &
               KB2, c, pc, m_0, kb, mloss0, tem0, tem1, rm, fg,      &
               nmod, ipr, A0, gridly, fra, frs, pha, pa, rhoa,       &
               xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,    &
               boundPosX, boundPosY, boundPosZ, grav, physic, curth, &
               fxcortex, nxcortex, fycortex, nycortex, fzcortex, nzcortex )

               if (mod(nstep,50) == 0  .and. rkin==1 .and. mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1) then
                  ! print*,'------- MAIN: after sweeps, old parameters --------'
                  ! print*,'Time step: ',nstep
                  ! print*,'Rkin:', rkin
                  ! print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell 324 217 217'                                                        
                  ! print*,'Conserved density: ',physic(324,217,217)%odenstw                                                        
                  ! print*,'Leptonic conserved density: ',physic(324,217,217)%odenstew
                  ! print*,'Pressure: ',physic(324,217,217)%opres
                  ! print*,'Energy and conserved energy: ',physic(324,217,217)%oenergy
                  ! print*,'Moment: ',physic(324,217,217)%omomenx,physic(324,217,217)%omomeny,physic(324,217,217)%omomenz
                  print*,'---------- MAIN: after sweeps, new ones ------------'
                  print*,'Time step: ',nstep
                  print*,'Rkin:', rkin
                  print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell 324 217 217'                                                        
                  print*,'Density and conserved density: ',physic(324,217,217)%densty,physic(324,217,217)%denstw
                  print*,'Pressure: ',physic(324,217,217)%pres
                  print*,'Internal energy: ',physic(324,217,217)%eps
                  print*,'Velocity: ',physic(324,217,217)%velx,physic(324,217,217)%vely,physic(324,217,217)%velz
                  print*,'Energy and conserved energy: ',physic(324,217,217)%energy
                  print*,'Moment: ',physic(324,217,217)%momenx,physic(324,217,217)%momeny,physic(324,217,217)%momenz
                  print*,'----------------------------------------------------'
               endif

               ! Newton-Raphson for recovering the physical variables
               call getprfq3d( basenm, nx, ny, nz, mnx1, mnx5, mny5, mnz1, mnz5,     &
                        small, smallp, vwi, vk, timej,            &
			                vwp, rhowp, uwp, yc, me, mp,             &
			               rkin, xzn, yzn, zzn, physic )

               ! if (mpi_coords(1)==0 .and. mpi_coords(2)==0 .and. mpi_coords(3)==1) then
               !    print*,'------------- MAIN: after getprfq3d ----------------'
               !    print*,'Time step: ',nstep
               !    print*,'Rkin:', rkin
               !    print*,'Nodes ',mpi_coords(1),mpi_coords(2),mpi_coords(3),'and cell 324 217 217'                                                        
               !    print*,'Density and conserved density: ',physic(324,217,217)%densty,physic(324,217,217)%denstw
               !    print*,'Pressure: ',physic(324,217,217)%pres
               !    print*,'Internal energy: ',physic(324,217,217)%eps
               !    print*,'Velocity: ',physic(324,217,217)%velx,physic(324,217,217)%vely,physic(324,217,217)%velz
               !    print*,'Energy and conserved energy: ',physic(324,217,217)%energy
               !    print*,'Moment: ',physic(324,217,217)%momenx,physic(324,217,217)%momeny,physic(324,217,217)%momenz
               !    print*,'----------------------------------------------------'
               ! endif

            ENDDO

            ! FINAL OF THE TIME EVOLUTION
            if (final == 1) then
#ifdef HDF
               call restrt( basenm, 0, mn6, mn5, mn1, mnx1, mnx5, &
               mny1, mny5, mnz1, mnz5, freqpl, dtini,             &
               H_openr, H_opend, rst_file, rst_rec, rst_filed,    &
               datfil, rstfil, outfil, tdlfil, suffixd, suffix,   &
               nrst, nout1,                                       & 
               nx, nxo, ny0, ny, nyold, nyh, nz, nzo, nyt,        &
               nxtot, nytot, nztot,                               &
               nsdim, igodu, igeomx, igeomy, igeomz, nstep,       &
               bndmxz, bndmny, bndmxy, bndmnz, bndmxx, bndmnx,    &
               small, epsiln, omg1, omg2, boundPosY,              &
               dt, dk1, dk2, psh1, psh2, temb, machb, csb2,       &
               nn, mm1, mm2, tem0, tem1,                          &          
               oang, ipr, ms, nmod, A0, fra, frs, pha,            &
               dely, delz, tout1, trst, gamma_ad,                 &
               gridlx, gridly, gridlz, pb, rhob, rhoeb, velb,     &
               timej, grav, rhoa, pa, dpady, uwi, vwi, vk,        &
               radius, cdenst, shift, fg, rm, brt, rorb,          &
               t0, t1, vptbd, rptbd, pptbd, alpn, betn, gamn,     &
               rbn, m_0, R_b, c, pc, rho_0a, q0, me, mp,          &
               xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, &
               delx, iextrax, iextraz, gridlx2, gridlz2,          &
               nx2, nz2, ixsym, izsym, ncl, ic, jc, kc, mclouds,  &
               physic, xCPUs, yCPUs, zCPUs  )
#endif 
               NoSalir = .FALSE.
            endif

            ! TIME STEPS FOR TAKEN A SNAPSHOT (trstrt)
            if (nrst >= nrstrt .OR. trst >= trstrt) then
#ifdef HDF
               call restrt( basenm, 0, mn6, mn5, mn1, mnx1, mnx5, &
               mny1, mny5, mnz1, mnz5, freqpl, dtini,             &
               H_openr, H_opend, rst_file, rst_rec, rst_filed,    &
               datfil, rstfil, outfil, tdlfil, suffixd, suffix,   &
               nrst, nout1,                                       & 
               nx, nxo, ny0, ny, nyold, nyh, nz, nzo, nyt,        &
               nxtot, nytot, nztot,                               &
               nsdim, igodu, igeomx, igeomy, igeomz, nstep,       &
               bndmxz, bndmny, bndmxy, bndmnz, bndmxx, bndmnx,    &
               small, epsiln, omg1, omg2, boundPosY,              &
               dt, dk1, dk2, psh1, psh2, temb, machb, csb2,       &
               nn, mm1, mm2, tem0, tem1,                          &          
               oang, ipr, ms, nmod, A0, fra, frs, pha,            &
               dely, delz, tout1, trst, gamma_ad,                 &
               gridlx, gridly, gridlz, pb, rhob, rhoeb, velb,     &
               timej, grav, rhoa, pa, dpady, uwi, vwi, vk,        &
               radius, cdenst, shift, fg, rm, brt, rorb,          &
               t0, t1, vptbd, rptbd, pptbd, alpn, betn, gamn,     &
               rbn, m_0, R_b, c, pc, rho_0a, q0, me, mp,          &
               xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, &
               delx, iextrax, iextraz, gridlx2, gridlz2,          &
               nx2, nz2, ixsym, izsym, ncl, ic, jc, kc, mclouds,  &
               physic, xCPUs, yCPUs, zCPUs  )          
#endif
            endif

            ! Write test output file (O)
            if ((nsdim == 1 .OR. basenm == 'RSCH' .OR. basenm == 'GSCH' .OR.   &
            basenm == 'GSSR' .OR. basenm == 'ESSR' .OR. basenm == 'RSSR') &
            .AND. (nout1 >= nout .OR. tout1 >= tout)) then
               call pltout( basenm, nsdim, nx, ny, nz, mnx5, mny5, mnz5, &
               igeomx, nstep, gamma_ad, dt, timej,                 &
               xznl, xzn, xznr, yznl, yzn, yznr, zzn, physic,      &
               rstfil, outfil, tdlfil, datfil, suffix, suffixd,    &
               nout1, tout1 )
            endif

            ! Write time delay files
            if (itdelay == 1 .AND. iwrtd == 1) then
               ny = ny0
               call tdelay( basenm, igeomx, igeomy, igeomz, nx, ny, nz, &
               mnx1, mnx5, mny1, mny5, mnz1, mnz5, nstep, timej, gamma_ad,    &
               gridlx, gridly, gridlz, xznl, yznl, zznl, physic,  &
               datfil, rstfil, outfil, tdlfil, suffixd, suffix )

               if (cartrank == 0) then
                  write(6,*) 'MAIN: File for time delay', & 
                  'printed at time ', & 
                  timej
               endif
               timed = timed + dtimed
               iwrtd = 0
            endif

            ! NEW TIME STEP VALUE
            call tstep( basenm, nx, ny, nz, mnx5, mny5, mnz5, nsdim,   &
            igeomy, itstp, nstep, timej, cfl,                          &
            meff, R_b, rho_0a, kb, me, mp, c, m_0,                     &
            xznl, xznr, xzn, yznl, yznr, zznl, zznr, physic, dt)

            ! Check whether we need to add any time delay
            if (itdelay == 1) then
               dt = MIN(dt, dtimed)
               if ((cartrank == 0) .AND. (dt == dtimed)) then
                  write(6,*) 'MAIN: Time step reset to ', dt
               endif

               if (timej+dt >= timed) then
                  dt = timed - timej
                  if (cartrank == 0) then
                     write(6,*) 'MAIN: Time step reset to ', dt
                  endif
                  iwrtd = 1
               endif
            endif

            ! NEXT STEP (iterative variable controlling the loop)
            nstep = nstep + 1
            
         ENDIF
      ENDDO

      call par_barrier( )
      if (cartrank .ne. 0) then
         call freeWAMemorySpace( )
         
         deallocate( xznl, xzn, xznr, &
         yznl, yzn, yznr,        &
         zznl, zzn, zznr,        &
         rhoa, pa, dpady,        &
         grav )
         deallocate( physic )
         deallocate( boundPosX, boundPosY, boundPosZ ) 
      endif

      call par_end()
      print*, 'MAIN, END', cartrank
      
!----------------------------------------------------------------------------END

      END PROGRAM RATPENAT
