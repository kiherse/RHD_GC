!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!     This subroutine reads the initial parameters.
!     Adapted from Prometheus.
!     MODIFIED 16-JUL-1993 to include test 'RCSR'
!     MODIFIED 15-JUL-1994 to include test 'REFS' 
!     MODIFIED 07-MAY-1998 to change input format specifiers
!===============================================================================
      subroutine input( suffix, suffixd, basenm, rst_file, rst_rec,  &
      nend, irstrt, nrstrt, nout, itstp, iapprs,        &
      irk, igodu, niter, nsdim, ixsym, izsym,           &
      iextrax, iextraz, nx, nx2, nxo, ny, nyold, nyh,   &
      nz, nz2, nzo, mnx1, mnx5, mnx6, mnz1, mnz5, mnz6, &
      igeomx, igeomy, igeomz,                           &
      xd, alpha, tmax, trstrt, tout, cfl, dtini,        &
      small, smlrho, smallp, smallu,                    &
      gridlx, gridlx2, gridlz, gridlz2, gridly,         &
      epsiln, omg1, omg2, freqpl)
      USE parallel
      USE constantes
      USE workarrays
      implicit none
# include "types.h"
!------------------------------------------------------Input  Variables
      
!------------------------------------------------------Output Variables
      character(2):: suffix
      character(3):: suffixd
      character(4):: basenm

      integer(is) :: nend, irstrt, nrstrt, nout, itstp, iapprs,   &
      irk, igodu, niter, nsdim, ixsym, izsym,      &
      iextrax, iextraz, nx, nx2, nxo, ny, nyold, nyh,   &
      nz, nz2, nzo,  mnx1, mnx5, mnx6, mnz1, mnz5, mnz6, &
      igeomx, igeomy, igeomz, freqpl, rst_file, rst_rec

      real(rs)    :: xd, alpha, tmax, trstrt, tout, cfl, dtini, &
      small, smlrho, smallp, smallu,             &
      gridlx, gridlx2, gridlz, gridlz2, gridly,  &
      epsiln, omg1, omg2

!-------------------------------------------------------Local Variables
      character(72) :: text
      character(5)  :: txtxt
      character(8)  :: label, label1
      character(256):: fich

      integer(is)   :: ierr

      real(rs)      :: aux
      
!-------------------------------------------------------------------------BEGIN
      ierr  = 0
      aux   = 1.0/SQRT(2.0)
      txtxt = '.....'

      fich= TRIM(path_wkd)//'/'//'inpt.dat'
      
      open( 1, file=fich, form='formatted', status='old', iostat=ierr )
      if (ierr/=0) then
         call ERROR( 'input: Error openning ', fich, __FILE__, __LINE__ )
      endif

      read (1,'(a72)') text
# ifdef DEBUG
      write(utrace,*) text
# endif
      read (1,'(a72)') text
# ifdef DEBUG
      write(utrace,*) text
# endif
      read (1,'(a72)') text
# ifdef DEBUG
      write(utrace,*) text
      write(utrace,*)'------------------------------------------------'
# endif

!BASENM is the root for the output, plot and restart file names (roots
!RST_, RBW_, RSR_, RBWI, RRST, RRBW, RCSR, REFS, RSCH stand for special tests 
      read (1,'(a8,1x,a4)') label, basenm

# ifdef DEBUG
      write(utrace,*) label, txtxt, basenm
# endif
      label1 = 'basenm'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!RST_FILE
      read (1,'(a8,1x,i10)') label, rst_file
# ifdef DEBUG
      write(utrace,*) label, txtxt, rst_file
# endif
      label1 = 'rst_file'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!RST_REC
      read (1,'(a8,1x,i10)') label, rst_rec
# ifdef DEBUG
      write(utrace,*) label, txtxt, rst_rec
# endif
      label1 = 'rst_rec'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif


!XD gives the abcise of the position of the initial discontinuity at a 
!high y=gridly/2 (only for rotated one-dimensional tests).
      read (1,'(a8,1x,d23.2)') label, xd

# ifdef DEBUG
      write(utrace,*) label, txtxt, xd
# endif

      label1 = 'xd'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!ALPHA is the rotated angle (in degrees) of the initial discontinuity 
!with respect to the x=constant axis (counterclockwise) (only for
!rotated one-dimensional tests).
      read (1,'(a8,1x,d23.2)') label, alpha

# ifdef DEBUG
      write(utrace,*) label, txtxt, alpha
# endif

      alpha = 2.0*pi*alpha/360.0
      label1 = 'alpha'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!NEND is the total number of timesteps.
      read (1,'(a8,1x,i10)') label, nend
# ifdef DEBUG
      write(utrace,*) label, txtxt, nend
# endif
      label1 = 'nend'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!The program stops when time is >= TMAX. 
      read (1,'(a8,1x,d23.2)') label, tmax
# ifdef DEBUG
      write(utrace,*) label, txtxt, tmax
# endif
      label1 = 'tmax'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!SUFFIX is the suffix for the output and restart file names.
      read (1,'(a8,1x,a2)') label, suffix
# ifdef DEBUG
      write(utrace,*) label, txtxt, suffix
# endif
      label1 = 'suffix'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!SUFFIXD is the suffix for the 3D output file names.
      read (1,'(a8,1x,a3)') label,suffixd
# ifdef DEBUG
      write(utrace,*) label, txtxt, suffixd
# endif
      label1 = 'suffixd'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!If IRSTRT == 1, the program restarts from the file BASENM//'R'//SUFFIX 
      read (1,'(a8,1x,i10)') label, irstrt
# ifdef DEBUG
      write(utrace,*) label, txtxt, irstrt
# endif
      label1 = 'irstrt'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!A restart file is written every NRSTRT timesteps.
      read (1,'(a8,1x,i10)') label, nrstrt
# ifdef DEBUG
      write(utrace,*) label, txtxt, nrstrt
# endif
      label1 = 'nrstrt'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!A restart file is written every TRSTRT time units.
      read (1,'(a8,1x,d23.2)') label, trstrt
# ifdef DEBUG
      write(utrace,*) label, txtxt, trstrt
# endif
      label1 = 'trstrt'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!An output file is written every NOUT timesteps. (Only for 1D problems).
      read (1,'(a8,1x,i10)') label, nout
# ifdef DEBUG
      write(utrace,*) label, txtxt, nout
# endif
      label1 = 'nout'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!An output file is written every TOUT time units. (Only for 1D problems).
      read (1,'(a8,1x,d23.2)') label, tout
# ifdef DEBUG
      write(utrace,*) label, txtxt, tout
# endif
      label1 = 'tout'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!A transversal plane to the axis is written every FREQPL processors. 
      read (1,'(a8,1x,i10)') label, freqpl
# ifdef DEBUG
      write(utrace,*) label, txtxt, freqpl
# endif
      label1 = 'freqpl'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!DT is written on the screen every ITSTP timesteps.
      read (1,'(a8,1x,i10)') label, itstp
# ifdef DEBUG
      write(utrace,*) label, txtxt, itstp
# endif
      label1 = 'itstp'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!CFL is the time-limiting factor (< 1).
      read (1,'(a8,1x,d23.2)') label, cfl
# ifdef DEBUG
      write(utrace,*) label, txtxt, cfl
# endif
      label1 = 'cfl'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!If IAPPRS is >= 1 then an approximate Riemann solver is used
      read (1,'(a8,1x,i10)') label, iapprs
# ifdef DEBUG
      write(utrace,*) label, txtxt, iapprs
# endif
      label1 = 'iapprs'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!IRK is the order of the Runge-Kutta used for time integration.
!If it is <= 0 subroutine STAT1D is used.
      read (1,'(a8,1x,i10)') label, irk
# ifdef DEBUG
      write(utrace,*) label, txtxt, irk
# endif
      label1 = 'irk'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!IGODU sets the spatial reconstruction technique (if it is == 1, then
!first order method; else, PPM)
      read (1,'(a8,1x,i10)') label, igodu
# ifdef DEBUG
      write(utrace,*) label, txtxt, igodu
# endif
      label1 = 'igodu'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!DTINI is the initial DT.
      read (1,'(a8,1x,d23.2)') label, dtini
# ifdef DEBUG
      write(utrace,*) label, txtxt, dtini
# endif
      label1 = 'dtini'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!NITER is the number of iterations to get the pressure.
      read (1,'(a8,1x,i10)') label, niter
# ifdef DEBUG
      write(utrace,*) label, txtxt, niter
# endif
      label1 = 'niter'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!SMALL is the threshold for nondimensional numbers (i.e. velocity).
      read (1,'(a8,1x,d23.2)') label, small
# ifdef DEBUG
      write(utrace,*) label, txtxt, small
# endif
      label1 = 'small'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!SMLRHO is the threshold for densities (rho, r).
      read (1,'(a8,1x,d23.2)') label, smlrho
# ifdef DEBUG
      write(utrace,*) label, txtxt ,smlrho
# endif
      label1 = 'smlrho'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!SMALLP is the threshold for pressure.
      read (1,'(a8,1x,d23.2)') label, smallp
# ifdef DEBUG
      write(utrace,*) label, txtxt, smallp
# endif
      label1 = 'smallp'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!SMALLU is the threshold for internal energy.
      read (1,'(a8,1x,d23.2)') label, smallu
# ifdef DEBUG
      write(utrace,*) label, txtxt, smallu
# endif
      label1 = 'smallu'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!NSDIM is the number of spatial dimensions (1, 2 or 3).
      read (1,'(a8,1x,i10)') label, nsdim
# ifdef DEBUG
      write(utrace,*) label, txtxt, nsdim
# endif
      label1 = 'nsdim'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

      if ((basenm=='RST_'.OR.basenm=='RBW_'.OR.   &
      basenm=='RSR_'.OR.basenm=='RBWI'.OR.   &
      basenm=='RCSR').AND. (nsdim/=1)) then
      nsdim = 1
# ifdef DEBUG
      write(utrace,*) '       NSDIM reset to 1' 
# endif
      endif

      if ((basenm=='RRST'.OR.basenm=='RRBW'.OR.   &
      basenm=='RRSR'.OR.basenm=='REFS'.OR.   &
      basenm=='RSCH'.OR.basenm=='PRST'.OR.   &
      basenm=='CRST') .AND. (nsdim/=2)) then
      nsdim = 2
# ifdef DEBUG
      write(utrace,*) '       NSDIM reset to 2'
# endif
      endif

      if ((nsdim==2.OR.nsdim==3).AND. (iapprs==0))then 
# ifdef DEBUG
         write(utrace,*) 'INPUT: The exact Riemann solver cannot be used in'
         write(utrace,*) '       two-dimensional calculations'
# endif
         ierr = ierr + 1
      endif

      if ((nsdim==2..OR.nsdim==3).AND. (irk.le.0)) then 
# ifdef DEBUG
         write(utrace,*) 'INPUT: Subroutine STAT1D can only be used in'
         write(utrace,*) '       one-dimensional calculations'
# endif
         ierr = ierr + 1
      endif

      if (nsdim.gt.3) then 
# ifdef DEBUG
         write(utrace,*) 'INPUT: Unsuitable number of spatial dimensions'
# endif
         ierr = ierr + 1
      endif

!IXSYM indicates whether the model is symmetric in the x-direction
!(only used in jet simulations and with planar coordinates)
      read (1,'(a8,1x,i10)') label, ixsym
# ifdef DEBUG
      write(utrace,*) label, txtxt, ixsym
# endif
      label1 = 'ixsym'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif
      if (basenm(1:1)/='J' .AND. basenm(1:1)/='H' .AND. &
      basenm(1:1)/='3' .AND. basenm(1:1)/='M' .AND. ixsym/=1) then 
      ixsym = 1
# ifdef DEBUG
      write(utrace,*) '       IXSYM reset to 1'
# endif
      endif

!IZSYM indicates whether the model is symmetric in the z-direction
!(only used in jet simulations and with planar coordinates)
      read (1,'(a8,1x,i10)') label, izsym
# ifdef DEBUG
      write(utrace,*) label, txtxt, izsym
# endif
      label1 = 'izsym'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!IEXTRAX == 1 forces the generation of an extended grid in the 
!x-direction with geometrically expanded zones (only used in jet
!simulations)
      read (1,'(a8,1x,i10)') label, iextrax
# ifdef DEBUG
      write(utrace,*) label, txtxt, iextrax
# endif
      label1 = 'iextrax'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif
      if (basenm(1:1)/='J' .AND. basenm(1:1)/='H' .AND. &
      basenm(1:1)/='3' .AND. basenm(1:1)/='M' .AND. iextrax/=0) then
      iextrax = 0  
# ifdef DEBUG
      write(utrace,*) '       IEXTRAX reset to 0'
# endif
      endif

!IEXTRAZ == 1 forces the generation of an extended grid in the 
!z-direction with geometrically expanded zones (only used in jet
!simulations)
      read (1,'(a8,1x,i10)') label, iextraz
# ifdef DEBUG
      write(utrace,*) label, txtxt, iextraz
# endif
      label1 = 'iextraz'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif
      if (basenm(1:1)/='J'.AND.basenm(1:1)/='H'.AND.  &
      basenm(1:1)/='3'.AND. basenm(1:1)/='M'.AND.iextraz/=0) then
      iextraz = 0  
# ifdef DEBUG
      write(utrace,*) '       IEXTRAZ reset to 0'
# endif
      endif


!GRIDLX is the length of the grid in the x-direction (in the case
!IEXTRAX == 1, is the length of the uniform innest x-grid; in the 
!case of IXSYM == 0 is the length of one half of the uniform innest x-grid)
      read (1,'(a8,1x,d23.2)') label, gridlx
# ifdef DEBUG
      write(utrace,*) label, txtxt, gridlx
# endif 
      label1 = 'gridlx'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif
      if ((basenm=='RST_'.OR.basenm=='RBW_'.OR.    &
      basenm=='RSR_'.OR.basenm=='RBWI'.OR.    &
      basenm=='RCSR'.OR.basenm=='RRST'.OR.    &
      basenm=='RRBW'.OR.basenm=='RRSR'.OR.    &
      basenm=='REFS'.OR.basenm=='CRST'.OR.    &
      basenm=='PRST') .AND. (gridlx/=1.0)) then
      gridlx = 1.0
# ifdef DEBUG
      write(utrace,*) '       GRIDLX reset to 1.0'
# endif
      endif

      if (basenm=='RSCH' .AND. gridlx/=aux) then
         gridlx = aux
# ifdef DEBUG
         write(utrace,*) '       GRIDLX reset to 1.0/sqrt(2.0)'
# endif
      endif

!GRIDLX2 is the length of the extended grid in the x-direction
!(only used when IEXTRAX == 1)
      read (1,'(a8,1x,d23.2)') label, gridlx2
# ifdef DEBUG
      write(utrace,*) label, txtxt, gridlx2
# endif
      label1 = 'gridlx2'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!     -GHANGES BY 3D PARALLELIZATION: Gridlx ir recalculated if IEXTRAX=1
      if (iextrax==1) then
         gridlx = gridlx + gridlx2
      endif
!     gridlx ir recalculated if IXSYM/=1
!      if (ixsym/=1) then
!         gridlx = 2.*gridlx 
!      endif

!     -END CHANGES

!GRIDLZ is the length of the grid in the z-direction (in the case
!IEXTRAZ == 1, is the length of the uniform innest z-grid; in the 
!case of IZSYM == 0 is the length of one half of the uniform innest z-grid)
      read (1,'(a8,1x,d23.2)') label, gridlz
# ifdef DEBUG
      write(utrace,*) label, txtxt, gridlz
# endif
      label1 = 'gridlz'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif


!GRIDLZ2 is the length of the extended grid in the z-direction
!(only used when IEXTRAZ == 1)
      read (1,'(a8,1x,d23.2)') label, gridlz2
# ifdef DEBUG
      write(utrace,*) label, txtxt, gridlz2
# endif
      label1 = 'gridlz2'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif


!     -GHANGES BY 3D PARALLELIZATION: Gridlz ir recalculated if IEXTRAZ=1
      if (iextraz==1) then
         gridlz = gridlz + gridlz2
      endif
!     gridlz ir recalculated if IZSYM/=1
!      if (izsym/=1) then
!         gridlz = 2.*gridlz 
!      endif
!     -END CHANGES

!GRIDLY is the length of the grid in the y-direction.
      read (1,'(a8,1x,d23.2)') label, gridly
# ifdef DEBUG
      write(utrace,*) label, txtxt, gridly
# endif
      label1 = 'gridly'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

      if ((basenm=='RRST'.OR.basenm=='RRBW'.OR.     &
      basenm=='PRST'.OR.basenm=='CRST'.OR.     &
      basenm=='RRSR'.OR.basenm=='RSCH').AND.   &
      gridly/=1.0) then
      gridly = 1.0
# ifdef DEBUG
      write(utrace,*) '       GRIDLY reset to 1.0'
# endif
      endif

      if (basenm=='RSCH' .AND. gridly/=aux) then
         gridly = aux
# ifdef DEBUG
         write(utrace,*) '       GRIDLY reset to 1.0/sqrt(2.0)'
# endif
      endif

      if (basenm=='REFS' .AND. gridly/=3.0) then
         gridly = 3.0
# ifdef DEBUG
         write(utrace,*) '       GRIDLY reset to 3.0'
# endif
      endif

!NX is the number of grid points in the x-direction (in the case
!IEXTRAX == 1, is the number of zones of the uniform innest x-grid; in 
!the case of IXSYM == 0 is the number of zones of one half of the 
!uniform innest x-grid)
      read (1,'(a8,1x,i10)') label, nx
!print*, 'INPUT: NX',nx
!call flush(6)
# ifdef DEBUG
      write(utrace,*) label, txtxt, nx
# endif
      label1 = 'nx'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

      if (nx<4) then
# ifdef DEBUG
         write(utrace,*) 'INPUT: Unsuitable nx'
# endif
         ierr = ierr + 1
      endif

!NX2 is the number of grid points of the extended grid in the x-direction
!(only used when IEXTRAX == 1)
      read (1,'(a8,1x,i10)') label, nx2
# ifdef DEBUG
      write(utrace,*) label, txtxt, nx2
# endif
      label1 = 'nx2'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif


!     -GHANGES BY 3D PARALLELIZATION
!     moved mnx_ to calculate_mn()
!     NX ir recalculated if IEXTRAX=1
      if (iextrax==1) then
         nx = nx + nx2
      endif
!     NX ir recalculated if IXSYM/=1
      ! Eliminado pq en inpt.dat ponemos los valores globales de nx y gridlx
!     if (ixsym/=1) then
         !!!nx = 2*nx   
      !endif
!     -END CHANGES

      
!NXO is the number of grid points of the grid in the x-direction in previous runs
      read (1,'(a8,1x,i10)') label, nxo
# ifdef DEBUG
      write(utrace,*) label, txtxt, nxo
# endif
      label1 = 'nxo'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif


!NY is the number of grid points in the y-direction.
      read (1,'(a8,1x,i10)') label, ny
# ifdef DEBUG
      write(utrace,*) label, txtxt, ny
# endif
      label1 = 'ny'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

      if (nsdim>1.AND.ny<4) then
# ifdef DEBUG
         write(utrace,*) 'INPUT: Unsuitable ny'
# endif
         ierr = ierr + 1
      endif

      if (nsdim==1) then
         ny = 1
      endif


!NYOLD is the OLD number of grid points in the y-direction.
      read (1,'(a8,1x,i10)') label, nyold
# ifdef DEBUG
      write(utrace,*) label, txtxt, nyold
# endif
      label1 = 'nyold'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif


!NZ is the number of grid points in the z-direction (in the case
!IEXTRAZ == 1, is the number of zones of the uniform innest z-grid; in 
!the case of IZSYM == 0 is the number of zones of one half of the 
!uniform innest z-grid)
      read (1,'(a8,1x,i10)') label, nz
# ifdef DEBUG
      write(utrace,*) label, txtxt, nz
# endif
      label1 = 'nz'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

      if (nsdim>2 .AND. nz<4) then
# ifdef DEBUG
         write(utrace,*) 'INPUT: Unsuitable nz'
# endif
         ierr = ierr + 1
      endif

!NZ2 is the number of grid points of the extended grid in the z-direction
!(only used when IEXTRAZ == 1)
      read (1,'(a8,1x,i10)') label, nz2
# ifdef DEBUG
      write(utrace,*) label, txtxt, nz2
# endif
      label1 = 'nz2'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

      

!     -GHANGES BY 3D PARALLELIZATION
!     moved mnx_ to calculate_mn()
!     NZ ir recalculated if IEXTRAZ=1
      if (iextraz==1) then
         nz = nz + nz2
      endif
!     NZ ir recalculated if IZSYM/=1
      ! Eliminado pq en inpt.dat ponemos los valores globales de nz y gridlz
      !if (izsym/=1) then
      !   nz = 2*nz 
      !endif
!     -END CHANGES

!NZO is the number of grid points of the grid in the x-direction in previous runs
      read (1,'(a8,1x,i10)') label, nzo
# ifdef DEBUG
      write(utrace,*) label, txtxt, nzo
# endif
      label1 = 'nzo'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!IGEOMX defines the type of the x coordinate. 
      read (1,'(a8,1x,i10)') label, igeomx
# ifdef DEBUG
      write(utrace,*) label, txtxt, igeomx
# endif
      label1 = 'igeomx'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

      if ((basenm=='RST_'.OR.basenm=='RBW_'.OR.   &
      basenm=='RSR_'.OR.basenm=='RBWI'.OR.   &
      basenm=='RRST'.OR.basenm=='RRBW'.OR.   &
      basenm=='SCHN'.OR.basenm=='REFS'.OR.   &
      basenm=='PRST'.OR.basenm=='CRST'.OR.   &
      basenm=='RSCH').AND.igeomx/=0) then
      igeomx = 0
# ifdef DEBUG
      write(utrace,*) '       IGEOMX reset to 0'
# endif
      endif

      if (basenm=='RCSR'.AND.igeomx/=1) then
         igeomx = 1
# ifdef DEBUG
         write(utrace,*) '       IGEOMX reset to 1'
# endif
      endif

      if (igeomx>2) then
# ifdef DEBUG
         write(utrace,*) 'INPUT: unsuitable igeomx'
# endif
         ierr = ierr + 1
      endif
      if (igeomx/=0.AND.ixsym==0) then 
# ifdef DEBUG
         write(utrace,*) 'INPUT: non-symmetric initial data in non'
         write(utrace,*) '       planar x-coordinate'
# endif
         ierr = ierr + 1
      endif
      if (igeomx/=0.AND.nsdim==3) then
# ifdef DEBUG
         write(utrace,*) 'INPUT: 3D simulation in non-planar coordinates'
# endif
         ierr = ierr + 1
      endif

!IGEOMY defines the type of the y coordinate.
      read (1,'(a8,1x,i10)') label, igeomy
# ifdef DEBUG
      write(utrace,*) label, txtxt, igeomy
# endif
      label1 = 'igeomy'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

      if ((basenm=='RRST'.OR.basenm=='RRBW'.OR.    &
      basenm=='RRSR'.OR.basenm=='REFS'.OR.    &
      basenm=='PRST'.OR.basenm=='CRST'.OR.    &
      basenm=='RSCH') .AND. igeomy/=0) then
      igeomy = 0
# ifdef DEBUG
      write(utrace,*) '       IGEOMY reset to 0'
# endif
      endif

      if (nsdim>1) then
         if ((igeomx==0.AND.igeomy/=0).OR.    &
         (igeomx==1.AND.igeomy/=0).OR.    &
         (igeomx==2.AND.igeomy/=4)) then
# ifdef DEBUG
         write(utrace,*) 'INPUT: unsuitable igeomy'
# endif
         ierr = ierr + 1
      endif
      endif

      if (igeomy/=0.AND.nsdim==3) then
# ifdef DEBUG
         write(utrace,*) 'INPUT:3D simulation in non-planar coordinates'
# endif
         ierr = ierr + 1
      endif

!IGEOMZ defines the type of the z coordinate.
      read (1,'(a8,1x,i10)') label, igeomz
# ifdef DEBUG
      write(utrace,*) label, txtxt, igeomz
# endif
      label1 = 'igeomz'

      if (nsdim>2) then
         if ((igeomx==0.AND.igeomy==0.AND.igeomz/=0).OR.    &
         (igeomx==1.AND.igeomy==0.AND.igeomz/=1).OR.    &
         (igeomx==2.AND.igeomy==4.AND.igeomz/=2)) then
# ifdef DEBUG
         write(utrace,*) 'INPUT: unsuitable igeomz'
# endif
         ierr = ierr + 1
      endif
      endif

      if (igeomz/=0.AND.nsdim==3) then
# ifdef DEBUG
         write(utrace,*) 'INPUT:3D simulation in non-planar coordinates'
# endif
         ierr = ierr + 1
      endif

!EPSILN is used in subroutine FLATEN.
      read (1,'(a8,1x,d23.2)') label, epsiln
# ifdef DEBUG
      write(utrace,*) label, txtxt, epsiln
# endif
      label1 = 'epsiln'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!OMG1 is used in subroutine FLATEN.
      read (1,'(a8,1x,d23.2)') label, omg1
# ifdef DEBUG
      write(utrace,*) label, txtxt, omg1
# endif
      label1 = 'omg1'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

!OMG2 is used in subroutine FLATEN (A value of 0.0 for OMG2 
!supresses the function of FLATEN).
      read (1,'(a8,1x,d23.2)') label, omg2
# ifdef DEBUG
      write(utrace,*) label, txtxt, omg2
# endif
      label1 = 'omg2'
      if (label/=label1) then
         call inputerror( 'INPUT', label, label1, ierr )
      endif

      close(1)

      if (ierr /= 0) then
         call ERROR( 'input', 'There are errors in the input file', &
         __FILE__, __LINE__ )
      endif

      nyh = ny

! Set ny homogeneous and ny to extended area
!      IF (basenm(1:1)=='J'.OR.basenm(1:1)=='3') THEN
!         ny = ny + (ny/5)
!      ENDIF

      RETURN
!----------------------------------------------------------------------------END
      END subroutine input


!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!     This subroutine provides beam, atmosphere, perturbation, boundary data
!     for steady jet simulations
!     CREATED  20-DEC-1996
!     MODIFIED 07-MAY-1998 to chnage input format specifiers
!===============================================================================
      subroutine calcule_mn( ixsym, izsym, iextrax, iextraz, nx, nx2, &
      nz, nz2, mnx1, mnx5, mnx6, mnz1, mnz5, mnz6)

      USE parallel
      USE constantes
      implicit none
# include "types.h"
!------------------------------------------------------Input  Variables
      integer(is) ::  ixsym, izsym, iextrax, iextraz, nx, nx2, nz, nz2        

!------------------------------------------------------Output Variables
      integer(is) :: mnx1, mnx5, mnx6, mnz1, mnz5, mnz6


      
!if (ixsym/=1) then

!if (iextrax==1) then
!mnx1 = 2*(nx + nx2) + 1
!mnx5 = 2*(nx + nx2) + 5
!mnx6 = 2*(nx + nx2) + 6
      
!else
!mnx1 = 2*nx + 1
!mnx5 = 2*nx + 5
!mnx6 = 2*nx + 6
      
!endif
!elseif (iextrax==1) then
!mnx1 = nx + nx2 + 1
!mnx5 = nx + nx2 + 5
!mnx6 = nx + nx2 + 6
!else
      mnx1 = nx + 1
      mnx5 = nx + 5
      mnx6 = nx + 6
!endif
      
!if (izsym/=1) then
!if (iextraz==1) then
!mnz1 = 2*(nz + nz2) + 1
!mnz5 = 2*(nz + nz2) + 5
!mnz6 = 2*(nz + nz2) + 6
!else
!mnz1 = 2*nz + 1
!mnz5 = 2*nz + 5
!mnz6 = 2*nz + 6
!endif
!elseif (iextraz==1) then
!mnz1 = nz + nz2 + 1
!mnz5 = nz + nz2 + 5
!mnz6 = nz + nz2 + 6
!else
      mnz1 = nz + 1
      mnz5 = nz + 5
      mnz6 = nz + 6
!endif

!print*,'CALCULE MN: mnx1,mnx5,mnx6,mnz1,mnz5,mnz6',mnx1,mnx5,mnx6,mnz1,mnz5,mnz6
!!!call flush(6)
      RETURN
!----------------------------------------------------------------------------END
      END subroutine calcule_mn



!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL Variables
! DESCRIPTION
! SEE ALSO
!     This subroutine provides beam, atmosphere, perturbation, boundary data
!     for steady jet simulations
!     CREATED  20-DEC-1996
!     MODIFIED 07-MAY-1998 to chnage input format specifiers
!===============================================================================
subroutine input_j( basenm,                                               &
                  suffixd, tdlfil, itdelay, gamma_ad, velb, rhob, rhoeb,  &
                  temb, machb, radius, cdenst, shift, q0, fg, ms,         &
                  dk1, dk2, rm, psh1, psh2, nn, mm1, mm2, ipr,            &
                  alpn, betn, gamn, rbn, brt, tem0, tem1, oang, csb2,     &
                  pb, mclouds, t0, t1, vptbd, rptbd, pptbd )

		
  USE parallel
  USE workarrays
  implicit none
# include "types.h"
!------------------------------------------------------ Input  Variables
  character, intent(in) :: basenm*(*)

!------------------------------------------------------ Output Variables
  character(3):: suffixd
  character(7):: tdlfil

  integer(is):: itdelay, ipr, ms, brt, mclouds

  real(rs)   :: gamma_ad, velb, rhob, rhoeb, machb, radius,   &
                temb, fg, dk1, dk2, cdenst, shift, q0,        &
                psh1, psh2, nn, mm1, mm2, oang, tem0, tem1,   & 
                csb2, pb, alpn, betn, gamn, rbn,              &
                t0, t1, vptbd, rptbd, pptbd, rm

!------------------------------------------------------- Local Variables
  character(72) :: text
  character(5)  :: txtxt
  character(8)  :: label, label1
  character(256):: fich

  integer(is)   :: ierr

  real(rs)      :: aux
            
!------------------------------------------------------------------------- BEGIN
  ierr  = 0
  txtxt = '.....'

  fich = TRIM(path_wkd)//'/'//'inpt_j.dat'

  open( 1, file=fich, form='formatted', status='old', iostat=ierr )
  if (ierr /= 0) then
    call ERROR( 'input_j: Error openning ', fich, __FILE__, __LINE__ )
  endif

  read (1,'(a72)') text
  read (1,'(a72)') text
  read (1,'(a72)') text

  !BEAM DATA
  !Adiabatic index of the EOS
  read (1,'(a8,1x,d23.2)') label, gamma_ad
  label1 = 'gamma'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Beam flow velocity
  read (1,'(a8,1x,d23.2)') label, velb
  label1 = 'velb'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Beam proper rest-mass density
  read (1,'(a8,1x,d23.2)') label, rhob
  label1 = 'rhob'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Beam proper leptonic rest-mass density
  read (1,'(a8,1x,d23.2)') label, rhoeb
  label1 = 'rhoeb'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Beam Temperature
  read (1,'(a8,1x,d23.2)') label, temb
  label1 = 'temb'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif


  !(Classical) beam Mach number
  read (1,'(a8,1x,d23.2)') label, machb
  label1 = 'machb'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  aux = gamma_ad - 1.0
  if (machb<velb/SQRT(aux)) then
# ifdef DEBUG
    write(utrace,*) 'INPUT_J: Unsuitable beam Mach number'
# endif
    ierr= ierr + 1
  endif

  csb2  = velb/machb
  csb2  = csb2*csb2

  pb    = aux*rhob*csb2/gamma_ad/(aux - csb2)

  !Jet radius in parsecs
  read (1,'(a8,1x,d23.2)') label, radius
  label1 = 'radius'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Core density in mp/cm^3
  read (1,'(a8,1x,d23.2)') label, cdenst
  label1 = 'cdenst'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Distance of the injection in the grid to the nucleus (in jet radii)
  read (1,'(a8,1x,d23.2)') label, shift
  label1 = 'shift'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Slope of the shear layer
  read (1,'(a8,1x,i10)') label, ms
  label1 = 'ms'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Mass loss by stars in the core in g/yr/pc^3
  read (1,'(a8,1x,d23.2)') label, q0
  label1 = 'q0'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

 !fg is the ratio of the galaxy core density to the nucleus density (unit density)
  read (1,'(a8,1x,d23.2)') label, fg
  label1 = 'fg'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !rm is the radius for temperarture profile (in R_b)
  read (1,'(a8,1x,d23.2)') label, rm
  label1 = 'rm'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Opening angle
  read (1,'(a8,1x,d23.2)') label, oang
  label1 = 'oang'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif 

  !EXTERNAL ATMOSPHERE
  !Beam-to-ambient pressure ratio at injection position
  read (1,'(a8,1x,d23.2)') label, dk1
  label1 = 'dk1'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  read (1,'(a8,1x,d23.2)') label, dk2
  label1 = 'dk2'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !External pressure scale heigth
  read (1,'(a8,1x,d23.2)') label, psh1
  label1 = 'psh1'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  read (1,'(a8,1x,d23.2)') label, psh2
  label1 = 'psh2'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Exponents
  read (1,'(a8,1x,d23.2)') label, nn
  label1 = 'nn'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  read (1,'(a8,1x,d23.2)') label, mm1
  label1 = 'mm1'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  read (1,'(a8,1x,d23.2)') label, mm2
  label1 = 'mm2'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !External atmosphere defined by profile (ipr=1) or opening angle (ipr=2)
  read (1,'(a8,1x,i10)') label, ipr
  label1 = 'ipr'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Temperature in the core
  read (1,'(a8,1x,d23.2)') label, tem0
  label1 = 'tem0'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Temperature in the outer component
  read (1,'(a8,1x,d23.2)') label, tem1
  label1 = 'tem1'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Clouds
  read (1,'(a8,1x,i10)') label, mclouds
  label1 = 'mclouds'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif


  !NUKER PROFILE
  !Exponents
  read (1,'(a8,1x,d23.2)') label, alpn
  label1 = 'alpn'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  read (1,'(a8,1x,d23.2)') label, betn
  label1 = 'betn'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  read (1,'(a8,1x,d23.2)') label, gamn
  label1 = 'gamn'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Core radius
  read (1,'(a8,1x,d23.2)') label, rbn
  label1 = 'rbn'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !BREMSSTRAHLUNG RADIATION
  read (1,'(a8,1x,i10)') label, brt
  label1 = 'brt'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !PERTURBATION (SQUARED)
  !Switch-on time
  read (1,'(a8,1x,d23.2)') label, t0
  label1 = 't0'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Switch-off time
  read (1,'(a8,1x,d23.2)') label, t1
  label1 = 't1'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Flow velocity
  read (1,'(a8,1x,d23.2)') label, vptbd
  label1 = 'vptbd'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Proper rest-mass density
  read (1,'(a8,1x,d23.2)') label, rptbd
  label1 = 'rptbd'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Pressure
  read (1,'(a8,1x,d23.2)') label, pptbd
  label1 = 'pptbd'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !TIME DELAY STUDIES
  !Storing data (y/n)
  read (1,'(a8,1x,i10)') label, itdelay
  label1 = 'itdelay'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  !Initial suffix for time delay files
  read (1,'(a8,1x,a3)') label, suffixd
  label1 = 'suffixd'
  if (label/=label1) then
    call inputerror( 'INPU_J', label, label1, ierr )
  endif

  tdlfil  = basenm//suffixd

  if (ierr /= 0) then
    call ERROR( 'input_j', 'There are errors in the input file', &
                __FILE__, __LINE__ )
  endif
!----------------------------------------------------------------------------END
END subroutine input_j


!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!===============================================================================
      subroutine inputerror( rutina, label, label1, ierr )
      implicit none
# include "types.h"
!-------------------------------------------------------Input Variables
      character(6)::  rutina
      character(8) ::  label, label1

!------------------------------------------------------Output Variables
      integer(is)  ::  ierr

!-------------------------------------------------------------------------BEGIN
      write(6,*) ' '
      write(6,*) rutina, ': Incorrect input deck'
      write(6,*) '       label = ',label, '  expected label = ', label1
      ierr = ierr+1
!----------------------------------------------------------------------------END
      END subroutine inputerror    
