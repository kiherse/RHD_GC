#include "types.h"

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================s
MODULE constantes
  real(rs), parameter :: pi = 3.14159265358979323844
END MODULE constantes

MODULE tipos
  type TboundPosX
    real(rs):: xzn
    real(rs):: xznl
    real(rs):: xznr
    real(rs):: dx
  end type TboundPosX

  type TboundPosY
    real(rs):: yzn
    real(rs):: yznl
    real(rs):: yznr
    real(rs):: dx
  end type TboundPosY

  type TboundPosZ
    real(rs):: zzn
    real(rs):: zznl
    real(rs):: zznr
    real(rs):: dx
  end type TboundPosZ

  type TboundPhysicMini
    real(rs):: pres
    real(rs):: densty
    real(rs):: denstye
  end type TboundPhysicMini

  type TboundPhysicShort
    real(rs):: velx
    real(rs):: vely
    real(rs):: velz
    real(rs):: densty
    real(rs):: denstye
    real(rs):: eps
  end type TboundPhysicShort

  type TboundPhysic
    real(rs):: velx
    real(rs):: vely
    real(rs):: velz
    real(rs):: densty
    real(rs):: denstye   
    real(rs):: eps
    real(rs):: pres
    real(rs):: tracer
  end type TboundPhysic

  type TplanePhysic
    real(rs):: velx
    real(rs):: vely
    real(rs):: velz
    real(rs):: densty
    real(rs):: denstye   
    real(rs):: eps
    real(rs):: tracer
  end type TplanePhysic

  type Tbounds
    real(rs):: rho
    real(rs):: pres
    real(rs):: tem
  end type

  type Tphysic
    real(rs):: velx
    real(rs):: vely
    real(rs):: velz
    real(rs):: momenx
    real(rs):: momeny
    real(rs):: momenz
    real(rs):: densty
    real(rs):: denstye   
    real(rs):: eps
    real(rs):: pres
    real(rs):: denstw
    real(rs):: denstew   
    real(rs):: energy
    real(rs):: tracer
    real(rs):: omomenx
    real(rs):: omomeny
    real(rs):: omomenz
    real(rs):: opres
    real(rs):: odenstw
    real(rs):: odenstew   
    real(rs):: oenergy
    real(rs):: otracer
    real(rs):: stems
    real(rs):: gravx, gravy, gravz
  end type Tphysic

  type TstepInfoReal
    real(rs):: dt
    real(rs):: sound
    real(rs):: v
  end type TstepInfoReal

  type TstepInfoInt
    integer(is):: ic
    integer(is):: jc
    integer(is):: kc
  end type TstepInfoInt

  type Twaspace
    ! -4:mn6 vectors
    real(rs), pointer :: wa100(:), wa101(:), wa102(:), wa103(:), wa104(:), &
                         wa105(:), wa106(:), wa107(:), wa108(:), wa109(:), &
                         wa110(:), wa111(:), wa112(:), wa113(:), wa114(:), &
                         wa115(:), wa116(:), wa117(:), wa118(:), wa119(:), &
                         wa120(:), wa121(:), wa122(:), wa123(:), wa124(:), &
                         wa125(:), wa126(:), wa127(:), wa128(:), wa129(:), &
                         wa130(:), wa131(:), wa132(:), wa133(:), wa134(:), &
                         wa135(:), wa136(:), wa137(:), wa138(:), wa139(:), &
                         wa140(:), wa141(:), wa142(:), wa143(:), wa144(:), &
                         wa145(:), wa146(:), wa147(:), wa148(:), wa149(:), &
                         wa150(:), wa151(:), wa152(:), wa153(:), wa154(:), &
                         wa155(:), wa156(:), wa157(:), wa158(:), wa159(:), &
                         wa160(:), wa161(:), wa162(:), wa163(:), wa164(:), &
                         wa165(:), wa166(:), wa167(:), wa168(:), wa169(:), &
                         wa170(:), wa171(:), wa172(:), wa173(:), wa174(:), &
                         wa175(:), wa176(:), wa177(:), wa178(:), wa179(:), &
                         wa180(:), wa181(:), wa182(:), wa183(:), wa184(:), &
                         wa185(:), wa186(:), wa187(:), wa188(:), wa189(:), &
                         wa190(:), wa191(:), wa192(:), wa193(:), wa194(:), &
                         wa195(:), wa196(:), wa197(:), wa198(:), wa199(:), &
                         wa1101(:), wa1102(:), wa1103(:), wa1104(:),       &
                         wa1105(:), wa1106(:) 

    ! -4:mn5 vectors
    real(rs), pointer :: wa200(:), wa201(:), wa202(:), wa203(:), wa204(:), &
                         wa205(:), wa206(:), wa207(:), wa208(:), wa209(:), &
                         wa210(:), wa211(:), wa212(:), wa213(:), wa214(:), &
                         wa215(:), wa216(:), wa217(:), wa218(:), wa219(:), &
                         wa220(:), wa221(:), wa222(:), wa223(:), wa224(:), &
                         wa225(:), wa226(:), wa227(:), wa228(:), wa229(:), &
                         wa230(:), wa231(:), wa232(:), wa233(:), wa234(:), &
                         wa235(:), wa236(:), wa237(:), wa238(:), wa239(:), &
                         wa240(:), wa241(:), wa242(:), wa243(:), wa244(:), &
                         wa245(:), wa246(:), wa247(:), wa248(:), wa249(:), &
                         wa250(:), wa251(:), wa252(:), wa253(:), wa254(:), &
                         wa255(:), wa256(:), wa257(:), wa258(:), wa259(:), &
                         wa260(:), wa261(:), wa262(:), wa263(:), wa264(:), &
                         wa265(:), wa266(:), wa267(:), wa268(:), wa269(:), &
                         wa270(:), wa271(:), wa272(:), wa273(:), wa274(:), &
                         wa275(:), wa276(:), wa277(:), wa278(:), wa279(:), &
                         wa280(:), wa281(:), wa282(:), wa283(:), wa284(:), &
                         wa285(:), wa286(:), wa287(:), wa288(:), wa289(:), &
                         wa290(:), wa291(:) 

    ! mn1 vectors
    real(rs), pointer :: wa300(:), wa301(:), wa302(:)

  end type Twaspace
END MODULE tipos

MODULE memoria
  integer(is), parameter :: maxVars = 200

  integer(is):: nbvars
  data          nbvars  /0/

  real(rs) :: maxMem, usedMem
  data        maxMem, usedMem  /0.0,0.0/

  type Tvar
    character(32):: vname
    integer(is)  :: esize
    integer(is)  :: bsize
  end type Tvar

  type(Tvar) :: vtable(maxVars)
 
  CONTAINS
  subroutine memPush( nele, elesize, var )
    implicit none
!------------------------------------------ Input Variables
    integer(is) :: nele, elesize
    character(*):: var

!------------------------------------------------------------------------- BEGIN
    usedMem = usedMem + REAL(nele) * REAL(elesize)
    maxMem = MAX( maxMem, usedMem )
    if (nbvars < maxVars) then
      nbvars = nbvars + 1
      vtable(nbvars)%vname = var
      vtable(nbvars)%esize = nele
      vtable(nbvars)%bsize = elesize
    else
      write(*,*) "Table of variables is full. Increase maxVars param"
    endif
!--------------------------------------------------------------------------- END
  end subroutine memPush

  subroutine memPop( var )
    implicit none
!------------------------------------------ Input Variables
    character(*):: var
!------------------------------------------ Local Variables
    integer(is) :: ii

!------------------------------------------------------------------------- BEGIN
    if (nbvars > 0) then
      ii= 1
      do while ((LGT(vtable(ii)%vname,var).OR.LLT(vtable(ii)%vname,var)) &
                .AND. ii<=nbvars)
        ii= ii+1
      enddo
      if (ii <= nbvars) then
        usedMem = usedMem - REAL(vtable(ii)%esize) * REAL(vtable(ii)%bsize)

        if (nbvars > 1) then
          vtable(ii) = vtable(nbvars)
        endif
        nbvars = nbvars - 1
      else
        write(*,*) 'Variable not found in table: ', var
      endif
    endif
!--------------------------------------------------------------------------- END
  end subroutine memPop

  subroutine showMem(  )
    implicit none
!------------------------------------------ Local Variables
    integer(is)   ::ii
    integer(is*2) ::kk

!------------------------------------------------------------------------- BEGIN
    if (nbvars > 0) then
      write(*,*) '------------------- Memory report ------------------------'
      write(*,*) ' There are allocated variables: ', nbvars
      do ii= 1, nbvars
        kk = (vtable(ii)%esize*vtable(ii)%bsize)
        write(*,*) ii, vtable(ii)%esize, vtable(ii)%bsize, kk, &
                   TRIM(vtable(ii)%vname)
      enddo
      write(*,*) '----------------------------------------------------------'
    endif
!--------------------------------------------------------------------------- END
  end subroutine showMem
END MODULE memoria


!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
MODULE parallel
  integer(is)              :: nuproc, nbproc, sizeTotal(3), utrace, traza
  integer(is), allocatable :: iniTramo(:,:), finTramo(:,:), sizeTramo(:,:)
END MODULE parallel

!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
MODULE workarrays
  USE tipos
  USE parallel

  ! command line arguments
  character(len=256) :: outRoot
  character(len=256) :: path_wkd

# ifdef OPENMP
  integer, parameter :: MAXTHREADS = 32

# else
  integer, parameter :: MAXTHREADS = 1
# endif
  ! number of OpenMP threads
  integer :: omp_th, curIter
  ! coords of OpenMP threads
  integer, dimension(:), allocatable :: ocurx, ocury, ocurz
  ! number of OpenMP threads per dimension
  integer :: onumx, onumy, onumz

  ! MPI dims, coords and periods
  integer :: mpi_dims(3), mpi_coords(3), cartrank, CART_WORLD
  logical :: periods(3)
  ! neighbor ranks
  integer :: xleft, xright, yleft, yright, zleft, zright

  type(Twaspace), save :: wa(MAXTHREADS)

  integer(is)          :: nele, maxele
!  integer(is), pointer :: walis(:,:)

  ! mnx1,mny1 vectors
  real(rs), pointer :: wa400(:,:), wa401(:,:), wa402(:,:)

  !mnx1,mny1,mnz1 array for tracer evolution
  real(rs), pointer :: wa409(:,:,:)

  !1024,1024,1024 array for data collection
!  real(rs), pointer :: wa4009(:,:,:)    

  ! Local:  {nx,nz}, {2,ny},  {6,nx,nz} vectors
  ! Global:          {2,nyh}, {6,nx,nz} vectors
  real(rs), pointer :: wa410(:,:), wa411(:,:), wa412(:,:,:), &
                       wa420(:,:), wa421(:,:,:)

  ! TboundPhysic nx,5,nz vector
  type(TboundPhysic), pointer :: wa500(:,:,:)

  ! TboundPhysicShort nx,nz vector
  type(TboundPhysicShort), pointer :: wa600(:,:)

  ! TstepInfoReal & TstepInfoInt vectors (size 0:cartrank)
  type(TstepInfoReal), pointer :: wa700(:)
  type(TstepInfoInt),  pointer :: wa701(:)

  ! TboundPhysic arrays used in intercambiaBoundPhysic
  type(TboundPhysic), pointer :: boundsXs1(:,:,:), boundsYs1(:,:,:),boundsZs1(:,:,:) 
  type(TboundPhysic), pointer :: boundsXr1(:,:,:), boundsYr1(:,:,:),boundsZr1(:,:,:)
  type(TboundPhysic), pointer :: boundsXs2(:,:,:), boundsYs2(:,:,:),boundsZs2(:,:,:) 
  type(TboundPhysic), pointer :: boundsXr2(:,:,:), boundsYr2(:,:,:),boundsZr2(:,:,:)

  CONTAINS

 subroutine getBoundsMemorySpace( nx, ny, nz )
    USE memoria
    implicit none
!------------------------------------------ Input Variables
    integer(is):: nx, ny, nz
!------------------------------------------ Local Variables
    integer(is) :: ierr
 
 allocate( boundsXr1(5,ny,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif

 !call memPush( 5*ny*nz, rs, 'boundsXr1' )
 
 allocate( boundsYr1(nx,5,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYr1' )
 
 allocate( boundsZr1(nx,ny,5), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZr1' )

  allocate( boundsXr2(5,ny,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( 5*ny*nz, rs, 'boundsXr2' )
 
 allocate( boundsYr2(nx,5,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYr2' )
 
 allocate( boundsZr2(nx,ny,5), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZr2' )
 

 allocate( boundsXs1(5,ny,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( 5*ny*nz, rs, 'boundsXs1' )
 
 allocate( boundsYs1(nx,5,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYs1' )
 
 allocate( boundsZs1(nx,ny,5), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZs1' )

  allocate( boundsXs2(5,ny,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( 5*ny*nz, rs, 'boundsXs2' )
 
 allocate( boundsYs2(nx,5,nz), stat= ierr )
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
  endif
  !call memPush( nx*5*nz, rs, 'boundsYs2' )
 
 allocate( boundsZs2(nx,ny,5), stat= ierr ) 
 if (ierr /= 0) then
      call ERROR( 'intercambiaBoundPhysic', 'Allocation of bounds', &
                  __FILE__, __LINE__  ) 
 endif
 !call memPush( nx*ny*5, rs, 'boundsZs2' )

!--------------------------------------------------------------------------- END
  end subroutine getBoundsMemorySpace

subroutine freeBoundsMemorySpace()
    USE memoria
    implicit none

!------------------------------------------------------------------------- BEGIN

    deallocate( boundsXs1 )
    !call memPop( 'boundsXs1' )
    deallocate( boundsYs1 )
    !call memPop( 'boundsYs1' )
    deallocate( boundsZs1 )
    !call memPop( 'boundsZs1' )
    deallocate( boundsXs2 )
    !call memPop( 'boundsXs2' )
    deallocate( boundsYs2 )
    !call memPop( 'boundsYs2' )
    deallocate( boundsZs2 )
    !call memPop( 'boundsZs2' )  
    deallocate( boundsXr1 )
    !call memPop( 'boundsXr1' )
    deallocate( boundsYr1 )
    !call memPop( 'boundsYr1' )
    deallocate( boundsZr1 )
    !call memPop( 'boundsZr1' )
    deallocate( boundsXr2 )
    !call memPop( 'boundsXr2' )
    deallocate( boundsYr2 )
    !call memPop( 'boundsYr2' )
    deallocate( boundsZr2 )
    !call memPop( 'boundsZr2' )

!--------------------------------------------------------------------------- END
  end subroutine freeBoundsMemorySpace

  subroutine getWAMemorySpace( mn6, mn5, mn1, mnx1, mny1, mnz1, nx, ny, nyh, nz, mpi_dims1, mpi_dims2, mpi_dims3)
    implicit none
!------------------------------------------ Input Variables
    integer(is):: mn6, mn5, mn1, mnx1, mny1, mnz1, nx, ny, nyh, nz, mpi_dims1, mpi_dims2, mpi_dims3

!------------------------------------------ Local Variables
    integer(is):: tid, nth, ierr

#   ifdef OPENMP
      integer(is):: omp_get_max_threads
#   endif

!------------------------------------------------------------------------- BEGIN
#   ifdef OPENMP

      nth = omp_get_max_threads()
      if (nth.GT.MAXTHREADS) then
        call ERROR( 'getMemorySpace', 'nth > MAXTHREADS', __FILE__, __LINE__  ) 
      endif
#   else
      nth = 1
#   endif

   !print*, 'getWAMemorySpace, nth',nth


    DO tid= 1, nth
      allocate( wa(tid)%wa100(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa101(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa102(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa103(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa104(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa105(-4:mn6), wa(tid)%wa106(-4:mn6), wa(tid)%wa107(-4:mn6), wa(tid)%wa108(-4:mn6), wa(tid)%wa109(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa110(-4:mn6), wa(tid)%wa111(-4:mn6), wa(tid)%wa112(-4:mn6), wa(tid)%wa113(-4:mn6), wa(tid)%wa114(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa115(-4:mn6), wa(tid)%wa116(-4:mn6), wa(tid)%wa117(-4:mn6), wa(tid)%wa118(-4:mn6), wa(tid)%wa119(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa120(-4:mn6), wa(tid)%wa121(-4:mn6), wa(tid)%wa122(-4:mn6), wa(tid)%wa123(-4:mn6), wa(tid)%wa124(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa125(-4:mn6), wa(tid)%wa126(-4:mn6), wa(tid)%wa127(-4:mn6), wa(tid)%wa128(-4:mn6), wa(tid)%wa129(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa130(-4:mn6), wa(tid)%wa131(-4:mn6), wa(tid)%wa132(-4:mn6), wa(tid)%wa133(-4:mn6), wa(tid)%wa134(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa135(-4:mn6), wa(tid)%wa136(-4:mn6), wa(tid)%wa137(-4:mn6), wa(tid)%wa138(-4:mn6), wa(tid)%wa139(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa140(-4:mn6), wa(tid)%wa141(-4:mn6), wa(tid)%wa142(-4:mn6), wa(tid)%wa143(-4:mn6), wa(tid)%wa144(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa145(-4:mn6), wa(tid)%wa146(-4:mn6), wa(tid)%wa147(-4:mn6), wa(tid)%wa148(-4:mn6), wa(tid)%wa149(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa150(-4:mn6), wa(tid)%wa151(-4:mn6), wa(tid)%wa152(-4:mn6), wa(tid)%wa153(-4:mn6), wa(tid)%wa154(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa155(-4:mn6), wa(tid)%wa156(-4:mn6), wa(tid)%wa157(-4:mn6), wa(tid)%wa158(-4:mn6), wa(tid)%wa159(-4:mn6),  stat= ierr)  
      allocate( wa(tid)%wa160(-4:mn6), wa(tid)%wa161(-4:mn6), wa(tid)%wa162(-4:mn6), wa(tid)%wa163(-4:mn6), wa(tid)%wa164(-4:mn6),  stat= ierr)   
      allocate( wa(tid)%wa165(-4:mn6), wa(tid)%wa166(-4:mn6), wa(tid)%wa167(-4:mn6), wa(tid)%wa168(-4:mn6), wa(tid)%wa169(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa170(-4:mn6), wa(tid)%wa171(-4:mn6), wa(tid)%wa172(-4:mn6), wa(tid)%wa173(-4:mn6), wa(tid)%wa174(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa175(-4:mn6), wa(tid)%wa176(-4:mn6), wa(tid)%wa177(-4:mn6), wa(tid)%wa178(-4:mn6), wa(tid)%wa179(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa180(-4:mn6), wa(tid)%wa181(-4:mn6), wa(tid)%wa182(-4:mn6), wa(tid)%wa183(-4:mn6), wa(tid)%wa184(-4:mn6),  stat= ierr)
      allocate( wa(tid)%wa185(-4:mn6), wa(tid)%wa186(-4:mn6), wa(tid)%wa187(-4:mn6), wa(tid)%wa188(-4:mn6), wa(tid)%wa189(-4:mn6), stat= ierr )
      allocate( wa(tid)%wa190(-4:mn6), wa(tid)%wa191(-4:mn6), wa(tid)%wa192(-4:mn6), wa(tid)%wa193(-4:mn6), wa(tid)%wa194(-4:mn6), stat= ierr )
      allocate( wa(tid)%wa195(-4:mn6), wa(tid)%wa196(-4:mn6), wa(tid)%wa197(-4:mn6), wa(tid)%wa198(-4:mn6), wa(tid)%wa199(-4:mn6), stat= ierr )
      allocate( wa(tid)%wa1101(-4:mn6), wa(tid)%wa1102(-4:mn6), wa(tid)%wa1103(-4:mn6), wa(tid)%wa1104(-4:mn6), stat= ierr )
      allocate( wa(tid)%wa1105(-4:mn6), wa(tid)%wa1106(-4:mn6), stat= ierr )



      if (ierr /= 0) then
        call ERROR( 'getMemorySpace', 'Allocation of wa1', &
                    __FILE__, __LINE__  ) 
      endif

      allocate( wa(tid)%wa200(-4:mn5), wa(tid)%wa201(-4:mn5), wa(tid)%wa202(-4:mn5), wa(tid)%wa203(-4:mn5), wa(tid)%wa204(-4:mn5),  stat= ierr)
      allocate( wa(tid)%wa205(-4:mn5), wa(tid)%wa206(-4:mn5), wa(tid)%wa207(-4:mn5), wa(tid)%wa208(-4:mn5), wa(tid)%wa209(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa210(-4:mn5), wa(tid)%wa211(-4:mn5), wa(tid)%wa212(-4:mn5), wa(tid)%wa213(-4:mn5), wa(tid)%wa214(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa215(-4:mn5), wa(tid)%wa216(-4:mn5), wa(tid)%wa217(-4:mn5), wa(tid)%wa218(-4:mn5), wa(tid)%wa219(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa220(-4:mn5), wa(tid)%wa221(-4:mn5), wa(tid)%wa222(-4:mn5), wa(tid)%wa223(-4:mn5), wa(tid)%wa224(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa225(-4:mn5), wa(tid)%wa226(-4:mn5), wa(tid)%wa227(-4:mn5), wa(tid)%wa228(-4:mn5), wa(tid)%wa229(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa230(-4:mn5), wa(tid)%wa231(-4:mn5), wa(tid)%wa232(-4:mn5), wa(tid)%wa233(-4:mn5), wa(tid)%wa234(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa235(-4:mn5), wa(tid)%wa236(-4:mn5), wa(tid)%wa237(-4:mn5), wa(tid)%wa238(-4:mn5), wa(tid)%wa239(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa240(-4:mn5), wa(tid)%wa241(-4:mn5), wa(tid)%wa242(-4:mn5), wa(tid)%wa243(-4:mn5), wa(tid)%wa244(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa245(-4:mn5), wa(tid)%wa246(-4:mn5), wa(tid)%wa247(-4:mn5), wa(tid)%wa248(-4:mn5), wa(tid)%wa249(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa250(-4:mn5), wa(tid)%wa251(-4:mn5), wa(tid)%wa252(-4:mn5), wa(tid)%wa253(-4:mn5), wa(tid)%wa254(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa255(-4:mn5), wa(tid)%wa256(-4:mn5), wa(tid)%wa257(-4:mn5), wa(tid)%wa258(-4:mn5), wa(tid)%wa259(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa260(-4:mn5), wa(tid)%wa261(-4:mn5), wa(tid)%wa262(-4:mn5), wa(tid)%wa263(-4:mn5), wa(tid)%wa264(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa265(-4:mn5), wa(tid)%wa266(-4:mn5), wa(tid)%wa267(-4:mn5), wa(tid)%wa268(-4:mn5), wa(tid)%wa269(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa270(-4:mn5), wa(tid)%wa271(-4:mn5), wa(tid)%wa272(-4:mn5), wa(tid)%wa273(-4:mn5), wa(tid)%wa274(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa275(-4:mn5), wa(tid)%wa276(-4:mn5), wa(tid)%wa277(-4:mn5), wa(tid)%wa278(-4:mn5), wa(tid)%wa279(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa280(-4:mn5), wa(tid)%wa281(-4:mn5), wa(tid)%wa282(-4:mn5), wa(tid)%wa283(-4:mn5), wa(tid)%wa284(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa285(-4:mn5), wa(tid)%wa286(-4:mn5), wa(tid)%wa287(-4:mn5), wa(tid)%wa288(-4:mn5), wa(tid)%wa289(-4:mn5), stat= ierr)
      allocate( wa(tid)%wa290(-4:mn5), wa(tid)%wa291(-4:mn5), stat= ierr)

      if (ierr /= 0) then
        call ERROR( 'getMemorySpace', 'Allocation of wa2', &
                    __FILE__, __LINE__  ) 
      endif

      allocate( wa(tid)%wa300(mn1),             &
                wa(tid)%wa301(mn1),             &
                wa(tid)%wa302(mn1), stat= ierr )
      if (ierr /= 0) then
        call ERROR( 'getMemorySpace', 'Allocation of wa3', &
                    __FILE__, __LINE__  ) 
      endif
    ENDDO


!    allocate( wa409(mnx1,mny1,mnz1), stat= ierr) 
!    if (ierr /= 0) then
!      call ERROR( 'getMemorySpace', 'Allocation of wa409', &
!                  __FILE__, __LINE__  )
!    endif    


    allocate( wa400(mnx1,mny1), wa401(mnx1,mny1), wa402(mnx1,mny1), stat= ierr )
    if (ierr /= 0) then
      call ERROR( 'getMemorySpace', 'Allocation of wa40', &
                  __FILE__, __LINE__  ) 
    endif

    allocate( wa410(nx,nz), wa411(2,ny), wa412(6,nx,nz), stat= ierr )
    if (ierr /= 0) then
      call ERROR( 'getMemorySpace', 'Allocation of wa41', &
                  __FILE__, __LINE__  ) 
    endif

    !allocate( wa420(2,nyh), wa421(6,nx,nz), stat= ierr )
    if (ierr /= 0) then
      call ERROR( 'getMemorySpace', 'Allocation of wa42', &
                  __FILE__, __LINE__  ) 
    endif

    allocate( wa500(nx,5,nz), stat= ierr )
    if (ierr /= 0) then
      call ERROR( 'getMemorySpace', 'Allocation of wa5', &
                  __FILE__, __LINE__  )     
    endif

    allocate( wa600(nx,nz), stat= ierr )
    if (ierr /= 0) then
      call ERROR( 'getMemorySpace', 'Allocation of wa6', &
                  __FILE__, __LINE__  ) 
    endif

    allocate( wa700(mpi_dims1*mpi_dims2*mpi_dims3), wa701(mpi_dims1*mpi_dims2*mpi_dims3), stat= ierr )
    if (ierr /= 0) then
      call ERROR( 'getMemorySpace', 'Allocation of wa7', &
                  __FILE__, __LINE__  ) 
    endif

    maxele = 1000
!    allocate( walis(3,maxele), stat= ierr )
!    if (ierr /= 0) then
!      call ERROR( 'getMemorySpace', 'Allocation of lista', &
!                  __FILE__, __LINE__  ) 
!    endif
!--------------------------------------------------------------------------- END
  end subroutine getWAMemorySpace

  subroutine freeWAMemorySpace()
    implicit none
!------------------------------------------ Local Variables
    integer(is):: tid, nth
#   ifdef OPENMP
      integer(is):: omp_get_max_threads
#   endif

!------------------------------------------------------------------------- BEGIN
#   ifdef OPENMP
      nth = omp_get_max_threads()
#   else
      nth = 1
#   endif

!    deallocate( walis )
    deallocate( wa700, wa701 )
    deallocate( wa600 )
    deallocate( wa500 )
    deallocate( wa410, wa411, wa412)
    !deallocate( wa420, wa421 )
    deallocate( wa400, wa401, wa402 )
!    deallocate( wa409 )

    DO tid= 1, nth
      deallocate( wa(tid)%wa300, wa(tid)%wa301, wa(tid)%wa302 )

      deallocate( wa(tid)%wa200, wa(tid)%wa201, wa(tid)%wa202, wa(tid)%wa203, wa(tid)%wa204, & 
                  wa(tid)%wa205, wa(tid)%wa206, wa(tid)%wa207, wa(tid)%wa208, wa(tid)%wa209, &
                  wa(tid)%wa210, wa(tid)%wa211, wa(tid)%wa212, wa(tid)%wa213, wa(tid)%wa214, &
                  wa(tid)%wa215, wa(tid)%wa216, wa(tid)%wa217, wa(tid)%wa218, wa(tid)%wa219, &
                  wa(tid)%wa220, wa(tid)%wa221, wa(tid)%wa222, wa(tid)%wa223, wa(tid)%wa224, &
                  wa(tid)%wa225, wa(tid)%wa226, wa(tid)%wa227, wa(tid)%wa228, wa(tid)%wa229, &
                  wa(tid)%wa230, wa(tid)%wa231, wa(tid)%wa232, wa(tid)%wa233, wa(tid)%wa234, &
                  wa(tid)%wa235, wa(tid)%wa236, wa(tid)%wa237, wa(tid)%wa238, wa(tid)%wa239, &
                  wa(tid)%wa240, wa(tid)%wa241, wa(tid)%wa242, wa(tid)%wa243, wa(tid)%wa244, &
                  wa(tid)%wa245, wa(tid)%wa246, wa(tid)%wa247, wa(tid)%wa248, wa(tid)%wa249, &
                  wa(tid)%wa250, wa(tid)%wa251, wa(tid)%wa252, wa(tid)%wa253, wa(tid)%wa254, &
                  wa(tid)%wa255, wa(tid)%wa256, wa(tid)%wa257, wa(tid)%wa258, wa(tid)%wa259, &
                  wa(tid)%wa260, wa(tid)%wa261, wa(tid)%wa262, wa(tid)%wa263, wa(tid)%wa264, &
                  wa(tid)%wa265, wa(tid)%wa266, wa(tid)%wa267, wa(tid)%wa268, wa(tid)%wa269, &
                  wa(tid)%wa270, wa(tid)%wa271, wa(tid)%wa272, wa(tid)%wa273, wa(tid)%wa274, &
                  wa(tid)%wa275, wa(tid)%wa276, wa(tid)%wa277, wa(tid)%wa278, wa(tid)%wa279, &
                  wa(tid)%wa280, wa(tid)%wa281, wa(tid)%wa282, wa(tid)%wa283, wa(tid)%wa284, &
                  wa(tid)%wa285, wa(tid)%wa286, wa(tid)%wa287, wa(tid)%wa288, wa(tid)%wa289, &
                  wa(tid)%wa290, wa(tid)%wa291 )

      deallocate( wa(tid)%wa100, wa(tid)%wa101, wa(tid)%wa102, wa(tid)%wa103, wa(tid)%wa104, &
                  wa(tid)%wa105, wa(tid)%wa106, wa(tid)%wa107, wa(tid)%wa108, wa(tid)%wa109, &
                  wa(tid)%wa110, wa(tid)%wa111, wa(tid)%wa112, wa(tid)%wa113, wa(tid)%wa114, &
                  wa(tid)%wa115, wa(tid)%wa116, wa(tid)%wa117, wa(tid)%wa118, wa(tid)%wa119, &
                  wa(tid)%wa120, wa(tid)%wa121, wa(tid)%wa122, wa(tid)%wa123, wa(tid)%wa124, &
                  wa(tid)%wa125, wa(tid)%wa126, wa(tid)%wa127, wa(tid)%wa128, wa(tid)%wa129, &
                  wa(tid)%wa130, wa(tid)%wa131, wa(tid)%wa132, wa(tid)%wa133, wa(tid)%wa134, &
                  wa(tid)%wa135, wa(tid)%wa136, wa(tid)%wa137, wa(tid)%wa138, wa(tid)%wa139, &
                  wa(tid)%wa140, wa(tid)%wa141, wa(tid)%wa142, wa(tid)%wa143, wa(tid)%wa144, &
                  wa(tid)%wa145, wa(tid)%wa146, wa(tid)%wa147, wa(tid)%wa148, wa(tid)%wa149, &
                  wa(tid)%wa150, wa(tid)%wa151, wa(tid)%wa152, wa(tid)%wa153, wa(tid)%wa154, &
                  wa(tid)%wa155, wa(tid)%wa156, wa(tid)%wa157, wa(tid)%wa158, wa(tid)%wa159, &
                  wa(tid)%wa160, wa(tid)%wa161, wa(tid)%wa162, wa(tid)%wa163, wa(tid)%wa164, &
                  wa(tid)%wa165, wa(tid)%wa166, wa(tid)%wa167, wa(tid)%wa168, wa(tid)%wa169, &
                  wa(tid)%wa170, wa(tid)%wa171, wa(tid)%wa172, wa(tid)%wa173, wa(tid)%wa174, &
                  wa(tid)%wa175, wa(tid)%wa176, wa(tid)%wa177, wa(tid)%wa178, wa(tid)%wa179, &
                  wa(tid)%wa180, wa(tid)%wa181, wa(tid)%wa182, wa(tid)%wa183, wa(tid)%wa184, &
                  wa(tid)%wa185, wa(tid)%wa186, wa(tid)%wa187, wa(tid)%wa188, wa(tid)%wa189, &
                  wa(tid)%wa190, wa(tid)%wa191, wa(tid)%wa192, wa(tid)%wa193, wa(tid)%wa194, &
                  wa(tid)%wa195, wa(tid)%wa196, wa(tid)%wa197, wa(tid)%wa198, wa(tid)%wa199, &
                  wa(tid)%wa1101, wa(tid)%wa1102, wa(tid)%wa1103, wa(tid)%wa1104,            &
                  wa(tid)%wa1105, wa(tid)%wa1106 )
    ENDDO
!--------------------------------------------------------------------------- END
  end subroutine freeWAMemorySpace
END MODULE workarrays




!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL VARIABLES
! DESCRIPTION
! SEE ALSO
!===============================================================================
MODULE tiempo
  integer(is), parameter:: maxcont = 10

  integer(is) :: timer(2,maxcont)
  real(rs)    :: abs_time(maxcont), neventos(maxcont)

  CONTAINS

  subroutine initTiempo( uf )
    implicit none
!------------------------------------------ Input Variables
    integer(is) :: uf

!------------------------------------------ Local Variables
    integer(is) :: ii, counter, rate, maxCounter
    real(rs)    :: tt

!------------------------------------------------------------------------- BEGIN
    do ii= 1, maxcont
      timer(1,ii)  = 0
      timer(2,ii)  = 0
      neventos(ii) = 0.0
      abs_time(ii) = 0.0
    enddo

    call system_clock( counter, rate, maxCounter )
    tt = REAL(maxCounter)/REAL(rate)

    write(uf,*) ' '
    write(uf,'(A,F14.6,A,F9.6,A)' ) 'Max. time that can be measured = ', &
                                        tt, ' (s) = ', tt/3600.0, ' (h)'
    write(uf,'(A,E16.9,A)') 'Time precision = ', 1.0/REAL(rate), ' (s)'
    write(uf,*) ' '
!--------------------------------------------------------------------------- END
  end subroutine initTiempo

  subroutine medirTiempo( job, cont )
    implicit none
!------------------------------------------ Input Variables
    integer(is) :: job, cont

!------------------------------------------ Local Variables
    integer(is):: rate, contMax, itime

!------------------------------------------------------------------------- BEGIN
    if (cont < 1 .OR. cont > maxcont) then
      call ERROR( 'medirTiempo', 'Wrong cont value', __FILE__, __LINE__ )
    endif

    IF (job == 1) THEN
      call system_clock( timer(1,cont), rate, contMax )

    ELSE IF (job == 2) THEN
      call system_clock( timer(2,cont), rate, contMax )

      if (timer(2,cont) >= timer(1,cont)) then
        itime = timer(2,cont) - timer(1,cont)
      else
        itime = (contMax - timer(1,cont)) + (timer(2,cont) + 1)
      endif

      abs_time(cont) = abs_time(cont) + (REAL(itime) / REAL(rate))
      neventos(cont) = neventos(cont) + 1.0
    ELSE
      call ERROR( 'medirTiempo', 'Wrong job value', __FILE__, __LINE__ )
    ENDIF
!--------------------------------------------------------------------------- END
  end subroutine medirTiempo

  subroutine showTiempo( uf, cont, message )
    implicit none
!------------------------------------------ Input Variables
    integer(is)  :: uf, cont
    character(*) :: message

!------------------------------------------ Local Variables
    real(rs)   :: tt

!------------------------------------------------------------------------- BEGIN
    tt = abs_time(cont)

    if (neventos(cont) <= 1.0) then
      write(uf,'(A,A,A,F14.6,A,F9.6,A)') '   ', message, '     time = ', &
            tt, ' (s) = ', tt/3600.0E0, ' (h)'
    else
      write(uf,'(A,A,A,F14.6,A,F9.6,A,A,F14.6,A)') &
        '   ', message, '     time = ', tt, ' (s) = ', tt/3600.0E0, ' (h)', &
        '    Mean time =', tt/neventos(cont), ' (s)'
    endif
    !flush(uf)
!--------------------------------------------------------------------------- END
  end subroutine showTiempo
END MODULE tiempo
