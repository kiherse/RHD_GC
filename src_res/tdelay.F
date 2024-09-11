!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL Variables
! DESCRIPTION
! SEE ALSO
!     This subroutine writes files for time delay.
!===============================================================================
subroutine tdelay( basenm, igeomx, igeomy, igeomz, nx, ny, nz, &
                   mnx1, mnx5, mny1, mny5, mnz1, mnz5, nstep, timej, gamma_ad,   &
                   gridlx, gridly, gridlz, xznl, yznl, zznl, physic, &
                   datfil, rstfil, outfil, tdlfil, suffixd, suffix )
  USE tipos
  USE parallel
  USE workarrays
  implicit none
# include "types.h"
!------------------------------------------------------- Input Variables
  character(4)  :: basenm

  integer(is)   :: igeomx, igeomy, igeomz, nx, ny, nz, &
                   mnx1, mnx5, mny1, mny5, mnz1, mnz5, nstep

  real(rs)      :: timej, gamma_ad, gridlx, gridly, gridlz

  real(rs)      :: xznl(mnx1), yznl(mny1), zznl(mnz1)

  type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)

!------------------------------------------------------ Output Variables
  character(8)  :: datfil
  character(7)  :: rstfil, outfil, tdlfil
  character(3)  :: suffixd
  character(2)  :: suffix

!------------------------------------------------------- Local Variables
  character(4),external :: int2ch4

  character(256):: fich

  integer(is)   :: ierr, nxt, i, j, k

!------------------------------------------------------------------------- BEGIN
  fich = TRIM(outRoot)//tdlfil//'_'//int2ch4(cartrank)

  open( 2, file= fich, form= 'unformatted', iostat= ierr )
  if (ierr /=  0) then
    call ERROR( 'tdelay: Error openning file ', fich, __FILE__, __LINE__  )
  endif

  nxt = 1
  do j= 1, ny
    do k= 1, nz
      do i= 1, nx
        if (physic(i,j,k)%tracer > 1.0E-2) then
          nxt = nxt+1
        endif
      enddo
    enddo
  enddo


  write(2,*) timej, gamma_ad, xznl, yznl, zznl, gridlx, gridly, gridlz, &
             nx, ny, nz, nstep, igeomx, igeomy, igeomz
  write(2,*) nxt

  do j= 1, ny
    do k= 1, nz
      do i= 1, nx
        if (physic(i,j,k)%tracer > 1.0E-2) then
          write(2,*) i, j, k,                                   &
                     physic(i,j,k)%densty, physic(i,j,k)%denstye, &
                     physic(i,j,k)%velx,  &
                     physic(i,j,k)%vely, physic(i,j,k)%velz,    &
                     physic(i,j,k)%pres, physic(i,j,k)%tracer 
        endif
      enddo
    enddo
  enddo
         
  close(2)

  call filnam( 'd', basenm, datfil, rstfil, outfil, tdlfil, suffixd, suffix )
!----------------------------------------------------------------------------END
END subroutine tdelay
