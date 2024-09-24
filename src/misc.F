
!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL VARIABLES
!     DESCRIPTION
!     SEE ALSO
!===============================================================================
      subroutine ERROR( routine, message, fich, line )
      implicit none
# include "types.h"
!------------------------------------------Input Variables
      character(*) :: routine
      character(*) :: message
      character(*) :: fich
      integer(is)  :: line

!-------------------------------------------------------------------------BEGIN
      write(*,*) '*************************************************************'
      write(*,*) 'ERROR'
      write(*,*) ' '
      write(*,*) routine, ': ', message
      write(*,*) ' '
      write(*,*) '           File: ', fich
      write(*,*) '           Line: ', line
      write(*,*) ' '
      write(*,*) '*************************************************************'

      STOP 'ERROR'
!---------------------------------------------------------------------------END
      END subroutine ERROR

!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL VARIABLES
!     DESCRIPTION
!     SEE ALSO
!===============================================================================
      function int2ch4( number )
      implicit none
# include "types.h"
!------------------------------------------Input Variables
      integer(is):: number

!------------------------------------------Output Variables
      character(4):: int2ch4

!------------------------------------------Local Variables
      integer(is):: tt, cero, jj
      character  :: cc(4)

!-------------------------------------------------------------------------BEGIN
      if (number > 9999 .OR. number < 0) then
         call ERROR( 'int2ch4', 'The number is out of range', __FILE__, __LINE__  )
      else
         cero = ICHAR( '0' )    ! 48 is the ascii code of zero
         tt   = number
         do jj = 4, 1, -1
            cc(jj) = CHAR(cero + MOD(tt,10))
            tt = tt/10
         enddo
         int2ch4 = cc(1)//cc(2)//cc(3)//cc(4)
      endif

      RETURN
!---------------------------------------------------------------------------END
      END function int2ch4

!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!===============================================================================
      subroutine showParameters( suffixd, rstfil, itdelay,                 &
      gamma_ad, velb, rhob, machb, pb, dk1, dk2, psh1, psh2, &
      nn, mm1, mm2, t0, t1, vptbd, rptbd, pptbd )
      implicit none
# include "types.h"
!------------------------------------------Input  Variables
      character(3):: suffixd
      character(7):: rstfil

      integer(is) :: itdelay

      real(rs)    :: gamma_ad, velb, rhob, machb, pb, dk1, dk2, psh1, psh2, &
      nn, mm1, mm2, t0, t1, vptbd, rptbd, pptbd

!------------------------------------------Local  Variables
      character(5):: txtxt

!--------------------------------------------------------------------------BEGIN
      txtxt = '.....'

      write(6,*) ' '
      write(6,*) '---------------- MAIN (showParameters) ------------------'
      write(6,*) 'gamma   ', txtxt, gamma_ad
      write(6,*) 'velb    ', txtxt, velb
      write(6,*) 'rhob    ', txtxt, rhob
      write(6,*) 'machb   ', txtxt, machb
      write(6,*) 'pb      ', txtxt, pb
      write(6,*) 'dk1     ', txtxt, dk1
      write(6,*) 'dk2     ', txtxt, dk2
      write(6,*) 'psh1    ', txtxt, psh1
      write(6,*) 'psh2    ', txtxt, psh2
      write(6,*) 'nn      ', txtxt, nn
      write(6,*) 'mm1     ', txtxt, mm1
      write(6,*) 'mm2     ', txtxt, mm2
      write(6,*) 't0      ', txtxt, t0
      write(6,*) 't1      ', txtxt, t1
      write(6,*) 'vptbd   ', txtxt, vptbd
      write(6,*) 'rptbd   ', txtxt, rptbd
      write(6,*) 'pptbd   ', txtxt, pptbd
      write(6,*) 'itdelay ', txtxt, itdelay
      write(6,*) 'suffixd ', txtxt, suffixd
      write(6,*) '--------------------------------------------------------'
      write(6,*) ' '
!----------------------------------------------------------------------------END
      END subroutine showParameters
