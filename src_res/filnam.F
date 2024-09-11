!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL Variables
!     DESCRIPTION
!     SEE ALSO
!     This subroutine constructs new filenames for output and restart
!     files.
!     MODIFIED 22-DEC-1995 to construct names for time delay files
!===============================================================================
      subroutine filnam( chd, basenm, datfil, rstfil, outfil, tdlfil, &
      suffixd, suffix )
      implicit none
# include "types.h"
!-------------------------------------------------------Input Variables
      character(1):: chd
      character(4):: basenm

!------------------------------------------------------Output Variables
      character(8):: datfil
      character(7):: rstfil, outfil, tdlfil
      character(3):: suffixd
      character(2):: suffix

!-------------------------------------------------------Local Variables
      character(1):: sf1, sf2, sf3

      integer(is) :: isf1, isf2, isf3

!-------------------------------------------------------------------------BEGIN
      IF (chd/='d') THEN
         if (suffix(2:2)=='z'.or.suffix(2:2)=='Z') then
            sf1  = suffix(1:1)
            isf1 = ICHAR(sf1)
            sf2  = suffix(2:2)
            isf2 = ICHAR(sf2)

            isf1 = isf1 + 1
            isf2 = isf2 - 25
            suffix(1:1) = CHAR(isf1)
            suffix(2:2) = CHAR(isf2)
         else
            sf2  = suffix(2:2)
            isf2 = ICHAR(sf2)
            isf2 = isf2 + 1
            suffix(2:2) = CHAR(isf2)
         endif

         rstfil = basenm // 'R' // suffix
         outfil = basenm // 'O' // suffix

      ELSE
         if (suffixd(3:3)=='z'.or.suffixd(3:3)=='Z') then
            sf1  = suffixd(1:1)
            isf1 = ICHAR(sf1)
            sf2  = suffixd(2:2)
            isf2 = ICHAR(sf2)
            sf3  = suffixd(3:3)
            isf3 = ICHAR(sf3)

            if (suffixd(2:2)=='z'.or.suffixd(2:2)=='Z') then
               isf1 = isf1 + 1
               isf2 = isf2 - 25
               isf3 = isf3 - 25
               suffixd(1:1) = CHAR(isf1)
               suffixd(2:2) = CHAR(isf2)
               suffixd(3:3) = CHAR(isf3)
            else
               isf2 = isf2 + 1
               isf3 = isf3 - 25
               suffixd(2:2) = CHAR(isf2)
               suffixd(3:3) = CHAR(isf3)
            endif
         else
            sf3  = suffixd(3:3)
            isf3 = ICHAR(sf3)
            isf3 = isf3 + 1
            suffixd(3:3) = CHAR(isf3)
         endif

         tdlfil = basenm // suffixd
         datfil = basenm // 'D' // suffixd
      ENDIF
!---------------------------------------------------------------------------END
      END subroutine filnam
