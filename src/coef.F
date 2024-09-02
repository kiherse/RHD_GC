!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL VARIABLES
!     DESCRIPTION
!     SEE ALSO
!     ADAPTED FROM PROMETHEUS
!     MODIFIED 09-JUL-99: ONE BUG FIXED IN THE STENCILS OF do LOOPS
!===============================================================================

      subroutine coef( tid, MMxx, N, DX, COEFF1, COEFF2, COEFF3, COEFF4, COEFF5 )
      USE workarrays
      implicit none
# include "types.h"

!-------------------------------------------------------Input Variables

      integer(is):: tid, MMxx, N

! size -4:mn5
      real(rs)   :: DX(-4:*)

!------------------------------------------------------Output Variables
! size -4:mn5
      real(rs)   :: COEFF1(-4:*), COEFF2(-4:*), COEFF3(-4:*), &
      COEFF4(-4:*), COEFF5(-4:*)

!------------------------------------------------------Local Variables

      integer(is):: i

      real(rs), pointer :: SCRCH1(:), SCRCH2(:), SCRCH3(:), &
      SCRCH4(:), SCRCH5(:)

!--------------------------------------------------------------------------BEGIN

! -4:mn6 vectors
      SCRCH1 => wa(tid)%wa138
      SCRCH2 => wa(tid)%wa139
      SCRCH3 => wa(tid)%wa140
      SCRCH4 => wa(tid)%wa141
      SCRCH5 => wa(tid)%wa142

      do i= MMxx-2, N+4
         SCRCH1(i) = DX(i) + DX(i-1)
         SCRCH2(i) = SCRCH1(i) + DX(i)
         SCRCH3(i) = SCRCH1(i) + DX(i-1)
      enddo

      do i= MMxx-2, N+3
         SCRCH4(i) = DX(i)/(SCRCH1(i)+ DX(i+1))
         COEFF1(i) = SCRCH4(i)*SCRCH3(i)/SCRCH1(i+1)
         COEFF2(i) = SCRCH4(i)*SCRCH2(i+1)/SCRCH1(i)
      enddo

      do i= MMxx-2, N+2
         SCRCH5(i) = 0.5*((SCRCH1(i)+SCRCH1(i+2))+(SCRCH1(i+2)+SCRCH1(i)))
         SCRCH4(i) = 1.0/SCRCH5(i)
         COEFF3(i) = -SCRCH4(i)*DX(i)*SCRCH1(i)/SCRCH3(i+1)
         COEFF4(i) = SCRCH4(i)*DX(i+1)*SCRCH1(i+2)/SCRCH2(i+1)

         SCRCH5(i) = 0.5*((DX(i+1)*COEFF3(i)+DX(i)*COEFF4(i))+ &
         (DX(i)*COEFF4(i)+DX(i+1)*COEFF3(i)))
         COEFF5(i) = DX(i)-2.0*SCRCH5(i)
         COEFF5(i) = COEFF5(i)/SCRCH1(i+1)
      enddo
      
!--------------------------------------------------------------------END

      END subroutine coef

!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL VARIABLES
!     DESCRIPTION
!     SEE ALSO
!     ADAPTED FROM PROMETHEUS
!     MODIFIED 09-JUL-99: One bug fixed in the stencils of do loops
!===============================================================================

      subroutine interp( tid, MMxx, N, small, COEFF1, COEFF2, COEFF3, COEFF4, COEFF5, &
      U, UM, UP, DELU )
      USE workarrays
      implicit none
# include "types.h"

!-------------------------------------------------------Input Variables

      integer(is):: tid, MMxx, N

      real(rs)   :: small

! size -4:mn5
      real(rs)   :: COEFF1(-4:*), COEFF2(-4:*), COEFF3(-4:*),  &
      COEFF4(-4:*), COEFF5(-4:*)

! size -4:mn5
      real(rs)   :: U(-4:*)

!------------------------------------------------------Output Variables

! size -4:mn5
      real(rs)   :: UM(-4:*), UP(-4:*), DELU(-4:*)

!-------------------------------------------------------Local Variables

      integer(is):: i

      real(rs)   :: SDELU, small2

      real(rs), pointer :: SCRCH1(:), SCRCH2(:), SCRCH3(:),  &
      SCRCH4(:), SCRCH5(:)

!---------------------------------------------------------------------------BEGIN

! -4:mn6 vectors
      SCRCH1 => wa(tid)%wa138
      SCRCH2 => wa(tid)%wa139
      SCRCH3 => wa(tid)%wa140
      SCRCH4 => wa(tid)%wa141
      SCRCH5 => wa(tid)%wa142


      do i= MMxx-3, N+4
         SCRCH1(i) = U(i)-U(i-1)
      enddo
!      print*,'interp 1',cartrank

      do i=MMxx-2,N+3
         DELU(i)  = COEFF1(i)*SCRCH1(i+1)+COEFF2(i)*SCRCH1(i)
      enddo

!      print*,'interp 2',cartrank


      do i= MMxx-2, N+3
         if (SCRCH1(i+1)*SCRCH1(i) > 0.0) then
            SCRCH2(i)= MIN(ABS(SCRCH1(i)), ABS(SCRCH1(i+1)))
            SDELU    = DELU(i)/ABS(DELU(i))
            DELU(i)  = MIN(ABS(DELU(i)),2.0*SCRCH2(i))*SDELU
         else
            DELU(i)  = 0.0
         endif
      enddo

!      print*,'interp 3',cartrank


!      small2 = SMALL*SMALL
      do i= MMxx-2, N+2
         UP(i)= 0.5*((U(i) + COEFF5(i)*SCRCH1(i+1)) + (COEFF5(i)*SCRCH1(i+1)+U(i)))

         SCRCH3(i)= COEFF3(i)*DELU(i+1) + COEFF4(i)*DELU(i)
         SCRCH4(i)= COEFF4(i)*DELU(i)   + COEFF3(i)*DELU(i+1)
         SCRCH5(i)= 0.5*(SCRCH3(i)+SCRCH4(i))

         UP(i)    = UP(i)+SCRCH5(i)

!         if (ABS(UP(i)) <= small2) then
!UP(i) = 0.0
!         endif

         UM(i+1) = UP(i)
      enddo

!            print*,'interp 4',cartrank

!-----------------------------------------------------------------------------END

      END subroutine interp

!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL VARIABLES
!     DESCRIPTION
!     SEE ALSO
!     ADAPTED FROM PROMETHEUS
!===============================================================================

      subroutine detect( tid, MMxx, N, smlrho, small, &
      U, DELU, X, DX, P, RHO, UM, UP )
      USE workarrays
      implicit none
# include "types.h"

!-------------------------------------------------------Input Variables

      integer(is):: tid, MMxx, N

      real(rs)   :: gamma_ad, smlrho, small

! size -4:mn5
      real(rs)   :: U(-4:*), DELU(-4:*), &
      X(-4:*), DX(-4:*),   &
      P(-4:*), RHO(-4:*)

!------------------------------------------------------Output Variables

! size -4:mn5
      real(rs)   :: UM(-4:*), UP(-4:*)

!-------------------------------------------------------Local Variables
      integer(is):: i

      real(rs)   :: ETA1,ETA2,EPSLN,AK0

      real(rs), pointer :: SCRCH1(:), SCRCH2(:), SCRCH3(:),   &
      SCRCH4(:)

      real(rs), pointer :: D2U(:), ETA(:), ETATIL(:)

!---------------------------------------------------------------------------BEGIN

! -4:mn6 vectors
      SCRCH1 => wa(tid)%wa138
      SCRCH2 => wa(tid)%wa139
      SCRCH3 => wa(tid)%wa140
      SCRCH4 => wa(tid)%wa141

      D2U    => wa(tid)%wa142
      ETA    => wa(tid)%wa143
      ETATIL => wa(tid)%wa144

!     The following parameters are set as in colella & woodward
!     (jcp, 54, 174 (1984))
      ETA1  = 5.0
      ETA2  = 0.05
      EPSLN = 0.1
      AK0   = 1.0

      do i= MMxx-2, N+3
         SCRCH1(i) = DX(i) + DX(i-1)
         SCRCH2(i) = SCRCH1(i) + DX(i+1)
         SCRCH3(i) = U(i) - U(i-1)
         SCRCH1(i) = SCRCH3(i)/SCRCH1(i)
      enddo

      do i= MMxx-2, N+2
         D2U(i)    = (SCRCH1(i+1) - SCRCH1(i))/SCRCH2(i)
         SCRCH4(i) = X(i)-X(i-1)
         SCRCH4(i) = SCRCH4(i)*SCRCH4(i)*SCRCH4(i)
      enddo

      do i= MMxx-1, N+1
         SCRCH1(i) = D2U(i+1)*D2U(i-1)
         SCRCH3(i) = ABS(U(i+1)-U(i-1))
         SCRCH3(i) = SCRCH3(i)-EPSLN*MIN(ABS(U(i+1)),ABS(U(i-1)))
      enddo

      do i= MMxx-1, N+1
         if ((U(i+1)-U(i-1))==0.0) then
            SCRCH2(i) = SMALL*SMLRHO
         else
            SCRCH2(i) = U(i+1) - U(i-1)
         endif

         if ((SCRCH1(i)>0.0) .OR. (SCRCH3(i)<0.0)) then
            ETATIL(i) = 0.0
         else
            ETATIL(i) = (D2U(i-1)-D2U(i+1))*(SCRCH4(i)+SCRCH4(i+1))
            ETATIL(i) = ETATIL(i)/(X(i+1)-X(i-1))/SCRCH2(i)
         endif
      enddo

      do i= MMxx-1, N+1
         ETA(i)    = MAX(ETA1-ETA1, MIN( ETA1*(ETATIL(i)-ETA2), ETA1/ETA1))
         SCRCH1(i) = ABS(P(i+1)-P(i-1))/MIN(P(i+1),P(i-1))
         SCRCH2(i) = ABS(RHO(i+1)-RHO(i-1))/MIN(RHO(i+1),RHO(i-1))
      enddo
! USE GAMMA=4./3. because it is too expensive to compute with Synge EOS
      do i= MMxx-1, N+1
         if (4./3.*AK0*SCRCH2(i) < SCRCH1(i)) then
            ETA(i) = 0.0
         endif
      enddo

      do i= MMxx-1, N+1, 1
         SCRCH1(i) = U(i-1) + 0.5*DELU(i-1)
         SCRCH2(i) = U(i+1) - 0.5*DELU(i+1)
         UM(i)     = UM(i)*(1.0-ETA(i)) + SCRCH1(i)*ETA(i)
         UP(i)     = UP(i)*(1.0-ETA(i)) + SCRCH2(i)*ETA(i)
      enddo

!-----------------------------------------------------------------------------END

      END subroutine detect

!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL VARIABLES
!     DESCRIPTION
!     SEE ALSO
!     ADAPTED FROM PROMETHEUS
!===============================================================================

      subroutine monot( tid, MMxx, N, small, U, UM, UP, DU, U6 )
      USE workarrays
      implicit none
# include "types.h"

!-------------------------------------------------------Input Variables

      integer(is):: tid, MMxx, N

      real(rs)   :: small

! size -4:mn5
      real(rs)   :: U(-4:*)

!------------------------------------------------------Output Variables

! size -4:mn5
      real(rs)   :: UM(-4:*), UP(-4:*), DU(-4:*), U6(-4:*)

!-------------------------------------------------------Local Variables

      integer(is):: i

      real(rs), pointer :: SCRCH1(:), SCRCH2(:), SCRCH3(:)

!---------------------------------------------------------------------------BEGIN

! -4:mn6 vectors
      SCRCH1 => wa(tid)%wa138
      SCRCH2 => wa(tid)%wa139
      SCRCH3 => wa(tid)%wa140


      do i= MMxx-1, N+1
         DU(i)     = UP(i)-UM(i)
         SCRCH1(i) = UP(i)-U(i)
         SCRCH1(i) = SCRCH1(i)*(UM(i)-U(i))
      enddo

      do i= MMxx-1, N+1
         if (SCRCH1(i)>=0.0) then
            UM(i) = U(i)
            UP(i) = U(i)
         endif
      enddo

      do i= MMxx-1, N+1
         DU(i)     = UP(i)-UM(i)
         SCRCH1(i) = (UP(i)-U(i))*(UM(i)-U(i))
         if (SCRCH1(i)==0.0) then
            SCRCH2(i) = UM(i)
            SCRCH3(i) = UP(i)
         else
            SCRCH2(i) = 3.0*U(i)-2.0*UP(i)
            SCRCH3(i) = 3.0*U(i)-2.0*UM(i)
         endif
      enddo

      do i= MMxx-1, N+1
         if (DU(i)*(UM(i)-SCRCH2(i))<0.0) then
            UM(i) = SCRCH2(i)
         endif

         if (DU(i)*(SCRCH3(i)-UP(i))<0.0) then
            UP(i) = SCRCH3(i)
         endif

         if (ABS(UM(i))<=SMALL*SMALL) then
            UM(i) = 0.0
         endif

         if (ABS(UP(i))<=SMALL*SMALL) then
            UP(i) = 0.0
         endif
      enddo

      do i= MMxx-1, N+1
         DU(i) = UP(i)-UM(i)
         U6(i) = 6.0*U(i)-3.0*(UM(i)+UP(i))
      enddo

!-----------------------------------------------------------------------------END

      END subroutine monot

!===============================================================================
!     NAME
!     F90 SPECIFICATION
!     ARGUMENTS
!     GLOBAL VARIABLES
!     DESCRIPTION
!     SEE ALSO
!     ADAPTED FROM PROMETHEUS
!     THIS SUBROUTINE FLATTENS ZONE STRUCTURE IN REGIONS WHERE SHOCKS
!     ARE TOO THIN.
!     THIS VERSION ONLY USES THE SIMPLEST FORM OF DISSIPATION AS
!     DESCRIBED IN THE APPENDIX OF COLELLA AND WOODWARD (JCP, VOL. 54,
!     174 (1984)).
!===============================================================================

      subroutine flaten( tid, MMxx, N, igodu, epsiln, omg1, omg2, P, VEL, FLATN, FLATN1 )
      USE workarrays
      implicit none
# include "types.h"

!-------------------------------------------------------Input Variables

      integer(is):: tid, MMxx, N

      integer(is):: igodu

      real(rs)   :: epsiln, omg1, omg2

! size -4:mn5
      real(rs)   :: P(-4:*), VEL(-4:*)

!------------------------------------------------------Output Variables

! size -4:mn5
      real(rs)   :: FLATN(-4:*), FLATN1(-4:*)

!-------------------------------------------------------Local Variables

      integer(is):: i               

      real(rs)   :: DP2

      real(rs), pointer :: SCRCH1(:), SCRCH2(:), SCRCH3(:)

      real(rs), pointer :: DP(:), DVEL(:)

!---------------------------------------------------------------------------BEGIN

! -4:mn6 vectors
      SCRCH1 => wa(tid)%wa138
      SCRCH2 => wa(tid)%wa139
      SCRCH3 => wa(tid)%wa140

! -4:mn5 vectors
      DP   => wa(tid)%wa275
      DVEL => wa(tid)%wa276


      do i= MMxx-3, N+3
         DP(i)     = P(i+1)-P(i-1)
         DVEL(i)   = VEL(i+1)-VEL(i-1)
         SCRCH1(i) = EPSILN*MIN(P(i+1),P(i-1))-ABS(DP(i))

         if (SCRCH1(i)<0.0 .AND. DVEL(i)<0.0) then
            SCRCH1(i) = 1.0
         else
            SCRCH1(i) = 0.0
         endif
      enddo

      do i= MMxx-2, N+2
         DP2 = P(i+2)-P(i-2)
         if (DP2==0.0) then
            if (DP(i)==0.0) then
               SCRCH2(i) = -OMG1
            else
               SCRCH2(i) = 1.0-OMG1
            endif
         else
            SCRCH2(i) = DP(i)/DP2 - OMG1
         endif
         SCRCH3(i) = SCRCH1(i)*MAX(SCRCH2(i)-SCRCH2(i), SCRCH2(i)*OMG2)
      enddo

      do i= MMxx-1, N+1
         if (DP(i)<0.0) then
            SCRCH2(i) = SCRCH3(i+1)
         else
            SCRCH2(i) = SCRCH3(i-1)
         endif
      enddo

      do i= MMxx-1, N+1
         FLATN(i) = MAX(SCRCH2(i),SCRCH3(i))
         IF (FLATN(i) == 0.) THEN
            FLATN(i) = MAX(FLATN(i)-FLATN(i), 0.)                                                                            
         ELSE
            FLATN(i) = MAX(FLATN(i)-FLATN(i), MIN(FLATN(i)/FLATN(i) , FLATN(i)))
         ENDIF
         IF (IGODU==1) then
            FLATN(i) = 1.0
         endif
         FLATN1(i)= 1.0-FLATN(i)
      enddo

!-----------------------------------------------------------------------------END

      END subroutine flaten
