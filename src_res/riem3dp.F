!===============================================================================
! NAME
! F90 SPECIFICATION
! ARGUMENTS
! GLOBAL Variables
! DESCRIPTION
! SEE ALSO
!     THIS SUBROUTINE GIVES A LINEARIZED SOLUTION OF THE RIEMANN PROBLEM
!     IN RELATIVISTIC HYDRODYNAMICS IN THE PRESENCE OF A TRANVERSAL
!     VELOCITY FIELD, ACCORDING TO MARQUINA'S SOLVER (JCP, SUBMITTED)
!     MODIFIED 7/12/01 IN ORDER TO AVOID INVERTING MATRICES
!     FOLLOWING DONAT ET AL. (1998). PERUCHO.
!===============================================================================
subroutine riem3dp( tid, M, N, RHOL, RHOR, UL, UR, VELL, VELR,           &
                    VELTL, VELTR, VELTTL, VELTTR, WL, WR,                &
                    CSL, CSR, PL, PR, HL, HR, DPDRHL, DPDRHR, DPDUL,     &
                    DPDUR, RL, RR, ML, MR, MTL, MTR, MTTL, MTTR, EL, ER, &
                    FR, FM, FMT, FMTT, FE )


  USE workarrays
  implicit none
# include "types.h"
!------------------------------------------------------- Input Variables
  integer(is):: tid, M, N

  real(rs):: gamma_ad

  ! size -4:mn6
  real(rs):: RHOL(-4:*),   RHOR(-4:*),   UL(-4:*),    UR(-4:*),    &
             VELL(-4:*),   VELR(-4:*),   VELTL(-4:*), VELTR(-4:*), &
             VELTTL(-4:*), VELTTR(-4:*), WL(-4:*),    WR(-4:*),    &
             CSL(-4:*),    CSR(-4:*),    PL(-4:*),    PR(-4:*),    &
             HL(-4:*),     HR(-4:*),     DPDRHL(-4:*),             &
             DPDRHR(-4:*), DPDUL(-4:*),  DPDUR(-4:*),              &
             RL(-4:*),     RR(-4:*),     ML(-4:*),    MR(-4:*),    &
             MTL(-4:*),    MTR(-4:*),    MTTL(-4:*),  MTTR(-4:*),  &
             EL(-4:*),     ER(-4:*)

!------------------------------------------------------ Output Variables
  ! size -4:mn6
  real(rs):: FR(-4:*), FM(-4:*), FMT(-4:*), FMTT(-4:*), FE(-4:*)

!------------------------------------------------------- Local Variables
  integer(is):: i

  real(rs):: KAPPAL, KAPPAR, DELTAL, DELTAR, S1L, S5L, S1R, S5R, &
             SCRCH6, SCRCH7, SCRCH8, SCRCH9,                     &
             SCRCH10, SCRCH11, SCRCH12, SCRCH13,                 &
             aux1, aux2, aux3, aux4
  
  real(rs), pointer :: A1L(:), A5L(:), A1R(:), A5R(:), &
                       PSI12(:), PSI52(:)

  real(rs), pointer :: LAMB1L(:), LAMB1R(:),                   &
                       LAMB2L(:), LAMB2R(:),                   &
                       LAMB5L(:), LAMB5R(:),                   &
                       R11(:), R15(:),                         &
                       R21(:), R25(:),                         &
                       R31(:), R35(:),                         &
                       R41(:), R45(:),                         &
                       R51(:), R55(:),                         &
                       PSI1(:), PSI2(:), PSI5(:),              &
                       Q1(:),  Q2(:),  Q3(:),  Q4(:),  Q5(:),  &
                       Q1L(:), Q2L(:), Q3L(:), Q4L(:), Q5L(:), &
                       Q1R(:), Q2R(:), Q3R(:), Q4R(:), Q5R(:), &
                       SCRCH1(:), SCRCH2(:), SCRCH3(:), SCRCH4(:)

!------------------------------------------------------------------------- BEGIN

  ! -4:mn6 vectors
  A1L   => wa(tid)%wa138
  A5L   => wa(tid)%wa139
  A1R   => wa(tid)%wa140
  A5R   => wa(tid)%wa141
  PSI12 => wa(tid)%wa142
  PSI52 => wa(tid)%wa143

  LAMB1L => wa(tid)%wa150
  LAMB1R => wa(tid)%wa151
  LAMB2L => wa(tid)%wa152
  LAMB2R => wa(tid)%wa153
  LAMB5L => wa(tid)%wa154
  LAMB5R => wa(tid)%wa155
  R11    => wa(tid)%wa156
  R15    => wa(tid)%wa157
  R21    => wa(tid)%wa158
  R25    => wa(tid)%wa159
  R31    => wa(tid)%wa160
  R35    => wa(tid)%wa161
  R41    => wa(tid)%wa162
  R45    => wa(tid)%wa163
  R51    => wa(tid)%wa164
  R55    => wa(tid)%wa165
  PSI1   => wa(tid)%wa166
  PSI2   => wa(tid)%wa167
  PSI5   => wa(tid)%wa168
  Q1     => wa(tid)%wa169
  Q2     => wa(tid)%wa170
  Q3     => wa(tid)%wa171
  Q4     => wa(tid)%wa172
  Q5     => wa(tid)%wa173
  Q1L    => wa(tid)%wa174
  Q2L    => wa(tid)%wa175
  Q3L    => wa(tid)%wa176
  Q4L    => wa(tid)%wa177
  Q5L    => wa(tid)%wa178
  Q1R    => wa(tid)%wa179
  Q2R    => wa(tid)%wa180
  Q3R    => wa(tid)%wa181
  Q4R    => wa(tid)%wa182
  Q5R    => wa(tid)%wa183
  SCRCH1 => wa(tid)%wa184
  SCRCH2 => wa(tid)%wa185
  SCRCH3 => wa(tid)%wa186
  SCRCH4 => wa(tid)%wa187

!     EIGENVALUES
  DO i= M, N
    SCRCH1(i) = 1.0 - 1.0/WL(i)**2

    SCRCH2(i) = 1.0-(VELL(i)*VELL(i)+(SCRCH1(i)-VELL(i)*VELL(i))*CSL(i)*CSL(i))

    SCRCH2(i) = SQRT((1.0-SCRCH1(i))*SCRCH2(i))
    SCRCH3(i) = 1.0-SCRCH1(i)*CSL(i)*CSL(i)

    LAMB1L(i) = VELL(i)*(1.0-CSL(i)*CSL(i))-CSL(i)*SCRCH2(i)
    LAMB1L(i) = LAMB1L(i)/SCRCH3(i)
    LAMB2L(i) = VELL(i)
    LAMB5L(i) = VELL(i)*(1.0-CSL(i)*CSL(i))+CSL(i)*SCRCH2(i)
    LAMB5L(i) = LAMB5L(i)/SCRCH3(i)

    SCRCH1(i) = 1.0 - 1.0/WR(i)**2

    SCRCH2(i) = 1.0-(VELR(i)*VELR(i)+(SCRCH1(i)-VELR(i)*VELR(i))*CSR(i)*CSR(i))

    SCRCH2(i) = SQRT((1.0-SCRCH1(i))*SCRCH2(i))
    SCRCH3(i) = 1.0-SCRCH1(i)*CSR(i)*CSR(i)

    LAMB1R(i) = VELR(i)*(1.0-CSR(i)*CSR(i))-CSR(i)*SCRCH2(i)
    LAMB1R(i) = LAMB1R(i)/SCRCH3(i)
    LAMB2R(i) = VELR(i)
    LAMB5R(i) = VELR(i)*(1.0-CSR(i)*CSR(i))+CSR(i)*SCRCH2(i)
    LAMB5R(i) = LAMB5R(i)/SCRCH3(i)

    PSI1(i)   = MAX(ABS(LAMB1L(i)),ABS(LAMB1R(i)))
    PSI2(i)   = MAX(ABS(LAMB2L(i)),ABS(LAMB2R(i)))
    PSI5(i)   = MAX(ABS(LAMB5L(i)),ABS(LAMB5R(i)))

!     10/12/01...SAVING CALCULATIONS
    PSI12(i) = PSI1(i) - PSI2(i)
    PSI52(i) = PSI5(i) - PSI2(i)
  ENDDO


!     LEFT NUMERICAL VISCOSITY
!     RIGHT EIGENVECTORS
!     RIJ STANDS FOR THE COMPONENT i OF THE RIGHT EIGENVECTOR J
!     RIGHT EIGENVECTORS CHANGED TO DONAT ET AL.'S 98 EXPRESSIONS
!     13/12/01
  DO i= M, N
    SCRCH2(i) = 1.0 - VELL(i)*VELL(i)

    A1L(i) =  SCRCH2(i) / (1.0 - VELL(i)*LAMB1L(i))
    A5L(i) =  SCRCH2(i) / (1.0 - VELL(i)*LAMB5L(i))

    SCRCH3(i) = HL(i)*WL(i)
    SCRCH4(i) = SCRCH3(i)*A1L(i)

    R11(i) = 1.0
    R21(i) = SCRCH4(i)*LAMB1L(i)
    R31(i) = SCRCH3(i)*VELTL(i)
    R41(i) = SCRCH3(i)*VELTTL(i)
    R51(i) = SCRCH4(i)-1.0

    SCRCH4(i) = SCRCH3(i)*A5L(i)

    R15(i) = 1.0
    R25(i) = SCRCH4(i)*LAMB5L(i)
    R35(i) = R31(i)
    R45(i) = R41(i)
    R55(i) = SCRCH4(i)-1.0
  ENDDO


!     7/12/01
!     COMPUTATION OF FACTORS WHICH ARE USED IN THE SUMATION OF LEFT
!     EIGENVECTORS TIMES THE CHARACTERISTIC Variables (S1L, S5L)
  DO i= M, N
    KAPPAL = DPDUL(i) / (DPDUL(i)-RHOL(i)*CSL(i)*CSL(i))
    SCRCH13 = 1.0 - KAPPAL
    DELTAL = (HL(i)**3)*WL(i)*SCRCH13*SCRCH2(i)*(A1L(i)*LAMB1L(i)-A5L(i)*LAMB5L(i))

    SCRCH6  = VELL(i)-LAMB5L(i)
    SCRCH7  = VELL(i)-LAMB1L(i)
    SCRCH8  = RHOL(i)*HL(i)*SCRCH13
    SCRCH9  = PL(i)*(2.0*KAPPAL-1.0)
    SCRCH10 = PL(i)*SCRCH13
    SCRCH11 = HL(i)**2/DELTAL
    SCRCH12 = WL(i)**2

    S1L =  SCRCH11*( SCRCH12*A5L(i)*SCRCH6*(SCRCH8 +SCRCH9)  +          &
                    SCRCH10*(VELL(i) +SCRCH6/(1.0-VELL(i)*LAMB5L(i))) )

    S5L = -SCRCH11* ( SCRCH12*A1L(i)*SCRCH7*(SCRCH8 + SCRCH9)  +    &
                      SCRCH10*(VELL(i) +SCRCH7/(1.0-VELL(i)*LAMB1L(i))) )

    Q1L(i) = PSI2(i)*RL(i) + (PSI12(i)*R11(i)*S1L + PSI52(i)*R15(i)*S5L)

    Q2L(i) = PSI2(i)*ML(i) + (PSI12(i)*R21(i)*S1L + PSI52(i)*R25(i)*S5L)

    Q3L(i) = PSI2(i)*MTL(i) + (PSI12(i)*R31(i)*S1L + PSI52(i)*R35(i)*S5L)

    Q4L(i) = PSI2(i)*MTTL(i) + (PSI12(i)*R41(i)*S1L + PSI52(i)*R45(i)*S5L)

    Q5L(i) = PSI2(i)*EL(i) + (PSI12(i)*R51(i)*S1L + PSI52(i)*R55(i)*S5L)
  ENDDO


!     RIGHT NUMERICAL VISCOSITY
!     RIGHT EIGENVECTORS
!     RIJ STANDS FOR THE COMPONENT i OF THE RIGHT EIGENVECTOR J
  DO i= M, N
    SCRCH2(i) = 1.0 - VELR(i)*VELR(i)

    A1R(i) =  SCRCH2(i) / (1.0 - VELR(i)*LAMB1R(i))
    A5R(i) =  SCRCH2(i) / (1.0 - VELR(i)*LAMB5R(i))

    SCRCH3(i) = HR(i)*WR(i)
    SCRCH4(i) = SCRCH3(i)*A1R(i)

    R11(i) = 1.0
    R21(i) = SCRCH4(i)*LAMB1R(i)
    R31(i) = SCRCH3(i)*VELTR(i)
    R41(i) = SCRCH3(i)*VELTTR(i)
    R51(i) = SCRCH4(i)-1.0

    SCRCH4(i) = SCRCH3(i)*A5R(i)

    R15(i) = 1.0
    R25(i) = SCRCH4(i)*LAMB5R(i)
    R35(i) = R31(i)
    R45(i) = R41(i)
    R55(i) = SCRCH4(i)-1.0
  ENDDO

!     7/12/01
!     COMPUTATION OF FACTORS WHICH ARE USED IN THE SUMATION OF LEFT
!     EIGENVECTORS TIMES THE CHARACTERISTIC Variables (S1R, S4R)
  DO i= M, N
    KAPPAR = DPDUR(i) / (DPDUR(i)-RHOR(i)*CSR(i)*CSR(i))

    SCRCH13 = 1.0 - KAPPAR

    DELTAR = (HR(i)**3)*WR(i)*SCRCH13*SCRCH2(i)*  &
             (A1R(i)*LAMB1R(i)-A5R(i)*LAMB5R(i))


    SCRCH6 = VELR(i)-LAMB5R(i)
    SCRCH7 = VELR(i)-LAMB1R(i)
    SCRCH8 = RHOR(i)*HR(i)*SCRCH13
    SCRCH9 = PR(i)*(2.0*KAPPAR-1.0)
    SCRCH10 = PR(i)*SCRCH13
    SCRCH11 = HR(i)**2/DELTAR
    SCRCH12 = WR(i)**2

    S1R =  SCRCH11* ( SCRCH12*A5R(i)*SCRCH6*(SCRCH8 +SCRCH9)  +  &
                      SCRCH10* (VELR(i) +SCRCH6/(1.0-VELR(i)*LAMB5R(i))) )

    S5R = -  SCRCH11*( SCRCH12*A1R(i)*SCRCH7*(SCRCH8 + SCRCH9)  +  &
                       SCRCH10*(VELR(i) +SCRCH7/(1.0-VELR(i)*LAMB1R(i))) )

    Q1R(i) = PSI2(i)*RR(i)   + (PSI12(i)*R11(i)*S1R + PSI52(i)*R15(i)*S5R)
    Q2R(i) = PSI2(i)*MR(i)   + (PSI12(i)*R21(i)*S1R + PSI52(i)*R25(i)*S5R)
    Q3R(i) = PSI2(i)*MTR(i)  + (PSI12(i)*R31(i)*S1R + PSI52(i)*R35(i)*S5R)
    Q4R(i) = PSI2(i)*MTTR(i) + (PSI12(i)*R41(i)*S1R + PSI52(i)*R45(i)*S5R)
    Q5R(i) = PSI2(i)*ER(i)   + (PSI12(i)*R51(i)*S1R + PSI52(i)*R55(i)*S5R)

!       TOTAL NUMERICAL VISCOSITY
    Q1(i)= Q1R(i)-Q1L(i)
    Q2(i)= Q2R(i)-Q2L(i)
    Q3(i)= Q3R(i)-Q3L(i)
    Q4(i)= Q4R(i)-Q4L(i)
    Q5(i)= Q5R(i)-Q5L(i)
  ENDDO


  DO i= M, N
    aux1 = RL(i) *VELL(i) +  RR(i)*VELR(i) - Q1(i)
    aux2 = RR(i) *VELR(i) +  RL(i)*VELL(i) - Q1(i)
    FR(i)= 0.25*(aux1+aux2)

    aux3 = ML(i) *VELL(i)+PL(i)
    aux4 = MR(i) *VELR(i)+PR(i)
    aux1 = aux3+aux4-Q2(i)
    aux2 = aux4+aux3-Q2(i)
    FM(i)= 0.25*(aux1+aux2)

    aux1 = MTL(i)*VELL(i)      +MTR(i)*VELR(i)      - Q3(i)
    aux2 = MTR(i)*VELR(i)      +MTL(i)*VELL(i)      - Q3(i)
    FMT(i)= 0.25*(aux1+aux2)

    aux1 = MTTL(i)*VELL(i)     +MTTR(i)*VELR(i)     - Q4(i)
    aux2 = MTTR(i)*VELR(i)     +MTTL(i)*VELL(i)     - Q4(i)
    FMTT(i)= 0.25*(aux1+aux2)

    aux3 = ML(i)-RL(i)*VELL(i)
    aux4 = MR(i)-RR(i)*VELR(i)
    aux1 = aux3+aux4  -Q5(i)
    aux2 = aux4+aux3  -Q5(i)
    FE(i)= 0.25*(aux1+aux2)
  ENDDO
!----------------------------------------------------------------------------END
END subroutine riem3dp
