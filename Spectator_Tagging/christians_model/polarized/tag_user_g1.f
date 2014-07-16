*DECK TAGGN
      SUBROUTINE TAGGN(G1, X, QQ, IPN)
C
C     PACKAGE TAG -- DEUTERON DIS WITH SPECTATOR TAGGING
C     AUTHOR C. WEISS (WEISS.AT.JLAB.ORG)
C
C     USER-DEFINED ROUTINE
C     NUCLEON STRUCTURE FUNCTION G1
C
C     INPUT:
C     X     X VARIABLE FOR NUCLEON
C     QQ    Q2 (GEV**2)
C     IPN   SWITCH: 1  PROTON, 2  NEUTRON
C
C     OUTPUT:
C     G1    NUCLEON STRUCTURE FUNCTION G1
C
C     V2 30JUN14
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
C     ...GRSV2000, STANDARD SCENARIO, LO
C
      ISET = 3
      CALL PARPOL(ISET, X, QQ, U, D, UB, DB, ST, GL, G1P, G1N)
C
      IF      (IPN.EQ.1) THEN
         G1 = G1P
      ELSE IF (IPN.EQ.2) THEN
         G1 = G1N
      ENDIF
CC
CC     ...CHECK: EXPLICIT CALCULATION OF STRUCTURE FUNCTION
CC        FROM QUARK DENSITIES (LO ONLY)
CC
CC     ...QUARK CHARGES
CC
C      EU   =  0.6666666666666D0
C      ED   = -0.3333333333333D0
C      ES   = -0.3333333333333D0
CC
C      IF      (IPN.EQ.1) THEN
C         G1 = (EU**2*(U + UB)
C     *      +  ED**2*(D + DB)
C     *      +  ES**2*(  2*ST))/X/2.D0
C      ELSE IF (IPN.EQ.2) THEN
C         G1 = (EU**2*(D + DB)
C     *      +  ED**2*(U + UB)
C     *      +  ES**2*(  2*ST))/X/2.D0
C      ENDIF
CC
      END
C
C---------------------------------------------------------------------
C
*DECK GRSV00
      BLOCK DATA GRSV00
C
C     INITIALIZATION OF GRSV2000 PDF PARAMETRIZATION
C
C     NOTE: IF MULTIPLE PDF SETS ARE USED IN THE SAME RUN,
C     THE INITIALIZATION MUST BE RESET BY THE CALLING PROGRAM.
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      COMMON /INTINI/  IINI
      DATA IINI /0/
C
      END
C
C---------------------------------------------------------------------
C
