*DECK TAGFN
      SUBROUTINE TAGFN(F2, X, QQ, IPN, SCLING)
C
C     PACKAGE TAG -- DEUTERON DIS WITH SPECTATOR TAGGING
C     AUTHOR C. WEISS (WEISS.AT.JLAB.ORG)
C
C     USER-DEFINED ROUTINE
C     NUCLEON STRUCTURE FUNCTION F2
C
C     INPUT:
C     X     X VARIABLE FOR NUCLEON
C     QQ    Q2 (GEV**2)
C     IPN   SWITCH: 1  PROTON, 2  NEUTRON
C     SCLING Scaling factor to vary F2N. 1 if no scaling desired.
C
C     OUTPUT:
C     F2    NUCLEON STRUCTURE FUNCTION F2
C
C     V1 13MAY14
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
C     ...QUARK CHARGES
C
      EU   =  0.6666666666666D0
      ED   = -0.3333333333333D0
      ES   = -0.3333333333333D0
C
C     ...GRV98 PDF PARAMETRIZATION
C
      ISET = 1
      CALL GRV98PA(ISET, X, QQ, UV, DV, US, DS, SS, GL)
C
      IF (IPN.EQ.1)      THEN
         F2 = SCLING * 
     *       (EU**2*(UV + 2*US)
     *      + ED**2*(DV + 2*DS)
     *      + ES**2*(     2*SS))
      ELSE IF (IPN.EQ.2) THEN
         F2 = SCLING * 
     *       (EU**2*(DV + 2*DS)
     *      + ED**2*(UV + 2*US)
     *      + ES**2*(     2*SS))
      ENDIF
C
C     PRINT *,"F2 returned as ", F2
      END
C
C---------------------------------------------------------------------
C
*DECK GRV98
      BLOCK DATA GRV98
C
C     INITIALIZATION OF GRV98 PDF PARAMETRIZATION
C
C     NOTE: IF MULTIPLE PDF SETS ARE USED IN THE SAME RUN,
C     THE INITIALIZATION MUST BE RESET BY THE CALLING PROGRAM.
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      COMMON /INTINIP/ IINIP
      COMMON /INTINIF/ IINIF
      DATA IINIP /0/
      DATA IINIF /0/
C
      END
C
C---------------------------------------------------------------------
C
