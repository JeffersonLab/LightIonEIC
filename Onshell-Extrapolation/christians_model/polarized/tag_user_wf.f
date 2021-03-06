*DECK TAGWF
      SUBROUTINE TAGWF(PSI, P)
C
C     PACKAGE TAG -- DEUTERON DIS WITH SPECTATOR TAGGING
C     AUTHOR C. WEISS (WEISS.AT.JLAB.ORG)
C
C     USER-DEFINED ROUTINE
C     DEUTERON NON-RELATIVISTIC WAVE FUNCTION PSI(P)
C     NORMALIZED SUCH THAT   INTEGRAL (D**3 P) PSI(P)**2 = 1
C
C     INPUT:
C     P    MODULUS OF 3-MOMENTUM (GEV)
C
C     OUTPUT:
C     PSI  DEUTERON NON-RELATIVISTIC WAVE FUNCTION PSI(P)
C
C     ALL DIMENSIONFUL QUANTITIES IN GEV UNITS
C
C     V1 13MAY14
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      PARAMETER (PI = 0.31415 92653 58979 D+01)
C
C     ...HULTHEN WAVE FUNCTION
C
      A = 0.045647
      B = 0.2719
C
      C = PI**2*(A - B)**2/A/B/(A + B)
C
      PSI = (  1.D0/(P**2 + A**2)
     *       - 1.D0/(P**2 + B**2))/SQRT(C)
C
      END
C
C---------------------------------------------------------------------
C
