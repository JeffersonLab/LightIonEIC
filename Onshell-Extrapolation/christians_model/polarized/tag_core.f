C     PACKAGE TAG -- DEUTERON DIS WITH SPECTATOR TAGGING
C     AUTHOR C. WEISS (WEISS.AT.JLAB.ORG)
C
C     CORE ROUTINES
C     TAGXP -- DEUTERON TAGGED CROSS SECTION, POLARIZED
C     TAGGD -- DEUTERON TAGGED STRUCTURE FUNCTIONS, POLARIZED
C     TAGX  -- DEUTERON TAGGED CROSS SECTION, UNPOLARIZED
C     TAGFD -- DEUTERON TAGGED STRUCTURE FUNCTIONS, UNPOLARIZED
C     TAGSP -- DEUTERON SPECTRAL FUNCTION
C
C     USER-DEFINED ROUTINES CALLED
C     TAGGN -- NUCLEON STRUCTURE FUNCTION, POLARIZED
C     TAGFN -- NUCLEON STRUCTURE FUNCTION, UNPOLARIZED
C     TAGWF -- DEUTERON NON-RELATIVISTIC WAVE FUNCTION
C
C     ADDITIONAL USER-DEFINED ROUTINES (OPTIONAL)
C     TAGRES -- RESIDUE OF SPECTRAL FUNCTION AT TPRIME = 0
C
C     V1 13MAY14
C     V2 30JUN14   ADDED POLARIZED E-D (TAGXP, TAGGN)
C
C---------------------------------------------------------------------
C
*DECK TAGXP
      SUBROUTINE TAGXP(FSIG, Y, EPS, D1, SED, X, QQ, ALR, PTR, 
     *   IPN, IPOL)
C
C     DEUTERON TAGGED ELECTROPRODUCTION CROSS SECTION
C
C     OUTPUT:
C     FSIG  TAGGED ELECTROPRODUCTION CROSS SECTION
C           DSIGMA = FSIG *DX *DQ**2 *(D**3 PR/ER)
C
C           DIMENSION OF DSIGMA IS NANOBARN
C               ''       FSIG      NANOBARN/GEV**4
C
C     Y     ELECTRON FRACTIONAL ENERGY LOSS
C     EPS   VIRTUAL PHOTON POLARIZATION PARAMETER
C     D1    DEPOLARIZATION FACTOR (COEFFICIENT OF G1)
C
C     INPUT:
C     SED   ELECTRON-DEUTERON S-INVARIANT (SQUARED CM ENERGY)
C     X     NUCLEON SCALING VARIABLE X = 2*XD
C     QQ    4-MOMENTUM TRANSFER Q**2
C     ALR   RECOIL MOMENTUM, LC FRACTION 
C     PTR   RECOIL MOMENTUM, TRANSVERSE
C
C     IPN   TAGGED NUCLEON TYPE
C           IPN = 1   PROTON  ACTIVE, NEUTRON TAGGED
C                 2   NEUTRON ACTIVE, PROTON  TAGGED
C
C     IPOL   RELATIVE ELECTRON-DEUTERON POLARIZATION
C            IPOL  =  (2*HELICITY-E) * HELICITY-D  =  (1, -1)
C
C            IPOL = 1  CROSS SECTION SIGMA(++) = SIGMA(--)
C                  -1  CROSS SECTION SIGMA(+-) = SIGMA(-+)
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      PARAMETER (PI = 0.31415 92653 58979 D+01)
C
C     ...FINE STRUCTURE CONSTANT
C
      PARAMETER (ALEM = 1.D0/137.D0)
C
C     ...NUCLEON AND DEUTERON MASS
C        UN   AVERAGE NUCLEON MASS (M(PROTON) + M(NEUTRON))/2   
C        ED   DEUTERON BINDING ENERGY
C        DEUTERON MASS CALCULATED AS  2*UN - ED
C
      PARAMETER (UN = 0.939D0, ED = 0.00222, UD = 2*UN - ED)
C
C     ...CONVERSION FACTOR FOR CROSS SECTION UNITS
C        1/GEV**2  =  0.398 *10**6  NANOBARN
C
      PARAMETER (CONVNB = 0.389D6)
C
      XD   = X/2
      Y    = QQ/XD/(SED - UD**2)
C
      CORR = (Y*XD*UD)**2/QQ 
      EPS  = (1 - Y - CORR)/(1 - Y + Y**2/2 + CORR) 
C
C     ...FACTOR FOR ELECTROPRODUCTION CROSS SECTION
C        EXTRA FACTOR 1/2 FROM DIFFERENTIAL DXD = DX/2 
C        CONVERSION GEV**(-2) TO NANOBARN
C
      FAC  = 2*PI *ALEM**2 *Y**2 /QQ**2 /(1 - EPS)
      FACD = FAC /2.D0 *CONVNB
C
C     ...COEFFICIENT OF POLARIZED STRUCTURE FUNCTION G1
C        (DEPOLARIZATION FACTOR)
C
      GAMS = 4 *XD**2 *UD**2/ QQ
      D1 = (1 - EPS)/Y*(2 - Y*(1 + GAMS*Y/2))
C
C     ...STRUCTURE FUNCTIONS (UNPOLARIZED, POLARIZED)
C
      CALL TAGFD(F2, X, QQ, ALR, PTR, IPN)
C
      CALL TAGGD(G1, X, QQ, ALR, PTR, IPN)
C
C     ...CROSS SECTION
C
      FSIG = FACD*(F2/XD - IPOL*2*D1*G1)
C
      END
C
C---------------------------------------------------------------------
C
*DECK TAGGD
      SUBROUTINE TAGGD(G1, X, QQ, ALR, PTR, IPN)
C
C     DEUTERON TAGGED STRUCTURE FUNCTION IN IMPULSE APPROXIMATION
C     POLARIZED STRUCTURE FUNCTION G1
C
C     G1    TAGGED STRUCTURE FUNCTION G1
C
C     X     NUCLEON SCALING VARIABLE X = 2*XD
C     QQ    4-MOMENTUM TRANSFER Q**2
C     ALR   RECOIL MOMENTUM, LC FRACTION 
C     PTR   RECOIL MOMENTUM, TRANSVERSE
C
C     IPN = 1   PROTON  ACTIVE, NEUTRON TAGGED
C           2   NEUTRON ACTIVE, PROTON  TAGGED
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
C     ...SPECTRAL FUNCTION
C
      CALL TAGSP(S, ALR, PTR)
C
C     ...STRUCTURE FUNCTION
C
      XEFF = X/(2 - ALR)
C
      CALL TAGGN(GN1, XEFF, QQ, IPN)
C
      G1 = 2.D0/(2 - ALR) *S *GN1
C
      END
C
C---------------------------------------------------------------------
C
*DECK TAGX
      SUBROUTINE TAGX(FSIG, SED, X, QQ, ALR, PTR, IPN)
C
C     DEUTERON TAGGED ELECTROPRODUCTION CROSS SECTION
C
C     FSIG  TAGGED ELECTROPRODUCTION CROSS SECTION
C           DSIGMA = FSIG *DX *DQ**2 *(D**3 PR/ER)
C
C           DIMENSION OF DSIGMA IS NANOBARN
C               ''       FSIG      NANOBARN/GEV**4
C
C     SED   ELECTRON-DEUTERON S-INVARIANT (SQUARED CM ENERGY)
C     X     NUCLEON SCALING VARIABLE X = 2*XD
C     QQ    4-MOMENTUM TRANSFER Q**2
C     ALR   RECOIL MOMENTUM, LC FRACTION 
C     PTR   RECOIL MOMENTUM, TRANSVERSE
C
C     IPN = 1   PROTON  ACTIVE, NEUTRON TAGGED
C           2   NEUTRON ACTIVE, PROTON  TAGGED
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      PARAMETER (PI = 0.31415 92653 58979 D+01)
C
C     ...FINE STRUCTURE CONSTANT
C
      PARAMETER (ALEM = 1.D0/137.D0)
C
C     ...NUCLEON AND DEUTERON MASS
C        UN   AVERAGE NUCLEON MASS (M(PROTON) + M(NEUTRON))/2   
C        ED   DEUTERON BINDING ENERGY
C        DEUTERON MASS CALCULATED AS  2*UN - ED
C
      PARAMETER (UN = 0.939D0, ED = 0.00222, UD = 2*UN - ED)
C
C     ...CONVERSION FACTOR FOR CROSS SECTION UNITS
C        1/GEV**2  =  0.398 *10**6  NANOBARN
C
      PARAMETER (CONVNB = 0.389D6)
C
      XD   = X/2
C
      Y    = QQ/XD/(SED - UD**2)
C
      CORR = (Y*XD*UD)**2/QQ 
      EPS  = (1 - Y - CORR)/(1 - Y + Y**2/2 + CORR) 
C
C     ...FACTOR FOR ELECTROPRODUCTION CROSS SECTION
C        EXTRA FACTOR 1/2 FROM DIFFERENTIAL DXD = DX/2 
C        CONVERSION GEV**(-2) TO NANOBARN
C
      FAC  = 2*PI *ALEM**2 *Y**2 /QQ**2 /(1 - EPS)
      FACD = FAC /2.D0 *CONVNB
C
C     ...STRUCTURE FUNCTIONS
C
      CALL TAGFD(F2, X, QQ, ALR, PTR, IPN)
C
C     ...CROSS SECTION
C
      FSIG = FACD*(F2/XD)
C
      END
C
C---------------------------------------------------------------------
C
*DECK TAGFD
      SUBROUTINE TAGFD(F2, X, QQ, ALR, PTR, IPN)
C
C     DEUTERON TAGGED STRUCTURE FUNCTION IN IMPULSE APPROXIMATION
C
C     F2    TAGGED STRUCTURE FUNCTION
C
C     X     NUCLEON SCALING VARIABLE X = 2*XD
C     QQ    4-MOMENTUM TRANSFER Q**2
C     ALR   RECOIL MOMENTUM, LC FRACTION 
C     PTR   RECOIL MOMENTUM, TRANSVERSE
C
C     IPN = 1   PROTON  ACTIVE, NEUTRON TAGGED
C           2   NEUTRON ACTIVE, PROTON  TAGGED
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
C     ...SPECTRAL FUNCTION
C
      CALL TAGSP(S, ALR, PTR)
C
C     ...STRUCTURE FUNCTION
C
      XEFF = X/(2 - ALR)
C
      CALL TAGFN(FN2, XEFF, QQ, IPN)
C
      F2 = S*FN2
C
      END
C
C---------------------------------------------------------------------
C
*DECK TAGSP
      SUBROUTINE TAGSP(S, ALR, PTR)
C
C     DEUTERON LIGHT-CONE SPECTRAL FUNCTION
C     NORMALIZED SUCH THAT  INTEGRAL (D ALR/ALR) (D**2 PTR) S = 1
C
C     ALR   RECOIL MOMENTUM, LC FRACTION 
C     PTR   RECOIL MOMENTUM, TRANSVERSE
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      PARAMETER (UN = 0.939D0)
C
C     ...EFFECTIVE CM MOMENTUM
C
      P = SQRT((PTR**2 + UN**2)/ALR/(2 - ALR) - UN**2)
C
C     ...SPECTRAL FUNCTION
C
      CALL TAGWF(PSI, P)
C
      S = SQRT(P**2 + UN**2) *PSI**2 /(2 - ALR)
C
      END
C
C---------------------------------------------------------------------
C
