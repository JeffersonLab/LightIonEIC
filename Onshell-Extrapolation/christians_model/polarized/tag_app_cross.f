*DECK CROSS
      PROGRAM CROSS
C
C     PACKAGE TAG -- DEUTERON DIS WITH SPECTATOR TAGGING
C     AUTHOR C. WEISS (WEISS.AT.JLAB.ORG)
C
C     APPLICATION PROGRAM
C
C     TABULATE SPIN-DEPENDENT CROSS SECTION IN X
C
C     V2 30JUN14
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      PARAMETER (UN = 0.939D0, ED = 0.00222, UD = 2*UN - ED)
C
      OPEN(1, FILE = 'CROSS.IN',   STATUS = 'OLD')
      OPEN(2, FILE = 'CROSS.OUT',  STATUS = 'OLD')
      REWIND(1)
      REWIND(2)
C
C     ...READ INPUT PARAMETERS (DOCUMENTATION SEE FILE)
C
      READ(1, *) SED
      READ(1, *) QQ
      READ(1, *) ALR
      READ(1, *) X
      READ(1, *) TPMIN
      READ(1, *) TPMAX
      READ(1, *) NSTEP

C
C     SET UP TP BINNING
C
      TPBIN = (TPMAX - TPMIN)/NSTEP
      CALL TPLOW(TP0, 1.D0, 0)
C
C
C     ...WRITE OUTPUT FILE HEADER
C
      WRITE(2, 91000) 'OUTPUT CROSS'
      WRITE(2, 91000) 'TABULATE SPIN-DEP CROSS SECTION IN X'
      WRITE(2, 91000)
      WRITE(2, 91001) SED, 'SED'
      WRITE(2, 91001) QQ,  'QQ'
      WRITE(2, 91001) ALR, 'ALR'
      WRITE(2, 91001) X, 'X'
      WRITE(2, 91001) TPMIN, 'TPMIN'
      WRITE(2, 91001) TPMAX, 'TPMAX'
      WRITE(2, 91000)
      WRITE(2, 91000) 'COL1 TP'
C     WRITE(2, 91000) 'COL2  Y'
C     WRITE(2, 91000) 'COL3  EPS (VIRTUAL PHOTON POLARIZATION)'
C     WRITE(2, 91000) 'COL4  D1 (DEPOLARIZATION FACTOR)'
C     WRITE(2, 91000) 'COL5  FSIGMA, PROTON  ACTIVE, IPOL =  1'
C     WRITE(2, 91000) 'COL6  FSIGMA, PROTON  ACTIVE, IPOL = -1'
      WRITE(2, 91000) 'COL2  FSIGMA, NEUTRON ACTIVE, IPOL =  1'
      WRITE(2, 91000) 'COL3  FSIGMA, NEUTRON ACTIVE, IPOL = -1'
      WRITE(2, 91000) 'COL4  ASYM, CALCULATED ASYMETRY'
91000 FORMAT('#*',                A)
91001 FORMAT('#', T3, F12.6, T20, A)
91002 FORMAT('#', T3, I12,   T20, A)
C
      DO 10000 I = 0, NSTEP

            TP = TPMIN + I*TPBIN 

            PR2 = (TP + TP0)/2
            ER  = SQRT(PR2 + UN**2)
            PRZ = ALR*UD/2 - ER
            PTR = SQRT(PR2 - PRZ**2)

C
C
C        ...CROSS SECTION
C
         IPN = 1
         CALL TAGXP(FPP, Y, EPS, D1, SED, X, QQ, ALR, PTR, IPN,  1)
         CALL TAGXP(FPM, Y, EPS, D1, SED, X, QQ, ALR, PTR, IPN, -1)
         IPN = 2
         CALL TAGXP(FNP, Y, EPS, D1, SED, X, QQ, ALR, PTR, IPN,  1)
         CALL TAGXP(FNM, Y, EPS, D1, SED, X, QQ, ALR, PTR, IPN, -1)

C
C        PRINT 90000,    X, FPP, FPM, FNP, FNM
C         WRITE(2, 90000) X, FPP, FPM, FNP, FNM
         ASYM = (FNP - FNM) / (FNP + FNM)
C        WRITE(2, 90000) TP, Y, EPS, D1, FPP, FPM, FNP, FNM, ASYM
         WRITE(2, 80000) TP, FNP, FNM, ASYM
90000    FORMAT(4(1X, F8.5), 5(1X, E10.5))
80000    FORMAT(1(1X, F8.5), 3(1X, E10.5))
10000    CONTINUE
C
C
      END
C
C---------------------------------------------------------------------
C

*DECK TPMIN
      SUBROUTINE TPLOW(TP, ALR, ILIM)
C
C     KINEMATIC LIMIT OF TPRIME
C
C     ILIM = 0:  ABSOLUTE LIMIT, CORRESPONDS TO ALPHAR = 2*MN/MD
C            1:  LIMIT FOR FIXED ALPHAR
C
C     TP IS THE TRUE (NEGATIVE) VALUE OF TPRIME
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      PARAMETER (UN = 0.939D0, ED = 0.00222, UD = 2*UN - ED)
C
      TP = -2*ED*UN + ED**2/2.
      IF (ILIM.EQ.1) TP = TP - 2*(ALR*UD/4. - UN**2/ALR/UD)**2
C
      END
C
