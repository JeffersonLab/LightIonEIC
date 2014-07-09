*DECK EVTP
      PROGRAM EVTP
C
C     PACKAGE TAG -- DEUTERON DIS WITH SPECTATOR TAGGING
C     AUTHOR C. WEISS (WEISS.AT.JLAB.ORG)
C
C     APPLICATION PROGRAM
C     TABULATE EVENT NUMBERS IN MEASUREMENT OF TPRIME-DISTRIBUTION
C
C     V1 13MAY14
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      PARAMETER (PI = 0.31415 92653 58979 D+01)
      PARAMETER (UN = 0.939D0, ED = 0.00222, UD = 2*UN - ED)
C

C     Create 3 different output files for the 3 different values used of F2N
      OPEN(101, FILE = 'EVTP.IN',   STATUS = 'OLD')
      OPEN(1, FILE = 'EVTP1.OUT',  STATUS = 'OLD')
      OPEN(2, FILE = 'EVTP2.OUT', STATUS = 'OLD')
      OPEN(3, FILE = 'EVTP4.OUT', STATUS = 'OLD')
      REWIND(101)
C
C     ...READ INPUT PARAMETERS (DOCUMENTATION SEE FILE)
C
      READ(101, *) ULUMI
      READ(101, *) IPN
      READ(101, *) SED
      READ(101, *) X
      READ(101, *) XBIN
      READ(101, *) QQ
      READ(101, *) QQBIN
      READ(101, *) ALR
      READ(101, *) ALRBIN
      READ(101, *) TPBIN
      READ(101, *) NBIN
C
C     ...KINEMATIC LIMIT IN TPRIME, ABSOLUTE AND FOR GIVEN ALPHAR
C
      CALL TPMIN(TP0, 1.D0, 0)
      CALL TPMIN(TP1,  ALR, 1)
C
C     ...RESIDUE OF SPECTRAL FUNCTION
C
      CALL TAGRES(RES, ALR)
C
C     Run the simulation 3 times while scaling F2N by different amounts

      DO 88 I = 1, 2

C        Scaling F2N by 0.9, 1.0 then 1.1
         SCLRTO = 0.1
         SCLING = 1.0 !+ (3.0 - I)*SCLRTO

C        ...FREE NUCLEON STRUCTURE FUNCTION (INPUT MODEL)
C
          CALL TAGFN(F2N, X, QQ, IPN, SCLING)
C
C         ...WRITE OUTPUT FILE HEADER
C

          
          WRITE(I, 91000) 'OUTPUT EVTP'
          WRITE(I, 91000) 
     * 'EVENT NRS IN MEASUREMENT OF TPRIME DISTRIBUTION'
          WRITE(I, 91003) ULUMI,  'ULUMI (1/NB)'
          WRITE(I, 91002) IPN,    'IPN'
          WRITE(I, 91001) SED,    'SED'
          WRITE(I, 91001) X,      'X'
          WRITE(I, 91001) XBIN,   'XBIN'
          WRITE(I, 91001) QQ,     'QQ'
          WRITE(I, 91001) QQBIN,  'QQBIN'
          WRITE(I, 91001) ALR,    'ALR'
          WRITE(I, 91001) ALRBIN, 'ALRBIN'
          WRITE(I, 91001) TPBIN,  'TPBIN'
          WRITE(I, 91002) NBIN,   'NBIN'
          WRITE(I, 91001) TP0,    'TP0 -- TMPMIN(ABS)'
          WRITE(I, 91001) TP1,    'TP1 -- TMPMIN(ALR)'
          WRITE(I, 91001) RES,    'RES'
          WRITE(I, 91001) F2N,    'F2N'
          WRITE(I, 91000)
          WRITE(I, 91000) 'COL1  TP     TPRIME (GEV**2)'
          WRITE(I, 91000) 
     * 'COL2  PR2    SQUARED RECOIL 3-MOM RF (GEV**2)'
          WRITE(I, 91000) 
     * 'COL3  PTR    RECOIL TRANSVERSE MOMENTUM (GEV)'
          WRITE(I, 91000) 
     * 'COL4  FSIG   TAGGED DIFF CROSS SECN (NB/GEV**4)'
          WRITE(I, 91000) 'COL5  UNUM   EVENT NUMBER IN BIN'
          WRITE(I, 91000) 'COL6  F2/SPOL   TAGGED F2/POLE FACTOR'
91000     FORMAT('#*',                A)
91001     FORMAT('#', T3, F12.6, T20, A)
91002     FORMAT('#', T3, I12,   T20, A)
91003     FORMAT('#', T3, E12.5, T20, A)
C
          DO 10001 IBIN = 1, NBIN
C         PRINT *, F2
C
C    ...TPRIME VALUES AT BIN END AND CENTER
C
             TPLO = -(IBIN - 1)*TPBIN
             TP   = -(IBIN)*TPBIN
C
C             IF (TPLO.GT.TP1) THEN
C
C     ...BIN PARTLY OUTSIDE PHYSICAL REGION: SKIP
C
C               PR2  = 0.D0
C               PTR  = 0.D0
C               FSIG = 0.D0
C               UNUM = 0.D0
C               S    = 0.D0
C
c        ELSE
C
C    ...BIN INSIDE PHYSICAL REGION: CALCULATE CROSS SECTION
C       AND NUMBER OF EVENTS
C
C    ...CALCULATE PTR
C
                PR2 = (TP0 - TP)/2
                ER  = SQRT(PR2 + UN**2)
                PRZ = ALR*UD/2 - ER
                PTR = SQRT(PR2 - PRZ**2)
C
                CALL TAGX(FSIG, SED, X, QQ, ALR, PTR, IPN, SCLING)
C
                PHASR = PI*UD/4./ER*TPBIN*ALRBIN
C     
                UNUM = ULUMI *XBIN*QQBIN*FSIG *PHASR
C
C    ...TAGGED STRUCTURE FUNCTION
C
                CALL TAGFD(F2, X, QQ, ALR, PTR, IPN, SCLING)
C
C     ...SPECTRAL FUNCTION AND POLE PART
C
                CALL TAGSP(S, ALR, PTR)
C
C            ENDIF
C
             SPOL = RES/(-TP)**2 

             IF(I.EQ.1) THEN
                DENOM = 1.0
             ELSE IF (I.EQ.2) THEN
                DENOM = SPOL
             ENDIF
             
C            ANORM = FSIG/SPOL
C            PRINT *, SPOL, ANORM, FSIG/ANORM
C
C        PRINT 90000,    TP, PR2, PTR, FSIG, UNUM, F2/SPOL
             WRITE(I, 90000) TP, PR2, PTR, FSIG/DENOM, UNUM, F2/SPOL
90000        FORMAT(1(1X, F7.4), 6(1X, E10.4))
C
10001 CONTINUE
88    CONTINUE

      END
C
C---------------------------------------------------------------------
C
*DECK TPMIN
      SUBROUTINE TPMIN(TP, ALR, ILIM)
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
C---------------------------------------------------------------------
C
