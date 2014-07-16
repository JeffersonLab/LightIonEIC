*DECK G1X
      PROGRAM G1X
C
C     PACKAGE TAG -- DEUTERON DIS WITH SPECTATOR TAGGING
C     AUTHOR C. WEISS (WEISS.AT.JLAB.ORG)
C
C     APPLICATION PROGRAM
C
C     TABULATE NUCLEON STRUCTURE FUNCTION F2 IN X
C
C     V2 30JUN14
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      OPEN(1, FILE = 'G1X.IN',   STATUS = 'OLD')
      OPEN(2, FILE = 'G1X.OUT',  STATUS = 'OLD')
      REWIND(1)
      REWIND(2)
C
C     ...READ INPUT PARAMETERS (DOCUMENTATION SEE FILE)
C
      READ(1, *) QQ
      READ(1, *) DL1
      READ(1, *) DL2
      READ(1, *) NSTEP
C
C     ...WRITE OUTPUT FILE HEADER
C
      WRITE(2, 91000) 'OUTPUT G1X'
      WRITE(2, 91000) 'TABULATE STRUCTURE FUNCTION G1 IN X'
      WRITE(2, 91000)
      WRITE(2, 91001) QQ, 'QQ'
      WRITE(2, 91000)
      WRITE(2, 91000) 'COL1  X'
      WRITE(2, 91000) 'COL2  G1, PROTON'
      WRITE(2, 91000) 'COL3  G2, NEUTRON'
91000 FORMAT('#*',                A)
91001 FORMAT('#', T3, F12.6, T20, A)
91002 FORMAT('#', T3, I12,   T20, A)
C
      DO 10000 I = 0, NSTEP
C
         DL = DL1 + (I*(DL2 - DL1))/NSTEP
         X = 10**DL
C
C        ...STRUCTURE FUNCTION G1
C
         CALL TAGGN(G1P, X, QQ, 1)
         CALL TAGGN(G1N, X, QQ, 2)
C
C        PRINT 90000,    X, G1P, G1N
         WRITE(2, 90000) X, G1P, G1N
90000    FORMAT(F10.5, 4(1X, E12.6))
C
10000 CONTINUE
C
      END
C
C---------------------------------------------------------------------
C