!! CERNLIB E104 modified to be used with (G)JR GRIDS:
!! Name changed from fint to dfint.
!! Real variables changed to double precision.
!! External references to CERNLIB (error handling) routines removed.
          DOUBLE PRECISION FUNCTION DFINT(NARG,ARG,NENT,ENT,TABLE)
          INTEGER   NENT(9), INDEX(32)
          DOUBLE PRECISION ARG(9),   ENT(9),   TABLE(9), WEIGHT(32)
          DFINT  =  0d0
          IF(NARG .LT. 1  .OR.  NARG .GT. 5)  GOTO 300
          LMAX      =  0
          ISTEP     =  1
          KNOTS     =  1
          INDEX(1)  =  1
          WEIGHT(1) =  1d0
          DO 100    N  =  1, NARG
             X     =  ARG(N)
             NDIM  =  NENT(N)
             LOCA  =  LMAX
             LMIN  =  LMAX + 1
             LMAX  =  LMAX + NDIM
             IF(NDIM .GT. 2)  GOTO 10
             IF(NDIM .EQ. 1)  GOTO 100
             H  =  X - ENT(LMIN)
             IF(H .EQ. 0.)  GOTO 90
             ISHIFT  =  ISTEP
             IF(X-ENT(LMIN+1) .EQ. 0d0)  GOTO 21
             ISHIFT  =  0
             ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
             GOTO 30
  10         LOCB  =  LMAX + 1
  11         LOCC  =  (LOCA+LOCB) / 2
             IF(X-ENT(LOCC))  12, 20, 13
  12         LOCB  =  LOCC
             GOTO 14
  13         LOCA  =  LOCC
  14         IF(LOCB-LOCA .GT. 1)  GOTO 11
             LOCA    =  MIN0( MAX0(LOCA,LMIN), LMAX-1 )
             ISHIFT  =  (LOCA - LMIN) * ISTEP
             ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
             GOTO 30
  20         ISHIFT  =  (LOCC - LMIN) * ISTEP
  21         DO 22  K  =  1, KNOTS
                INDEX(K)  =  INDEX(K) + ISHIFT
  22            CONTINUE
             GOTO 90
  30         DO 31  K  =  1, KNOTS
                INDEX(K)         =  INDEX(K) + ISHIFT
                INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
                WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
                WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
  31            CONTINUE
             KNOTS  =  2*KNOTS
  90         ISTEP  =  ISTEP * NDIM
 100         CONTINUE
          DO 200    K  =  1, KNOTS
             I  =  INDEX(K)
             DFINT  =  DFINT + WEIGHT(K) * TABLE(I)
 200         CONTINUE
          RETURN
 300      WRITE(*,1000) NARG
          STOP
1000      FORMAT( 7X, 24HFUNCTION DFINT... NARG =,I6,
     +              17H NOT WITHIN RANGE)
          END
