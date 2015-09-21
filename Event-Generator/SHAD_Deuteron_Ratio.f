
      subroutine SHAD_Deuteron_Ratio(x,Q2,pt,pz,F2D_ratio)
      implicit double precision (A-H,O-Z)
      dimension x_ar(11),Q2_ar(13),pt_ar(7),pz_ar(13),shad4(13,11,7,13)
      dimension shad_temp_pz(13),shad_temp_pt(7),shad_temp_x(11),
     >shad_temp_Q2(13),B(13),C(13),D(13)
      
C Number of grid points in each variable
        Nx=11
        NQ2=13
        Npt=7
        Npz=13

C Reading the date file
        open(1,file='Deuteron_shadowing_tagged_F2ratio_2014_grid.dat',
     >status='unknown')
           
        do iq=1,NQ2     
        read(1,*) Q2_ar(iq)   
        
        do ix=1,Nx
        do ipt=1,Npt
        do iz=1,Npz

        read(1,*) x_ar(ix),pt_ar(ipt),pz_ar(iz),shad4(iq,ix,ipt,iz)
        
        enddo
        enddo
        enddo
        enddo

        close(1)

C Interpolation in four variables
        do iq=1,NQ2
        do ix=1,Nx
        do ipt=1,Npt  
        do iz=1,Npz  
           
        shad_temp_pz(iz)=shad4(iq,ix,ipt,iz)   
        enddo

        call spline(Npz,pz_ar,shad_temp_pz,B,C,D)
        shad_temp_pt(ipt)=seval(Npz,pz,pz_ar,shad_temp_pz,B,C,D)
        enddo

        call spline(Npt,pt_ar,shad_temp_pt,B,C,D)
        shad_temp_x(ix)=seval(Npt,pt,pt_ar,shad_temp_pt,B,C,D)
        enddo
        
        call spline(Nx,x_ar,shad_temp_x,B,C,D)
        shad_temp_Q2(iq)=seval(Nx,x,x_ar,shad_temp_x,B,C,D)
        enddo

        call spline(NQ2,Q2_ar,shad_temp_Q2,B,C,D)
        F2D_ratio=seval(NQ2,Q2,Q2_ar,shad_temp_Q2,B,C,D)

c      print *,'Hey print from inside function'
c      print *,x,Q2,pt,pz,F2D_ratio

        return
        end



C ---------------------------------------------------------------------
      SUBROUTINE SPLINE(N,X,Y,B,C,D)
C ---------------------------------------------------------------------
c***************************************************************************
C CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
C INTERPOLATION SUBROUTINES ARE TAKEN FROM
C G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
C COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
c      IMPLICIT REAL*8(A-H,O-Z)
      implicit double precision (A-H,O-Z)

      DIMENSION X(13),Y(13),B(13),C(13),D(13)
      NM1=N-1
      IF(N.LT.2) RETURN
      IF(N.LT.3) GO TO 250
      D(1)=X(2)-X(1)
      C(2)=(Y(2)-Y(1))/D(1)
      DO 210 I=2,NM1
        D(I)=X(I+1)-X(I)
        B(I)=2.0D0*(D(I-1)+D(I))
        C(I+1)=(Y(I+1)-Y(I))/D(I)
        C(I)=C(I+1)-C(I)
 210             CONTINUE
      B(1)=-D(1)
      B(N)=-D(N-1)
      C(1)=0.0D0
      C(N)=0.0D0
      IF(N.EQ.3) GO TO 215
      C(1)=C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N)=C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1)=C(1)*D(1)**2.0D0/(X(4)-X(1))
      C(N)=-C(N)*D(N-1)**2.0D0/(X(N)-X(N-3))
 215       CONTINUE
      DO 220 I=2,N
        T=D(I-1)/B(I-1)
        B(I)=B(I)-T*D(I-1)
        C(I)=C(I)-T*C(I-1)
 220             CONTINUE
      C(N)=C(N)/B(N)
      DO 230 IB=1,NM1
        I=N-IB
        C(I)=(C(I)-D(I)*C(I+1))/B(I)
 230             CONTINUE
      B(N)=(Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.0D0*C(N))
      DO 240 I=1,NM1
        B(I)=(Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.0D0*C(I))
        D(I)=(C(I+1)-C(I))/D(I)
        C(I)=3.0D0*C(I)
 240             CONTINUE
      C(N)=3.0D0*C(N)
      D(N)=D(N-1)
      RETURN
 250       CONTINUE
      B(1)=(Y(2)-Y(1))/(X(2)-X(1))
      C(1)=0.0D0
      D(1)=0.0D0
      B(2)=B(1)
      C(2)=0.0D0
      D(2)=0.0D0
      RETURN
      END
c
c***************************************************************************
C ---------------------------------------------------------------------
      double precision FUNCTION SEVAL(N,XX,X,Y,B,C,D)
C ---------------------------------------------------------------------
c***************************************************************************
C CALCULATE THE DISTRIBUTION AT XX BY CUBIC SPLINE INTERPOLATION.
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION X(13),Y(13),B(13),C(13),D(13)
      DATA I/1/
      IF(I.GE.N) I=1
      IF(XX.LT.X(I)) GO TO 310
      IF(XX.LE.X(I+1)) GO TO 330
 310       CONTINUE
      I=1
      J=N+1
 320       CONTINUE
      K=(I+J)/2
      IF(XX.LT.X(K)) J=K
      IF(XX.GE.X(K)) I=K
      IF(J.GT.I+1) GO TO 320
 330       CONTINUE
      DX=XX-X(I)
      SEVAL=Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
      RETURN
      END
c

 
