      REAL FUNCTION mct_eic_report(ff)
*******************************************************
      INCLUDE ?
      INTEGER ievent,first
      INTEGER i,j,k,l
      INTEGER il,jl,kl

      DATA first/0/

      REAL Q2AVG(10)
      DATA Q2AVG/0.6,1.2,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5/
      REAL Q2AVGmin(10)
      DATA Q2AVGmin/0.5,0.7,1.5,3.0,4.,5.,6.,7.,8.,9./
      REAL Q2AVGmax(10)
      DATA Q2AVGmax/0.7,1.5,3.0,4.0,5.,6.,7.,8.,9.,10./
      
      REAL Q2AVG2(10)
      DATA Q2AVG2/0.7,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5/
      REAL Q2AVGmin2(10)
      DATA Q2AVGmin2/0.5,1.,2.,3.,4.,5.,6.,7.,8.,9./
      REAL Q2AVGmax2(10)
      DATA Q2AVGmax2/1.0,2.,3.,4.,5.,6.,7.,8.,9.,10./

      REAL XAVGmin(6)
c      DATA XAVGmin/0.0,0.2,0.4,0.6,0.8/
      DATA XAVGmin/0.05,0.1,0.2,0.4,0.6,0.8/
      REAL XAVGmax(6)
c      DATA XAVGmax/0.2,0.4,0.6,0.8,1.0/
      DATA XAVGmax/0.1,0.2,0.4,0.6,0.8,1.0/
      REAL XAVG(6)
      DATA XAVG/0.07,0.15,0.3,0.5,0.7,0.9/

      REAL weight
      Real pi,alpha_em,mn,mp
      DATA pi,alpha_em,mn,mp/3.141592,0.0073,0.93955,0.93835/
********************************************************
C     Histos Booking
********************************************************
      IF(first.eq.0.and.ff.eq.1)THEN
         first=1

c     event is weighted by cross-section
         call hbook2(1000,' p?RT! vs. [a]?R!',
     &        100,0.7,1.3,100,0.0,0.3,0.)
         call hbook2(1001,' t^,! vs. [a]?R!',
     &        100,0.7,1.3,100,0.0001,0.18,0.)
         
c     event is un-weighted 
         call hbook2(2000,' p?RT! vs. [a]?R!',
     &        100,0.7,1.3,100,0.0,0.3,0.)
         call hbook2(2001,' t^,! vs. [a]?R!',
     &        100,0.7,1.3,100,0.0001,0.18,0.)


 

      ENDIF
*************************************************************
C     Calculation parts
*************************************************************      
      if(crs.eq.0)return
c$$$      collfreq = 7.485*1.E8
      collfreq = 1
      fnanobarn2cm2 = 1.E-33
      flum = 1.E33
c$$$      weight = (crs1*fnanobarn2cm2*flum/collfreq)
      weight = crs1
c     weight =1.0

      ievent = ievent + 1
      r2d = 180/3.14159
      call hf2(1000,alpha_r,p_rt,weight)
      call hf2(1001,alpha_r,tprim,weight)

      call hf2(2000,alpha_r,p_rt,1.)
      call hf2(2001,alpha_r,tprim,1.)




C     **************************************************  
      END

