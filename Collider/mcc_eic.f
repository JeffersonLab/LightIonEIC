      REAL FUNCTION mcc_eic(ff)
*******************************************************
      INCLUDE ?
      INTEGER ievent,first
      INTEGER i,j,k,l
      INTEGER il,jl,kl

      DATA first/0/

     
      DOUBLE Q2AVG(1)
      DATA Q2AVG/17.5/
      DOUBLE Q2AVGmin(1)
      DATA Q2AVGmin/15.0/
      DOUBLE Q2AVGmax(1)
      DATA Q2AVGmax/20.0/
      
      DOUBLE XAVGmin(3)
c$$$      DATA XAVGmin/0.02,0.04,0.06/
c$$$      DOUBLE XAVGmax(3)
c$$$      DATA XAVGmax/0.04,0.06,0.08/
c$$$      DOUBLE XAVG(3)
c$$$      DATA XAVG/0.035,0.055,0.075/

      DATA XAVGmin/0.04,0.06,0.08/
      DOUBLE XAVGmax(3)
      DATA XAVGmax/0.06,0.08,0.10/
      DOUBLE XAVG(3)
      DATA XAVG/0.05,0.07,0.09/


c$$$
c$$$      DOUBLE XAVGmin(3)
c$$$      DATA XAVGmin/0.0251,0.0398,0.0631/
c$$$      DOUBLE XAVGmax(3)
c$$$      DATA XAVGmax/0.0398,0.0631,0.1/
c$$$      DOUBLE XAVG(3)
c$$$      DATA XAVG/0.03245,0.05145,0.08155/

      DOUBLE WTPRIM,tprimx,tprim0
      DOUBLE weighty,wx
      DOUBLE pi,alpha_em,mn,mp,Edeu,xE_recol
      DATA pi,alpha_em,mn,mp,Edeu
     &/3.141592,0.0073,0.93955,0.93827,0.00222/
      REAL phibin,rphi
      DOUBLE mdeu

  
c     recall 
      DOUBLE WCRS1,WCRS2,WCRS0,WXJACOB
      DOUBLE spectral,WALPHA_R
      DOUBLE CVAL_TminMIN,CVAL_TminMAX

********************************************************
C     Histos Booking
********************************************************
      IF(first.eq.0.and.ff.eq.1)THEN
         first=1

         DO il=1,1
c            DO jl=1,5
            DO jl=1,1

               call hbook2(777,' tprime vs. alphaR',
     &                 100,0.7,1.4,100,0.0,0.2,0.) 
               call hbook2(776,' tprime vs. alphaR',
     &                 100,0.95,1.05,100,0.0,0.01,0.) 
 
              call hbook1(10000+1000*il+100*jl,'crs(si) vs. tprime',
     &                 20,0.0,0.1,0.)  
               call hbook1(20000+1000*il+100*jl,'F2D vs tprime',
     &                 20,0.0,0.1,0.)  
               call hbook1(50000+1000*il+100*jl,'F2D vs tprime',
     &                 20,0.0,0.1,0.)  
               call hbook1(30000+1000*il+100*jl,'crs weight by 1',
     &                 20,0.0,0.1,0.)  

               DO kl=1,19
               call hbook1(90000+1000*il+100*jl+kl,'crs weight by CRS0',
     &                 300,-0.1,0.1,0.)
               call hbook1(40000+1000*il+100*jl+kl,'crs weight by CRS2',
     &                 300,-0.1,0.1,0.)
               call hbook1(70000+1000*il+100*jl+kl,'crs weight by CRS1',
     &                 300,-0.1,0.1,0.)


               ENDDO

            ENDDO
         ENDDO

      ENDIF
*************************************************************
C     Calculation parts
*************************************************************      

      WCRS0 = CRS0
      WCRS1 = CRS1
      WCRS2 = CRS2
      WXJACOB = XJACOB
      WALPHA_R = alpha_r
      WTPRIM = TPRIM

c$$$      collfreq = 7.485*1.E8
      collfreq = 1
      fnanobarn2cm2 = 1.E-33
      flum = 1.E33
      weight = (crs*fnanobarn2cm2*flum/collfreq)
c     weight =1.0
      mdeu = 2*mn-Edeu

      ievent = ievent + 1
      r2d = 180/3.14159
      call hf2(1000,x,q2,weight)
c     1/2 is took into account in Event-Generator
c     I should remove this later, Once I have new MC event from EG.
c$$$      xE_recol = x/WXJACOB/2.
      xE_recol = x/WXJACOB/2.
      tprim0 = 2*(xE_recol**2-mp**2)-wtprim
c$$$      tprimx = wtprim-tprim0
      tprimx = wtprim-tprim0

c      print *,xE_recol,2*(xE_recol**2-mp**2),wtprim,tprim0
      CALL HF2(777,alpha_r,tprimx,1.)
      CALL HF2(776,alpha_r,tprimx,1.)

      DO il=1,1
         VAL_X= XAVG(il)
         VAL_XMIN=XAVGmin(il)
         VAL_XMAX=XAVGmax(il)
         

c         DO jl=1,5
         DO jl=1,1

            VAL_Q2 = Q2AVG(jl)
            VAL_Q2MIN = Q2AVGmin(jl)
            VAL_Q2MAX = Q2AVGmax(jl)

            p_r_mag = sqrt(p_rx**2+p_ry**2+p_rz**2)
            rphi = atan2(p_ry/p_r_mag,p_rx/p_r_mag)


            if(x.ge.VAL_XMIN.AND.x.lt.VAL_XMAX)then

               if(q2.ge.VAL_Q2MIN.AND.q2.lt.VAL_Q2MAX)then

                   theta_rr=acos(p_rz/p_r_mag)

		   theta_rd_lab=acos(p_rz/p_r_mag)*r2d
                      

                   IF(walpha_r.ge.0.9800.and.walpha_r.le.1.0200)THEN
                       CALL HF1(10000+1000*il+100*jl,tprim-tprim0,crs1)
                       CALL HF1(20000+1000*il+100*jl,tprim-tprim0,crs2)
                       CALL HF1(30000+1000*il+100*jl,tprim-tprim0,1.)
                       CALL HF1(50000+1000*il+100*jl,tprim-tprim0,crs0)
                   ENDIF
               
c tprime _min should be D_binding energy * M_D = 2*10^-3 GeV * 2 GeV = 4*10^-3 GeV^2
                  DO kl=1,19
                     CVAL_TminMIN = (0.005+0.005*(kl-1))
                     CVAL_TminMAX = (0.005+0.005*(kl))
c$$$                     if(wtprim.ge.CVAL_TminMIN.AND.
c$$$     &                    wtprim.lt.CVAL_TminMAX)then

                     if(tprimx.ge.CVAL_TminMIN.AND.
     &                    tprimx.lt.CVAL_TminMAX)then

                        e_r = sqrt(mp**2+(p_rx**2+p_ry**2+p_rz**2))
                        spectral = wcrs/wcrs1

c                      If(walpha_r.ge.0.9800.and.walpha_r.le.1.0200)then
cc                        If(walpha_r.ge.0.9800.and.walpha_r.le.1.000)then
                      If(walpha_r.ge.1.0000.and.walpha_r.le.1.0200)then
                           call hf1(90000+1000*il+100*jl+kl,
     &                          rphi,crs0)
                           call hf1(40000+1000*il+100*jl+kl,
     &                          rphi,crs2)
                           call hf1(70000+1000*il+100*jl+kl,
     &                          rphi,crs1)

                        ENDIF


                     endif
                    
                  ENDDO 
c  alpha cut (end)
c               ENDIF  

               endif

            endif

         ENDDO
      ENDDO




C     **************************************************  
      END
