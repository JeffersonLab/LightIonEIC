     REAL FUNCTION mcc_eic2(ff)
*******************************************************
      INCLUDE ?
      INTEGER ievent,first
      INTEGER i,j,k,l
      INTEGER il,jl,kl

      DATA first/0/

c$$$      REAL Q2AVG(7)
c$$$      DATA Q2AVG/5.5,6.5,8.0,10.0,12.0,14.0,17.5/
c$$$      REAL Q2AVGmin(7)
c$$$      DATA Q2AVGmin/5.0,6.0,7.1,8.0,9.0,11.0,13.0/
c$$$      REAL Q2AVGmax(7)
c$$$      DATA Q2AVGmax/6.0,7.0,8.0,8.9,11.0,13.0,15.0/
      REAL Q2AVG(1)
      DATA Q2AVG/17.5/
      REAL Q2AVGmin(1)
      DATA Q2AVGmin/15.0/
      REAL Q2AVGmax(1)
      DATA Q2AVGmax/20.0/

      DATA XAVGmin/0.04,0.06,0.08/
      DOUBLE XAVGmax(3)
      DATA XAVGmax/0.06,0.08,0.10/
      DOUBLE XAVG(3)
      DATA XAVG/0.05,0.07,0.09/
c$$$
c$$$      REAL XAVGmin(4)
c$$$      DATA XAVGmin/0.02,0.04,0.06,0.08/
c$$$      REAL XAVGmax(4)
c$$$      DATA XAVGmax/0.04,0.06,0.08,0.10/
c$$$      REAL XAVG(4)
c$$$      DATA XAVG/0.035,0.055,0.075,0.095/


      DOUBLE weighty
      Real pi,alpha_em,mn
      DATA pi,alpha_em,mn/3.141592,0.0073,0.93955/
********************************************************
C     Histos Booking
********************************************************
      IF(first.eq.0.and.ff.eq.1)THEN
         first=1

         DO il=1,4
c$$$            DO jl=1,5
            DO jl=1,1

               DO kl=1,19
            call hbook1(80000+1000*il+100*jl+kl,'crs FGIS',
     &                 300,-0.1,0.1,0.)
            call hbook1(50000+1000*il+100*jl+kl,'crs FSIG*Jacob',
     &                 300,-0.1,0.1,0.)
            call hbook1(60000+1000*il+100*jl+kl,'crs FSIG*Jacob1',
     &                 300,-0.1,0.1,0.)

               ENDDO

            ENDDO
         ENDDO

      ENDIF
*************************************************************
C     Calculation parts
*************************************************************      
      collfreq = 1
      fnanobarn2cm2 = 1.E-33
      flum = 1.E33
      weight = (crs*fnanobarn2cm2*flum/collfreq)
c     weight =1.0

      ievent = ievent + 1
      r2d = 180/3.14159
      call hf2(1000,x,q2,weight)
      

      DO il=1,4
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


                  DO kl=1,19
                     CVAL_TminMIN = (0.005+0.005*(kl-1))
                     CVAL_TminMAX = (0.005+0.005*(kl))
                     if(tprim.ge.CVAL_TminMIN.AND.
     &                    tprim.lt.CVAL_TminMAX)then

                        If(alpha_r.gt.0.98.and.alpha_r.lt.1.02)then
                           call hf1(80000+1000*il+100*jl+kl,
     &                          rphi,crs0)
                           call hf1(50000+1000*il+100*jl+kl,
     &                          rphi,crs2)
                            call hf1(60000+1000*il+100*jl+kl,
     &                          rphi,crs1)
                        ENDIF
c     
                     endif
                    
                  ENDDO 

               endif
            endif

         ENDDO
      ENDDO




C     **************************************************  
      END
