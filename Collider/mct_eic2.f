      REAL FUNCTION mct_eic2(ff)
*******************************************************
      INCLUDE ?
      INTEGER ievent,first
      INTEGER i,j,k,l
      INTEGER il,jl,kl

      DATA first/0/

      REAL Q2AVG(10)
      DATA Q2AVG/0.6,1.2,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5/
      REAL Q2AVGmin(10)
      DATA Q2AVGmin/0.5,0.7,1.5,3.,4.,5.,6.,7.,8.,9./
      REAL Q2AVGmax(10)
      DATA Q2AVGmax/0.7,1.5,3.,4.,5.,6.,7.,8.,9.,10./
      
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

      DOUBLE weighty
      Real pi,alpha_em,mn
      DATA pi,alpha_em,mn/3.141592,0.0073,0.93955/
********************************************************
C     Histos Booking
********************************************************
      IF(first.eq.0.and.ff.eq.1)THEN
         first=1
c$$$
c$$$               call hbook2(1000,'Q^2 vs. x?b!',
c$$$     &        100,0.0,2.0,100,1.0,15.0,0.)

         DO il=1,5
            DO jl=1,9

               DO kl=1,19
            call hbook1(80000+1000*il+100*jl+kl,'crs si0',
c     &                 100,0.0,1.e3,0.)
     &                  120,-3.2,3.2,0.)



               ENDDO

            ENDDO
         ENDDO

      ENDIF
*************************************************************
C     Calculation parts
*************************************************************      
      if(crs.eq.0)return
c$$$      collfreq = 7.485*1.E8
      collfreq = 1
      fnanobarn2cm2 = 1.E-33
      flum = 1.E33
      weight = (crs*fnanobarn2cm2*flum/collfreq)
c     weight =1.0

      ievent = ievent + 1
      r2d = 180/3.14159
      call hf2(1000,x,q2,weight)
      

      DO il=1,5
         VAL_X= XAVG(il)
         VAL_XMIN=XAVGmin(il)
         VAL_XMAX=XAVGmax(il)

         DO jl=1,9

         if(x.lt.0.2)then         
            VAL_Q2 = Q2AVG(jl)
            VAL_Q2MIN = Q2AVGmin(jl)
            VAL_Q2MAX = Q2AVGmax(jl)
         else
            VAL_Q2 = Q2AVG2(jl)
            VAL_Q2MIN = Q2AVGmin2(jl)
            VAL_Q2MAX = Q2AVGmax2(jl)
         endif
            
            p_r_mag = sqrt(p_rx**2+p_ry**2+p_rz**2)
            rphi = atan2(p_ry/p_r_mag,p_rx/p_r_mag)

            if(x.ge.VAL_XMIN.AND.x.lt.VAL_XMAX)then
               
               if(q2.ge.VAL_Q2MIN.AND.q2.lt.VAL_Q2MAX)then
                   theta_rr=acos(p_rz/p_r_mag)
                   theta_rd_lab=acos(p_rz/p_r_mag)*r2d

                   If(alpha_r.gt.0.97.and.alpha_r.lt.1.03)then

                  DO kl=1,19
                     CVAL_TminMIN = (0.001+0.005*(kl-1))
                     CVAL_TminMAX = (0.001+0.005*(kl))
                     if(tprim.ge.CVAL_TminMIN.AND.
     &                    tprim.lt.CVAL_TminMAX)then

                           call hf1(80000+1000*il+100*jl+kl,
     &                          rphi,1.)
c     
                        endif
                    
                  ENDDO 
               endif

               endif
            endif

         ENDDO
      ENDDO




C     **************************************************  
      END

