      Subroutine fsgen_phot(ierr,nevent)
c  ---------------
c
c  Notes:
c  ------
c
c  Author:    S.Stepanyan     Created:  Tue Mar 25 14:32:25 EST 1997
c  -------
c
c  Major revisions:
c  ----------------
c     
c
      implicit none
      save
#include "bcs.inc"
#include "wcs.inc"
#include "leptou.inc"
#include "ludat.inc"
#include "lujet.inc"
#include "parm1.inc"
#include "parm2.inc"
#include "parm3.inc"
#include "gbd.inc"
c
      integer mlink,nw,lucomp
      integer ierr,nevent,err,i,ne_direct,ne_file
      REAL sinthetae1,theta_e1,phi_e1,brem
      REAL getrand,ULMASS,PI,r2d
      data pi,r2d/3.141592654,57.29578/
      integer jmcvx,jmctk,ntrk,id
      logical eventin,first,user_fc,wtfact
      data first/.true./
      real e_electron,e_photon,xb,fmg,rm2,phiminel,phimaxel
      data fmg,rm2/5.07,0.5881/
      REAL Q2Space(200)
      real etrg, phi,theta
      Real Fq2
      External Fq2,user_fc,lucomp
c
      ne_file=0
      ierr=0
      if(first)then
         wtfact=.true.
         nw=0
         first=.false.
         if(nq2.gt.0..and.q2max.gt.0..and.q2min.gt.0.)call FunLXP(Fq2,Q2Space,q2min,q2max)
      endif
 91   continue
      ne_direct=0
      PHMASS=0.
      TOTMASS=0.
      LINE=0
      ne_file=ne_file+1
      np=npart
c target definition
      do i=1,3
         ptarg(i)=0.
      enddo
      ptarget=0
      if(IDtargetN.ne.0)then
 92      continue
         CALL fermi_bonn(ptarg)
         ptarget=0
         do i=1,3
            ptarget=ptarget+ptarg(i)*ptarg(i)
         enddo
         ptarget=sqrt(ptarget)
         if(IDspect.gt.0)then
            spmass=Pmas(lucomp(IDspect),1)
            pspect=ptarget
            espect=sqrt(pspect*pspect+spmass*spmass)
         endif
         etarget=tmass_all-espect
         tmass=etarget*etarget-ptarget*ptarget
         if(tmass.lt.0.)go to 92
         tmass=sqrt(tmass)
      endif
c
      DO I=1,NP
 55      continue
         if(kp(i).eq.71)then
            amass(i)=pmas(kp(i),1)+getrand(-1.)*pmas(kp(i),2)
            call hf1(301,amass(i),1.)
         else
            AMASS(I) = ULMASS (KP(I),1)
          endif
          if(AMASS(I).lt.minmass(i))goto55
          TOTMASS=TOTMASS+AMASS(I)
      ENDDO
      IF(bin.eq.0)THEN
        if(ne_file.eq.1)then
           call bdrop(jw,'E')
           call frbos(jw,jin,'E',err)
          if(err.eq.-1)then
            write(lout,*)'End of file ',inam
            ierr=-1
            go to 99
          elseif(err.eq.1)then
            write(lout,*)'Error reading the file ',inam
            ierr=1
            go to 99
          endif
          jmcvx=mlink(jw,'MCVX',0)
          jmctk=mlink(jw,'MCTK',0)
          eventin=jmctk.gt.0
          if(eventin)then
            ntrk=jw(jmctk-4)
            if(ntrk.gt.1)then
               write(lout,*)'To many electrons in the input'
               ierr=10
               goto99
            endif
            e1=rjw(jmctk+4)
            theta_e1=acos(rjw(jmctk+3))
            phi_e1=atan2(rjw(jmctk+2),rjw(jmctk+1))
            id=jw(jmctk+7)
            q2=2.*ebeam*e1*(1.-cos(theta_e1))
            u=ebeam-e1
            w2=pmass**2+2.*pmass*u-q2
            TECM=sqrt(w2)
            TECM=sqrt(Tmass**2+2*u*Tmass-Q2)
          endif
        endif
        if(ne_file.gt.100)ne_file=0
 86     call massgen()
        IF(TECM.le.TOTMASS+0.001)then
           GOTO86
        endif
        CALL lu1ent(1,id,e1,theta_e1,phi_e1)
        line=1
        call lulist(1)
        phmass=q2
        u=ebeam-e1
        qv=SQRT(Ebeam**2+E1**2-2.*Ebeam*E1*COS(THETA_E1))
      ELSEIF(Q2MAX.NE.0..or.Q2Min.ne.0.)THEN
         ID=11
 98      ne_direct=0
C     
         if(WMAX.gt.0.)then
            W=WMIN+getrand(-1.)*(WMAX-WMIN)
            W2=W**2
C
            if(q2max.gt.0.)then
               If(nq2.gt.0)then
                  call FunLUX(Q2Space,q2,1)
               Else
                  Q2=Q2MIN+getrand(-1.)*(Q2MAX-Q2MIN)
               EndIf
               PHMASS=-Q2
            elseif(xmax.gt.0)then
               xb=xmin+getrand(-1.)/100*(xmax-xmin)
               if(cohl.eq.0.)then
                  q2=xb*(w2-pmass**2)/(1.-xb)
               else
                  q2=(w2-pmass**2-pmass*fmg*xb*rm2)/(pmass*fmg*xb-1)
               endif
            else
               print *,'Wrong definition Q2 or Xb must be defined'
               stop
            endif
         elseif(xmax.gt.0)then
            xb=xmin+getrand(-1.)/100*(xmax-xmin)
            if(q2max.gt.0.)then
               If(nq2.gt.0)then
 211              call FunLUX(Q2Space,q2,1)
                  if(q2.lt.q2min.or.q2.gt.q2max)goto211
               Else
                  Q2=Q2MIN+getrand(-1.)*(Q2MAX-Q2MIN)
               EndIf
               PHMASS=-Q2
               if(cohl.eq.0.)then
                  w2=pmass**2+q2/xb-q2
               else
                  w2=pmass**2+pmass*xb*fmg*(q2+rm2)-q2
               endif
               w=sqrt(w2)
            else
               print *,'Wrong definition, Q2 or W must be defined'
               stop
            endif
         else
            print *,'Wrong definition, Xb or W must be defined'
            stop
         endif
C
         U=(W2-PMASS**2+Q2)/2./PMASS
         IF(U.GT.Ebeam-0.02)GO TO 98
c         TECM=W
c         TECM=sqrt(Tmass**2+2*u*Tmass-Q2)
c         IF(TECM.le.TOTMASS)GOTO98
         E1=Ebeam-U
         SINTHETAE1=SQRT(Q2/4./Ebeam/E1)
         THETA_E1=2.*ASIN(SINTHETAE1)
         do i=1,10
            if(pidfc(i).eq.11)then
               if(ffc(i).gt.0)then
                  if(theta_e1*r2d.lt.thfcm(i).or.theta_e1*r2d.gt.thfcx(i)) go to 98
                  phimaxel=phfcx(i)/r2d
                  phiminel=phfcm(i)/r2d
                  PHI_E1=phiminel+(phimaxel-phiminel)*getrand(-1.)
                  go to 96
               elseif(ffc(i).eq.0)then
                  phimaxel=2*pi
                  phiminel=0.
                  PHI_E1=phiminel+(phimaxel-phiminel)*getrand(-1.)
                  go to 96
               elseif(ffc(i).lt.0)then
                  phimaxel=phfcx(i)/r2d
                  phiminel=phfcm(i)/r2d
                  PHI_E1=phiminel+(phimaxel-phiminel)*getrand(-1.)
                  if(user_fc(id,E1,theta_e1,phi_e1))go to 96
                  go to 98
               endif
            endif
         enddo
C     
 96      qv=SQRT(Ebeam**2+E1**2-2.*Ebeam*E1*COS(THETA_E1))
c         phi_e1=pi
         phi_q=phi_e1+pi
         if(phi_q.gt.pi)phi_q=phi_q-2.*pi
         qx=-e1*cos(phi_e1)*SIN(THETA_E1)
         qy=-e1*sin(phi_e1)*SIN(THETA_E1)
         qz=Ebeam-E1*COS(THETA_E1)
         theta_q=acos(qz/qv)
         TECM=sqrt((etarget+u)**2-
     -        (qx+ptarg(1))**2-(qy+ptarg(2))**2-(qz+ptarg(3))**2)
 87      call massgen()
         if(ne_direct.gt.10000)go to 98
         ne_direct=ne_direct+1
         IF(TECM.le.TOTMASS+0.00001)then
            GOTO87
         endif
         if(wtfact)then
            CALL GBD(1)
            if(wt.gt.wt0)then
               wt0=wt
               nw=0
            else
               nw=nw+1
            endif
            if(nw.gt.10000)wtfact=.false.
            go to 91
         endif
         line=1
         CALL LU1ENT(line,ID,E1,THETA_E1,PHI_E1)
         if(IDtargetN.ne.0.and.IDspect.ne.0)then
            line=2
            phi=atan2(-ptarg(2),-ptarg(1))
            theta=acos(-ptarg(3)/pspect)
            CALL LU1ENT(line,IDspect,Espect,theta,phi)
            k(line,1)=1
         endif
      ELSE
 97      continue
         ne_direct=0
         if(wtfact)then
            qv=WMAX
         else
            qv=WMIN+getrand(-1.)*(WMAX-WMIN)
            brem=getrand(-1.)/Wmin
            if(1./qv.lt.brem)go to 97
         endif
         W2=PMASS**2+2.*PMASS*qv
         Q2=0.0
         U=qv
         e_electron=ebeam
         e_photon=U
         qx=0.
         qy=0.
         qz=qv
         TECM=sqrt((etarget+u)**2-
     -        (qx+ptarg(1))**2-(qy+ptarg(2))**2-(qz+ptarg(3))**2)
 88      call massgen()
c         print *,ne_direct,u,tecm,totmass,wmax,ptarg
         if(ne_direct.gt.10000)go to 91
         ne_direct=ne_direct+1
         IF(TECM.le.TOTMASS+0.00001)then
            GOTO88
         endif
         if(wtfact)then
            CALL GBD(1)
            if(wt.gt.wt0)then
               wt0=wt
               nw=0
            else
               nw=nw+1
            endif
            if(nw.gt.10000)wtfact=.false.
            go to 91
         else
            id=22
            ierr = 0
            line=1
            CALL LU1ENT(line,ID,u,0.,0.)
            k(1,1)=1
            if(IDtargetN.ne.0.and.IDspect.ne.0)then
             line=2
             phi=atan2(-ptarg(2),-ptarg(1))
             theta=acos(-ptarg(3)/pspect)
             CALL LU1ENT(line,IDspect,Espect,theta,phi)
             k(line,1)=1
           endif
         endif
      ENDIF
 99   return
      end
c
      subroutine massgen()
      implicit none
#include "leptou.inc"
#include "ludat.inc"
#include "lujet.inc"
#include "parm1.inc"
#include "parm3.inc"
#include "gbd.inc"
      integer i
      real ulmass, getrand, dm
      external ulmass
c
      totmass=0.
      DO I=1,NP
 55      continue
         if(kp(i).eq.71)then
            dm=min(tecm-totmass,pmas(kp(i),2))
            amass(i)=pmas(kp(i),1)+getrand(-1.)*(dm-pmas(kp(i),1))
c            print *,i,kp(i),pmas(kp(i),1),pmas(kp(i),2),amass(i)
            call hf1(301,amass(i),1.)
         else
            AMASS(I) = ULMASS (KP(I),1)
          endif
          if(AMASS(I).lt.minmass(i))goto55
          TOTMASS=TOTMASS+AMASS(I)
      ENDDO
c
      return
      end
