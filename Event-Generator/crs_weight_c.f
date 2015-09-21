      program crs_weight_c
c   test version-1.0
        common/par/pi,pm,pmp,pmn,dm,eb
        character (len=99) :: filename
        integer :: i, j, lines,IPN
      
c        real, dimension(999,20) :: aa
        real hvec0(14),h(14)

        real ee(11),evec1(8)
        real sp(11),hvec1(8)
        real zz(11),hvec2(8)
        real ww(11),hvec3(8)
        REAL qv(3),qv_mag
        REAL pfs(3),pfs_q(3) 
        real pfscm(3)
        real mprot,mneut
        real R2D,alpha_em

        double precision :: f2nn,f2nn3,RES,FSIG,F2D,FLD,S,E_recoil
        double precision :: xkjfac,xjacobian,crs0,crs1,crs2
        double precision :: PI
        double precision :: E_ebeam,E_ibeam,Pol_ebeam,Pol_ibeam
        integer nptl,pid,pid1,pid2
        integer egen,sgen,xgen
        integer echarge,scharge,xcharge
        R2D = 180/3.14
        alpha_em = 0.0073

        ntgt =2
        nproton = 1



c     deutron mass and momentum
        dm = 1.875
        mprot =0.93827
        mneut =0.93955

c     Christian's definition
        ED = 0.00222
        UD = 2*mneut - ED
        PI = 0.314159265358979D+01

c     selection of model (M. Sargsian): 1,  (C. Weiss): 2
        imodel = 2
c     Call F2 parameterization for : (proton):1, (neutron):2, (P.Jimenez-Delgad): 3
        IPN =2
        
*********************************
* Initialization
*********************************	

C for GEMC format 
        iu = 77
        write (filename, '( "ed_semi_eic",I2,".dat" )' )  iu
        OPEN(unit=iu,file=filename, status='new')
          
C for ANALYSIS format 
       iw = 79
        write (filename, '( "ed_semi_eic",I2,".dat" )' )  iw
        OPEN(unit=iw,file=filename, status='new')
          

c  fsgen 100K event generated.....

        OPEN(unit=99,file='gn_ed_eic.input', status='old',
     &       ACTION='READ')

c      Misak's Cross-section model (IA) ............. LC/VN
c      initialization edenx function.........
        IF(imodel.eq.1)THEN
           icase = 0
           call edenx(se,sq,q2,q0,ktr,alpha_r,p_rt,x,s,si,
     &          Fd_L,Fd_T,Fd_TL,Fd_TT,f2d,f1n,f2n,sun,icase,lc,lbj)
           x = 0.0
        ELSE IF(imodel.eq.2)THEN
C     ...KINEMATIC LIMIT IN TPRIME, ABSOLUTE AND FOR GIVEN ALPHAR
           CALL TPMIN(TP0, 1.D0, 0)
           CALL TPMIN(TP1,  alpha_r, 1)
           IF(IPN.eq.3)THEN
              call JR14NLO08SFinit
           ENDIF
      ENDIF


c if more than 99999 causes the memory problem...
         DO ll=1,99999
C additional head information for beam particles (E_e, E_ion, Pol_e, Pol_ion)
           read (99,*) h1,h2
     &           ,h(1),h(2),h(3),h(4),h(5),h(6),h(7),h(8),h(9),h(10)
     &           ,h(11),h(12),h(13),h(14)
           read (99,*) ee(1),ee(2),ee(3),ee(4),ee(5),ee(6)
     &          ,ee(7),ee(8),ee(9),ee(10),ee(11)
           read (99,*) sp(1),sp(2),sp(3),sp(4),sp(5),sp(6)
     &          ,sp(7),sp(8),sp(9),sp(10),sp(11)
c           read (99,*) zz(1),zz(2),zz(3),zz(4),zz(5),zz(6)
c     &          ,zz(7),zz(8),zz(9),zz(10),zz(11)
  
           nptl =2
        
           nptlevnt = h1
           ievntnum  = h2

c      igen = daughter particle (2nd generation)

        DO iv=1,14
           hvec0(iv) =  h(iv)
        ENDDO

        xinv_xbj    = hvec0(1) 
        xinv_q2     = hvec0(2) 
        xinv_se     = hvec0(3) 
        xinv_sq     = hvec0(4) 
        xinv_alpha  = hvec0(5) 
        xinv_pRT    = hvec0(6) 
        xinv_tprim  = hvec0(7) 
        xinv_jacob  = hvec0(8) 
c     p_RT2 = tempVar - pSpectator_Rest_z**2
        xinv_pRT0   = hvec0(9)
c     sqrt(invts.tPrime/2)
        xinv_pRTs   = hvec0(10)

        E_ebeam   = hvec0(11)
        E_ibeam   = hvec0(12)
        xebeam_pol = hvec0(13)
        xibeam_pol = hvec0(14)

c     electron --------------------- trigger
        pid = ee(2)
        egen = 1
        echarge = -1
        DO iev=4,11
           ieev = iev -3
           evec1(ieev)= ee(iev)
        ENDDO
c     hadron 1 --------------------- spectator proton
        pid1 = sp(2)
        sgen = 2
        scharge = 1
        DO ihv1=4,11
           ihhv1 = ihv1 -3
           hvec1(ihhv1)= sp(ihv1)
        ENDDO
c     hadron 2 --------------------- anything dummy hadron X
c        pid2 = zz(2)
c        xgen = 3
c        xcharge = 0
c        DO ihv2=4,11
c           ihhv2 = ihv2 -3
c           hvec2(ihhv2)= zz(ihv2)
c        ENDDO




        pe = sqrt(evec1(1)**2+evec1(2)**2+evec1(3)**2)
        pe_x = evec1(1)
        pe_y = evec1(2)
        pe_z = evec1(3)
        pe_E = evec1(4)
        pe_M = evec1(5)
        ve_x = evec1(6)
        ve_y = evec1(7)
        ve_z = evec1(8)


        q0 = xinv_q2/2.0/mprot/xinv_xbj  
       epr  = E_ebeam - pe
c theta_e unit : radian
c        th_e = acos(evec1(3)/pe)
c        th_e = 2.0*asin(sqrt(xinv_q2/4.0/E_ebeam/epr))        
       th_e = acos(1-xinv_q2/(2*E_ebeam*epr))
c fsgen : Q2 range set 0. - 10.GeV2
c        q2 = 4.*E_ebeam*pe*sin(th_e/2.)**2
        qv_mag   = sqrt(xinv_q2 + q0**2)

        cs_the_qe = (E_ebeam - pe*cos(th_e))/qv_mag
   
        yy = q0/E_ebeam

c virtual photon component
        qx=-pe*cos(atan2(pe_y,pe_x))*sin(th_e)
        qy=-pe*sin(atan2(pe_y,pe_x))*sin(th_e)
        qz=E_ebeam-pe*cos(th_e)
        qv(1)  = qx
        qv(2)  = qy
        qv(3)  = qz
        qv_mag = sqrt(dot(qv,qv))

        call extractang(the_q,the_qd,phi_q,phi_qd,qv)
        beta  = sqrt(dot(qv,qv))/(q0+dm)

        ed =  sqrt(dm**2 + E_ibeam**2)

c invariant energy square of e + d system (k_e + p_d)**2
        se = 2.0*E_ebeam*(ed + E_ibeam) + dm**2
        sq = -q2 + 2.0*q0*ed + 2.0*qv_mag*cs_the_qe*E_ibeam + dm**2

c spectator proton 
        p_r = sqrt(hvec1(1)**2+hvec1(2)**2+hvec1(3)**2)
c transverse momentum of recoil nucleon vs q ?? <--------
        p_rt= sqrt(hvec1(1)**2+hvec1(2)**2)
        p_rx= hvec1(1)
        p_ry= hvec1(2)
        p_rz= hvec1(3)
        v_rx= hvec1(6)
        v_ry= hvec1(7)
        v_rz= hvec1(8)
        xnthetap_r = acos(p_rz/p_r)
c        E_r = sqrt(pm**2+p_r**2) : same as below
        E_r = hvec1(4)
        pm = hvec1(5)
       
        pfs(1) = p_rx
        pfs(2) = p_ry
        pfs(3) = p_rz
        pfs_mag = sqrt(dot(pfs,pfs))
c---------------------------------------------------------------
c     Calculate C.M varibles
c---------------------------------------------------------------
        call rotate(-1,the_q,phi_q,pfs,pfs_q)
        call extractang(the_hq,the_hqd,phi_hq,phi_hqd,pfs_q)
        call boost(-beta,pfs_q,E_r,pfscm,E_h_cm)
        call extractang(the_h_cm,the_hd_cm,phi_h_cm,phi_hd_cm,pfscm)
        
        pf_cm_mag=sqrt(dot(pfscm,pfscm))

        
c
c p_n = (P_D - p_S)  as a four-vector, then compute  p_n^2  = t
c compare this event-by-event with
c M_D^2 + M^2 - 2*P_D \cdot p_S = t


c dm**2-2.0*dm*e_r+pm**2+2.0*(dm-e_r)*q0 + 2.0*pr*qv*cos(thr)-q2
c invariant mass square definition
       	wn22 = dm**2-2.0*dm*E_r+pm**2
     &       +2.0*(dm-E_r)*q0 +2.0*p_r*qv_mag*cos(xnthetap_r)-q2

        p2q = (E_r*q0)-(p_rx*qx+p_ry*qy+p_rz*qz)
        xstar = q2/(2*p2q)

        wn2 = (dm-E_r)**2-p_r**2 + q2*(1-xstar)/xstar
        wn = sqrt(wn2)

        qv2pr = (p_rx*qx+p_ry*qy+p_rz*qz)
        theta_qp = the_h_cm
        phi_qp = phi_h_cm



        alpha_r = xinv_alpha 


C missing mass2 distribution
         e_x=E_ebeam+dm-pe-sqrt(p_r**2+pm**2)

         px_x=qx
     &        -p_r*cos(atan2(p_ry,
     &        p_rx))*sin(xnthetap_r)

         py_x=qy
     &        -p_r*sin(atan2(p_ry,
     &        p_rx))*sin(xnthetap_r)

         pz_x=qz-p_r*cos(xnthetap_r)



         p_x2=sqrt(px_x**2+py_x**2+pz_x**2)
         xmx2=sqrt(e_x**2-p_x2**2)
c Misak's definition
         xmx22 = sqrt( dm**2 - 2.0*E_r*dm + pm**2 +
     &         2.0*q0*(dm-E_r) - 2.0*qv_mag*p_rz - q2)



         tprime_KJ =  xinv_tprim
         t_chn_MK = tprime_KJ + pm**2



        IF(imodel.eq.1)THEN

c option nucleon -1: neutron, +1: proton
        ktr  = -1
c selection of bindidng model lc=1:light cone, 0: virtual nucleon
        lc = 1
c set of bjorken limit 0: no limit, 1: set value
        lbj = 1
c initialization 0: calculation, 1: PWIA
        icase = 1

          call edenx(xinv_se,xinv_sq,xinv_q2,q0,ktr,xinv_alpha,xinv_pRT,
     &  x,s,si,Fd_L,Fd_T,Fd_TL,Fd_TT,f2d,f1eff,f2eff,sun,icase,lc,lbj)



        icon =1
        thetr =  theta_qp
        phir =  phi_qp
CCC----	subroutine  spectral_un(q0,qv,pr,thr,phir,ktr,tr,sun,lc,icon)
       call spectral_un(q0,qv_mag,p_r,thetr,phir,1,E_r,sun,lc,1)        
c*******************************************************************
*         Spectral Function
*         q0 - virtual photon energy
*         qv - virtual photon momentum
*         pr - recoil nucleon momentum
*         thr - recoil nucleon polar angle vs q
*         phir - recoil nucleon azimuthal angle
*         ktr - type of the recoil nucleon (1):proton, (-1):neutron 
* ------------- (IPN for Weiss code (1):proton, (2):neutron )
*         tr  - recoil nucleons kinetic energy
*         sun  - spectral function 
*         lc  - 0 virtual nucleon, 1 light cone  approximation
*         icon- 0-initialization, 1 PWIA, 2 - FSI, 12- PWIA+FSI
********************************************************************


c     wave function with Paris potential
        CC =0.3939
c     wave function with Bonn potential
c        CC =0.3930
        
        if(lc.eq.1)then
           xnn = dm/(2*(dm-E_r))
        else
           xnn = 1/(2-alpha_r)
        endif

        ResPot = (CC/(1.414*pi*2*mprot))**2


        xstar2 = x/alpha_r


        if(sqrt(xinv_se).gt.0.0)then
           if(si.gt.0 .and. tprime_KJ.gt.0)then


              fac1 = (4*pi*alpha_em**2)/(x*q2**2)
              fac2 = (1-yy-((x**2)*(yy**2)*(mneut**2))/q2)
              fac = fac1*fac2

              spect = tprime_KJ**2/(E_r*ResPot*xnn)
            
              xjacobian = xinv_jacob

              crs0 = si
              crs1 = si*xjacobian
              crs2  = si*spect/fac
              dum2 =0.

    
c$$$  INPUT for the GEMC format........................
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
c     xinv_se,xinv_sq,xinv_q2,q0,ktr,xinv_alpha,xinv_pRT,
              write(iu,101) nptl,
     &  dble(xinv_xbj),dble(xinv_q2),dble(xinv_se),dble(xebeam_pol),
     &  dble(xinv_alpha),dble(xinv_pRT),dble(xinv_tprim),dble(FSIG),
     &  dble(xibeam_pol)
              write(iu,11) egen,echarge,pid,pe_x,pe_y,pe_z,pe_E,pe_M,
     &             ve_x,ve_y,ve_z
              write(iu,11) sgen,scharge,pid1,p_rx,p_ry,p_rz,E_r,pm,
     &             v_rx,v_ry,v_rz 
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
c

c     dummy variable re-assign
              wn=1
              xmx2 =1
c SHOULD BE SAME FORMAT AS MODEL==2
              write(iw,10) xinv_xbj,xinv_q2,xinv_pRT,xinv_alpha,
     &             crs0,crs1,e_x,xjacobian,crs2,
     &             xinv_tprim,t_chn_MK,theta_qp,xinv_pRTs,xinv_pRT0,
     &             p_rx,p_ry,p_rz,crs0

          endif
        endif








        ELSE IF(imodel.eq.2)THEN

C     ...RESIDUE OF SPECTRAL FUNCTION
           CALL TAGRES(RES,dble(xinv_alpha))

C     ...FREE NUCLEON STRUCTURE FUNCTION (INPUT MODEL)
* ------------- (IPN for Weiss code (1):proton, (2):neutron )
 
           CALL TAGF2N(f2nn,dble(xinv_xbj),dble(xinv_q2),2)


C    ...CALCULATE PTR
C
c$$$            PR2 = (TP0 - TP)/2
c$$$            ER  = SQRT(PR2 + UN**2)
c$$$            PRZ = ALR*UD/2 - ER
c$$$            PTR = SQRT(PR2 - PRZ**2)
C


C     -------------------------------------------- reference 
c$$$         CALL TAGCD(FSIG,
c$$$     &          dble(xinv_se),dble(xinv_xbj),dble(xinv_q2),
c$$$     &          dble(xinv_alpha),dble(xinv_pRT0),IPN=2)
c$$$
C     -------------------------------------------- No Smearing 
c$$$         CALL TAGCD(FSIG,
c$$$     &          dble(xinv_se),dble(xinv_xbj),dble(xinv_q2),
c$$$     &          dble(xinv_alpha),dble(xinv_pRT),2)
C

C     -----------------------    On Smearing p_RT from tPrime
         CALL TAGCD(FSIG,
     &          dble(xinv_se),dble(xinv_xbj),dble(xinv_q2),
     &          dble(xinv_alpha),dble(xinv_pRT),2)



            E_recoil= sqrt(dble(xinv_tprim)/2+mprot**2)
       PHASR = PI*UD/4./dble(E_recoil)*dble(xinv_tprim)*dble(xinv_alpha)
C     
            ULUMI = 10E-6
            UNUM = ULUMI*dble(xinv_xbj)*dble(xinv_q2)*FSIG *PHASR
C
C    ...TAGGED STRUCTURE FUNCTION
C
            CALL TAGFD(F2D,FLD,
     &           dble(xinv_xbj),dble(xinv_q2),
     &           dble(xinv_alpha),dble(xinv_pRT),2)
C
C     ...SPECTRAL FUNCTION AND POLE PART
C
            CALL TAGSD(S,dble(xinv_alpha),dble(xinv_pRT))
C
         SPOL = RES/(dble(xinv_tprim))**2 
         F2DExt =F2D/SPOL

         si = FSIG
c         print *,xinv_se,xinv_tprim,si

       if(sqrt(xinv_se).gt.0.0)then
           if(si.gt.0 .and. xinv_tprim.gt.0)then

c              xjacobian = xinv_xbj/E_r
              xjacobian = xinv_jacob

              xkjfac = xinv_xbj/E_r/(xinv_tprim**(1/4.5))/2.62


 
c     factor of 2 should be remove later because I corrected the original C++ source code

              crs0 = FSIG
              crs1 = FSIG*dble(xinv_xbj)/E_recoil
     &             /(dble(xinv_tprim)**(1/4.5))/2.62
              crs2 = FSIG*dble(xinv_xbj)
     &             /sqrt(dble(xinv_tprim)/2+mprot**2)

              dum2 =0.


c$$$  INPUT for the GEMC format........................
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
c     xinv_se,xinv_sq,xinv_q2,q0,ktr,xinv_alpha,xinv_pRT,
              write(iu,101) nptl,
     &  dble(xinv_xbj),dble(xinv_q2),dble(xinv_se),dble(xebeam_pol),
     &  dble(xinv_alpha),dble(xinv_pRT),dble(xinv_tprim),dble(FSIG),
     &  dble(xibeam_pol)
              write(iu,11) egen,echarge,pid,pe_x,pe_y,pe_z,pe_E,pe_M,
     &             ve_x,ve_y,ve_z
              write(iu,11) sgen,scharge,pid1,p_rx,p_ry,p_rz,E_r,pm,
     &             v_rx,v_ry,v_rz 
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
c

c     dummy variable re-assign
              wn=1
              xmx2 =1

c$$$  INPUT for the ANALYSIS format........................
              write(iw,10) xinv_xbj,xinv_q2,xinv_pRT,xinv_alpha,
     &             FSIG,F2DExt,e_x,xjacobian,xkjfac,
     &             xinv_tprim,t_chn_MK,theta_qp,xinv_pRTs,xinv_pRT0,
     &             p_rx,p_ry,p_rz,crs0

          endif
        endif



         ENDIF
C

C     This is data format for batch mode analysis 
 10     format(18(e16.8,','))

C     This is data format for GEMC

 101    format('      ',I2,'  ',e10.4,'  ',
     &       e10.4,'  ',e10.4,'  ',e10.4,'  ',e10.4,'  ',
     &       e10.4,'  ',e10.4,'  ',e10.4,'  ',e10.4,'  ',
     &       e10.4)

 
 11    format('   ',I2,'   ',I2,'  1   ','   ',I4,'   0   ','   0  '
     &       ,e10.3,'   ',e10.3,'   ',e10.3,
     &       '  ',e10.3,'  ',e10.3,'  ',
     &       e10.3,'  ',e10.3,'  ',e10.3)

        ENDDO

        close(99)
        close(iu)
        close(iw)

        end
      








C======================================================================
C
C    LOCAL SUBROUTINE....MODULATED CROSS-SECTION ROUTINE
C
C======================================================================





C======================================================================
      function dot(v1,v2)
C----------------------------------------------------------------------
C-
C-   Purpose : Make a DOT product
C-
C-   Inputs  : v1(3) and v2(3) are vector (vx,vy,vz) 
C-
C-   Outputs : dot
C----------------------------------------------------------------------
      implicit none

	  real 	dot
	  real 	v1(3)
	  real 	v2(3)

	  integer i

	  dot = 0
	  do i = 1,3
	    dot = dot + v1(i)*v2(i)
	  enddo

      return
	  end
C ---------------------------------------------------------------------




C======================================================================
      subroutine extractang(th,thd,ph,phd,vec)
C----------------------------------------------------------------------
C-
C-   Purpose : Extract theta and phi from the vector
C-
C-   Inputs  : vec(3) are vector (vx,vy,vz) 
C-
C-   Outputs : th, thd, ph, phd
C----------------------------------------------------------------------
      implicit none

      real 	dot
      real 	vec(3)
      real 	th, thd, ph, phd
      real    vec_mag

      vec_mag   = sqrt(dot(vec,vec))

      th  = acos(vec(3)/vec_mag)
      thd = th*90./acos(0.)
      ph    = atan2(vec(2),vec(1))
      phd   = ph*90./acos(0.)
	  return
	  end
C ---------------------------------------------------------------------



C======================================================================
      SUBROUTINE rotate(idir,the,phi,p1,p2)    
C----------------------------------------------------------------------
C-
C-   Purpose : Rotation of a 3-vector: {x,y,z} <---> {x',y',z'}
C-
C-   Inputs  : the, phi  are angles of new z' axis defined
C-                       in master reference frame (x,y,z).
C-             idir = 1  Rotation from (x',y',z') to (x,y,z) reference frame.
C-             idir =-1  Rotation from (x,y,z) to (x',y',z') reference frame.
C-             p1(3)      is 3-vector {Px,Py,Pz} or any 3-vector {x,y,z} 
C- 
C-   Outputs : p2(3)      is 3-vector {Px,Py,Pz} or any 3-vector {x,y,z} 
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------
C
C Input variables
      integer idir
      REAL    the, phi, p1(3), p2(3)
C local variables
      INTEGER j
      REAL    ROT(3,3), pv(3)  
C
      IF(the**2+phi**2 .GT. 1E-20) THEN   
        ROT(1,1) =  COS(the)*COS(phi)  
        ROT(1,2) = -SIN(phi)  
        ROT(1,3) =  SIN(the)*COS(phi)  
        ROT(2,1) =  COS(the)*SIN(phi)  
        ROT(2,2) =  COS(phi)   
        ROT(2,3) =  SIN(the)*SIN(phi)  
        ROT(3,1) = -SIN(the)  
        ROT(3,2) =  0. 
        ROT(3,3) =  COS(the)
C   
        DO j = 1,3    
          pv(j) = p1(j)
        ENDDO    
        DO j = 1,3    
          IF(idir.GE.0) THEN
            p2(j) = ROT(j,1)*pv(1) + ROT(j,2)*pv(2) + ROT(j,3)*pv(3)
          ELSE
            p2(j) = ROT(1,j)*pv(1) + ROT(2,j)*pv(2) + ROT(3,j)*pv(3)
          ENDIF 
        ENDDO
      ENDIF
C
      RETURN    
      END   
C ---------------------------------------------------------------------



C======================================================================
      SUBROUTINE boost(BETA,P1,E1,P2,E2)    
C----------------------------------------------------------------------
C-
C-   Purpose : Makes Lorentz Boost (Lab <--> CM)
C-
C-   Inputs  : P1(3),E1 is four vector (Px,Py,Pz,E) 
C-             BETA is CM velocity
C-                  If BETA < 0. Boost from Lab to center mass system.
C-                  If BETA > 0. Boost from center mass to Lab system.
C- 
C-   Outputs : P2(3),E2
C----------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------
C
C Input variables
      real BETA, P1(3), E1, P2(3), E2
C
C Local variables  
      real GAMMA
C
      P2(1) = P1(1)
      P2(2) = P1(2)
      IF(BETA**2.GT.1E-20) THEN    
        GAMMA = 1.0 / SQRT(1.0 - BETA**2)  
        P2(3) = GAMMA*(P1(3) + E1*BETA)    
        E2 = GAMMA*(E1 + P1(3)*BETA) 
      ENDIF
C 
      RETURN    
      END   
C ---------------------------------------------------------------------













C======================================================================
*DECK TPMIN
      SUBROUTINE TPMIN(TP, ALR, ILIM)
C ---------------------------------------------------------------------
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
C---------------------------------------------------------------------










