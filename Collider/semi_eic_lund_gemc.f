      program crs_weight
c   test version-1.0
        common/par/pi,pm,pmp,pmn,dm,eb
        character (len=99) :: filename
        integer :: i, j, lines,IPN

c        real, dimension(999,20) :: aa
        real hvec0(8),h(8)

        real ee(11),evec1(8)
        real sp(11),hvec1(8)
        real zz(11),hvec2(8)
        real ww(11),hvec3(8)
        REAL qv(3),qv_mag
        REAL pfs(3),pfs_q(3) 
        real pfscm(3)
        real mprot,mneut,pd
        real R2D,alpha_em

        double precision :: f2nn,f2nn3,RES,FSIG,F2D,S,E_recoil
       double precision :: xkjfac,xjacobian,crs0,crs1,crs2
        double precision :: PI

        integer nptl,pid,pid1,pid2
        integer egen,sgen,xgen
        integer echarge,scharge,xcharge
        R2D = 180/3.14
        alpha_em = 0.0073

        ntgt =2
        nproton = 1

c     Some variable definition Christian's model
        ED = 0.00222
        UD = 2*mneut - ED
        PI = 0.314159265358979D+01


C     ###########################################################
C     USER CHOICE                                          BEGIN
C     ###########################################################

c     Initial electron beam energy
        ei =5.
c     Deuteron mass 
        dm = 1.875
c     Fixed target case  deuteron momentum = 0.0      \\        pd = 0.0
        pd = 100.0
        mprot =0.93827
        mneut =0.93955
c     Selection of model (M. Sargsian): 1,  (C. Weiss): 2
        imodel = 2
c     Call F2 parameterization for IMODEL=2,  (C.Weiss) : 1,  (P.Jimenez-Delgad): 3
        IPN =1
C     ###########################################################
C     USER CHOICE                                            END
C     ###########################################################


*********************************
* Initialization
*********************************	

        iu = 77
        write (filename, '( "ed_semi_eic",I2,".dat" )' )  iu
        OPEN(unit=iu,file=filename, status='new')
          

c  fsgen 100K event generated.....fora moment

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
c         DO ll=10,20

           read (99,*) h1,h2,
     &           h(1),h(2),h(3),h(4),h(5),h(6),h(7),h(8)

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

        DO iv=1,8
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
c$$$        pid2 = zz(2)
c$$$        xgen = 3
c$$$        xcharge = 0
c$$$        DO ihv2=4,11
c$$$           ihhv2 = ihv2 -3
c$$$           hvec2(ihhv2)= zz(ihv2)
c$$$        ENDDO
c$$$



        pe = sqrt(evec1(1)**2+evec1(2)**2+evec1(3)**2)
        pe_x = evec1(1)
        pe_y = evec1(2)
        pe_z = evec1(3)
        pe_E = evec1(4)
        pe_M = evec1(5)
        ve_x = evec1(6)
        ve_y = evec1(7)
        ve_z = evec1(8)
c        q0   = ei - pe
        q0  = 2*xinv_q2/2/dm/xinv_xbj


c theta_e unit : radian
        th_e = acos(evec1(3)/pe)
c fsgen : Q2 range set 0. - 10.GeV2
        q2 = 4.*ei*pe*sin(th_e/2.)**2
        qv_mag   = sqrt(q2 + q0**2)
        cs_the_qe = (ei - pe*cos(th_e))/qv_mag
   
        yy = q0/ei

c virtual photon component
        qx=-pe*cos(atan2(pe_y,pe_x))*sin(th_e)
        qy=-pe*sin(atan2(pe_y,pe_x))*sin(th_e)
        qz=ei-pe*cos(th_e)
        qv(1)  = qx
        qv(2)  = qy
        qv(3)  = qz
        qv_mag = sqrt(dot(qv,qv))

        call extractang(the_q,the_qd,phi_q,phi_qd,qv)
        beta  = sqrt(dot(qv,qv))/(q0+dm)

        ed =  sqrt(dm**2 + pd**2)

c invariant energy square of e + d system (k_e + p_d)**2
        se = 2.0*ei*(ed + pd) + dm**2
        sq = -q2 + 2.0*q0*ed + 2.0*qv_mag*cs_the_qe*pd + dm**2

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


c option nucleon -1: neutron, +1: proton
        ktr  = -1
c selection of bindidng model lc=1:light cone, 0: virtual nucleon
        lc = 1
c set of bjorken limit 0: no limit, 1: set value
        lbj = 1
c initialization 0: calculation, 1: PWIA
        icase = 1

c        alpha_r = 2*(E_r - p_rz)/(ed - pd) 
        alpha_r = xinv_alpha 


C missing mass2 distribution
         e_x=ei+dm-pe-sqrt(p_r**2+pm**2)

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


         P_v_x = hvec2(1)
         P_v_y = hvec2(2)
         P_v_z = hvec2(3)
         xE_v  = hvec2(4)
         xM_v  = hvec2(5)

         V_v_x = hvec2(6)
         V_v_y = hvec2(7)
         V_v_z = hvec2(8)

         P_v = sqrt(P_v_x**2+P_v_y**2+P_v_z**2)
         vmdum2 = E_v**2-P_v**2
         vmdum  = sqrt(vmdum2)


         tprime_KJ =  xinv_tprim
         t_chn_MK = tprime_KJ + pm**2

        IF(imodel.eq.1)THEN

          call edenx(xinv_se,xinv_sq,xinv_q2,q0,ktr,xinv_alpha,xinv_pRT,
     &    x,s,si,Fd_L,Fd_T,Fd_TL,Fd_TT,f2d,f1eff,f2eff,sun,icase,lc,lbj)


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
* ------------- (IPN for Weiss code (2):proton, (1):neutron )
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
            
c              crs1 = si/fac
c              xjacobian = xinv_xbj/E_r
              xjacobian = xinv_jacob

              crs0 = si
              crs1 = si*xjacobian
              crs2  = si*spect/fac
              dum2 =0.



c$$$  INPUT for the GEMC format........................
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
c$$$c     xinv_se,xinv_sq,xinv_q2,q0,ktr,xinv_alpha,xinv_pRT,

              write(iu,10) nptl,ntgt,nproton,
     &             x,xinv_q2,xinv_alpha,xinv_pRT,tprime_KJ,crs1,dum2
              write(iu,11) egen,echarge,pid,pe_x,pe_y,pe_z,pe_E,pe_M,
     &             ve_x,ve_y,ve_z
              write(iu,11) sgen,scharge,pid1,p_rx,p_ry,p_rz,E_r,pm,
     &             v_rx,v_ry,v_rz
c$$$c     remove for the GEMC simulation input
c$$$             write(iu,11) xgen,xcharge,pid2,P_v_x,P_v_y,P_v_z,xE_v,
c$$$     &             xM_v,V_v_x,V_v_y,V_v_z
c$$$ 
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
c

c     dummy variable re-assign
              wn=1
              xmx2 =1
c$$$
c$$$              write(iu,10) xinv_xbj,xinv_q2,xinv_pRT,xinv_alpha,
c$$$     &             si,p_x2,e_x,xmx2,wn,
c$$$     &             xinv_tprim,t_chn_MK,theta_qp,crs1,crs2,
c$$$     &             p_rx,p_ry,p_rz,yy
c$$$ 
          endif
        endif


c$$$  INPUT for the GEMC format........................
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
c$$$
c$$$ 10     format('      ',I2,'  ',I2,'  ',I2,'  ',f10.4,'  ',
c$$$     &       f12.6,'  ',f12.6,'  ',f12.6,'  ',e14.8,'  ',e14.8,'  ',e14.8)
c$$$ 
c$$$ 11    format('   ',I2,'   ',I2,'  1   ','   ',I4,'   0   ','   0  '
c$$$     &       ,e14.8,'   ',e14.8,'   ',e14.8,
c$$$     &       '  ',e14.8,'  ',e14.8,'  ',
c$$$     &       e14.8,'  ',e14.8,'  ',e14.8)
c$$$
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  








        ELSE IF(imodel.eq.2)THEN

         IF(xinv_xbj.gt.1. .OR. xinv_xbj.lt.0.)RETURN
    
 

C     ...RESIDUE OF SPECTRAL FUNCTION
           CALL TAGRES(RES,(xinv_alpha))


C     ...FREE NUCLEON STRUCTURE FUNCTION (INPUT MODEL)
* ------------- (IPN for Weiss code (1):proton, (2):neutron )
 
           CALL TAGFN(f2nn,dble(xinv_xbj),dble(xinv_q2),IPN)




C    ...CALCULATE PTR
C
c$$$            PR2 = (TP0 - TP)/2
c$$$            ER  = SQRT(PR2 + UN**2)
c$$$            PRZ = ALR*UD/2 - ER
c$$$            PTR = SQRT(PR2 - PRZ**2)
C

         CALL TAGX(FSIG,
     &          (xinv_se),dble(xinv_xbj),dble(xinv_q2),
     &          (xinv_alpha),(xinv_pRT),IPN)
C
            E_recoil= sqrt((xinv_tprim)/2+mprot**2)
       PHASR = PI*UD/4./(E_recoil)*(xinv_tprim)*(xinv_alpha)
C     
            ULUMI = 10E-6
            UNUM = ULUMI*(xinv_xbj)*(xinv_q2)*FSIG *PHASR
C
C    ...TAGGED STRUCTURE FUNCTION
C
            CALL TAGFD(F2D,
     &           dble(xinv_xbj),dble(xinv_q2),
     &           (xinv_alpha),(xinv_pRT),2)
C
C     ...SPECTRAL FUNCTION AND POLE PART
C
            CALL TAGSP(S,(xinv_alpha),(xinv_pRT))
C
         SPOL = RES/((xinv_tprim))**2 
         F2DExt =F2D/SPOL

         si = FSIG

       if(sqrt(xinv_se).gt.0.0)then
           if(si.gt.0 .and. xinv_tprim.gt.0)then

c              xjacobian = xinv_xbj/E_r
              xjacobian = xinv_jacob

              xkjfac = xinv_xbj/E_r/(xinv_tprim**(1/4.5))/2.62



              CALL TPMIN(xinv_tpmin,(xinv_alpha),1)

c     factor of 2 should be remove later because I corrected the original C++ source code

              crs0 = FSIG
              crs1 = FSIG*dble(xinv_xbj)/E_recoil
     &             /(dble(xinv_tprim)**(1/4.5))/2.62
              crs2 = FSIG*dble(xinv_xbj)
     &             /sqrt(dble(xinv_tprim)/2+mprot**2)
              dum2 =0.


c              print *,xjacobian,xinv_jacob

c$$$  INPUT for the GEMC format........................
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
c$$$c     xinv_se,xinv_sq,xinv_q2,q0,ktr,xinv_alpha,xinv_pRT,

              write(iu,10) nptl,ntgt,nproton,
     &             x,xinv_q2,xinv_alpha,xinv_pRT,tprime_KJ,crs1,dum2
              write(iu,11) egen,echarge,pid,pe_x,pe_y,pe_z,pe_E,pe_M,
     &             ve_x,ve_y,ve_z
              write(iu,11) sgen,scharge,pid1,p_rx,p_ry,p_rz,E_r,pm,
     &             v_rx,v_ry,v_rz
c$$$c     remove for the GEMC simulation input
c$$$             write(iu,11) xgen,xcharge,pid2,P_v_x,P_v_y,P_v_z,xE_v,
c$$$     &             xM_v,V_v_x,V_v_y,V_v_z
c$$$ 
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
c

c     dummy variable re-assign
              wn=1
              xmx2 =1
c$$$
c$$$              write(iu,10) xinv_xbj,xinv_q2,xinv_pRT,xinv_alpha,
c$$$     &             FSIG,p_x2,e_x,xjacobian,xkjfac,
c$$$     &             xinv_tprim,t_chn_MK,theta_qp,crs1,crs2,
c$$$     &             p_rx,p_ry,p_rz,crs0


          endif
        endif



         ENDIF
C

c$$$
c$$$ 10     format(18(e16.8,','))


c$$$  INPUT for the GEMC format........................
c$$$c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  

 10     format('      ',I2,'  ',I2,'  ',I2,'  ',f10.4,'  ',
     &    f12.6,'  ',f12.6,'  ',f12.6,'  ',e14.8,'  ',e14.8,'  ',e14.8)
 
 11    format('   ',I2,'   ',I2,'  1   ','   ',I4,'   0   ','   0  '
     &       ,e14.8,'   ',e14.8,'   ',e14.8,
     &       '  ',e14.8,'  ',e14.8,'  ',
     &       e14.8,'  ',e14.8,'  ',e14.8)

c XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  





        ENDDO

        close(99)
        close(iu)

        end
      


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




**************************************************************
      subroutine edenx(se,sq,q2,q0,ktr,alpha_r,p_rt,x,s,si,
     &       Fd_L,Fd_T,Fd_TL,Fd_TT,f2d,f1eff,f2eff,sun,icase,lc,lbj)
************************************************************************
* edenX v.3
*        edenx calculates the differential cross section of 
*        e+d -> e+n+X reaction, in which N is detected in the target 
*        fragmentation region
*
*        INPUTS
*        Se    - invariant energy square of e and d system (k_e + p_d)^2
*        Sq    - invariant energy square of q and d system  (q + p_d)^2
*        q2    - Q2
*        q0    - transfered energy
*        alpha_r =  2*(E_r - p_rz)/(E_d - p_dz) 
*        E_r - energy of recoil particle
*        p_rz - z component of the recoil nucleon
*        E_d  - energy of the deuteron
*        p_dz - momentum of the deuteron
*        z - axis is defined by the direction of the virtual photon
*        p_rt - transverse momentum  of recoil nucleon vs q
*
*        ktr   - type of n (-1 - neutron, 1 - proton)
*        icase - initialization (0) - calculation 1 - PWIA, 
*        2-FSI, 12 - IA+FSI (not active yet)
*        lc    - 0 - virtual nucleon , 1 - light cone
*        lbj   - 0 - no Bjorken limit, Bjorken Limit 
*
*        OUTPUT
*        s   - cross secetion in dsigma/dx/dQ2 d^3p_s/E_s  - in nb/GeV^4  (with azimuthal dependence not active yet)
*        si  - integrated cross section by recoil nucleon's azimuthal angle /(2pi)
*                dsigma/dx/dQ2 dalpha_r/alpha_r, d^2p_t 
*        x     - bjorken x  (output) = trivial x

*        F2d - F2 structure function of the Deuteron (not active yet)


* Dec    2013 -  updated for colloder kinematics
* July - 2004
* June - 2005
* FIU, Miami
* 
***************************************************************************
	implicit none
	real aem,pi,pm,pmp,pmn,dm,eb
	real ei,er,q0,qv,q2,x,y,sin_the_2,the
	real pr,thr,phir,Tr,e_r,al_r,pt_r,wn2,pr_z,p_rz,p_rt
	real term0,term,a2,a3,asi,s_a,s_b
	integer icase,lc,ktr,lbj
	real sun,s,si,f2d,Fd_L,FD_T,FD_TL,FD_TT,F2d_si,F1d_si
        real f1eff,f2eff
        real se,sq
        real alpha_r
	common/aem/aem
        common/par/pi,pm,pmp,pmn,dm,eb
	s   = 0.0
	si  = 0.0
	f2d = 0.0
***********************************************
*         Initialization
************************************************
	if(icase.eq.0)then
	aem = 1.0/137.0
        pi  = acos(-1.0)       
        pm  = 0.938272         
        dm  = 1.875628
        pmp = pm
        pmn = 0.939565  
	eb  = 0.0022245
	call spectral_un(q0,qv,pr,thr,phir,ktr,Tr,sun,lc,icase)
	return
	endif
***************************************************
*        Start of calculations
***************************************************
*	q0 = q2/2.0/pm/x
	qv = sqrt(q2 + q0**2)
*	y  = q0/ei
        y  = sq/se
        x  = 2.0*q2/(sq-dm**2+q2)
*        print *,"x",x
*	er = ei-q0
*	print *,er
*	if(er.le.0.0)return
*	sin_the_2 = sqrt(q2/4.0/ei/er)		
*	if(sin_the_2**2.gt.1.0)return
*	the = 2.0*asin(sin_the_2)
	
*        print *,"yx",y,x,sq,se,dm,q2

*         print *, se,sq,q2,q0,ktr,alpha_r,p_rt

*************************************************
*        Recoil nucleon kinematics
*************************************************
        p_rz = (pm**2*(1.0-alpha_r**2) + p_rt**2)/2.0/pm/alpha_r
        pr = sqrt(p_rt**2 + p_rz**2)
	e_r  = sqrt(pm**2 + pr**2)

        Tr = er - pm
* 	al_r = (e_r - pr*cos(thr))/pm
*	pt_r = pr*sin(thr)
*        print *,"yx",y,x,sq,se,p_rz,pr,e_r


*****************************************************
*    Checking the produced mass
*****************************************************
*	wn2 = 
*     &  dm**2-2.0*dm*e_r+pm**2+2.0*(dm-e_r)*q0 + 2.0*pr*qv*cos(thr)-q2
*	if(wn2.lt.pm**2)return



******************************************************************
*   {4\pi \alpha^2\over (x Q^4)} (1-y - {x^2 y^2 m^2\over Q^2})
******************************************************************
	term0 = 2.0*pi*2.0*aem**2/x/q2**2 * (1-y-(x*y*pm)**2/q2)
	term  = term0*0.389379*1000.*1000.   !nb/GeV^2

*	a2  = (q2/2.0/qv**2+tan(the/2.0)**2)*q0/pm
*	a3  = sqrt(q2/qv**2+tan(the/2.0)**2)
*	asi = 2.0*tan(the/2.0)**2*q0/pm
****************************************************************************
*     Structure Functions calculations
****************************************************************************	
        phir = 0.0
	call 
     &d_structure_func(ktr,q2,x,alpha_r,p_rt,phir,Tr,
     &Fd_L,Fd_T,Fd_TL,Fd_TT,
     &F2d_si,F1d_si,f1eff,f2eff,sun,icase,lc,lbj)



*	s_a  = Fd_L+ a2*Fd_T + a3*cos(phir)*Fd_TL + cos(2.0*phir)*Fd_TT
	s_b  = F2d_si !nor + asi*F1d_si
*	s_b  = Fd_L + a2*Fd_T
	f2d  = F2d_si

        s = 0.0
*	s    = term* s_a
	si   = term* s_b

*        print *,"si",term,s_b,si
*        print *, "fd_l",fd_tl
	return
	end


	subroutine 
     &  d_structure_func(its,q2,x,als,pt,phir,Tr,
     &                  Fd_L,Fd_T,Fd_TL,Fd_TT,
     &                  F2d_si,F1d_si,f1eff,f2eff,sun,icon,lc,lbj)
********************************************************************************************
*  Subroutine calculates deuteron DIS structure functions for semiinclusive scattering
*
*  its    - type of the recoil nucleon -1 - neutron 1 - proton
*  q2     - Q2
*  x      - Bjorken x
*  als    - alpha of the recoil nucleon
*  pt     - transverse momentum of the recoil nucleon
*  phit   - azimuthal angle of the recoil nucleon
*  Tr     - Kinetic energy of the recoil nucleon
*  Fd_L   - FL  structure function of the deuteron
*  Fd_T   - FT  structure function of the deuteron
*  Fd_TL  - FTL structure function of the deuteron
*  FD_TT  - FTT structure function of the duteron
*  F2d_si - F2  structure function for azimuthal averaged cross section
*  F1d_si - F1  structure function for azimuthal averaged cross section
*  icon   - caseses 0 - initialization, 1- PWIA, 2 - FSI, 12-IA+FSI 
*  lc     - 0 - virtual nucleon, 1 - light cone
*  lbj    - 0 - no Bj limit, 1 - Bjorken limit
*
*
* Misak Sargsian
* July 2004
* June 2005
* FIU, TUN, Miami
**********************************************************************************************
	implicit none
	real aem,pi,pm,pmp,pmn,dm,eb
	real q2,q0,qv,x,sin_del,cos_del,alq
	real als,pt,phir,Tr,psz,ps,es,al,pr,thr
	real pq,pplus,qplus,xtil,starm2,w2n,q0_off
	integer icon,lc,its,lbj,itn,ktr
	real FD_L_N,FD_T_N,FD_TL_N,FD_TT_N,F2d_si_N,F1d_si_N
	real Fd_L,FD_T,FD_TL,FD_TT,F2d_si,F1d_si,F1eff,F2eff
	real f_n,sun,sunflux
        common/par/pi,pm,pmp,pmn,dm,eb
	Fd_L   = 0.0
	Fd_T   = 0.0
	Fd_TL  = 0.0
	Fd_TT  = 0.0
	F2d_si = 0.0
	F1d_si = 0.0

	if(als.le.0.0.or.als.ge.2.0)return
	q0 = q2/2.0/pm/x
	qv = sqrt(q2+q0**2)

*	sin_del = sqrt(q2)/qv
*	cos_del =       q0/qv
*	alq = (q0-qv)/pm    !\alpha_q
	psz = (pm**2 + pt**2 - als**2*pm**2)/(2.0*als*pm)
	ps  = sqrt(pt**2 + psz**2)
	es = sqrt(pm**2 + ps**2)
	al = 2.0-als  ! momentum fraction of struck nucleon


*        print *,"**pt,psz**",pt,psz
***************************************************************
*        definining some off-shell and light cone quantities
****************************************************************
	if(lc.eq.0)then     ! virtual nucleon approximation
	pq   = (dm-es)*q0 + psz*qv
	elseif(lc.eq.1)then ! Light cone approximation
	pplus = dm - (pm**2+pt**2)/(pm*als)
	qplus = -Q2/(pm*alq)
	pq    = (1.0/2.0)*(pplus*pm*alq + qplus*pm*al) 
	endif



*	xtil = x/(2.0-als)   !\tilde x in Bjorken limit
*	if(pq.le.0.0)return
*        xtil = q2/(2.0*pq)
*	starm2 = (dm-es)**2 - ps**2
*	w2n   = -q2 + 2.0*pq + starm2  ! final produced mass^2
*        write(21,*)"wn**",sqrt(w2n)
*	diff = w2n - starm2
*********************************************************************
*    Since Bodek parameterization is defined by the final masses
*    we define \tilde x using the final masses 
*********************************************************************	
*	q0_off = (w2n-pm**2+q2)/(2.0*pm)
*	            xtil = q2/(2.0*pm*q0_off)
	if(lbj.eq.1)xtil = x/al
	if(xtil.le.0.0.or.xtil.ge.1.0)return

*___________________________________________________________
* calculation of effective structure functions of nucleons
*___________________________________________________________
	itn   = -its
	call struct_func_nucleon(itn,xtil,q2,al,pt,F1eff,F2eff)
*	f2eff = 1.0
*	print *,"xof",xtil,F2eff
*===========================================================
	if(lbj.eq.0)then !non Bjorken limit
*===========================================================
************************************************************
*   four structure functions for semiinclusive scattering
************************************************************
*	FD_L_N  = -sin_del**2*(q0/pm)*F1eff + 
*     &            (1.0+cos_del)**2*(al+pq/q2*alq)**2*(q0*pm)/pq*F2eff
*	FD_T_N  = 2.0*F1eff + pt**2/pq*F2eff
*	FD_TL_N = 2.0*(1.0+cos_del)*pt/pm*(al+pq/q2*alq)*(q0*pm)/pq*F2eff
*	FD_TT_N = sin_del**2/2.0 *(pt**2/pm**2)* q0*pm/pq *F2eff

*	print *,"aha",xtil,Fd_L_N,f1eff,f2eff,itn


**************************************************************
* Two structure functions - averaged over phi
**************************************************************
*	F2d_si_N = ((1.0+cos_del)**2*(al+pq/q2*alq)**2 
*     &              + sin_del**2/2.0* (pt**2/pm**2))
*     &             *(q0*pm)/pq * F2eff
*	F1d_si_N = F1eff + pt**2/(2.0*pq)*F2eff
*==========================================================
	elseif(lbj.eq.1)then ! Bjorken limit
*===========================================================
************************************************************
* four structure functions for semiinclusive scattering
************************************************************
	FD_L_N  = 0.0 !-2.0*x*F1eff + al*F2eff
	FD_T_N  = 0.0 ! 2.0*F1eff
	FD_TL_N = 0.0! 2.0*pt/pm*F2eff
	FD_TT_N =  0.0
**************************************************************
* Two structure functions - averaged over phi
**************************************************************
	F2d_si_N = al* F2eff
	F1d_si_N = F1eff
	endif
*______________________________________________________________
*   Calculating the spectral function
*______________________________________________________________
	pr = ps
	             thr = 0.0
	if(pr.ne.0.0)thr = acos(psz/ps)
	ktr = its
	call spectral_un(q0,qv,pr,thr,phir,ktr,Tr,sun,lc,icon)	
****************************************************************
*    Calculates 1/n factrod from polext paper
****************************************************************
	           f_n = al
	if(lc.eq.0)f_n = 2.*(dm-es)/dm
*	print *, "mis",es,als,ps
	sunflux = sun/f_n

	Fd_L   = FD_L_N   * sunflux
	Fd_T   = FD_T_N   * sunflux
	Fd_TL  = FD_TL_N  * sunflux
	Fd_TT  = FD_TT_N  * sunflux
	F2d_si = F2d_si_N * sunflux
	F1d_si = F1d_si_N * sunflux
	return
	end



	subroutine struct_func_nucleon(itn,x,q2,al,pt,F1,F2)
	implicit none
	real pi,pm,pmp,pmn,dm,eb
	real q0,q2,x,w2
	real al,pt
	integer itn
	real F2a,F1a,F1,F2
        common/par/pi,pm,pmp,pmn,dm,eb
	q0 = q2/(2.0*pm*x)
	w2 = pm**2 + 2.0*pm*q0 - q2 
	if(itn.eq.-1)then
	call struct_bodek_n(q0,q2,w2,F2a,F1a)
*	f2a = 1.
	elseif(itn.eq.1)then
	call struct_bodek_p(q0,q2,w2,F2a,F1a)
	endif
	F1 = F1a
	F2 = F2a
	return
	end

	


	function spectral(als,pt,icase)
************************************************************************
*         Structure function of the Deuteron in GeV^-3
*         als - Light Cone Momentum Fraction of Spectator Nucleon
*         pt  - magnitude of the transverse momentum of the spectator
*         icase - (1 - virtual Nucleon, 2 - LC)
*************************************************************************
	common/param/pm,pi,dm
	spectral = 0.0
	if(als.le.0.0.or.als.ge.2.0)return

	if(icase.eq.1)then    ! VN
*	psz = (pm**2 + pt**2 - als**2*(dm/2.0)**2)/(2.0*als*(dm/2.0))
	psz = (pm**2 + pt**2 - als**2*pm**2)/(2.0*als*pm)
	ps  = sqrt(pt**2 + psz**2)
	es = sqrt(pm**2 + ps**2)
	if(es.ge.dm)return
*	factor = (2.0-als)*es*dm/(2.0*(dm-es))
	factor = es/als*dm/(2.0*(dm-es)) 
	spectral = factor *fd(ps,0)
	elseif(icase.eq.2)then    !LC
	pk = sqrt((pm**2 + pt**2)/(als*(2.0-als))-pm**2)
	ek = sqrt(pm**2 + pk**2)
*	spectral = ek/(2.0-als)*fd(pk,0)
	spectral = ek/als*fd(pk,0)
	endif
	return
	end
	


	subroutine  spectral_un(q0,qv,pr,thr,phir,ktr,tr,sun,lc,icon)
*******************************************************************
*         Spectral Function
*         q0 - virtual photon energy
*         qv - virtual photon momentum
*         pr - recoil nucleon momentum
*         thr - recoil nucleon polar angle vs q
*         phir - recoil nucleon azimuthal angle
*         ktr - type of the recoil nucleon 1- proton -1 - neutron
*         tr  - recoil nucleons kinetic energy
*         sun  - spectral function 
*         lc  - 0 virtual nucleon, 1 light cone  approximation
*         icon- 0-initialization, 1 PWIA, 2 - FSI, 12- PWIA+FSI
*
* Misak Sargsian
* June 2005
* FIU, TUN, Miami
* 
* Sep 24 - truncated delta
********************************************************************
	implicit none
	real pi,pm,pmp,pmn,dm,eb
	real q0,qv
	real pr,thr,phir,tr,ts
	real s,spec,sun
	integer kj,ktm,ktr,ifsi,lc,icon
        common/par/pi,pm,pmp,pmn,dm,eb
	common/kinetic/ts,ifsi
	sun   = 0.0
	if(icon.eq.0)then
        call dmatrix_x(q0,qv,kj,ktm,ktr,pr,thr,phir,spec,lc,icon)
        return
        endif
************************************** 
*       initial  nucleon's isospin
**************************************
        ktm =-ktr 
	ts  = tr
        s   = 0.0
        do kj = -1,1,1
        call dmatrix_x(q0,qv,kj,ktm,ktr,pr,thr,phir,spec,lc,icon)
        s = s + spec
        enddo

        sun  = s/3.0

        return
        end




       subroutine dmatrix_x(q0,qv,kj,ktm,ktr,pr,thr,phir,spec,lc,icon)
**********************************************************************
* dmatrix_x - v.1.00
* Density matrix of the deuteron for d(e,e'Nr)X reaction
*
* q0   - virtual photon energy
* qv   - virtual photon momentum
* kj   - total spin projectsion of d2 target 1, 0, -1             - input
* ktm  - isospin projection of   m  nucleon 1->1/2,-1 -1/2        - input
* ktr  - isospin projection of   r  nucleon 1->1/2,-1 -1/2        - input  
* pr   - momentum of             r  nucleon  GeV/c                - input
* thr  - polar angle             r  nucleon                       - input
* phir - azimuthal angle         r  nucleon                       - input
*
* m - knocked out
* r - recoil 
*
* spec - the value of the density matrix                          - output
* icon - initialization (0),  PWIA (-1), FSI (2), PWIA+FSI (12)
*
* The azimuthal angles are defined in the reference frame where 
* q || z and xz is the (eq), e is the initial electron
*
*
* Misak Sargsian
* FIU, December, 2003
*
*
********************************************************************** 
	implicit none
	real pi,pm,pmp,pmn,dm,eb
	real q0,qv
	real pr,thr,phir,ts,p_r,p_r_x,p_r_y,p_r_z,p_r_t,e_r,pvec_r
	real p_m,thm,phim
	real als,pk,ek,pkz,prk,thrk,phirk
	real wf_re,wf_im,wf0_re,wf0_im,wf1_re,wf1_im
	real xxx,fd
	real fct0,fctm
	real spec0,spec
	integer kj,ktm,ktr,lc,lcone,icon,ifsi
	integer ksm,ksr
	

        common/par/pi,pm,pmp,pmn,dm,eb
        dimension pvec_r(3)
	common/vnlc/lcone
	common/kinetic/ts,ifsi
        spec  = 0.0
 	if(icon.eq.0)then
****************************************************************
*       initialization 
****************************************************************
        pi  = acos(-1.0)       
        pm  = 0.938279         
        dm  = 1.875628
        pmp = pm
        pmn = 0.939565  
	eb  = 0.0022245
        xxx = fd(0.,1)                    ! initialize the Paris WF
****************************************************************
*       end initialization 
****************************************************************
        return
        endif

****************************************************************
*       INPUT      
****************************************************************
	lcone = lc
************************************** 
*       spectator/recoil nucleon
*       angles are  with respect to q
**************************************
	p_r    =  pr
 	p_r_x  =  p_r * sin(thr)*cos(phir) 
        p_r_y  =  p_r * sin(thr)*sin(phir)
        p_r_z  =  p_r * cos(thr)
	p_r_t  = sqrt(p_r_x**2 + p_r_y**2)
	e_r    = sqrt(pm**2 + p_r**2)
        pvec_r(1) =  p_r_x
        pvec_r(2) =  p_r_y
        pvec_r(3) =  p_r_z

        p_m  = p_r
        thm  = pi - thr
        phim = pi + phir
     

	fct0 = e_r !* dm/(2.0*(dm-e_r))
	fctm = 1.0 !sqrt((dm-e_r)/e_r)

	if(lc.eq.1)then
****************************************************
*  Light Cone approximation
****************************************************
	als = (e_r-p_r_z)/pm
	pk  = sqrt((pm**2+p_r_t**2)/(als*(2.0-als))-pm**2)
	ek  = sqrt(pm**2 + pk**2)
	pkz = ek*(1.0-als)
  	call polar(p_r_x,p_r_y,pkz,prk,thrk,phirk)
        p_m  = prk
        thm  = pi - thrk
        phim = pi + phirk

	fct0 = ek/(2.0-als)
	fctm = sqrt(pm/ek)
	endif


      spec0 = 0.0
      do ksm = -1,1,2        
        do ksr = -1,1,2
         wf0_re  = 0.0
         wf0_im  = 0.0
         wf1_re  = 0.0
         wf1_im  = 0.0
         wf_re   = 0.0
         wf_im   = 0.0
	 
      if(icon.eq.1.or.icon.ge.12)then    ! A0  - contribution
      ifsi = 0
      call dfunction(kj,ktm,ksm,ksr,p_m,thm,phim,wf0_re,wf0_im)
      endif

      if(icon.eq.2.or.icon.eq.12)then    ! A1 - contribution
      ifsi = 1
      call A1_fct(q0,qv,ktm,ksm,ktr,ksr,pvec_r,kj,wf1_re,wf1_im)
*	print *, "wfs",wf1_re,wf1_im
      endif
	
      wf_re = wf0_re + wf1_re  
      wf_im = wf0_im + wf1_im  
      spec0 = spec0 +  wf_re**2 + wf_im**2
      enddo
      enddo
*	print *,p_r,wf1_re,wf1_im
      spec =  fct0*spec0
      return
      end


      subroutine A1_fct(q0,qv,ktm,ksm,ktr,ksr,p_rv,kj,wf1_re,wf1_im)
******************************************************************
* inputs 
* q0  - virtual photon energy
* qv  - virtual photon momentum
* ktm - isospin projection of m nucleon 1 - proton, -1 neutron
* ksm - spin projection of m  nucleon 1 -> (1/2)  -1 ->(-1/2)
* ktr - isospin  projection of recoil nucleon 1 - proton, -1 neutron
* ksr - spin projection of recoil nucleon   1 -> (1/2)  -1 ->(-1/2)
* p_r - x,y,z components of recoil nucleon momentum
* kj  - spin projection of d2 -1,0,1 
*
* wf1_re
* wf1_im
*******************************************************************
	implicit none
	real pi,pm,pmp,pmn,dm,eb
	real phl,phu,un_1_re,un_1_im
	real q0,qv,Q2,q0s,qvs,qpl,qmn
	real p_m_x,p_m_y,p_m_z
	real p_rv, p_r,E_r,p_r_x,p_r_y,p_r_z,p_r_t,T_r
	real als,pk,ek,pkz
	real eps,qa,qb,sum_re,sum_im
	real wf1_re,wf1_im
	real s,delta_0,starm2,w_n_2,w_n0_2,del_w2
	integer ktm,ksm,ktr,ksr,kj
	integer itm,ism,itr,isr,ij,lc
	
        external phl, phu, un_1_re, un_1_im
        dimension p_rv(3)
        common/par/pi,pm,pmp,pmn,dm,eb
	common/vphoton/Q2,q0s,qvs
        common/missing/  p_m_x,  p_m_y,  p_m_z
        common/recoil_r/ p_r_x,  p_r_y,  p_r_z
        common/conf/itm,ism,itr,isr,ij
        common/enginv/s
        common/delta0/delta_0
	common/vnlc/lc
	wf1_re = 0.0
	wf1_im = 0.0
	
        itm    = ktm
        ism    = ksm
        itr    = ktr
        isr    = ksr
        ij     = kj
        p_r_x =  p_rv(1)
        p_r_y =  p_rv(2)
        p_r_z =  p_rv(3)

	q2  = qv**2 - q0**2
	q0s = q0
	qvs = qv


* W_n  - produced on missing nucleon with momentum pm
**********************************************************
        p_r  = sqrt(p_r_x**2 + p_r_y**2 + p_r_z**2)
        E_r = sqrt(pm**2+p_r**2)

        starm2 = dm**2 - 2.0*E_r*dm + pm**2 
	w_n_2   = starm2 + 2.0*q0*(dm-E_r) + 2.0*qv*p_r_z - Q2
	w_n0_2  = pm**2  + 2.0*q0*pm - Q2
*	w_n0_2  = 2.0*q0*pm - Q2
        T_r = E_r - pm
	del_w2 = (w_n_2 - w_n0_2)
*	del_w2 =  w_n_2 - pm**2
	if(del_w2.lt.0.0)del_w2=0.0
        delta_0 =  T_r*(dm+q0)/qv  + del_w2/(2.0*qv)
*	delta_0 = 0.0
*        delta_0 =  T_r*(dm+q0)/qv  + (w_n_2 -pm**2)/(2.0*qv)
*	print *,"delta_0",delta_0,w_n_2,w_n0_2,qv,p_r_z
	if(lc.eq.1)then
****************************************************
*  Light Cone approximation
****************************************************
	p_r_t = sqrt(p_r_x**2 + p_r_y**2)
	als = (E_r-p_r_z)/pm
	pk  = sqrt((pm**2+p_r_t**2)/(als*(2.0-als))-pm**2)
	ek  = sqrt(pm**2 + pk**2)
	pkz = ek*(1.0-als)
	p_r_z = pkz
	qpl = q0+qv
	qmn = q0-qv
	delta_0 = (dm+qmn)/qpl*(p_r_t**2/(pm*als)) !nor  + 
!     & + (w_n_2 - w_n0_2)/(qpl)
	endif


        p_m_x =  -p_r_x
        p_m_y =  -p_r_y
        p_m_z =  -p_r_z

********************************************************
*	print *,delta_0
**************************************************
* This block calculates the s of the final NN system
* We calculate it using the fact that s = (m_d+q)^2 = m_d^2 + 2.0*m_d*q0 - Q2
**************************************************
        s = dm**2 + 2.0*dm*q0 - Q2
***************************************************************

        
        eps = 0.0001
        qa  = 0.0
        qb  = 1.0
        sum_re = 0.0
        call gadap2(qa,qb,phl,phu,un_1_re,eps,sum_re) 
        eps = 0.0001
        qa  = 0.0
        qb  = 1.0
        sum_im = 0.0
  	call gadap2(qa,qb,phl,phu,un_1_im,eps,sum_im)                

        wf1_re = sum_re
        wf1_im = sum_im         

*        write(6,*)"sums",sum_re,sum_im
        return
        end

        function phl(q)
	phl = 0.0
	return
	end
	function phu(q)
	phu = 2.0*acos(-1.0)
	return
	end

	function  un_1_re(q,phiq)
	implicit none
	real q,phiq
	real pi,pm,pmp,pmn,dm,eb
	real Q2,q0,qv
	real  p_m_x, p_m_y, p_m_z
	real  p_r_x, p_r_y, p_r_z,p_r,E_r
	real  p_1_x, p_1_y, p_1_z,p1,th1,phi1
	real  p_r1_x, p_r1_y, p_r1_z,p_r1,E_r1
	real s,t,delta_0,pr1pr
	real q_x,q_y
	real fct,un_1_re,un_1_re_pole,un_1_re_pv
	real wfre,wfim,awfre,awfim
	real f_XN_on_re, f_XN_on_im,F_XN_off_re,F_XN_off_im
	integer ktm,ksm,ktr,ksr,kj,lc
	real mx,mx2
        real xi,tnorm
	real wn2,E_r_pr,del_w2,delta_pr

        common/par/pi,pm,pmp,pmn,dm,eb
	common/vphoton/Q2,q0,qv
        common/missing/  p_m_x,  p_m_y,  p_m_z
        common/recoil_r/ p_r_x,  p_r_y,  p_r_z
        common/conf/ktm,ksm,ktr,ksr,kj
        common/enginv/s
        common/delta0/delta_0
	common/vnlc/lc


        un_1_re = 0.0
	q_x      = q*cos(phiq)
	q_y      = q*sin(phiq)

        p_r  = sqrt(p_r_x**2 + p_r_y**2 + p_r_z**2)
        E_r = sqrt(pm**2+p_r**2)
	fct = 1.0
	if(lc.eq.1)fct = pm/E_r




 	p_1_z   = fct*p_m_z !+ delta_0
	p_1_y   = p_m_y + q_y
	p_1_x   = p_m_x + q_x

*  	call polar(p_1_x,p_1_y,0.0,p1,th1,phi1)

**************************************************
* calculating delta_pr - delta that depends on pr
**************************************************
	p1 = sqrt(p_1_x**2 + p_1_y**2)
	p_1_z = 0.0
	E_r_pr = sqrt(pm**2 + p1**2)	
*	E_r_pr = pm
	mx2 =  dm**2 - 2.0*E_r_pr*dm + pm**2 +
     &         2.0*q0*(dm-E_r_pr) - 2.0*qv*p_1_z - Q2

*	mx2 = pm**2  + 2.0*q0*pm - Q2
        wn2 =  dm**2 - 2.0*E_r*dm    + pm**2 +  
     &	       2.0*q0*(dm-E_r)    + 2.0*qv*p_r_z - Q2

	if(mx2.le.0.0.or.wn2.lt.pm**2)return
	del_w2 = wn2-mx2
	if(del_w2.lt.0.0)del_w2=0.0

        delta_pr =  (E_r-E_r_pr)*(dm+q0)/qv  + del_w2/(2.0*qv)
	
*	print *,"deltas",delta_pr,delta_0
	p_1_z = p_m_z + delta_pr


*	print *,p_m_x,p_m_y
***************************************************
*  Recoil Nucleon's momentum and energy in 
*  the intermediate state
***************************************************

	p_r1_z  = - p_1_z
	p_r1_y  = - p_1_y
	p_r1_x  = - p_1_x	

  	call polar(p_1_x,p_1_y,p_1_z,p1,th1,phi1)

       
        wfre  = 0.0
        wfim  = 0.0
        awfre = 0.0
        awfim = 0.0

        call   dfunction(kj,ktm,ksm,ksr,p1,th1,phi1,wfre,wfim)
        call dfunction_a(kj,ktm,ksm,ksr,p1,th1,phi1,awfre,awfim)

*******************************************************************
*         calculating t that enters in the FSI through (p_r-p_r1)^2
********************************************************************
        p_r1 = sqrt(p_r1_x**2 + p_r1_y**2 + p_r1_z**2)
        E_r1 = sqrt(pm**2+p_r1**2)
*        p_r  = sqrt(p_r_x**2 + p_r_y**2 + p_r_z**2)
*        E_r  = sqrt(pm**2+p_r**2)


*****************************************************************************
*        Calculation of xi factor
*****************************************************************************
	mx2 =  dm**2 - 2.0*E_r1*dm + pm**2 +
     &         2.0*q0*(dm-E_r1) + 2.0*qv*p_r1_z - Q2
	             mx = 0.0
        if(mx2.lt.0.0)return
	if(mx2.gt.0.0)mx = sqrt(mx2)
	xi = sqrt((s-(pm-mx)**2)*(s-(pm+mx)**2))/
     &        (2.0*sqrt(E_r*E_r1)*qv)
* 	xi = 1.0
*****************************************************************************
*****************************************************************************
*        Calculation of norm term
*****************************************************************************
	tnorm = sqrt((dm-E_r)/(dm-E_r1))
*	tnorm = 1.0
*****************************************************************************




	pr1pr = p_r1_x*p_r_x +  p_r1_y*p_r_y +  p_r1_z*p_r_z 

	t    = 2.0*pm**2 - 2.0*E_r*E_r1 + 2.0*pr1pr
      
*	t = 0.0
	if(t.gt.0.0)return


        f_XN_on_re = 0.0
        f_XN_on_im = 0.0
        call  f_XN_on(s,t,f_XN_on_re,f_XN_on_im)

        f_XN_off_re = 0.0
        f_XN_off_im = 0.0
        call f_XN_off(s,t,f_XN_off_re,f_XN_off_im)
*	f_XN_off_re = 0.0
*	f_XN_off_im = 0.0


              
        un_1_re_pole = -(f_XN_on_re * wfim + f_XN_on_im * wfre )
        un_1_re_pv   = -(f_XN_off_re*awfre - f_XN_off_im*awfim)
        
*        write(6,*)"pv",un_1_re_pole,un_1_re_pv,awfre,awfim
        un_1_re = xi*(1.0/4.0)*q*(un_1_re_pole + un_1_re_pv)/(2.0*pi)**2 
     &            *tnorm

        return
	end


        function  un_1_im(q,phiq)
	implicit none
	real q,phiq
	real pi,pm,pmp,pmn,dm,eb
	real Q2,q0,qv
	real  p_m_x, p_m_y, p_m_z
	real  p_r_x, p_r_y, p_r_z,p_r,E_r
	real  p_1_x, p_1_y, p_1_z,p1,th1,phi1
	real  p_r1_x, p_r1_y, p_r1_z,p_r1,E_r1
	real s,t,delta_0,pr1pr
	real q_x,q_y
	real fct,un_1_im,un_1_im_pole,un_1_im_pv
	real wfre,wfim,awfre,awfim
	real f_XN_on_re, f_XN_on_im,F_XN_off_re,F_XN_off_im
	integer ktm,ksm,ktr,ksr,kj,lc
	real mx,mx2
        real xi,tnorm
	real wn2,E_r_pr,del_w2,delta_pr


        common/par/pi,pm,pmp,pmn,dm,eb
	common/vphoton/Q2,q0,qv
        common/missing/  p_m_x,  p_m_y,  p_m_z
        common/recoil_r/ p_r_x,  p_r_y,  p_r_z
        common/conf/ktm,ksm,ktr,ksr,kj
        common/enginv/s
        common/delta0/delta_0
	common/vnlc/lc

        un_1_im = 0.0
        
	q_x      = q*cos(phiq)
	q_y      = q*sin(phiq)

        p_r  = sqrt(p_r_x**2 + p_r_y**2 + p_r_z**2)
        E_r = sqrt(pm**2+p_r**2)
	fct = 1.0
	if(lc.eq.1)fct = pm/E_r

 	p_1_z   = fct*p_m_z !+ delta_0
	p_1_y   = p_m_y + q_y
	p_1_x   = p_m_x + q_x


*  	call polar(p_1_x,p_1_y,0.0,p1,th1,phi1)


**************************************************
* calculating delta_pr - delta that depends on pr
**************************************************
	p1 = sqrt(p_1_x**2+p_1_y**2)
	p_1_z = 0.0

	E_r_pr = sqrt(pm**2 + p1**2)	
*	E_r_pr = pm
	mx2 =  dm**2 - 2.0*E_r_pr*dm + pm**2 +
     &         2.0*q0*(dm-E_r_pr) - 2.0*qv*p_1_z - Q2

*	mx2 = pm**2  + 2.0*q0*pm - Q2
        wn2 =  dm**2 - 2.0*E_r*dm    + pm**2 +  
     &	       2.0*q0*(dm-E_r)    + 2.0*qv*p_r_z - Q2

	if(mx2.le.0.0.or.wn2.lt.pm**2)return

	del_w2 = wn2-mx2
	if(del_w2.le.0.0)del_w2=0.0
        delta_pr =  (E_r-E_r_pr)*(dm+q0)/qv  + del_w2/(2.0*qv)
	
	p_1_z = p_m_z + delta_pr

***************************************************
*  Recoil Nucleon's momentum and energy in 
*  the intermediate state
***************************************************

	p_r1_z  = - p_1_z
	p_r1_y  = - p_1_y
	p_r1_x  = - p_1_x	


  	call polar(p_1_x,p_1_y,p_1_z,p1,th1,phi1)
       
        wfre = 0.0
        wfim = 0.0
        awfre = 0.0
        awfim = 0.0
           
        call   dfunction(kj,ktm,ksm,ksr,p1,th1,phi1,wfre,wfim)
        call dfunction_a(kj,ktm,ksm,ksr,p1,th1,phi1,awfre,awfim)

*******************************************************************
*         calculating t that enters in the FSI through (p_r-p_r1)^2
********************************************************************
        p_r1 = sqrt(p_r1_x**2 + p_r1_y**2 + p_r1_z**2)
        E_r1 = sqrt(pm**2+p_r1**2)
*        p_r  = sqrt(p_r_x**2 + p_r_y**2 + p_r_z**2)
*        E_r  = sqrt(pm**2+p_r**2)

*****************************************************************************
*        Calculation of xi factor
*****************************************************************************
	mx2 =  dm**2 - 2.0*E_r1*dm + pm**2 +
     &         2.0*q0*(dm-E_r1) + 2.0*qv*p_r1_z - Q2
	             mx = 0.0
        if(mx2.lt.0.0)return
	if(mx2.gt.0.0)mx = sqrt(mx2)
	xi = sqrt((s-(pm-mx)**2)*(s-(pm+mx)**2))/
     &        (2.0*sqrt(E_r*E_r1)*qv)
*	xi = 1.0
*****************************************************************************
*****************************************************************************
*        Calculation of norm term
*****************************************************************************
	tnorm = sqrt((dm-E_r)/(dm-E_r1))
*	tnorm = 1.0
*****************************************************************************




	pr1pr = p_r1_x*p_r_x +  p_r1_y*p_r_y +  p_r1_z*p_r_z 

	t    = 2.0*pm**2 - 2.0*E_r*E_r1 + 2.0*pr1pr
*	t = 0.0
	if(t.gt.0.0)return

        f_XN_on_re = 0.0
        f_XN_on_im = 0.0       
        call  f_XN_on(s,t,f_XN_on_re,f_XN_on_im)

        f_XN_off_re = 0.0
        f_XN_off_im = 0.0
        call f_XN_off(s,t,f_XN_off_re,f_XN_off_im)
*	f_XN_off_re = 0.0
*	f_XN_off_im = 0.0


        un_1_im_pole =  (f_XN_on_re * wfre - f_XN_on_im * wfim)
        un_1_im_pv   = -(f_XN_off_re*awfim + f_XN_off_im*awfre)

*	write(6,*)"pvim",un_1_im_pv
        un_1_im = xi*(1.0/4.0)*q*(un_1_im_pole + un_1_im_pv)/(2.0*pi)**2 
     &            *tnorm
        return
	end





        SUBROUTINE SQUAR_EQ(A,B,C,X1,X2,I_TEST)
        I_TEST = 0
        DSKM  = B**2 - 4.0*A*C
        IF(DSKM.LT.0.0)THEN
                IF(DSKM.GT.-1.0E-4)THEN
                DSKM = 0.0
                GOTO 111
                ELSE
                I_TEST = 1
                ENDIF
        RETURN
        ENDIF
111     X1 = ( -B + SQRT(DSKM) ) / 2.0 / A
        X2 = ( -B - SQRT(DSKM) ) / 2.0 / A
        RETURN
        END



 	subroutine polar(p_x,p_y,p_z,p,th,phi)
	implicit none
	real p_x,p_y,p_z,p,th,phi
	real arg,sn_th
        p = sqrt(p_x**2 + p_y**2 + p_z**2)
        th = 0.0
        phi = 0.0
        if(p.eq.0.0)return
	arg = p_z/p
	if(arg.gt. 1.0)arg =  1.0
	if(arg.lt.-1.0)arg = -1.0
	th  = acos(arg)
	sn_th = sin(th)
	phi = 0.0
        if(p_y.eq.0.0.and.p_x.lt.0.0)then
        phi = acos(-1.0)
        return
        elseif(p_y.eq.0.0.and.p_x.gt.0.0)then
        phi = 0.0
        return
        endif
	if(sn_th.ne.0.0)then
	arg = p_x/p/sn_th
	if(arg.gt. 1.0)arg =  1.0
	if(arg.lt.-1.0)arg = -1.0
	phi = acos(arg)
	if((p_y/sn_th).lt.0.0)phi = 2.0*acos(-1.0) - phi
	endif
	return
	end


    	subroutine f_XN_on(s,t,f_XN_onre,f_XN_onim)
        sigma_tot = sigma_xn_tot(s) / 0.389385   ! Gev^{-2}
        f_XN_onre = sigma_tot*a_xn(s)*exp(bn(s)*t/2.0)
        f_XN_onim = sigma_tot        *exp(bn(s)*t/2.0)
	
        return
        end
       
        subroutine f_XN_off(s,t,f_XN_offre,f_XN_offim)
        sigma_tot  = sigma_xn_tot(s) / 0.389385   ! Gev^{-2}
        f_XN_offre = sigma_tot*a_xn(s)*exp(bn(s)*t/2.0)
        f_XN_offim = sigma_tot        *exp(bn(s)*t/2.0)
        return
        end



       function sigma_xn_tot(s)
***************************************************************************
*      Total Cross section DISX-n scattering
***************************************************************************
       common/par/pi,pm,pmp,pmn,dm,eb ! pm -> mass of the proton
*       stot = 47.3 + 0.513*alog(p)**2 - 4.27*alog(p)
       sigma_xn_tot = 40.0 !60. !40. !80.0 !60.0
       return
       end

       function bn(s)
***************************************************************************
*      DISX_N scattering Slope parameter B
***************************************************************************
       common/par/pi,pm,pmp,pmn,dm,eb ! pm -> mass of the proton
       bn  = 8.0
       return 
       end

       function a_xn(s)
****************************************************************************
*       Re/Im  part of DISX - N scattering (alpha_xn)
****************************************************************************
       common/par/pi,pm,pmp,pmn,dm,eb  ! pm -> mass of the proton
*        a_pn=-0.56207 +0.24223E-01*p -0.50362E-03*p**2 
*     &       +0.48408E-05*p**3 -0.17331E-07*p**4
	a_xn = -0.2
        return
       end

        





****************************************
* version 2
*
* Misak Sargsian
*
* June 2005
* FIU, TUN Miami
*****************************************
      subroutine struct_bodek_p(gp0,gp2,fm2,f2,f1)        
*******************************************************
* Structure function f2, f2 of protons within Bodek   *
* parameterization                                    *
*******************************************************
      PM  = 0.938279   
      R   = 0.18
*      gp0 = gp2/2.0/pm/x                                              
*      FM2 = PM**2+2.*PM*GP0-GP2                                           
      W2H = 0.                                                            
      IF(FM2.LT.PM**2) RETURN                                           
      WI  = SQRT(FM2)                                                      
      W2H = GP_H(GP0,GP2)*B(WI,GP2)/GP0  
      W1H = (1.+GP0**2/GP2)/(1.+R)*W2H
      f2  = gp0 * w2h                           
      f1  = pm  * w1h                                   
      RETURN                                                            
      END                                                               
C     .........................                                         
      FUNCTION GP_H(Q0,Q2)                                                
      PM = 0.938279                                        
      XX = Q2/(2.*PM*Q0)                                         
      GI = 2.*PM*Q0                                              
      WW = (GI+1.642)/(Q2+0.376)                                 
      T  = (1.-1./WW)                                             
      WP = 0.256*T**3+2.178*T**4+0.898*T**5-6.716*T**6+3.756*T**7
      GP_H=WP*WW*Q2/(2.*PM*Q0) 
      RETURN                                                            
      END                                                               
C---------------------------------                                      





C..................................             
      subroutine struct_bodek_n(gp0,gp2,fm2,f2,f1)         
*******************************************************
* Structure function f2, f2 of neutrons within Bodek  *
* parameterization                                    *
*******************************************************
      common/ths_ths/ths0
      PM  = 0.938279     
      R   = 0.18            
*      gp0 = gp2/2.0/pm/x                                             
*      FM2=PM**2+2.*PM*GP0-GP2                                           
      W2NT=0.                                                            
      IF(FM2.LT.PM**2) RETURN                                           
      WI=SQRT(FM2)                                                      
      W2NT =GP_N(GP0,GP2)*B(WI,GP2)/GP0       
      W1NT = (1.+GP0**2/GP2)/(1.+R)*W2NT
      f2  = gp0 * w2NT                           
      f1  = pm  * w1NT       
*      write(15,*)"bodek",ths0,f1,f2,gp0,gp2,wi
      RETURN                                                            
      END                                                               
C     .........................                                         
      FUNCTION GP_N(Q0,Q2)                                                
      PM = 0.938279                                        
      XX = Q2/(2.*PM*Q0)                                         
      GI = 2.*PM*Q0                                              
      WW = (GI+1.642)/(Q2+0.376)                                 
      T  = (1.-1./WW)                                             
      WN = 0.064*T**3+0.225*T**4+4.106*T**5-7.079*T**6+3.055*T**7    
      GP_N=WN*WW*Q2/(2.*PM*Q0) 
      RETURN                                                            
      END                                                               
C---------------------------------                                      

 
C  *******************************                                       
C  *    BODEK PARAMETERIZATION    *                                       
C  *******************************                                       
      FUNCTION B(WM,QSQ)                                                 
      DIMENSION C(24)                                                    
      INTEGER LSPIN(4)                                                   
      DATA PMSQ/0.880324/,PM2/1.876512/,PM/0.938256/                     
      DATA NRES/4/,NBKG/5/                                               
      DATA LSPIN/1,2,3,2/                                                
      DATA C/1.0741163,0.75531124,3.3506491,1.7447015,3.5102405,1.040004 
     *,1.2299128,0.10625394,0.48132786,1.5101467,0.081661975,0.65587179, 
     *1.7176216,0.12551987,0.7473379,1.953819,0.19891522,-0.17498537,    
     *0.0096701919,-0.035256748,3.5185207,-0.59993696,4.7615828,0.411675  
     *89/                                                                
      B=0.                                                               
      IF(WM.LE.0.939)RETURN                                               
      WSQ=WM**2                                                          
      OMEGA=1.+(WSQ-PMSQ)/QSQ                                            
      X=1./OMEGA                                                         
      XPX=C(22)+C(23)*(X-C(24))**2                                       
      PIEMSQ=(C(1)-PM)**2                                                
********************************************************
*     added part
********************************************************
      B1 = 0.0
      IF(WM.EQ.C(1))GOTO 11 
******************************************************** 
      B1=AMAX1(0.,(WM-C(1)))/(WM-C(1))*C(2)       ! 0/0                  
********************************************************
11    EB1=C(3)*(WM-C(1))                                                 
      IF(EB1.GT.25.)GO TO 1                                              
      B1=B1*(1.0-EXP(-EB1))                                              
*********************************************************
*     added part
*********************************************************  
      B2 = 0.0
      IF(WM.EQ.C(4))GOTO 12  
*********************************************************
1     B2=AMAX1(0.,(WM-C(4)))/(WM-C(4))*(1.-C(2))   ! 0/0                 
*********************************************************
12    EB2=C(5)*(WSQ-C(4)**2)                                             
      IF(EB2.GT.25.) GO TO 2                                             
      B2=B2*(1.-EXP(-EB2))                                               
2     CONTINUE                                                           
      BBKG=B1+B2                                                         
      BRES=C(2)+B2                                                       
      RESSUM=0.                                                          
      DO 30 I=1,NRES                                                     
      INDEX=(I-1)*3+1+NBKG                                               
      RAM=C(INDEX)                                                       
      IF(I.EQ.1)RAM=C(INDEX)+C(18)*QSQ+C(19)*QSQ**2                      
      RMA=C(INDEX+1)                                                     
      IF(I.EQ.3)RMA=RMA*(1.+C(20)/(1.+C(21)*QSQ))                        
      RWD=C(INDEX+2)                                                     
      QSTARN=SQRT(AMAX1(0.,((WSQ+PMSQ-PIEMSQ)/(2.*WM))**2-PMSQ))         
      QSTARO=SQRT(AMAX1(0.,((RMA**2-PMSQ+PIEMSQ)/(2.*RMA))**2-PIEMSQ))   
      IF(QSTARO.LE.1.E-10)GO TO 40                                       
      TERM=6.08974*QSTARN                                                
      TERMO=6.08974*QSTARO                                               
      J=2*LSPIN(I)                                                       
      K=J+1                                                              
      GAMRES=RWD*(TERM/TERMO)**K*(1.+TERMO**J)/(1.+TERM**J)              
      GAMRES=GAMRES/2.                                                   
      BRWIG=GAMRES/((WM-RMA)**2+GAMRES**2)/3.1415926                     
      RES=RAM*BRWIG/PM2                                                  
      GO TO 30                                                           
40    RES=0.                                                             
30    RESSUM=RESSUM+RES                                                  
      B=BBKG*(1.+(1.-BBKG)*XPX)+RESSUM*(1.-BRES)                         
      RETURN                                                             
      END                                                                


      subroutine dfunction(kj,kt,ksp,ksn,p,theta,phi,wf_rea,wf_ima)
****************************************************************************************
*       Subroutine calculates the deuteron wave function by spin and isospin 
*       components for given deuteron spin projection
*       kj  = 1,0,-1 - deuteron spin projection
*       kt  = 1,-1 - struck out nucleon isospin 1-proton -1 neutron
*       ksp = 1,-1 - spin projection of proton
*       ksn = 1,-1 - spin projection of neutron
*       p, theta, phi - momentum(GeV/c), polar (rad) and azimuthal(rad) angles of 
*                       relative p-n momentum
*       wf_rea - real part of the wave function
*       wf_ima - imaginary part of the wave function
*
*      20-July-2003  
*      Miami
*
*****************************************************************************************



                  fis =  1.0 ! struck proton
      if(kt.eq.-1)fis = -1.0 ! struck neutron


      if(kj.eq.1)then ! deuteron spin projection 1
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      wf_re = uu(p) + wd(p)/sqrt(8.0)*(3.0*cos(theta)**2-1.0)
      wf_im = 0.0
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      wf_re = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi)   
      wf_im = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      wf_re = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi)   
      wf_im = wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      wf_re = wd(p)/sqrt(8.0)*3.0*sin(theta)**2*cos(2.0*phi)
      wf_im = wd(p)/sqrt(8.0)*3.0*sin(theta)**2*sin(2.0*phi)
      endif

      elseif(kj.eq.0)then ! deuteron spin projection 0
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      wf_re =  wd(p)*3.0/2.0*cos(theta)*sin(theta)*cos(phi)
      wf_im = -wd(p)*3.0/2.0*cos(theta)*sin(theta)*sin(phi)
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      wf_re = uu(p)/sqrt(2.0) - wd(p)/2.0 * (3.0*cos(theta)**2-1.0)
      wf_im = 0.0
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      wf_re = uu(p)/sqrt(2.0) - wd(p)/2.0 * (3.0*cos(theta)**2-1.0)
      wf_im = 0.0
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      wf_re = -wd(p)*3.0/2.0*cos(theta)*sin(theta)*cos(phi)
      wf_im = -wd(p)*3.0/2.0*cos(theta)*sin(theta)*sin(phi)
      endif

      elseif(kj.eq.-1)then 
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      wf_re =  wd(p)/sqrt(8.0)*3.0*sin(theta)**2*cos(2.0*phi)
      wf_im = -wd(p)/sqrt(8.0)*3.0*sin(theta)**2*sin(2.0*phi)
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      wf_re =  wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi)   
      wf_im = -wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi)        
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      wf_re =  wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi)   
      wf_im = -wd(p)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      wf_re = uu(p) + wd(p)/sqrt(8.0)*(3.0*cos(theta)**2-1.0)
      wf_im = 0.0
      endif
      endif
      wf_rea = fis*wf_re 
      wf_ima = fis*wf_im
      return
      end

      subroutine dfunctionc(kj,kt,ksp,ksn,p,theta,phi,wf_rea,wf_ima)
****************************************************************************************
*       Subroutine calculates the deuteron wave function by spin and isospin 
*       components for given deuteron spin projection
*       kj  = 1,0,-1 - deuteron spin projection
*       kt  = 1,-1 - struck out nucleon isospin 1-proton -1 neutron
*       ksp = 1,-1 - spin projection of proton
*       ksn = 1,-1 - spin projection of neutron
*       p, theta, phi - momentum(GeV/c), polar (rad) and azimuthal(rad) angles of 
*                       relative p-n momentum
*       wf_rea - real part of the wave function
*       wf_ima - imaginary part of the wave function
*
*      1-December-2003  
*      Miami
*
*****************************************************************************************
                  fis =  1.0 ! struck proton
      if(kt.eq.-1)fis = -1.0 ! struck neutron
      call zeta(kj,ksp,ksn,ze_re,ze_im)
      call stensor(kj,ksp,ksn,theta,phi,s_re,s_im)
      s8 = sqrt(8.0)

      wf_re = uu(p)*ze_re + wd(p)/s8 * s_re 
      wf_im =             + wd(p)/s8 * s_im 
      wf_rea = fis*wf_re 
      wf_ima = fis*wf_im
      return
      end


      subroutine dfunction_a(kj,kt,ksp,ksn,p,theta,phi,awf_rea,awf_ima)
****************************************************************************************
*       Subroutine calculates the a-deuteron wave function by spin and isospin 
*       components for given deuteron spin projection
*       kj  = 1,0,-1 - deuteron spin projection
*       kt  = 1,-1 - knocked out nucleon isospin 1-proton -1 neutron
*       ksp = 1,-1 - spin projection of proton
*       ksn = 1,-1 - spin projection of neutron
*       p, theta, phi - momentum(GeV/c), polar (rad) and azimuthal(rad) angles of 
*                       relative p-n momentum
*       awf_rea - real part of the wave function
*       awf_ima - imaginary part of the wave function
*
*      20-July-2003  
*      Miami
*
*****************************************************************************************

      pz = p*cos(theta) 
      pmin = 0.005
                  fis =  1.0 ! knocked-out proton
      if(kt.eq.-1)fis = -1.0 ! knocked-out neutron

      call zeta(kj,ksp,ksn,ze_re,ze_im)

      call stensor(kj,ksp,ksn,theta,phi,s_re,s_im)
      theta0 = acos(-1.0)/2.0
      call stensor(kj,ksp,ksn,theta0,phi,s0_re,s0_im)
                           tn2 = 0.0
      if(cos(theta).ne.0.0)tn2 = (sin(theta)/cos(theta))**2
      
      pt = abs(p*sin(theta))
      s8 = sqrt(8.0)

      awf_re= uu1(p,pt)*ze_re + 
     &        wd1(p,pt)/s8*s_re + tn2*wd2(pt)/s8*(s_re-s0_re) 
      awf_im= wd1(p,pt)/s8*s_im + tn2*wd2(pt)/s8*(s_im-s0_im) 

      awf_rea = pz*(fis*awf_re)
      awf_ima = pz*(fis*awf_im)
      return
      end


      subroutine zeta(kj,ksp,ksn,ze_re,ze_im)
****************************************************************************************
*       Subroutine calculates the spin wave function of the s component of the deuteron 
*       for given deuteron spin projection
*       kj  = 1,0,-1 - deuteron spin projection
*       ksp = 1,-1 - spin projection of proton
*       ksn = 1,-1 - spin projection of neutron
*       ze_re - real part of the tensor 
*       ze_im - imaginary part of the tensor wave 
*
*      1-December-2003  
*      FIU
*      Miami
*
*****************************************************************************************
      ze_re = 0.0
      ze_im = 0.0
      if(kj.eq.1)then      ! deuteron spin projection 1
      if(ksp.eq.1.and.ksn.eq.1)ze_re = 1             ! spins up
      elseif(kj.eq.0)then  ! deuteron spin projection 0
      if(ksp.eq.1.and.ksn.eq.-1)ze_re = 1.0/sqrt(2.0)! proton up neutron down
      if(ksp.eq.-1.and.ksn.eq.1)ze_re = 1.0/sqrt(2.0)! proton down neutron up
      elseif(kj.eq.-1)then ! deuteron spin projection -1
      if(ksp.eq.-1.and.ksn.eq.-1)ze_re = 1.0         ! spins down
      endif
      return
      end



      subroutine stensor(kj,ksp,ksn,theta,phi,s_re,s_im)
****************************************************************************************
*       Subroutine calculates the tensor S by spin and isospin 
*       components for given deuteron spin projection
*       S = [3 (\sigma_p p)(\sigma_n p)/p^2 - \sigma_p\sigma_n]\zeta
*       kj  = 1,0,-1 - deuteron spin projection
*       kt  = 1,-1 - struck out nucleon isospin 1-proton -1 neutron
*       ksp = 1,-1 - spin projection of proton
*       ksn = 1,-1 - spin projection of neutron
*       p, theta, phi - momentum(GeV/c), polar (rad) and azimuthal(rad) angles of 
*                       relative p-n momentum
*       s_re - real part of the tensor 
*       s_im - imaginary part of the tensor wave 
*
*      1-December-2003  
*      FIU
*      Miami
*
*****************************************************************************************

      if(kj.eq.1)then ! deuteron spin projection 1
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      s_re = 3.0*cos(theta)**2-1.0
      s_im = 0.0
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      s_re = 3.0*cos(theta)*sin(theta)*cos(phi)   
      s_im = 3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      s_re = 3.0*cos(theta)*sin(theta)*cos(phi)   
      s_im = 3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      s_re = 3.0*sin(theta)**2*cos(2.0*phi)
      s_im = 3.0*sin(theta)**2*sin(2.0*phi)
      endif

      elseif(kj.eq.0)then ! deuteron spin projection 0
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      s_re =  sqrt(2.0)*3.0*cos(theta)*sin(theta)*cos(phi)
      s_im = -sqrt(2.0)*3.0*cos(theta)*sin(theta)*sin(phi)
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      s_re = -sqrt(2.0)*(3.0*cos(theta)**2-1.0)
      s_im =  0.0
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      s_re = -sqrt(2.0)*(3.0*cos(theta)**2-1.0)
      s_im = 0.0
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      s_re = -sqrt(2.0)*3.0*cos(theta)*sin(theta)*cos(phi)
      s_im = -sqrt(2.0)*3.0*cos(theta)*sin(theta)*sin(phi)
      endif

      elseif(kj.eq.-1)then 
          if(ksp.eq.1.and.ksn.eq.1)then ! spins up
      s_re =  3.0*sin(theta)**2*cos(2.0*phi)
      s_im = -3.0*sin(theta)**2*sin(2.0*phi)
      elseif(ksp.eq.1.and.ksn.eq.-1)then ! proton up neutron down
      s_re =  3.0*cos(theta)*sin(theta)*cos(phi)   
      s_im = -3.0*cos(theta)*sin(theta)*sin(phi)        
      elseif(ksp.eq.-1.and.ksn.eq.1)then ! proton down neutron up
      s_re =  3.0*cos(theta)*sin(theta)*cos(phi)   
      s_im = -3.0*cos(theta)*sin(theta)*sin(phi)   
      elseif(ksp.eq.-1.and.ksn.eq.-1)then ! proton down neutron down
      s_re = 3.0*cos(theta)**2-1.0
      s_im = 0.0
      endif
      endif
      return
      end



      function uu(p)
      uu = U(p/0.197328)/sqrt(0.197328**3)
      return
      end
      function uu1(p,pt)
      x  = p /0.197328
      xt = pt/0.197328
      uu1 = U_A(x,xt)/sqrt(0.197328**3)/0.197328
*      uu1 = U_A(x,xt)/0.197328**(5.0/2.0)
      return
      end

      function wd(p)
      wd = W(p/0.197328)/sqrt(0.197328**3)
      return
      end
      function wd1(p,pt)
      x  = p /0.197328
      xt = pt/0.197328
      wd1 = W_A(x,xt)/sqrt(0.197328**3)/0.197328
      return
      end

      function wd2(pt)
      xt = pt/0.197328
      wd2 = W_AA(xt)/sqrt(0.197328**3)/0.197328
      return
      end


      FUNCTION FD(X,i)                                                  
C ************************************************                      
C *  DEUTRON WAVE FUNCTION WITH PARIS POTENTIAL  *                      
C ************************************************                      
      COMMON/PARIS/C(13),D(13),BM(13)                                   
      real *8 c,d,bm
      if(i.eq.1)then
      C(1)=0.88688076                                                   
      C(2)=-0.34717093                                                  
      C(3)=-3.050238                                                    
      C(4)=56.207766                                                    
      C(5)=-749.57334                                                   
      C(6)=5336.5279                                                    
      C(7)=-22706.863                                                   A1507090
      C(8)=60434.4690                                                   A1507100
      C(9)=-102920.58                                                   A1507110
      C(10)=112233.57                                                   A1507120
      C(11)=-75925.226                                                  A1507130
      C(12)=29059.715                                                   A1507140
      A=0.                                                              A1507150
      DO 401 J=1,12                                                     A1507160
401   A=A+C(J)                                                          A1507170
      C(13)=-A                                                          A1507180
      D(1)=0.023135193                                                  A1507190
      D(2)=-0.85604572                                                  A1507200
      D(3)=5.6068193                                                    A1507210
      D(4)=-69.462922                                                   A1507220
      D(5)=416.31118                                                    A1507230
      D(6)=-1254.6621                                                   A1507240
      D(7)=1238.783                                                     A1507250
      D(8)=3373.9172                                                    A1507260
      D(9)=-13041.151                                                   A1507270
      D(10)=19512.524                                                   A1507280
      DO 402 J=1,13                                                     A1507290
402   BM(J)=0.23162461+(J-1)                                            A1507300
      A=0.                                                              A1507310
      B=0.                                                              A1507320
      CC=0.                                                             A1507330
      DO 3 J=1,10                                                       A1507340
      A=A+D(J)/BM(J)**2                                                 A1507350
      B=B+D(J)                                                          A1507360
3     CC=CC+D(J)*BM(J)**2                                               A1507370
      D(11)=BM(11)**2/(BM(13)**2-BM(11)**2)/(BM(12)**2-BM(11)           A1507380
     ***2)*(-BM(12)**2*BM(13)**2*A+(BM(12)**2+BM(13)**2)*B-CC)          A1507390
      D(12)=BM(12)**2/(BM(11)**2-BM(12)**2)/(BM(13)**2-BM(12)           A1507400
     ***2)*(-BM(13)**2*BM(11)**2*A+(BM(13)**2+BM(11)**2)*B-CC)          A1507410
      D(13)=BM(13)**2/(BM(12)**2-BM(13)**2)/(BM(11)**2-BM(13)           A1507420
     ***2)*(-BM(11)**2*BM(12)**2*A+(BM(11)**2+BM(12)**2)*B-CC)          A1507430
      endif
      FD=(U(X/0.197328)**2+W(X/0.197328)**2)/0.197328**3                A1507440
      RETURN                                                            A1507450
      END                                                               A1507460
C ***** S PARTIAL WAVE ******                                           A1507470
      FUNCTION U(X)                                                     A1507480
      COMMON/PARIS/C(13),D(13),BM(13)                                   A1507490
      real *8 c,d,bm
      A=0.                                                              A1507500
      DO 1 J=1,13                                                       A1507510
1     A=C(J)/(X*X+BM(J)**2)+A                                           A1507520
      F=0.79788456                                                      A1507530
      U=A*F/SQRT(4.*3.14159265)                                         A1507540
      RETURN                                                            A1507550
      END                                                               A1507560
C  **** D PARTIAL WAVE *****                                            A1507570
      FUNCTION W(X)                                                     A1507580
      COMMON/PARIS/C(13),D(13),BM(13)                                   A1507590
      real *8 c,d,bm
      A=0.                                                              A1507600
      DO 1 J=1,13                                                       A1507610
1     A=D(J)/(X*X+BM(J)**2)+A                                           A1507620
      F=0.79788456                                                      A1507630
      W=A*F/SQRT(4.*3.14159265)                                         A1507640
      RETURN                                                            A1507650
      END                                                               A1507660




C ***** S' PARTIAL WAVE ******                                           
      FUNCTION U_A(X,XT)                                              
      COMMON/PARIS/C(13),D(13),BM(13)                                    
      real *8 c,d,bm

      U_A = 0.0
      A=0.                                                               
      DO 1 J=1,13                                                        
1     A=C(J)/(X*X+BM(J)**2)/sqrt(XT**2+BM(J)**2)+A                       
      F=0.79788456                                                       
      U_A=A*F/SQRT(4.*3.14159265)                                        
      RETURN                                                             
      END                                                                


C  **** D' 1  PARTIAL WAVE *****                                            
      FUNCTION W_A(X,XT)                                                     
      COMMON/PARIS/C(13),D(13),BM(13)                                   
      real *8 c,d,bm

      A=0.                                                              
      DO 1 J=1,13                                                       
1     A=D(J)/(X*X+BM(J)**2)/sqrt(XT**2+BM(J)**2)+ A                   
      F=0.79788456                                                      
      W_A=A*F/SQRT(4.*3.14159265)                                         
      RETURN                                                            
      END                                                               

C  **** D' 2  PARTIAL WAVE *****                                            
      FUNCTION W_AA(XT)                                                     
      COMMON/PARIS/C(13),D(13),BM(13)                                   
      real *8 c,d,bm
      A=0.                                                              
      DO 1 J=1,13                                                       
1     A=D(J)/(BM(J)**2)/sqrt(XT**2+BM(J)**2)+ A                   
      F=0.79788456                                                      
      W_AA = A*F/SQRT(4.*3.14159265)                                         
      RETURN                                                            
      END                                                               

C  **** D PARTIAL WAVE *****                                             
      FUNCTION W_Askzbi(X,XT)                                                 
      COMMON/PARIS/C(13),D(13),BM(13)                                    
      real *8 c,d,bm
      W_Askzbi = 0.0
      A=0.     
      IF(X.GT.0.0)THEN
      DO 1 J=1,13    
      TERM1 = (3.0*XT**2+2.0*BM(J)**2)/(X*X+BM(J)**2)
     &        /sqrt(XT**2+BM(J)**2)
      TERM2 = 0.0 !3.0*XT/X**2
      TERM  = TERM1   !- TERM2
      A   = D(J)/BM(J)**2 * TERM + A
1     CONTINUE
      ENDIF
      F=0.79788456                                                       
      W_Askzbi=A*F/SQRT(4.*3.14159265)                                        
      RETURN                                                             
      END                                                                






*********************************                                       
*   Subroutine for integration                                          
********************************                                        
      SUBROUTINE GADAP(A0,B0,F,EPS,SUM)                                 
      COMMON/GADAP1/ NUM,IFU                                            
      EXTERNAL F                                                        
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)     
    1 FORMAT(16H GADAP:I TOO BIG)                                       
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)          
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.3   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)    
      F2(1)=F(0.5*(A0+B0))  
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)    
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(A(I)+B(I)-W1)   
      F2(I+1)=F3(I) 
      F3(I+1)=F(B(I)-A(I+2)+W1) 
      F1(I+2)=F(U2) 
      F2(I+2)=F2(I) 
      F3(I+2)=F(B(I+2)+A(I+2)-U2)   
      F1(I+3)=F(A(I)+A(I+2)-W1) 
      F2(I+3)=F1(I) 
      F3(I+3)=F(W1) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   


      SUBROUTINE GADAPu(A0,B0,F,EPS,SUM)                                 
      COMMON/GADAPu1/ NUM,IFU                                            
      EXTERNAL F                                            
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)     
*1     FORMAT(16H GADAPu:I TOO BIG)              
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)          
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.3   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)    
      F2(1)=F(0.5*(A0+B0))  
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)    
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(A(I)+B(I)-W1)   
      F2(I+1)=F3(I) 
      F3(I+1)=F(B(I)-A(I+2)+W1) 
      F1(I+2)=F(U2) 
      F2(I+2)=F2(I) 
      F3(I+2)=F(B(I+2)+A(I+2)-U2)   
      F1(I+3)=F(A(I)+A(I+2)-W1) 
      F2(I+3)=F1(I) 
      F3(I+3)=F(W1) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   

      SUBROUTINE GADAP2(A0,B0,FL,FU,F,EPS,SUM)  
      COMMON/GADAP_2/ NUM,IFU    
      EXTERNAL F,FL,FU  
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      X=0.5*(1+C)*A0+0.5*(1-C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F1(1)=FGADAP(X,AY,BY,F,EPS)   
      X=0.5*(A0+B0) 
      AY=FL(X)  
      BY=FU(X)  
      F2(1)=FGADAP(X,AY,BY,F,EPS)   
      X=0.5*(1-C)*A0+0.5*(1+C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F3(1)=FGADAP(X,AY,BY,F,EPS)   
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      X=A(I)+B(I)-W1    
      AY=FL(X)  
      BY=FU(X)  
      F1(I+1)=FGADAP(X,AY,BY,F,EPS) 
      F2(I+1)=F3(I) 
      X=B(I)-A(I+2)+W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+1)=FGADAP(X,AY,BY,F,EPS) 
      X=U2  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+2)=FGADAP(X,AY,BY,F,EPS) 
      F2(I+2)=F2(I) 
      X=B(I+2)+A(I+2)-U2    
      AY=FL(X)  
      BY=FU(X)  
      F3(I+2)=FGADAP(X,AY,BY,F,EPS) 
      X=A(I)+A(I+2)-W1  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+3)=FGADAP(X,AY,BY,F,EPS) 
      F2(I+3)=F1(I) 
      X=W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+3)=FGADAP(X,AY,BY,F,EPS) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   
      FUNCTION FGADAP(X,A0,B0,F,EPS)    
      COMMON/GADAP_2/ NUM,IFU    
      EXTERNAL F    
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(X,0.5*(1+C)*A0+0.5*(1-C)*B0)  
      F2(1)=F(X,0.5*(A0+B0))    
      F3(1)=F(X,0.5*(1-C)*A0+0.5*(1+C)*B0)  
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(X,A(I)+B(I)-W1) 
      F2(I+1)=F3(I) 
      F3(I+1)=F(X,B(I)-A(I+2)+W1)   
      F1(I+2)=F(X,U2)   
      F2(I+2)=F2(I) 
      F3(I+2)=F(X,B(I+2)+A(I+2)-U2) 
      F1(I+3)=F(X,A(I)+A(I+2)-W1)   
      F2(I+3)=F1(I) 
      F3(I+3)=F(X,W1)   
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  FGADAP=SUM    
      EPS=EPS/RED   
      RETURN    
      END   




      SUBROUTINE GADAPS2(A0,B0,FL,FU,F,EPS,SUM)  
      COMMON/GADAPS_2/ NUM,IFU    
      EXTERNAL F,FL,FU  
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      X=0.5*(1+C)*A0+0.5*(1-C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F1(1)=FGADAPS(X,AY,BY,F,EPS)   
      X=0.5*(A0+B0) 
      AY=FL(X)  
      BY=FU(X)  
      F2(1)=FGADAPS(X,AY,BY,F,EPS)   
      X=0.5*(1-C)*A0+0.5*(1+C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F3(1)=FGADAPS(X,AY,BY,F,EPS)   
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      X=A(I)+B(I)-W1    
      AY=FL(X)  
      BY=FU(X)  
      F1(I+1)=FGADAPS(X,AY,BY,F,EPS) 
      F2(I+1)=F3(I) 
      X=B(I)-A(I+2)+W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+1)=FGADAPS(X,AY,BY,F,EPS) 
      X=U2  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+2)=FGADAPS(X,AY,BY,F,EPS) 
      F2(I+2)=F2(I) 
      X=B(I+2)+A(I+2)-U2    
      AY=FL(X)  
      BY=FU(X)  
      F3(I+2)=FGADAPS(X,AY,BY,F,EPS) 
      X=A(I)+A(I+2)-W1  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+3)=FGADAPS(X,AY,BY,F,EPS) 
      F2(I+3)=F1(I) 
      X=W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+3)=FGADAPS(X,AY,BY,F,EPS) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   
      FUNCTION FGADAPS(X,A0,B0,F,EPS)    
      COMMON/GADAPS_2/ NUM,IFU    
      EXTERNAL F    
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(X,0.5*(1+C)*A0+0.5*(1-C)*B0)  
      F2(1)=F(X,0.5*(A0+B0))    
      F3(1)=F(X,0.5*(1-C)*A0+0.5*(1+C)*B0)  
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(X,A(I)+B(I)-W1) 
      F2(I+1)=F3(I) 
      F3(I+1)=F(X,B(I)-A(I+2)+W1)   
      F1(I+2)=F(X,U2)   
      F2(I+2)=F2(I) 
      F3(I+2)=F(X,B(I+2)+A(I+2)-U2) 
      F1(I+3)=F(X,A(I)+A(I+2)-W1)   
      F2(I+3)=F1(I) 
      F3(I+3)=F(X,W1)   
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  FGADAPS=SUM    
      EPS=EPS/RED   
      RETURN    
      END   





      SUBROUTINE simpson(A0,B0,F,n,SUM) 
C ********************************************************************** 
C PURPOSE - INTEGRATE A FUNCTION F(X) 
C METHOD - STUPID
C USAGE - CALL simpson(a0,b0,f,n,sum)
C PARAMETERS A0 - LOWER LIMIT (INPUT,REAL) 
C B0 - UPPER LIMIT (INPUT,REAL) 
C F - FUNCTION F(X) TO BE INTEGRATED. MUST BE 
C SUPPLIED BY THE USER. (INPUT,REAL FUNCTION) 
C n - NUMBER OF DIVISIONS 
C SUM - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL) 
C PRECISION - SINGLE 
C REQ'D PROG'S - F 
C
* FIU 
C....................................................................... 
*      implicit real*8(a-h,o-x)
*      integer*8 n
      external f
*      tn = n
      width = (b0-a0)/n
      sumi = 0.0
      do i = 1,n
      x = a0 + float(i)*width 
      sumi = sumi + f(x)* width
      enddo
     
      sum = sumi
      return
      end

***************************************************************************
*
* $Id
*
* $Log
*
* #include "gen/pilot.h"
*#if defined(CERNLIB_DOUBLE)
*      SUBROUTINE RADMUL
*     1 (F,N,A,B,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,NFNEVL,IFAIL)
*      CHARACTER NAME*(*)
*      PARAMETER (NAME = 'RADMUL')
*      CALL MTLPRT(NAME,'D120',
*     +'not available on this machine - see documentation')
*      RETURN
*      END

*      SUBROUTINE DADMUL
*     1 (F,N,A,B,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,NFNEVL,IFAIL)
*#include "gen/imp64.inc"

*#else
*      SUBROUTINE DADMUL
*     1 (F,N,A,B,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,NFNEVL,IFAIL)
* #include "gen/imp128.inc"
*     CHARACTER NAME*(*)
*      PARAMETER (NAME = 'DADMUL')
*      CALL MTLPRT(NAME,'D120',
*     +'not available on this machine - see documentation')
*      RETURN
*      END

      SUBROUTINE RADMUL
     1 (F,N,A,B,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,NFNEVL,IFAIL)
*#endif
 
      LOGICAL LDV
 
      DIMENSION A(*),B(*),WK(*)
      DIMENSION CTR(15),WTH(15),WTHL(15),Z(15)
      DIMENSION W(2:15,5),WP(2:15,3)
 
      PARAMETER (R1 = 1, HF = R1/2)
 
      PARAMETER (XL2 =  0.35856 85828 00318 073D0)
      PARAMETER (XL4 =  0.94868 32980 50513 796D0)
      PARAMETER (XL5 =  0.68824 72016 11685 289D0)
 
      PARAMETER (W2 =  980*R1/6561, W4 = 200*R1/19683)
      PARAMETER (WP2 =  245*R1/486, WP4 = 25*R1/729)
 
      DATA (W(N,1),W(N,3),N=2,15)
     1/-0.193872885230909911D+00,  0.518213686937966768D-01,
     2 -0.555606360818980835D+00,  0.314992633236803330D-01,
     3 -0.876695625666819078D+00,  0.111771579535639891D-01,
     4 -0.115714067977442459D+01, -0.914494741655235473D-02,
     5 -0.139694152314179743D+01, -0.294670527866686986D-01,
     6 -0.159609815576893754D+01, -0.497891581567850424D-01,
     7 -0.175461057765584494D+01, -0.701112635269013768D-01,
     8 -0.187247878880251983D+01, -0.904333688970177241D-01,
     9 -0.194970278920896201D+01, -0.110755474267134071D+00,
     A -0.198628257887517146D+01, -0.131077579637250419D+00,
     B -0.198221815780114818D+01, -0.151399685007366752D+00,
     C -0.193750952598689219D+01, -0.171721790377483099D+00,
     D -0.185215668343240347D+01, -0.192043895747599447D+00,
     E -0.172615963013768225D+01, -0.212366001117715794D+00/
 
      DATA (W(N,5),W(N+1,5),N=2,14,2)
     1/ 0.871183254585174982D-01,  0.435591627292587508D-01,
     2  0.217795813646293754D-01,  0.108897906823146873D-01,
     3  0.544489534115734364D-02,  0.272244767057867193D-02,
     4  0.136122383528933596D-02,  0.680611917644667955D-03,
     5  0.340305958822333977D-03,  0.170152979411166995D-03,
     6  0.850764897055834977D-04,  0.425382448527917472D-04,
     7  0.212691224263958736D-04,  0.106345612131979372D-04/
 
      DATA (WP(N,1),WP(N,3),N=2,15)
     1/-0.133196159122085045D+01,  0.445816186556927292D-01,
     2 -0.229218106995884763D+01, -0.240054869684499309D-01,
     3 -0.311522633744855959D+01, -0.925925925925925875D-01,
     4 -0.380109739368998611D+01, -0.161179698216735251D+00,
     5 -0.434979423868312742D+01, -0.229766803840877915D+00,
     6 -0.476131687242798352D+01, -0.298353909465020564D+00,
     7 -0.503566529492455417D+01, -0.366941015089163228D+00,
     8 -0.517283950617283939D+01, -0.435528120713305891D+00,
     9 -0.517283950617283939D+01, -0.504115226337448555D+00,
     A -0.503566529492455417D+01, -0.572702331961591218D+00,
     B -0.476131687242798352D+01, -0.641289437585733882D+00,
     C -0.434979423868312742D+01, -0.709876543209876532D+00,
     D -0.380109739368998611D+01, -0.778463648834019195D+00,
     E -0.311522633744855959D+01, -0.847050754458161859D+00/
 
      RESULT=0
      ABSERR=0
      IFAIL=3
      IF(N .LT. 2 .OR. N .GT. 15) RETURN
      IF(MINPTS .GT. MAXPTS) RETURN
 
      IFNCLS=0
      LDV=.FALSE.
      TWONDM=2**N
      IRGNST=2*N+3
      IRLCLS=2**N+2*N*(N+1)+1
      ISBRGN=IRGNST
      ISBRGS=IRGNST
      IF(MAXPTS .LT. IRLCLS) RETURN
      DO 10 J = 1,N
      CTR(J)=(B(J)+A(J))*HF
   10 WTH(J)=(B(J)-A(J))*HF
 
   20 RGNVOL=TWONDM
      DO 30 J = 1,N
      RGNVOL=RGNVOL*WTH(J)
   30 Z(J)=CTR(J)
      SUM1=F(N,Z)
 
      DIFMAX=0
      SUM2=0
      SUM3=0
      DO 40 J = 1,N
      Z(J)=CTR(J)-XL2*WTH(J)
      F2=F(N,Z)
      Z(J)=CTR(J)+XL2*WTH(J)
      F2=F2+F(N,Z)
      WTHL(J)=XL4*WTH(J)
      Z(J)=CTR(J)-WTHL(J)
      F3=F(N,Z)
      Z(J)=CTR(J)+WTHL(J)
      F3=F3+F(N,Z)
      SUM2=SUM2+F2
      SUM3=SUM3+F3
      DIF=ABS(7*F2-F3-12*SUM1)
      DIFMAX=MAX(DIF,DIFMAX)
      IF(DIFMAX .EQ. DIF) IDVAXN=J
   40 Z(J)=CTR(J)
 
      SUM4=0
      DO 70 J = 2,N
      J1=J-1
      DO 60 K = J,N
      DO 50 L = 1,2
      WTHL(J1)=-WTHL(J1)
      Z(J1)=CTR(J1)+WTHL(J1)
      DO 50 M = 1,2
      WTHL(K)=-WTHL(K)
      Z(K)=CTR(K)+WTHL(K)
   50 SUM4=SUM4+F(N,Z)
   60 Z(K)=CTR(K)
   70 Z(J1)=CTR(J1)
 
      SUM5=0
      DO 80 J = 1,N
      WTHL(J)=-XL5*WTH(J)
   80 Z(J)=CTR(J)+WTHL(J)
   90 SUM5=SUM5+F(N,Z)
      DO 100 J = 1,N
      WTHL(J)=-WTHL(J)
      Z(J)=CTR(J)+WTHL(J)
      IF(WTHL(J) .GT. 0) GO TO 90
  100 CONTINUE
 
      RGNCMP=RGNVOL*(WP(N,1)*SUM1+WP2*SUM2+WP(N,3)*SUM3+WP4*SUM4)
      RGNVAL=W(N,1)*SUM1+W2*SUM2+W(N,3)*SUM3+W4*SUM4+W(N,5)*SUM5
      RGNVAL=RGNVOL*RGNVAL
      RGNERR=ABS(RGNVAL-RGNCMP)
      RESULT=RESULT+RGNVAL
      ABSERR=ABSERR+RGNERR
      IFNCLS=IFNCLS+IRLCLS
 
      IF(LDV) THEN
  110  ISBTMP=2*ISBRGN
       IF(ISBTMP .GT. ISBRGS) GO TO 160
       IF(ISBTMP .LT. ISBRGS) THEN
        ISBTPP=ISBTMP+IRGNST
        IF(WK(ISBTMP) .LT. WK(ISBTPP)) ISBTMP=ISBTPP
       ENDIF
       IF(RGNERR .GE. WK(ISBTMP)) GO TO 160
       DO 130 K = 0,IRGNST-1
  130  WK(ISBRGN-K)=WK(ISBTMP-K)
       ISBRGN=ISBTMP
       GO TO 110
      ENDIF
  140 ISBTMP=(ISBRGN/(2*IRGNST))*IRGNST
      IF(ISBTMP .GE. IRGNST .AND. RGNERR .GT. WK(ISBTMP)) THEN
       DO 150 K = 0,IRGNST-1
  150  WK(ISBRGN-K)=WK(ISBTMP-K)
       ISBRGN=ISBTMP
       GO TO 140
      ENDIF
 
  160 WK(ISBRGN)=RGNERR
      WK(ISBRGN-1)=RGNVAL
      WK(ISBRGN-2)=IDVAXN
      DO 170 J = 1,N
      ISBTMP=ISBRGN-2*J-2
      WK(ISBTMP+1)=CTR(J)
  170 WK(ISBTMP)=WTH(J)
      IF(LDV) THEN
       LDV=.FALSE.
       CTR(IDVAX0)=CTR(IDVAX0)+2*WTH(IDVAX0)
       ISBRGS=ISBRGS+IRGNST
       ISBRGN=ISBRGS
       GO TO 20
      ENDIF
      RELERR=ABSERR/ABS(RESULT)
      IF(ISBRGS+IRGNST .GT. IWK) IFAIL=2
      IF(IFNCLS+2*IRLCLS .GT. MAXPTS) IFAIL=1
      IF(RELERR .LT. EPS .AND. IFNCLS .GE. MINPTS) IFAIL=0
      IF(IFAIL .EQ. 3) THEN
       LDV=.TRUE.
       ISBRGN=IRGNST
       ABSERR=ABSERR-WK(ISBRGN)
       RESULT=RESULT-WK(ISBRGN-1)
       IDVAX0=WK(ISBRGN-2)
       DO 190 J = 1,N
       ISBTMP=ISBRGN-2*J-2
       CTR(J)=WK(ISBTMP+1)
  190  WTH(J)=WK(ISBTMP)
       WTH(IDVAX0)=HF*WTH(IDVAX0)
       CTR(IDVAX0)=CTR(IDVAX0)-WTH(IDVAX0)
       GO TO 20
      ENDIF
      NFNEVL=IFNCLS
      RETURN
      END


      subroutine qgaus(func,a,b,ss)
      real a,b,ss,func
      external func
      integer j
      real dx,xm,xr,w(5),x(5)
      save w,x
      data w/0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 
     &       0.0666713443/
      data x/0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666,
     &       0.9739065285/

      xm = 0.5*(b+a)
      xr = 0.5*(b-a)
      ss = 0.0
      do j =1,5
         dx = xr*x(j)
         ss = ss + w(j)*(func(xm+dx)+func(xm-dx))
      enddo
      ss = xr*ss
      return
      end





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
*DECK TAGX
      SUBROUTINE TAGX(FSIG, SED, XXX, QQ, ALR, PTR, IPN)
C
C     DEUTERON TAGGED ELECTROPRODUCTION CROSS SECTION
C
C     FSIG  TAGGED ELECTROPRODUCTION CROSS SECTION
C           DSIGMA = FSIG *DX *DQ**2 *(D**3 PR/ER)
C
C           DIMENSION OF DSIGMA IS NANOBARN
C               ''       FSIG      NANOBARN/GEV**4
C
C     SED   ELECTRON-DEUTERON S-INVARIANT (SQUARED CM ENERGY)
C     X     NUCLEON SCALING VARIABLE X = 2*XD
C     QQ    4-MOMENTUM TRANSFER Q**2
C     ALR   RECOIL MOMENTUM, LC FRACTION 
C     PTR   RECOIL MOMENTUM, TRANSVERSE
C
C     IPN = 1   PROTON  ACTIVE, NEUTRON TAGGED
C           2   NEUTRON ACTIVE, PROTON  TAGGED
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      PARAMETER (PI = 0.31415 92653 58979 D+01)
C
C     ...FINE STRUCTURE CONSTANT
C
      PARAMETER (ALEM = 1.D0/137.D0)
C
C     ...NUCLEON AND DEUTERON MASS
C        UN   AVERAGE NUCLEON MASS (M(PROTON) + M(NEUTRON))/2   
C        ED   DEUTERON BINDING ENERGY
C        DEUTERON MASS CALCULATED AS  2*UN - ED
C
      PARAMETER (UN = 0.939D0, ED = 0.00222, UD = 2*UN - ED)
C
C     ...CONVERSION FACTOR FOR CROSS SECTION UNITS
C        1/GEV**2  =  0.398 *10**6  NANOBARN
C
      PARAMETER (CONVNB = 0.389D6)
C
      XD   = XXX/2
C
      Y    = QQ/XD/(SED - UD**2)
C
      CORR = (Y*XD*UD)**2/QQ 
      EPS  = (1 - Y - CORR)/(1 - Y + Y**2/2 + CORR) 
C
C     ...FACTOR FOR ELECTROPRODUCTION CROSS SECTION
C        EXTRA FACTOR 1/2 FROM DIFFERENTIAL DXD = DX/2 
C        CONVERSION GEV**(-2) TO NANOBARN
C
      FAC  = 2*PI *ALEM**2 *Y**2 /QQ**2 /(1 - EPS)
      FACD = FAC /2.D0 *CONVNB
C
C     ...STRUCTURE FUNCTIONS
C
      CALL TAGFD(F2D, XXX, QQ, ALR, PTR, IPN)
C
C     ...CROSS SECTION
C
      FSIG = FACD*(F2D/XD)
C
      END
C
C---------------------------------------------------------------------
C
C
*DECK TAGSP
      SUBROUTINE TAGSP(S, ALR, PTR)
C
C     DEUTERON LIGHT-CONE SPECTRAL FUNCTION
C     NORMALIZED SUCH THAT  INTEGRAL (D ALR/ALR) (D**2 PTR) S = 1
C
C     ALR   RECOIL MOMENTUM, LC FRACTION 
C     PTR   RECOIL MOMENTUM, TRANSVERSE
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      PARAMETER (UN = 0.939D0)
C
C     ...EFFECTIVE CM MOMENTUM
C
      P = SQRT((PTR**2 + UN**2)/ALR/(2 - ALR) - UN**2)
C
C     ...SPECTRAL FUNCTION
C
      CALL TAGWF(PSI, P)
C
      S = SQRT(P**2 + UN**2) *PSI**2 /(2 - ALR)
C
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE TAGRES(RES, ALR)
C
C     PACKAGE TAG -- DEUTERON DIS WITH SPECTATOR TAGGING
C     AUTHOR C. WEISS (WEISS.AT.JLAB.ORG)
C
C     USER-DEFINED ROUTINE
C     DEUTERON LIGHT-CONE SPECTRAL FUNCTION
C     RESIDUE AT NUCLEON POLE TPRIME = 0
C     S = RES/(-TPRIME)**2 + LESS SINGULAR
C
C     INPUT:
C     ALR    RECOIL MOMENTUM LIGHT-CONE FRACTION ALPHAR
C
C     OUTPUT:
C     RES    RESIDUE OF SPECTRAL FUNCTION (GEV UNITS)
C
C     V1 13MAY14
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      PARAMETER (PI = 0.31415 92653 58979 D+01)
      PARAMETER (UN = 0.939D0)
C
C     ...HULTHEN WAVE FUNCTION
C
      A = 0.045647
      B = 0.2719
C
      C = PI**2*(A - B)**2/A/B/(A + B)
C
      EPOL = SQRT(-A**2 + UN**2)
C
C     ...RESIDUE
C
      RES = 4*EPOL/C/ALR
C
      END
C

      SUBROUTINE TAGFN(F2NN, XXX, QQ, IPN)
C
C     PACKAGE TAG -- DEUTERON DIS WITH SPECTATOR TAGGING
C     AUTHOR C. WEISS (WEISS.AT.JLAB.ORG)
C
C     USER-DEFINED ROUTINE
C     NUCLEON STRUCTURE FUNCTION F2
C
C     INPUT:
C     X     X VARIABLE FOR NUCLEON
C     QQ    Q2 (GEV**2)
C     IPN   SWITCH: 1  PROTON, 2  NEUTRON
C
C     OUTPUT:
C     F2    NUCLEON STRUCTURE FUNCTION F2
C
C     V1 13MAY14
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C

      double precision JR14NLO08SFF2n
      double precision F2n(0:38),eF2n
      integer set

C     ...QUARK CHARGES
C
      EU   =  0.6666666666666D0
      ED   = -0.3333333333333D0
      ES   = -0.3333333333333D0
C
C     ...GRV98 PDF PARAMETRIZATION
C


      ISET = 1
      CALL GRV98PA(ISET, XXX, QQ, UV, DV, US, DS, SS, GL)

C
      IF (IPN.EQ.1)      THEN
         F2NN = EU**2*(UV + 2*US)
     *      + ED**2*(DV + 2*DS)
     *      + ES**2*(     2*SS)
      ELSE IF (IPN.EQ.2) THEN
         F2NN = EU**2*(DV + 2*DS)
     *      + ED**2*(UV + 2*US)
     *      + ES**2*(     2*SS)
      ELSE IF (IPN.EQ.3) THEN
c     Pedro Jimenez-Delgado
         do set = 0,38
            F2n(set) = JR14NLO08SFF2n(XXX,QQ,set)
         enddo
         F2NN = F2n(0)
c$$$         eF2n = dsqrt(eF2n)
       ENDIF
c$$$       if(IPN.eq.3)then
c$$$          F2NN=F2NN3
c$$$       endif
       IF(IPN.lt.3)RETURN
c       print *,IPN,F2NN
C
      END
C
C---------------------------------------------------------------------
C

*DECK TAGFD
      SUBROUTINE TAGFD(F2, XXX, QQ, ALR, PTR, IPN)
C
C     DEUTERON TAGGED STRUCTURE FUNCTION IN IMPULSE APPROXIMATION
C
C     F2    TAGGED STRUCTURE FUNCTION
C
C     X     NUCLEON SCALING VARIABLE X = 2*XD
C     QQ    4-MOMENTUM TRANSFER Q**2
C     ALR   RECOIL MOMENTUM, LC FRACTION 
C     PTR   RECOIL MOMENTUM, TRANSVERSE
C
C     IPN = 1   PROTON  ACTIVE, NEUTRON TAGGED
C           2   NEUTRON ACTIVE, PROTON  TAGGED
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
C     ...SPECTRAL FUNCTION
C
      CALL TAGSP(S, ALR, PTR)
C
C     ...STRUCTURE FUNCTION
C
      XEFF = XXX/(2 - ALR)
C
      CALL TAGFN(FN2, XEFF, QQ, IPN)
C
      F2 = S*FN2
C
      END
C
C---------------------------------------------------------------------

C
C---------------------------------------------------------------------
C
*DECK TAGWF
      SUBROUTINE TAGWF(PSI, P)
C
C     PACKAGE TAG -- DEUTERON DIS WITH SPECTATOR TAGGING
C     AUTHOR C. WEISS (WEISS.AT.JLAB.ORG)
C
C     USER-DEFINED ROUTINE
C     DEUTERON NON-RELATIVISTIC WAVE FUNCTION PSI(P)
C     NORMALIZED SUCH THAT   INTEGRAL (D**3 P) PSI(P)**2 = 1
C
C     INPUT:
C     P    MODULUS OF 3-MOMENTUM (GEV)
C
C     OUTPUT:
C     PSI  DEUTERON NON-RELATIVISTIC WAVE FUNCTION PSI(P)
C
C     ALL DIMENSIONFUL QUANTITIES IN GEV UNITS
C
C     V1 13MAY14
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      PARAMETER (PI = 0.31415 92653 58979 D+01)
C
C     ...HULTHEN WAVE FUNCTION
C
      A = 0.045647
      B = 0.2719
C
      C = PI**2*(A - B)**2/A/B/(A + B)
C
      PSI = (  1.D0/(P**2 + A**2)
     *       - 1.D0/(P**2 + B**2))/SQRT(C)
C
      END
C


*DECK GRV98
      BLOCK DATA GRV98
C
C     INITIALIZATION OF GRV98 PDF PARAMETRIZATION
C
C     NOTE: IF MULTIPLE PDF SETS ARE USED IN THE SAME RUN,
C     THE INITIALIZATION MUST BE RESET BY THE CALLING PROGRAM.
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      COMMON /INTINIP/ IINIP
      COMMON /INTINIF/ IINIF
      DATA IINIP /0/
      DATA IINIF /0/
C
      END

C---------------------------------------------------------------------
C
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                   *
*     G R V  -  P R O T O N  - P A R A M E T R I Z A T I O N S      *
*                                                                   *
*                          1998 UPDATE                              *
*                                                                   *
*                  For a detailed explanation see                   *
*                   M. Glueck, E. Reya, A. Vogt :                   *
*        hep-ph/9806404  =  DO-TH 98/07  =  WUE-ITP-98-019          *
*                  (To appear in Eur. Phys. J. C)                   *
*                                                                   *
*   This package contains subroutines returning the light-parton    *
*   distributions in NLO (for the MSbar and DIS schemes) and LO;    * 
*   the respective light-parton, charm, and bottom contributions    *
*   to F2(electromagnetic); and the scale dependence of alpha_s.    *
*                                                                   *
*   The parton densities and F2 values are calculated from inter-   *
*   polation grids covering the regions                             *
*         Q^2/GeV^2  between   0.8   and  1.E6 ( 1.E4 for F2 )      *
*            x       between  1.E-9  and   1.                       *
*   Any call outside these regions stops the program execution.     *
*                                                                   *
*   At Q^2 = MZ^2, alpha_s reads  0.114 (0.125) in NLO (LO); the    *
*   heavy quark thresholds, QH^2 = mh^2, in the beta function are   *
*            mc = 1.4 GeV,  mb = 4.5 GeV,  mt = 175 GeV.            *
*   Note that the NLO alpha_s running is different from GRV(94).    * 
*                                                                   *
*    Questions, comments etc to:  avogt@physik.uni-wuerzburg.de     *
*                                                                   *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
*
*
      SUBROUTINE GRV98PA (ISET, XXX, Q2, UV, DV, US, DS, SS, GL)
*********************************************************************
*                                                                   *
*   THE PARTON ROUTINE.                                             *
*                                     __                            *
*   INPUT:   ISET =  1 (LO),  2 (NLO, MS), or  3 (NLO, DIS)         *
*            X  =  Bjorken-x        (between  1.E-9 and 1.)         *
*            Q2 =  scale in GeV**2  (between  0.8 and 1.E6)         *
*                                                                   *
*   OUTPUT:  UV = u - u(bar),  DV = d - d(bar),  US = u(bar),       *
*            DS = d(bar),  SS = s = s(bar),  GL = gluon.            *
*            Always x times the distribution is returned.           *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINIP / IINIP , and the     *
*            integer variable  IINIP  has always to be zero when    *
*            GRV98PA is called for the first time or when  ISET     *
*            has been changed.                                      *
*                                                                   *
*   GRIDS:   1. grv98lo.grid, 2. grv98nlm.grid, 3. grv98nld.grid,   *
*            (1+1809 lines with 6 columns, 4 significant figures)   *
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NPART=6, NX=68, NQ=27, NARG=2)
      DIMENSION XUVF(NX,NQ), XDVF(NX,NQ), XDEF(NX,NQ), XUDF(NX,NQ),
     1          XSF(NX,NQ), XGF(NX,NQ), PARTON (NPART,NQ,NX-1), 
     2          QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      CHARACTER*80 LINE
      COMMON / INTINIP / IINIP
      SAVE XUVF, XDVF, XDEF, XUDF, XSF, XGF, NA, ARRF
*
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 0.8E0, 
     1           1.0E0, 1.3E0, 1.8E0, 2.7E0, 4.0E0, 6.4E0,
     2           1.0E1, 1.6E1, 2.5E1, 4.0E1, 6.4E1,
     3           1.0E2, 1.8E2, 3.2E2, 5.7E2,
     4           1.0E3, 1.8E3, 3.2E3, 5.7E3,
     5           1.0E4, 2.2E4, 4.6E4,
     6           1.0E5, 2.2E5, 4.6E5, 
     7           1.E6 /
       DATA XB / 1.0E-9, 1.8E-9, 3.2E-9, 5.7E-9, 
     A           1.0E-8, 1.8E-8, 3.2E-8, 5.7E-8, 
     B           1.0E-7, 1.8E-7, 3.2E-7, 5.7E-7, 
     C           1.0E-6, 1.4E-6, 2.0E-6, 3.0E-6, 4.5E-6, 6.7E-6,
     1           1.0E-5, 1.4E-5, 2.0E-5, 3.0E-5, 4.5E-5, 6.7E-5,
     2           1.0E-4, 1.4E-4, 2.0E-4, 3.0E-4, 4.5E-4, 6.7E-4,
     3           1.0E-3, 1.4E-3, 2.0E-3, 3.0E-3, 4.5E-3, 6.7E-3,
     4           1.0E-2, 1.4E-2, 2.0E-2, 3.0E-2, 4.5E-2, 0.06, 0.08,
     5           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     6           0.3, 0.325, 0.35, 0.375, 0.4,  0.45, 0.5, 0.55,
     7           0.6, 0.65,  0.7,  0.75,  0.8,  0.85, 0.9, 0.95, 1. /
*
*...CHECK OF X AND Q2 VALUES : 



      IF ( (XXX.LT.0.99D-9) .OR. (XXX.GT.1.D0) ) THEN
         WRITE(6,91) 
  91     FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
         STOP
      ENDIF
      IF ( (Q2.LT.0.799) .OR. (Q2.GT.1.01E6) ) THEN
         WRITE(6,92) 
  92     FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
         STOP
      ENDIF
      IF (IINIP .NE. 0) GOTO 16
*
*...INITIALIZATION, IF REQUIRED :
*
*    SELECTION AND READING OF THE GRID : 
*    (COMMENT: FIRST NUMBER IN THE FIRST LINE OF THE GRID)
      IF (ISET .EQ. 1) THEN
        OPEN (11,FILE='grv98lo.grid',STATUS='old')   ! 7.332E-05
      ELSE IF (ISET .EQ. 2) THEN
        OPEN (11,FILE='grv98nlm.grid',STATUS='old')  ! 1.015E-04
      ELSE IF (ISET .EQ. 3) THEN
        OPEN (11,FILE='grv98nld.grid',STATUS='old')  ! 1.238E-04
      ELSE
        WRITE(6,93)
  93    FORMAT (2X,'NO OR INVALID PARTON SET CHOICE')
        STOP
      END IF
      IINIP = 1
      READ(11,89) LINE
  89  FORMAT(A80)
      DO 15 M = 1, NX-1 
      DO 15 N = 1, NQ
      READ(11,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1            PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M) 
  90  FORMAT (6(1PE10.3))
  15  CONTINUE
      CLOSE(11)
*
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0V = XB(IX)**0.5 
        XB0S = XB(IX)**(-0.2) 
        XB1 = 1.-XB(IX)
        XUVF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0V)
        XDVF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0V)
        XDEF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**7 * XB0V) 
        XUDF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**7 * XB0S)
        XSF(IX,IQ)  = PARTON(5,IQ,IX) / (XB1**7 * XB0S)
        XGF(IX,IQ)  = PARTON(6,IQ,IX) / (XB1**5 * XB0S)
  20  CONTINUE
        XUVF(NX,IQ) = 0.E0
        XDVF(NX,IQ) = 0.E0
        XDEF(NX,IQ) = 0.E0
        XUDF(NX,IQ) = 0.E0
        XSF(NX,IQ)  = 0.E0
        XGF(NX,IQ)  = 0.E0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
*
*...CONTINUATION, IF INITIALIZATION WAS DONE PREVIOUSLY.
*
  16  CONTINUE
*
*...INTERPOLATION :
      XT(1) = DLOG(XXX)
      XT(2) = DLOG(Q2)
      X1 = 1.- XXX
      XV = XXX**0.5
      XS = XXX**(-0.2)
      UV = FINT(NARG,XT,NA,ARRF,XUVF) * X1**3 * XV
      DV = FINT(NARG,XT,NA,ARRF,XDVF) * X1**4 * XV
      DE = FINT(NARG,XT,NA,ARRF,XDEF) * X1**7 * XV
      UD = FINT(NARG,XT,NA,ARRF,XUDF) * X1**7 * XS
      US = 0.5 * (UD - DE)
      DS = 0.5 * (UD + DE)
      SS = FINT(NARG,XT,NA,ARRF,XSF)  * X1**7 * XS
      GL = FINT(NARG,XT,NA,ARRF,XGF)  * X1**5 * XS 
*
 60   RETURN
      END
*
*
*
*
      SUBROUTINE GRV98F2 (ISET, XXX, Q2, F2L, F2C, F2B, F2)
*********************************************************************
*                                                                   *
*   THE F2 ROUTINE.                                                 *
*                                                                   *
*   INPUT :  ISET = 1 (LO),  2 (NLO).                               *
*            X  =  Bjorken-x        (between  1.E-9 and 1.)         *
*            Q2 =  scale in GeV**2  (between  0.8 and 1.E4)         *
*                                                                   *
*   OUTPUT:  F2L = F2(light), F2C = F2(charm), F2B = F2(bottom,)    *
*            F2  = sum, all given for electromagnetic proton DIS.   *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINIF / IINIF , and the     *
*            integer variable  IINIF  has always to be zero when    *
*            GRV98F2 is called for the first time or when  ISET     *
*            has been changed.                                      *
*                                                                   *
*   GRIDS:   1. grv98lof.grid, 2. grv98nlf.grid.                    *
*            (1+1407 lines with 3 columns, 4 significant figures)   *
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NSTRF=3, NX=68, NQ=21, NARG=2)
      DIMENSION XF2LF(NX,NQ), XF2CF(NX,NQ), XF2BF(NX,NQ), 
     1          STRFCT (NSTRF,NQ,NX-1), QS(NQ), XB(NX), 
     3          XT(NARG), NA(NARG), ARRF(NX+NQ) 
      CHARACTER*80 LINE
      COMMON / INTINIF / IINIF
      SAVE XF2LF, XF2CF, XF2BF, NA, ARRF
*
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 0.8E0, 
     1           1.0E0, 1.3E0, 1.8E0, 2.7E0, 4.0E0, 6.4E0,
     2           1.0E1, 1.6E1, 2.5E1, 4.0E1, 6.4E1,
     3           1.0E2, 1.8E2, 3.2E2, 5.7E2,
     4           1.0E3, 1.8E3, 3.2E3, 5.7E3,
     5           1.0E4/ 
       DATA XB / 1.0E-9, 1.8E-9, 3.2E-9, 5.7E-9, 
     A           1.0E-8, 1.8E-8, 3.2E-8, 5.7E-8, 
     B           1.0E-7, 1.8E-7, 3.2E-7, 5.7E-7, 
     C           1.0E-6, 1.4E-6, 2.0E-6, 3.0E-6, 4.5E-6, 6.7E-6,
     1           1.0E-5, 1.4E-5, 2.0E-5, 3.0E-5, 4.5E-5, 6.7E-5,
     2           1.0E-4, 1.4E-4, 2.0E-4, 3.0E-4, 4.5E-4, 6.7E-4,
     3           1.0E-3, 1.4E-3, 2.0E-3, 3.0E-3, 4.5E-3, 6.7E-3,
     4           1.0E-2, 1.4E-2, 2.0E-2, 3.0E-2, 4.5E-2, 0.06, 0.08,
     5           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     6           0.3, 0.325, 0.35, 0.375, 0.4,  0.45, 0.5, 0.55,
     7           0.6, 0.65,  0.7,  0.75,  0.8,  0.85, 0.9, 0.95, 1. /
*
*...CHECK OF X AND Q2 VALUES : 
      IF ( (XXX.LT.0.99D-9) .OR. (XXX.GT.1.D0) ) THEN
         WRITE(6,91) 
  91     FORMAT (2X,'STR.FCT. INTERPOLATION: X OUT OF RANGE')
         STOP
      ENDIF
      IF ( (Q2.LT.0.799) .OR. (Q2.GT.1.01E4) ) THEN
         WRITE(6,92) 
  92     FORMAT (2X,'STR.FCT. INTERPOLATION: Q2 OUT OF RANGE')
         STOP
      ENDIF
      IF (IINIF .NE. 0) GOTO 16
*
*...INITIALIZATION, IF REQUIRED :
*
*    SELECTION AND READING OF THE GRID : 
*    (COMMENT: FIRST NUMBER IN THE FIRST LINE OF THE GRID)
      IF (ISET .EQ. 1) THEN
        OPEN (11,FILE='grv98lof.grid',STATUS='old')  !  7.907E-01
      ELSE IF (ISET.EQ.2) THEN
        OPEN (11,FILE='grv98nlf.grid',STATUS='old')  !  9.368E-01
      ELSE
        WRITE(6,93)
  93    FORMAT (2X,'NO OR INVALID STR.FCT. SET CHOICE')
        STOP
      END IF
      IINIF = 1
      READ(11,89) LINE
  89  FORMAT(A80)
      DO 15 M = 1, NX-1 
      DO 15 N = 1, NQ
      READ(11,90) STRFCT(1,N,M), STRFCT(2,N,M), STRFCT(3,N,M) 
  90  FORMAT (3(1PE10.3))
  15  CONTINUE
      CLOSE(11)
*
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XBS = XB(IX)**0.2 
        XB1 = 1.-XB(IX)
        XF2LF(IX,IQ) = STRFCT(1,IQ,IX) / XB1**3 * XBS
        XF2CF(IX,IQ) = STRFCT(2,IQ,IX) / XB1**7 * XBS
        XF2BF(IX,IQ) = STRFCT(3,IQ,IX) / XB1**7 * XBS 
  20  CONTINUE
        XF2LF(NX,IQ) = 0.E0
        XF2CF(NX,IQ) = 0.E0
        XF2BF(NX,IQ) = 0.E0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
*
*...CONTINUATION, IF INITIALIZATION WAS DONE PREVIOUSLY.
*
  16  CONTINUE
*
*...INTERPOLATION :
      XT(1) = DLOG(XXX)
      XT(2) = DLOG(Q2)
      X1 = 1.- XXX
      XS = XXX**(-0.2)
      F2L = FINT(NARG,XT,NA,ARRF,XF2LF) * X1**3 * XS
      F2C = FINT(NARG,XT,NA,ARRF,XF2CF) * X1**7 * XS
      F2B = FINT(NARG,XT,NA,ARRF,XF2BF) * X1**7 * XS
      F2  = F2L + F2C + F2B
*
 60   RETURN
      END
*
*
*
*
      FUNCTION FINT(NARG,ARG,NENT,ENT,TABLE)
*********************************************************************
*                                                                   *
*   THE INTERPOLATION ROUTINE (CERN LIBRARY ROUTINE E104)           *
*                                                                   *
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ARG(5),NENT(5),ENT(10),TABLE(10)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FINT=0.
   10 FAC=1.
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FINT=FINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END
*
*
*
*
      FUNCTION ALPHAS (Q2, NAORD)
*********************************************************************
*                                                                   *
*   THE ALPHA_S ROUTINE.                                            *
*                                                                   *
*   INPUT :  Q2    =  scale in GeV**2  (not too low, of course);    *
*            NAORD =  1 (LO),  2 (NLO).                             *
*                                                                   *
*   OUTPUT:  alphas_s/(4 pi) for use with the GRV(98) partons.      *  
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A - Z)
      INTEGER NF, K, I, NAORD
      DIMENSION LAMBDAL (3:6),  LAMBDAN (3:6), Q2THR (3)
*
*...HEAVY QUARK THRESHOLDS AND LAMBDA VALUES :
      DATA Q2THR   /  1.960,  20.25,  30625. /
      DATA LAMBDAL / 0.2041, 0.1750, 0.1320, 0.0665 /
      DATA LAMBDAN / 0.2994, 0.2460, 0.1677, 0.0678 /
*
*...DETERMINATION OF THE APPROPRIATE NUMBER OF FLAVOURS :
      NF = 3
      DO 10 K = 1, 3
      IF (Q2 .GT. Q2THR (K)) THEN
         NF = NF + 1
      ELSE
          GO TO 20
       END IF
  10   CONTINUE
*
*...LO ALPHA_S AND BETA FUNCTION FOR NLO CALCULATION :
  20   B0 = 11.- 2./3.* NF
       B1 = 102.- 38./3.* NF
       B10 = B1 / (B0*B0)
       IF (NAORD .EQ. 1) THEN
         LAM2 = LAMBDAL (NF) * LAMBDAL (NF)
         ALP  = 1./(B0 * DLOG (Q2/LAM2))
         GO TO 1
       ELSE IF (NAORD .EQ. 2) then
         LAM2 = LAMBDAN (NF) * LAMBDAN (NF)
         B1 = 102.- 38./3.* NF
         B10 = B1 / (B0*B0)
       ELSE
         WRITE (6,91)
  91     FORMAT ('INVALID CHOICE FOR ORDER IN ALPHA_S')
         STOP
       END IF
*
*...START VALUE FOR NLO ITERATION :
       LQ2 = DLOG (Q2 / LAM2)
       ALP = 1./(B0*LQ2) * (1.- B10*DLOG(LQ2)/LQ2)
*
*...EXACT NLO VALUE, FOUND VIA NEWTON PROCEDURE :
       DO 2 I = 1, 6
       XL  = DLOG (1./(B0*ALP) + B10)
       XLP = DLOG (1./(B0*ALP*1.01) + B10)
       XLM = DLOG (1./(B0*ALP*0.99) + B10)
       Y  = LQ2 - 1./ (B0*ALP) + B10 * XL
       Y1 = (- 1./ (B0*ALP*1.01) + B10 * XLP
     1       + 1./ (B0*ALP*0.99) - B10 * XLP) / (0.02D0*ALP)
       ALP = ALP - Y/Y1
  2    CONTINUE
*
*...OUTPUT :
  1    ALPHAS = ALP
       RETURN
       END