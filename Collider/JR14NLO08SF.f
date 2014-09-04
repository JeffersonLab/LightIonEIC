!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! JR14NLO08SF (Phys. Rev. D89 (2014) 074049 [arXiv:1403.1852])
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This package contains the JR14 NLO(MS-bar) dynamical (Q02 = 0.8 GeV2)
!! predictions for (neutral current, photon) structure functions (SF) of
!! the proton and the neutron.
!!
!! The grids are generated for 10^-5 <= x <= 1 and 2 <= Q2 <= 1000 (GeV2)
!! Outside these ranges the output is obtained through extrapolation
!! (of the fits results to 1 GeV2 <= Q2 <= 2 GeV2, and numerical beyond that)
!! and should NOT be used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The routines use a modification of the standard multidimensional
!! linear interpolation routine FINT (CERNLIB E104) distributed as dfint.f
!!
!! The file './JR14NLO08SF.grd', where "./" means "path from the working
!! directory to the file", is read.
!!
!! For questions, comments, problems etc please contact: pedro@jlab.org
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! JR14NLO08SFinit:
!!  Initialization routine of the package to be called (only once)
!!  before using any of the other subroutines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! JR14NLO08SF'function'(x,Q2,set) with 'function' = F2p,F2n,FLp,FLn:
!!    x == Bjorken-x
!!    Q2 == Q2 (GeV2)
!!    set == set to be used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set == Index specifying the set to be used:
!!        0 central value
!!        1,2,...,38 set corresponding to a displacement +1 (Dchi2=1) 
!!        from set 0 in the direction of the ith eigenvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      block data JR14NLO08SF
       implicit none
       integer shape(2)
       double precision grid(145)
       common /JR14NLO08SFgrid/ grid,shape
       data shape /114,31/
       data grid
     &   /1d-5,1.1d-5,1.4d-5,1.6d-5,1.8d-5,2d-5,2.25d-5,2.5d-5,2.8d-5,
     &    3.2d-5,3.6d-5,4d-5,4.5d-5,5d-5,5.6d-5,6.3d-5,7d-5,8d-5,9d-5,
     &    1d-4,1.1d-4,1.4d-4,1.6d-4,1.8d-4,2d-4,2.25d-4,2.5d-4,2.8d-4,
     &    3.2d-4,3.6d-4,4d-4,4.5d-4,5d-4,5.6d-4,6.3d-4,7d-4,8d-4,9d-4,
     &    1d-3,1.1d-3,1.4d-3,1.6d-3,1.8d-3,2d-3,2.25d-3,2.5d-3,2.8d-3,
     &    3.2d-3,3.6d-3,4d-3,4.5d-3,5d-3,5.6d-3,6.3d-3,7d-3,8d-3,9d-3,
     &    1d-2,1.1d-2,1.4d-2,1.6d-2,1.8d-2,2d-2,2.25d-2,2.5d-2,2.8d-2,
     &    3.2d-2,3.6d-2,4d-2,4.5d-2,5d-2,5.6d-2,6.3d-2,7d-2,8d-2,9d-2,
     &    0.10d0,0.125d0,0.15d0,0.175d0,0.20d0,0.225d0,0.25d0,0.275d0,
     &    0.30d0,0.325d0,0.35d0,0.375d0,0.40d0,0.425d0,0.45d0,0.475d0,
     &    0.50d0,0.525d0,0.55d0,0.575d0,0.60d0,0.625d0,0.65d0,0.675d0,
     &    0.70d0,0.725d0,0.75d0,0.775d0,0.80d0,0.825d0,0.85d0,0.875d0,
     &    0.9d0,0.92d0,0.94d0,0.96d0,0.98d0,1d0,
     &    1d0,1.25d0,1.6d0,2d0,2.5d0,3.16d0,4d0,5d0,6.3d0,8d0,
     &    1d1,1.25d1,1.6d1,2d1,2.5d1,3.16d1,4d1,5d1,6.3d1,8d1,
     &    1d2,1.25d2,1.6d2,2d2,2.5d2,3.16d2,4d2,5d2,6.3d2,8d2,
     &    1d3/
      end block data JR14NLO08SF

      subroutine JR14NLO08SFinit
       implicit none
       integer shape(2),i,j,k
       double precision grid(145),
     &                  F2p(114,31,0:38),F2n(114,31,0:38),
     &                  FLp(114,31,0:38),FLn(114,31,0:38)
       common /JR14NLO08SFgrid/ grid,shape
       common /JR14NLO08SFF2pc/ F2p
       common /JR14NLO08SFF2nc/ F2n
       common /JR14NLO08SFFLpc/ FLp
       common /JR14NLO08SFFLnc/ FLn
       open(10,file='./JR14NLO08SF.grd')
       do 10 k=0,38
          do 11 j=1,31
             do 12 i=1,114
                read(10,*) F2p(i,j,k)
   12        continue
   11     continue
   10  continue
       do 20 k=0,38
          do 21 j=1,31
             do 22 i=1,114
                read(10,*) F2n(i,j,k)
   22        continue
   21     continue
   20  continue
       do 30 k=0,38
          do 31 j=1,31
             do 32 i=1,114
                read(10,*) FLp(i,j,k)
   32        continue
   31     continue
   30  continue
       do 40 k=0,38
          do 41 j=1,31
             do 42 i=1,114
                read(10,*) FLn(i,j,k)
   42        continue
   41     continue
   40  continue
       close(10)
       return
      end subroutine JR14NLO08SFinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR14NLO08SFF2p(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(145),F2p(114,31,0:38),arg(2),x,Q2,dfint
       common /JR14NLO08SFgrid/ grid,shape
       common /JR14NLO08SFF2pc/ F2p
       arg(1) = x
       arg(2) = Q2
       JR14NLO08SFF2p = dfint(2,arg,shape,grid,F2p(1,1,set))
       return
      end function JR14NLO08SFF2p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR14NLO08SFF2n(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(145),F2n(114,31,0:38),arg(2),x,Q2,dfint
       common /JR14NLO08SFgrid/ grid,shape
       common /JR14NLO08SFF2nc/ F2n
       arg(1) = x
       arg(2) = Q2
       JR14NLO08SFF2n = dfint(2,arg,shape,grid,F2n(1,1,set))
       return
      end function JR14NLO08SFF2n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR14NLO08SFFLp(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(145),FLp(114,31,0:38),arg(2),x,Q2,dfint
       common /JR14NLO08SFgrid/ grid,shape
       common /JR14NLO08SFFLpc/ FLp
       arg(1) = x
       arg(2) = Q2
       JR14NLO08SFFLp = dfint(2,arg,shape,grid,FLp(1,1,set))
       return
      end function JR14NLO08SFFLp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function JR14NLO08SFFLn(x,Q2,set)
       implicit none
       integer shape(2),set
       double precision grid(145),FLn(114,31,0:38),arg(2),x,Q2,dfint
       common /JR14NLO08SFgrid/ grid,shape
       common /JR14NLO08SFFLnc/ FLn
       arg(1) = x
       arg(2) = Q2
       JR14NLO08SFFLn = dfint(2,arg,shape,grid,FLn(1,1,set))
       return
      end function JR14NLO08SFFLn
