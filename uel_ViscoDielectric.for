************************************************************************
!
! User element (UEL) for coupled large visco-elastic deformation and
!  dielectric behavior in 2D and 3D.  This is for plane strain,
!  axisymetric, and 3D.
!
! This particular UEL corresponds to the paper:
!  S. Wang, M. Decker, D.L. Henann, and S.A. Chester.
!  Modeling of dielectric viscoelastomers with application to
!  electromechanical instabilities.  Journal of the Mechanics
!  and Physics of Solids.
!
!
! Solution variables (or nodal variables) are the displacements and the
!  electric potential.  Depending on the particular analysis DOF's 1,2,3
!  are the displacements, and DOF 11 is the electric potential.
!
!
! Material behavior is Neo-Hookean or Arruda-Boyce rubber elasticity 
!  with ``nVisco'' finite-deformation Maxwell elements in parallel
!  and ideal (constant permittivity) dielectric behavior. (To make 
!  the material Neo-Hookean set ``lamL'' less than zero.)
!
! Mechanical behavior:
!
!                        |
!             ------------------------
!             |      |     |         |    Maxwell elements:
!             |      |     |         |  
!  Rubber     |      |     |         |
!  Elastic    |      \     \         \    Springs: Gneq(nVisco)
!  spring:    \      /     /         /    
!   Gshear    /      |     |         |    
!   Kbulk     \      |     |   ...   |  
!             /      |     |         |    Linearly viscous
!             |    | - | | - |     | - |   dashpots: tau(nVisco)
!             |    |___| |___|     |___|
!             |      |     |         |  
!             ------------------------
!                        |
!
!
!
!  The specific forms for the surface tension follow from the paper:
!     1) Modeling of elasto-capillary phenomena. D.L. Henann, and 
!        K. Bertoldi, 2014. Soft Matter.
!
!  The basic element formulation follows from the paper:
!     1) A finite element implementation of a coupled diffusion-
!        deformation theory for elastomeric gels. S.A. Chester,
!        C.V. Di Leo, and L. Anand, 2015. International Journal
!        of Solids and Structures.
!
!  The specific forms for viscoelasticity follow from the paper:
!     1) A micromechanically motivated diffusion-based transient network
!        model and its incorporation into finite rubber viscoelasticty.
!        C. Linder, M. Tkachuk, and C. Miehe, 2011. Journal of the
!        Mechanics and Physics of Solids.
!
! 
! This subroutine is for the following element types
!  > two-dimensional 4 node isoparametric element as shown below
!       with 1pt (reduced) or 4pt (full) gauss integration.
!  > three-dimensional 8 node isoparametric element as shown below
!       with 1pt (reduced) or 8pt (full) gauss integration.
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
!  Mechanical, traction- and pressure-type boundary conditions 
!   may be applied to the dummy mesh using the Abaqus built-in 
!   commands *Dload or *Dsload.
!
! Surface charge density boundary conditions are supported in
!  these elements.  Based on our convention, the face on which
!  the charge is applied is the "label", i.e.,
!  - U1,U2,U3,U4,... refer to charge applied to faces 1,2,3,4,...
!
! Surface tension effects are also supported in the 2D elements.
!  Based on our convention, the face on which the surface tension
!  acts is the ``label - 10'', i.e.,
!  - U11,U12,U13,U14,... refer to surface tension acting
!    on faces 1,2,3,4,... respectively,
!
!
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
!
!
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
!     Face numbering follows:
!       Face 1 = nodes 1,2,3,4
!       Face 2 = nodes 5,8,7,6
!       Face 3 = nodes 1,5,6,2
!       Face 4 = nodes 2,6,7,3
!       Face 5 = nodes 3,7,8,4
!       Face 6 = nodes 4,8,5,1
!
! David L. Henann, Shuolun Wang, and Shawn A. Chester, August 2014
!
***********************************************************************
!
! User element statement in the input file (set ? values as needed):
!
!  2D elements
!  *User Element,Nodes=4,Type=U?,Iproperties=3,Properties=?,Coordinates=2,Variables=?,Unsymm
!  1,2,11
!
!  3D elements
!  *User Element,Nodes=8,Type=U3,Iproperties=3,Properties=?,Coordinates=3,Variables=?,Unsymm
!  1,2,3,11
!
!
!     State Variables
!     --------------------------------------------------------------
!     Local SDV's (used for the solution procedure)
!       j = 0
!       do k = 1,nIntPt
!        do n = 1,nVisco
!         m = (n-1)*6
!         svars(1+m+j) = A(n,1,1) -- symmetric stretch-like internal
!         svars(2+m+j) = A(n,2,2)    variable for viscous mechanism n
!         svars(3+m+j) = A(n,3,3)    at integ pt k
!         svars(4+m+j) = A(n,2,3)
!         svars(5+m+j) = A(n,1,3)
!         svars(6+m+j) = A(n,1,2)
!        end loop over n, the viscous mechanisms
!        j = j + nlSdv
!       end loop over k, the integration points
!
!     In the input file, set 'nlSdv'= number of local SDV's
!
!     In the input file, set 'nInt'= number of volume integ pts
!
!     In the input file, set 'varibles'=(nlSdv*nIntPt)
!
!     In the input file, set 'IProperties'=3
!
!     In the input file, set 'Properties'=(4+2*nVisco)
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     Gshear  = props(1) ! Shear modulus of the equilibrium mechanism
!     lamL    = props(2) ! Locking stretch (set less than zero for NeoHookean)
!     Kbulk   = props(3) ! Bulk modulus of the equilibrium mechanism
!     permit  = props(4) ! Permittivity
!     Gneq(n) = props((n-1)*2+5) ! Shear modulus of viscous mechanism n
!     tau(n)  = props((n-1)*2+6) ! Relaxation time for viscous mechanism n
!
!     nVisco = jprops(1) ! Number of viscous mechanisms
!     nlSdv  = jprops(2) ! Number of local sdv's per integ pt
!     nInt   = jprops(3) ! Number of volume integ points
!
!***********************************************************************
!
! User subroutine (UVARM) for plotting the effective distortional stretch
!  on the deformed body as a contour.
!
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
      real*8 Lam1,Lam2,Lam3,detF

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

      ! Obtain the pricipal stretches (``DGP'') from the dummy mesh
      !
      CALL GETVRM('DGP',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     +     MATLAYO,LACCFLA)
      !
      ! These are the three principal stretches
      !
      Lam1 = array(1)
      Lam2 = array(2)
      Lam3 = array(3)
      !
      ! Compute detF based on the principal values
      !
      detF = Lam1*Lam2*Lam3
      !
      ! Compute the effective distortional stretch for plotting
      !
      uvar(1) = dsqrt((Lam1*Lam1 + Lam2*Lam2 + Lam3*Lam3)/3.d0)
     +     *(detF**(-1.d0/3.d0))

      
      RETURN
      END
!
!***********************************************************************
!
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      integer lenJobName,lenOutDir,nDim,nInt,nIntS,nVisco
      character*256 jobName,outDir,fileName

      parameter(nIntS=1) ! number of surface integration points, not used here


      !----------------------------------------------------------------
      ! 
      ! Perform initial checks
      !
      !
      ! Open the debug/error message file
      !
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     +     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')


      ! Check the procedure type, this should be a coupled
      !  temperature displacement or pore pressure displacement
      !  which are any of the following (64,65,72,73)
      !
      if((lflags(1).eq.64).or.(lflags(1).eq.65).or.
     +     (lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! all is good
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and chekc the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         write(80,*) 'Abaqus does not have the right procedure'
         write(80,*) 'go back and chekc the procedure type'
         write(80,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem, I think this only matters
      !  when it comes to output in viewer
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a small displacement analysis'
         write(80,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a linear purturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.0.0) return
      !
      ! Done with initial checks
      !
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! 
      ! Call the paricular element to perform the analysis
      !
      !
      ! Obtain the number of viscous mechanisms
      !
      nVisco = jprops(1)
      !
      ! Obtain the number of volume integration points
      !
      nInt = jprops(3)
      !
      !
      if(jtype.eq.1) then
         !
         ! This is a plane strain analysis
         !
         nDim = 2
         call UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS,nVisco)
         !
         !
      elseif(jtype.eq.2) then
         !
         ! This is an axisymmetric analysis
         !
         nDim = 2
         call UAX4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS,nVisco)
         !
         !
      elseif(jtype.eq.3) then
         !
         ! This is a 3D analysis
         !
         nDim = 3
         call U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS,nVisco)
         !
         !
      else
         !
         ! We have a problem...
         !
         write(*,*) 'Element type not supported, jtype=',jtype
         write(80,*) 'Element type not supported, jtype=',jtype
         call xit
         !
      endif
      !
      ! Done with this element, RHS and AMATRX already returned
      !  as output from the specific element routine called
      !
      !----------------------------------------------------------------


      return
      end subroutine uel

************************************************************************

      subroutine UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS,nVisco)

*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 u(nNode,2),du(nNode,ndofel),phi(nNode),coordsC(mcrd,nNode)
      real*8 uNew(nNode,ndofel),uOld(nNode,ndofel),u_t(nNode,ndofel)
      real*8 dPhi(nNode)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,faceFlag,nIntS,nVisco

      real*8 Iden(3,3),Le,theta0,phi0,Ru(2*nNode,1),Rphi(nNode,1),detF
      real*8 Kuu(2*nNode,2*nNode),Kphiphi(nNode,nNode),sh0(nNode)
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C,Vmol
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),ds,flux
      real*8 sh(nNode),detMapJ,phi_t,dsh(nNode,2),detMapJC,Limit,umeror
      real*8 dshC(nNode,2),mu_tau,mu_t,dMUdX(2,1),dMUdt,F_tau(3,3),DmDJ
      real*8 F_t(3,3),detF_tau,xi(nInt,2),detF_t,TR_tau(3,3),T_tau(3,3)
      real*8 SpUUMod(3,3,3,3),phi_tau,dPdt,DphiDmu,DphidotDmu,Mfluid
      real*8 Smat(3,1),Bmat(3,2*nNode),BodyForceRes(2*nNode,1),Qmat(4,4)
      real*8 Gmat(4,2*nNode),G0mat(4,2*nNode),Amat(4,4),wS(nIntS),deltaA
      real*8 xLocal(nIntS),yLocal(nIntS),Kphiu(nNode,2*nNode),D_tau(3,1)
      real*8 Kuphi(2*nNode,nNode),Nvec(1,nNode),ResFac,AmatUPhi(4,2)
      real*8 SpUPhiMod(3,3,3),SpPhiUMod(3,3,3),SpPhiPhiMod(3,3)
      real*8 A_t(nVisco,3,3),A_tau(nVisco,3,3),ER(3,1),detMapJ0
      real*8 bodyForce(2),bodyCharge,AmatPhiU(2,4),Amatphi(2,2)
      real*8 Dvec(2,1),DdsDu(1,2*nNode)



      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Get element parameters
      !
      nlSdv = jprops(2) ! number of local sdv's per integ point



      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rphi  = zero
      Kuu = zero
      Kphiphi = zero
      Kuphi = zero
      Kphiu = zero
      Energy = zero



      ! Obtain nodal displacements and chemical potentials
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         k = k + 1
         phi(i) = Uall(k)
         dPhi(i) = DUall(k,1)
      enddo


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      ! Impose any time-stepping changes on the increments of
      !  electric potential or displacement if you wish
      !
      ! electric potential increment
      !
      do i=1,nNode
         if(dabs(dPhi(i)).gt.1.d6) then
            pnewdt = 0.5
            return
         endif
      enddo
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     +     ((coordsC(2,1)-coordsC(2,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.0*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo


      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the
      !  `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures, 33, 3277-3296.
      !
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif



      ! Calculate the deformation gradient at the element centriod
      !  at the the begining and end of the increment for use in 
      !  the `F-bar' method. `Tau' represents the end of the increment
      !  and `t' the previous increment.
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      ! modify for plane-strain
      !
      Fc_tau(3,3) = one
      Fc_t(3,3) = one
      !
      ! 2D plane-strain implementation detF
      !
      detFc_t = Fc_t(1,1)*Fc_t(2,2) - Fc_t(1,2)*Fc_t(2,1)
      detFc_tau = Fc_tau(1,1)*Fc_tau(2,2) - Fc_tau(1,2)*Fc_tau(2,1)
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over body integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step.
            !  Give initial conditions.
            !
            do n=1,nVisco
               A_t(n,:,:) = Iden
            enddo
            !
         else
            !
            ! this is not the first increment, read old values
            !
            do n=1,nVisco
               k = (n-1)*6
               A_t(n,1,1) = svars(1+k+jj)
               A_t(n,2,2) = svars(2+k+jj)
               A_t(n,3,3) = svars(3+k+jj)
               A_t(n,2,3) = svars(4+k+jj)
               A_t(n,3,2) = A_t(n,2,3)
               A_t(n,1,3) = svars(5+k+jj)
               A_t(n,3,1) = A_t(n,1,3)
               A_t(n,1,2) = svars(6+k+jj)
               A_t(n,2,1) = A_t(n,1,2)
            enddo
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Obtain the referential electric field at this integration point
         !
         ER = zero
         do k=1,nNode
            do i=1,nDim
               ER(i,1) = ER(i,1) + phi(k)*dsh(k,i)
            enddo
         enddo



         ! Obtain, and modify the deformation gradient at this integration
         !  point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.  Also, take care of plane-strain or axisymetric
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! modify F(3,3) for plane-strain 
         !
         F_tau(3,3) = one
         F_t(3,3) = one
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 4 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            !  2D plane-strain implementation
            !
            detF_t = F_t(1,1)*F_t(2,2) - F_t(1,2)*F_t(2,1)
            detF_tau = F_tau(1,1)*F_tau(2,2) - F_tau(1,2)*F_tau(2,1)
            do i=1,nDim
               do j=1,nDim
                  F_tau(i,j) =((detFc_tau/detF_tau)**half)*F_tau(i,j)
                  F_t(i,j) = ((detFc_t/detF_t)**half)*F_t(i,j)
               enddo
            enddo
         endif
         call mdet(F_tau,detF)


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point 
         !
         call integ(props,nprops,dtime,nVisco,
     +        F_tau,ER,A_t,
     +        T_tau,D_tau,A_tau,deltaA,
     +        SpUUMod,SpUPhiMod,SpPhiUMod,SpPhiPhiMod,
     +        stat)
         !
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         do n=1,nVisco
            k = (n-1)*6
            svars(1+k+jj) = A_tau(n,1,1)
            svars(2+k+jj) = A_tau(n,2,2)
            svars(3+k+jj) = A_tau(n,3,3)
            svars(4+k+jj) = A_tau(n,2,3)
            svars(5+k+jj) = A_tau(n,1,3)
            svars(6+k+jj) = A_tau(n,1,2)
         enddo
         jj = jj + nlSdv ! setup for the next intPt



         ! Time stepping algorithim based on the constitutive response.
         !  Here based on the change of the viscous stretch increment.
         !
         Limit = 0.1d0
         umeror = deltaA/Limit
         !
         if(umeror.le.half) then
            pnewdt = 1.5d0
         elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
            pnewdt = 1.25d0
         elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
            pnewdt = 0.75d0
         else
            pnewdt = half
         endif


         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(1,2)
         !
         Bmat = zero
         do kk=1,nNode
            Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,2+nDim*(kk-1)) = dshC(kk,1)
         enddo
         !
         bodyForce = zero ! the body force vector may be specified here
         !
         BodyForceRes = zero
         do kk=1,nNode
            BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*bodyForce(1)
            BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*bodyForce(2)
         enddo
         !
         Ru = Ru + detmapJC*w(intpt)*
     +        (
     +        -matmul(transpose(Bmat),Smat)
     +        + BodyForceRes
     +        )


         ! Compute/update the electric potential residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         !
         bodyCharge = zero ! the free charge density may be specified here
         !
         Dvec(1,1) = D_tau(1,1)
         Dvec(2,1) = D_tau(2,1)
         !
         Rphi = Rphi + detmapJC*w(intpt)*
     +        (
     +        matmul(dshC,Dvec)
     +        + transpose(Nvec)*bodyCharge
     +        )



         ! Compute/update the displacement tangent matrix
         !
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
            Gmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Gmat(4,2+nDim*(kk-1)) = dshC(kk,2)
         enddo

         G0mat = zero
         do kk=1,nNode
            G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(3,1+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(4,2+nDim*(kk-1)) = dshC0(kk,2)
         enddo

         Amat = zero
         Amat(1,1) = SpUUMod(1,1,1,1)
         Amat(1,2) = SpUUMod(1,1,2,1)
         Amat(1,3) = SpUUMod(1,1,1,2)
         Amat(1,4) = SpUUMod(1,1,2,2)
         Amat(2,1) = SpUUMod(2,1,1,1)
         Amat(2,2) = SpUUMod(2,1,2,1)
         Amat(2,3) = SpUUMod(2,1,1,2)
         Amat(2,4) = SpUUMod(2,1,2,2)
         Amat(3,1) = SpUUMod(1,2,1,1)
         Amat(3,2) = SpUUMod(1,2,2,1)
         Amat(3,3) = SpUUMod(1,2,1,2)
         Amat(3,4) = SpUUMod(1,2,2,2)
         Amat(4,1) = SpUUMod(2,2,1,1)
         Amat(4,2) = SpUUMod(2,2,2,1)
         Amat(4,3) = SpUUMod(2,2,1,2)
         Amat(4,4) = SpUUMod(2,2,2,2)


         Qmat = zero
         Qmat(1,1) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
         Qmat(2,1) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
         Qmat(3,1) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
         Qmat(4,1) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
         Qmat(1,4) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
         Qmat(2,4) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
         Qmat(3,4) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
         Qmat(4,4) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            

         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            ! This is the tangent using the F-bar method with the
            !  4 node fully integrated linear element
            !
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +            matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           + matmul(transpose(Gmat),matmul(Qmat,(G0mat-Gmat)))
     +           )
         else
            !
            ! This is the tangent not using the F-bar method with all
            !  other elements
            !
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           )
         endif



         ! Compute/update the electric potential tangent matrix
         !
         Amatphi(1,1) = SpPhiPhiMod(1,1)
         Amatphi(1,2) = SpPhiPhiMod(1,2)
         Amatphi(2,1) = SpPhiPhiMod(2,1)
         Amatphi(2,2) = SpPhiPhiMod(2,2)
         !
         Kphiphi = Kphiphi - detmapJC*w(intPt)*
     +        (
     +        matmul(matmul(dshC,Amatphi),transpose(dshC))
     +        )


         

         ! Compute/update the electric potential - displacement tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         !
         AmatPhiU = zero
         AmatPhiU(1,1) = SpPhiUMod(1,1,1)
         AmatPhiU(1,2) = SpPhiUMod(1,2,1)
         AmatPhiU(1,3) = SpPhiUMod(1,1,2)
         AmatPhiU(1,4) = SpPhiUMod(1,2,2)
         AmatPhiU(2,1) = SpPhiUMod(2,1,1)
         AmatPhiU(2,2) = SpPhiUMod(2,2,1)
         AmatPhiU(2,3) = SpPhiUMod(2,1,2)
         AmatPhiU(2,4) = SpPhiUMod(2,2,2)
         !
         Kphiu = Kphiu - detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(dshC,AmatPhiU),Gmat)
     +        )


         ! Compute/update the displacement - chemical potential tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         !
         AmatUPhi = zero
         AmatUPhi(1,1) = SpUPhiMod(1,1,1)
         AmatUPhi(2,1) = SpUPhiMod(2,1,1)
         AmatUPhi(3,1) = SpUPhiMod(1,2,1)
         AmatUPhi(4,1) = SpUPhiMod(2,2,1)
         AmatUPhi(1,2) = SpUPhiMod(1,1,2)
         AmatUPhi(2,2) = SpUPhiMod(2,1,2)
         AmatUPhi(3,2) = SpUPhiMod(1,2,2)
         AmatUPhi(4,2) = SpUPhiMod(2,2,2)
         !
         Kuphi = Kuphi + detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(transpose(Gmat),AmatUPhi),transpose(dshC))
     +        )


      enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface terms here
      !
      !
      if(ndload.gt.0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i=1,ndload
            !
            ! based on my convention the face which the boundary condition
            !  acts on is the ``label''
            !
            face = jdltyp(i,1) ! label
            flux = adlmag(i,1) ! BC magnitude
            !
            ! For surface charge, flux = omega
            ! For surface tension, flux = gamma
            !
            if((face.ge.1).and.(face.le.4)) then
               !
               ! surface charge applied
               !
               if(face.eq.1) then
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,1)-coordsC(1,2))**two) + 
     +                 ((coordsC(2,1)-coordsC(2,2))**two))
                  Rphi(1,1) = Rphi(1,1) + half*flux*Le
                  Rphi(2,1) = Rphi(2,1) + half*flux*Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kphiu(1,1) = Kphiu(1,1) - 
     +                 (coordsC(1,1)-coordsC(1,2))*flux/(two*Le)
                  Kphiu(1,2) = Kphiu(1,2) -
     +                 (coordsC(2,1)-coordsC(2,2))*flux/(two*Le)
                  Kphiu(1,3) = Kphiu(1,3) -
     +                 (coordsC(1,2)-coordsC(1,1))*flux/(two*Le)
                  Kphiu(1,4) = Kphiu(1,4) -
     +                 (coordsC(2,2)-coordsC(2,1))*flux/(two*Le)
                  Kphiu(2,1) = Kphiu(2,1) -
     +                 (coordsC(1,1)-coordsC(1,2))*flux/(two*Le)
                  Kphiu(2,2) = Kphiu(2,2) -
     +                 (coordsC(2,1)-coordsC(2,2))*flux/(two*Le)
                  Kphiu(2,3) = Kphiu(2,3) -
     +                 (coordsC(1,2)-coordsC(1,1))*flux/(two*Le)
                  Kphiu(2,4) = Kphiu(2,4) -
     +                 (coordsC(2,2)-coordsC(2,1))*flux/(two*Le)
                  !
               elseif(face.eq.2) then
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,2)-coordsC(1,3))**two) + 
     +                 ((coordsC(2,2)-coordsC(2,3))**two))
                  Rphi(2,1) = Rphi(2,1) + half*flux*Le
                  Rphi(3,1) = Rphi(3,1) + half*flux*Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kphiu(2,3) = Kphiu(2,3) - 
     +                 (coordsC(1,2)-coordsC(1,3))*flux/(two*Le)
                  Kphiu(2,4) = Kphiu(2,4) -
     +                 (coordsC(2,2)-coordsC(2,3))*flux/(two*Le)
                  Kphiu(2,5) = Kphiu(2,5) -
     +                 (coordsC(1,3)-coordsC(1,2))*flux/(two*Le)
                  Kphiu(2,6) = Kphiu(2,6) - 
     +                 (coordsC(2,3)-coordsC(2,2))*flux/(two*Le)
                  Kphiu(3,3) = Kphiu(3,3) -
     +                 (coordsC(1,2)-coordsC(1,3))*flux/(two*Le)
                  Kphiu(3,4) = Kphiu(3,4) - 
     +                 (coordsC(2,2)-coordsC(2,3))*flux/(two*Le)
                  Kphiu(3,5) = Kphiu(3,5) -
     +                 (coordsC(1,3)-coordsC(1,2))*flux/(two*Le)
                  Kphiu(3,6) = Kphiu(3,6) -
     +                 (coordsC(2,3)-coordsC(2,2))*flux/(two*Le)
                  !
               elseif(face.eq.3) then
                  !
                  ! Modify the displacement residual, loop over nodes
                  !
                  Le = dsqrt(((coordsC(1,3)-coordsC(1,4))**two) + 
     +                 ((coordsC(2,3)-coordsC(2,4))**two))
                  Rphi(3,1) = Rphi(3,1) + half*flux*Le
                  Rphi(4,1) = Rphi(4,1) + half*flux*Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kphiu(3,5) = Kphiu(3,5) -
     +                 (coordsC(1,3)-coordsC(1,4))*flux/(two*Le)
                  Kphiu(3,6) = Kphiu(3,6) - 
     +                 (coordsC(2,3)-coordsC(2,4))*flux/(two*Le)
                  Kphiu(3,7) = Kphiu(3,7) -
     +                 (coordsC(1,4)-coordsC(1,3))*flux/(two*Le)
                  Kphiu(3,8) = Kphiu(3,8) -
     +                 (coordsC(2,4)-coordsC(2,3))*flux/(two*Le)
                  Kphiu(4,5) = Kphiu(4,5) -
     +                 (coordsC(1,3)-coordsC(1,4))*flux/(two*Le)
                  Kphiu(4,6) = Kphiu(4,6) -
     +                 (coordsC(2,3)-coordsC(2,4))*flux/(two*Le)
                  Kphiu(4,7) = Kphiu(4,7) -
     +                 (coordsC(1,4)-coordsC(1,3))*flux/(two*Le)
                  Kphiu(4,8) = Kphiu(4,8) - 
     +                 (coordsC(2,4)-coordsC(2,3))*flux/(two*Le)
                  !
               elseif(face.eq.4) then
                  !
                  ! Modify the displacement residual, loop over nodes
                  !
                  Le = dsqrt(((coordsC(1,4)-coordsC(1,1))**two) + 
     +                 ((coordsC(2,4)-coordsC(2,1))**two))
                  Rphi(4,1) = Rphi(4,1) + half*flux*Le
                  Rphi(1,1) = Rphi(1,1) + half*flux*Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kphiu(4,7) = Kphiu(4,7) - 
     +                 (coordsC(1,4)-coordsC(1,1))*flux/(two*Le)
                  Kphiu(4,8) = Kphiu(4,8) - 
     +                 (coordsC(2,4)-coordsC(2,1))*flux/(two*Le)
                  Kphiu(4,1) = Kphiu(4,1) - 
     +                 (coordsC(1,1)-coordsC(1,4))*flux/(two*Le)
                  Kphiu(4,2) = Kphiu(4,2) - 
     +                 (coordsC(2,1)-coordsC(2,4))*flux/(two*Le)
                  Kphiu(1,7) = Kphiu(1,7) - 
     +                 (coordsC(1,4)-coordsC(1,1))*flux/(two*Le)
                  Kphiu(1,8) = Kphiu(1,8) - 
     +                 (coordsC(2,4)-coordsC(2,1))*flux/(two*Le)
                  Kphiu(1,1) = Kphiu(1,1) - 
     +                 (coordsC(1,1)-coordsC(1,4))*flux/(two*Le)
                  Kphiu(1,2) = Kphiu(1,2) - 
     +                 (coordsC(2,1)-coordsC(2,4))*flux/(two*Le)
                  !
               endif
               !
            elseif((face.ge.11).and.(face.le.14)) then
               !
               ! surface tension
               !
               if(face.eq.11) then
                  !
                  ! Surface tension on face 1
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,1)-coordsC(1,2))**two) + 
     +                 ((coordsC(2,1)-coordsC(2,2))**two))
                  Ru(1,1) = Ru(1,1)-flux*(coordsC(1,1)-coordsC(1,2))/Le
                  Ru(2,1) = Ru(2,1)-flux*(coordsC(2,1)-coordsC(2,2))/Le
                  Ru(3,1) = Ru(3,1)-flux*(coordsC(1,2)-coordsC(1,1))/Le
                  Ru(4,1) = Ru(4,1)-flux*(coordsC(2,2)-coordsC(2,1))/Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kuu(1,1) = Kuu(1,1) + flux/Le - 
     +                 (coordsC(1,1)-coordsC(1,2))*
     +                 (coordsC(1,1)-coordsC(1,2))*flux/(Le**three)
                  Kuu(1,2) = Kuu(1,2) - 
     +                 (coordsC(1,1)-coordsC(1,2))*
     +                 (coordsC(2,1)-coordsC(2,2))*flux/(Le**three)
                  Kuu(1,3) = Kuu(1,3) - flux/Le - 
     +                 (coordsC(1,1)-coordsC(1,2))*
     +                 (coordsC(1,2)-coordsC(1,1))*flux/(Le**three)
                  Kuu(1,4) = Kuu(1,4) - 
     +                 (coordsC(1,1)-coordsC(1,2))*
     +                 (coordsC(2,2)-coordsC(2,1))*flux/(Le**three)
                  Kuu(2,1) = Kuu(2,1) - 
     +                 (coordsC(2,1)-coordsC(2,2))*
     +                 (coordsC(1,1)-coordsC(1,2))*flux/(Le**three)
                  Kuu(2,2) = Kuu(2,2) + flux/Le -
     +                 (coordsC(2,1)-coordsC(2,2))*
     +                 (coordsC(2,1)-coordsC(2,2))*flux/(Le**three)
                  Kuu(2,3) = Kuu(2,3) -
     +                 (coordsC(2,1)-coordsC(2,2))*
     +                 (coordsC(1,2)-coordsC(1,1))*flux/(Le**three)
                  Kuu(2,4) = Kuu(2,4) - flux/Le -
     +                 (coordsC(2,1)-coordsC(2,2))*
     +                 (coordsC(2,2)-coordsC(2,1))*flux/(Le**three)
                  Kuu(3,1) = Kuu(3,1) - flux/Le - 
     +                 (coordsC(1,2)-coordsC(1,1))*
     +                 (coordsC(1,1)-coordsC(1,2))*flux/(Le**three)
                  Kuu(3,2) = Kuu(3,2) - 
     +                 (coordsC(1,2)-coordsC(1,1))*
     +                 (coordsC(2,1)-coordsC(2,2))*flux/(Le**three)
                  Kuu(3,3) = Kuu(3,3) + flux/Le - 
     +                 (coordsC(1,2)-coordsC(1,1))*
     +                 (coordsC(1,2)-coordsC(1,1))*flux/(Le**three)
                  Kuu(3,4) = Kuu(3,4) - 
     +                 (coordsC(1,2)-coordsC(1,1))*
     +                 (coordsC(2,2)-coordsC(2,1))*flux/(Le**three)
                  Kuu(4,1) = Kuu(4,1) -
     +                 (coordsC(2,2)-coordsC(2,1))*
     +                 (coordsC(1,1)-coordsC(1,2))*flux/(Le**three)
                  Kuu(4,2) = Kuu(4,2) - flux/Le -
     +                 (coordsC(2,2)-coordsC(2,1))*
     +                 (coordsC(2,1)-coordsC(2,2))*flux/(Le**three)
                  Kuu(4,3) = Kuu(4,3) -
     +                 (coordsC(2,2)-coordsC(2,1))*
     +                 (coordsC(1,2)-coordsC(1,1))*flux/(Le**three)
                  Kuu(4,4) = Kuu(4,4) + flux/Le -
     +                 (coordsC(2,2)-coordsC(2,1))*
     +                 (coordsC(2,2)-coordsC(2,1))*flux/(Le**three)
                  !
               elseif(face.eq.12) then
                  !
                  ! Surface tension on face 2
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,2)-coordsC(1,3))**two) + 
     +                 ((coordsC(2,2)-coordsC(2,3))**two))
                  Ru(3,1) = Ru(3,1)-flux*(coordsC(1,2)-coordsC(1,3))/Le
                  Ru(4,1) = Ru(4,1)-flux*(coordsC(2,2)-coordsC(2,3))/Le
                  Ru(5,1) = Ru(5,1)-flux*(coordsC(1,3)-coordsC(1,2))/Le
                  Ru(6,1) = Ru(6,1)-flux*(coordsC(2,3)-coordsC(2,2))/Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kuu(3,3) = Kuu(3,3) + flux/Le - 
     +                 (coordsC(1,2)-coordsC(1,3))*
     +                 (coordsC(1,2)-coordsC(1,3))*flux/(Le**three)
                  Kuu(3,4) = Kuu(3,4) - 
     +                 (coordsC(1,2)-coordsC(1,3))*
     +                 (coordsC(2,2)-coordsC(2,3))*flux/(Le**three)
                  Kuu(3,5) = Kuu(3,5) - flux/Le - 
     +                 (coordsC(1,2)-coordsC(1,3))*
     +                 (coordsC(1,3)-coordsC(1,2))*flux/(Le**three)
                  Kuu(3,6) = Kuu(3,6) - 
     +                 (coordsC(1,2)-coordsC(1,3))*
     +                 (coordsC(2,3)-coordsC(2,2))*flux/(Le**three)
                  Kuu(4,3) = Kuu(4,3) -
     +                 (coordsC(2,2)-coordsC(2,3))*
     +                 (coordsC(1,2)-coordsC(1,3))*flux/(Le**three)
                  Kuu(4,4) = Kuu(4,4) + flux/Le -
     +                 (coordsC(2,2)-coordsC(2,3))*
     +                 (coordsC(2,2)-coordsC(2,3))*flux/(Le**three)
                  Kuu(4,5) = Kuu(4,5) -
     +                 (coordsC(2,2)-coordsC(2,3))*
     +                 (coordsC(1,3)-coordsC(1,2))*flux/(Le**three)
                  Kuu(4,6) = Kuu(4,6) - flux/Le -
     +                 (coordsC(2,2)-coordsC(2,3))*
     +                 (coordsC(2,3)-coordsC(2,2))*flux/(Le**three)
                  Kuu(5,3) = Kuu(5,3) - flux/Le - 
     +                 (coordsC(1,3)-coordsC(1,2))*
     +                 (coordsC(1,2)-coordsC(1,3))*flux/(Le**three)
                  Kuu(5,4) = Kuu(5,4) - 
     +                 (coordsC(1,3)-coordsC(1,2))*
     +                 (coordsC(2,2)-coordsC(2,3))*flux/(Le**three)
                  Kuu(5,5) = Kuu(5,5) + flux/Le - 
     +                 (coordsC(1,3)-coordsC(1,2))*
     +                 (coordsC(1,3)-coordsC(1,2))*flux/(Le**three)
                  Kuu(5,6) = Kuu(5,6) - 
     +                 (coordsC(1,3)-coordsC(1,2))*
     +                 (coordsC(2,3)-coordsC(2,2))*flux/(Le**three)
                  Kuu(6,3) = Kuu(6,3) -
     +                 (coordsC(2,3)-coordsC(2,2))*
     +                 (coordsC(1,2)-coordsC(1,3))*flux/(Le**three)
                  Kuu(6,4) = Kuu(6,4) - flux/Le -
     +                 (coordsC(2,3)-coordsC(2,2))*
     +                 (coordsC(2,2)-coordsC(2,3))*flux/(Le**three)
                  Kuu(6,5) = Kuu(6,5) -
     +                 (coordsC(2,3)-coordsC(2,2))*
     +                 (coordsC(1,3)-coordsC(1,2))*flux/(Le**three)
                  Kuu(6,6) = Kuu(6,6) + flux/Le -
     +                 (coordsC(2,3)-coordsC(2,2))*
     +                 (coordsC(2,3)-coordsC(2,2))*flux/(Le**three)
                  !
               elseif(face.eq.13) then
                  !
                  ! Surface tension on face 3
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,3)-coordsC(1,4))**two) + 
     +                 ((coordsC(2,3)-coordsC(2,4))**two))
                  Ru(5,1) = Ru(5,1)-flux*(coordsC(1,3)-coordsC(1,4))/Le
                  Ru(6,1) = Ru(6,1)-flux*(coordsC(2,3)-coordsC(2,4))/Le
                  Ru(7,1) = Ru(7,1)-flux*(coordsC(1,4)-coordsC(1,3))/Le
                  Ru(8,1) = Ru(8,1)-flux*(coordsC(2,4)-coordsC(2,3))/Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kuu(5,5) = Kuu(5,5) + flux/Le - 
     +                 (coordsC(1,3)-coordsC(1,4))*
     +                 (coordsC(1,3)-coordsC(1,4))*flux/(Le**three)
                  Kuu(5,6) = Kuu(5,6) - 
     +                 (coordsC(1,3)-coordsC(1,4))*
     +                 (coordsC(2,3)-coordsC(2,4))*flux/(Le**three)
                  Kuu(5,7) = Kuu(5,7) - flux/Le - 
     +                 (coordsC(1,3)-coordsC(1,4))*
     +                 (coordsC(1,4)-coordsC(1,3))*flux/(Le**three)
                  Kuu(5,8) = Kuu(5,8) - 
     +                 (coordsC(1,3)-coordsC(1,4))*
     +                 (coordsC(2,4)-coordsC(2,3))*flux/(Le**three)
                  Kuu(6,5) = Kuu(6,5) -
     +                 (coordsC(2,3)-coordsC(2,4))*
     +                 (coordsC(1,3)-coordsC(1,4))*flux/(Le**three)
                  Kuu(6,6) = Kuu(6,6) + flux/Le -
     +                 (coordsC(2,3)-coordsC(2,4))*
     +                 (coordsC(2,3)-coordsC(2,4))*flux/(Le**three)
                  Kuu(6,7) = Kuu(6,7) -
     +                 (coordsC(2,3)-coordsC(2,4))*
     +                 (coordsC(1,4)-coordsC(1,3))*flux/(Le**three)
                  Kuu(6,8) = Kuu(6,8) - flux/Le -
     +                 (coordsC(2,3)-coordsC(2,4))*
     +                 (coordsC(2,4)-coordsC(2,3))*flux/(Le**three)
                  Kuu(7,5) = Kuu(7,5) - flux/Le - 
     +                 (coordsC(1,4)-coordsC(1,3))*
     +                 (coordsC(1,3)-coordsC(1,4))*flux/(Le**three)
                  Kuu(7,6) = Kuu(7,6) - 
     +                 (coordsC(1,4)-coordsC(1,3))*
     +                 (coordsC(2,3)-coordsC(2,4))*flux/(Le**three)
                  Kuu(7,7) = Kuu(7,7) + flux/Le - 
     +                 (coordsC(1,4)-coordsC(1,3))*
     +                 (coordsC(1,4)-coordsC(1,3))*flux/(Le**three)
                  Kuu(7,8) = Kuu(7,8) - 
     +                 (coordsC(1,4)-coordsC(1,3))*
     +                 (coordsC(2,4)-coordsC(2,3))*flux/(Le**three)
                  Kuu(8,5) = Kuu(8,5) -
     +                 (coordsC(2,4)-coordsC(2,3))*
     +                 (coordsC(1,3)-coordsC(1,4))*flux/(Le**three)
                  Kuu(8,6) = Kuu(8,6) - flux/Le -
     +                 (coordsC(2,4)-coordsC(2,3))*
     +                 (coordsC(2,3)-coordsC(2,4))*flux/(Le**three)
                  Kuu(8,7) = Kuu(8,7) -
     +                 (coordsC(2,4)-coordsC(2,3))*
     +                 (coordsC(1,4)-coordsC(1,3))*flux/(Le**three)
                  Kuu(8,8) = Kuu(8,8) + flux/Le -
     +                 (coordsC(2,4)-coordsC(2,3))*
     +                 (coordsC(2,4)-coordsC(2,3))*flux/(Le**three)
                  !
               elseif(face.eq.14) then
                  !
                  ! Surface tension on face 4
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,4)-coordsC(1,1))**two) + 
     +                 ((coordsC(2,4)-coordsC(2,1))**two))
                  Ru(7,1) = Ru(7,1)-flux*(coordsC(1,4)-coordsC(1,1))/Le
                  Ru(8,1) = Ru(8,1)-flux*(coordsC(2,4)-coordsC(2,1))/Le
                  Ru(1,1) = Ru(1,1)-flux*(coordsC(1,1)-coordsC(1,4))/Le
                  Ru(2,1) = Ru(2,1)-flux*(coordsC(2,1)-coordsC(2,4))/Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kuu(7,7) = Kuu(7,7) + flux/Le - 
     +                 (coordsC(1,4)-coordsC(1,1))*
     +                 (coordsC(1,4)-coordsC(1,1))*flux/(Le**three)
                  Kuu(7,8) = Kuu(7,8) - 
     +                 (coordsC(1,4)-coordsC(1,1))*
     +                 (coordsC(2,4)-coordsC(2,1))*flux/(Le**three)
                  Kuu(7,1) = Kuu(7,1) - flux/Le - 
     +                 (coordsC(1,4)-coordsC(1,1))*
     +                 (coordsC(1,1)-coordsC(1,4))*flux/(Le**three)
                  Kuu(7,2) = Kuu(7,2) - 
     +                 (coordsC(1,4)-coordsC(1,1))*
     +                 (coordsC(2,1)-coordsC(2,4))*flux/(Le**three)
                  Kuu(8,7) = Kuu(8,7) -
     +                 (coordsC(2,4)-coordsC(2,1))*
     +                 (coordsC(1,4)-coordsC(1,1))*flux/(Le**three)
                  Kuu(8,8) = Kuu(8,8) + flux/Le -
     +                 (coordsC(2,4)-coordsC(2,1))*
     +                 (coordsC(2,4)-coordsC(2,1))*flux/(Le**three)
                  Kuu(8,1) = Kuu(8,1) -
     +                 (coordsC(2,4)-coordsC(2,1))*
     +                 (coordsC(1,1)-coordsC(1,4))*flux/(Le**three)
                  Kuu(8,2) = Kuu(8,2) - flux/Le -
     +                 (coordsC(2,4)-coordsC(2,1))*
     +                 (coordsC(2,1)-coordsC(2,4))*flux/(Le**three)
                  Kuu(1,7) = Kuu(1,7) - flux/Le - 
     +                 (coordsC(1,1)-coordsC(1,4))*
     +                 (coordsC(1,4)-coordsC(1,1))*flux/(Le**three)
                  Kuu(1,8) = Kuu(1,8) - 
     +                 (coordsC(1,1)-coordsC(1,4))*
     +                 (coordsC(2,4)-coordsC(2,1))*flux/(Le**three)
                  Kuu(1,1) = Kuu(1,1) + flux/Le - 
     +                 (coordsC(1,1)-coordsC(1,4))*
     +                 (coordsC(1,1)-coordsC(1,4))*flux/(Le**three)
                  Kuu(1,2) = Kuu(1,2) - 
     +                 (coordsC(1,1)-coordsC(1,4))*
     +                 (coordsC(2,1)-coordsC(2,4))*flux/(Le**three)
                  Kuu(2,7) = Kuu(2,7) -
     +               (coordsC(2,1)-coordsC(2,4))*
     +                 (coordsC(1,4)-coordsC(1,1))*flux/(Le**three)
                  Kuu(2,8) = Kuu(2,8) - flux/Le -
     +                 (coordsC(2,1)-coordsC(2,4))*
     +                 (coordsC(2,4)-coordsC(2,1))*flux/(Le**three)
                  Kuu(2,1) = Kuu(2,1) -
     +                 (coordsC(2,1)-coordsC(2,4))*
     +                 (coordsC(1,1)-coordsC(1,4))*flux/(Le**three)
                  Kuu(2,2) = Kuu(2,2) + flux/Le -
     +                 (coordsC(2,1)-coordsC(2,4))*
     +                 (coordsC(2,1)-coordsC(2,4))*flux/(Le**three)
                  !
               endif
               !
            else
               write(*,*) 'Unknown face=',face
               write(80,*) 'Unknown face=',face
               call xit
            endif

         enddo ! loop over ndload
      endif ! ndload.gt.0 or not
      !
      ! End loop over surface terms
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Rphi,Kuu,Kuphi,Kphiu,Kphiphi,
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine UPE4

************************************************************************

      subroutine UAX4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS,nVisco)
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 u(nNode,2),du(nNode,ndofel),phi(nNode),coordsC(mcrd,nNode)
      real*8 uNew(nNode,ndofel),uOld(nNode,ndofel),u_t(nNode,ndofel)
      real*8 dPhi(nNode)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,nVisco,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,faceFlag,nIntS

      real*8 Iden(3,3),Le,theta0,phi0,Ru(2*nNode,1),Rphi(nNode,1),detF
      real*8 Kuu(2*nNode,2*nNode),Kphiphi(nNode,nNode),sh0(nNode),AR_t
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C,AR0
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),ds,flux,ARc
      real*8 sh(nNode),detMapJ,phi_t,dsh(nNode,2),detMapJC,Limit,umeror
      real*8 dshC(nNode,2),mu_tau,mu_t,dMUdX(2,1),dMUdt,F_tau(3,3),DmDJ
      real*8 F_t(3,3),detF_tau,xi(nInt,2),detF_t,TR_tau(3,3),T_tau(3,3)
      real*8 SpUUMod(3,3,3,3),phi_tau,dPdt,DphiDmu,DphidotDmu,AR,deltaA
      real*8 SmatAx(4,1),BmatAx(4,2*nNode),BodyForceResAx(2*nNode,1)
      real*8 GmatAx(5,2*nNode),G0matAx(5,2*nNode),AmatAx(5,5),wS(nIntS)
      real*8 xLocal(nIntS),yLocal(nIntS),Kphiu(nNode,2*nNode),D_tau(3,1)
      real*8 Kuphi(2*nNode,nNode),Nvec(1,nNode),ResFac,AmatUPhi(5,2)
      real*8 SpUPhiMod(3,3,3),SpPhiUMod(3,3,3),SpPhiPhiMod(3,3)
      real*8 A_t(nVisco,3,3),A_tau(nVisco,3,3),ER(3,1),detMapJ0
      real*8 bodyForce(2),bodyCharge,AmatPhiU(2,5),Amatphi(2,2)
      real*8 Dvec(2,1),QmatAx(5,5),DdsDu(1,2*nNode),dAR(1,2*nNode)


      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Get element parameters
      !
      nlSdv = jprops(2) ! number of local sdv's per integ point



      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rphi  = zero
      Kuu = zero
      Kphiphi = zero
      Kuphi = zero
      Kphiu = zero
      Energy = zero



      ! Obtain nodal displacements and chemical potentials
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         k = k + 1
         phi(i) = Uall(k)
         dPhi(i) = DUall(k,1)
      enddo



      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo



      ! Impose any time-stepping changes on the increments of
      !  electric potential or displacement if you wish
      !
      ! electric potential increment
      !
      do i=1,nNode
         if(dabs(dPhi(i)).gt.1.d6) then
            pnewdt = 0.5
            return
         endif
      enddo
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     +     ((coordsC(2,1)-coordsC(2,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.0*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo


      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the
      !  `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures, 33, 3277-3296.
      !
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! For an axisymmetric problem, find the ``r'' that
      !  shows up in the integrals for axisymmetric
      !  and in the F(3,3), the factors of 2Pi are for integrals
      !  i.e., dV = 2 pi r dr dz
      !
      AR0  = zero
      ARc  = zero
      AR_t = zero
      do i=1,nNode
         ! radial coord in ref config at centroid
         AR0  = AR0 + sh0(i)*coords(1,i)
         ! radial coord in current config at centroid
         ARc  = ARc + sh0(i)*(coords(1,i) + u(i,1))
         ! radial coord in current config at centroid in previous step
         AR_t = AR_t + sh0(i)*(coords(1,i) + uOld(i,1))
      enddo



      ! Calculate the deformation gradient at the element centriod
      !  at the the begining and end of the increment for use in 
      !  the `F-bar' method. `Tau' represents the end of the increment
      !  and `t' the previous increment.
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      ! modify for axisymmetric
      !
      Fc_tau(3,3) = ARc/AR0
      Fc_t(3,3) = AR_t/AR0
      !
      ! axisymmetric implementation detF
      !
      call mdet(Fc_tau,detFc_tau)
      call mdet(Fc_t,detFc_t)
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over body integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif



      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step.
            !  Give initial conditions.
            !
            do n=1,nVisco
               A_t(n,:,:) = Iden
            enddo
            !
         else
            !
            ! this is not the first increment, read old values
            !
            do n=1,nVisco
               k = (n-1)*6
               A_t(n,1,1) = svars(1+k+jj)
               A_t(n,2,2) = svars(2+k+jj)
               A_t(n,3,3) = svars(3+k+jj)
               A_t(n,2,3) = svars(4+k+jj)
               A_t(n,3,2) = A_t(n,2,3)
               A_t(n,1,3) = svars(5+k+jj)
               A_t(n,3,1) = A_t(n,1,3)
               A_t(n,1,2) = svars(6+k+jj)
               A_t(n,2,1) = A_t(n,1,2)
            enddo
            !
         endif



         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! For an axisymmetric problem, find the ``r'' that
         !  shows up in the integrals for axisymmetric
         !  and in the F(3,3), the factors of 2Pi are for integrals
         !  i.e., dV = 2 pi r dr dz
         !
         !
         AR0  = zero
         AR   = zero
         AR_t = zero
         do i=1,nNode
            AR0 = AR0 + sh(i)*coords(1,i)
            AR  = AR  + sh(i)*(coords(1,i) + u(i,1))
            AR_t = AR_t + sh(i)*(coords(1,i) + uOld(i,1))
         enddo
         AR0  = two*Pi*AR0
         AR   = two*Pi*AR
         AR_t = two*Pi*AR_t


         ! Obtain the referential electric field at this integration point
         !
         ER = zero
         do k=1,nNode
            do i=1,nDim
               ER(i,1) = ER(i,1) + phi(k)*dsh(k,i)
            enddo
         enddo


         ! Obtain, and modify the deformation gradient at this integration
         !  point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.  Also, take care of plane-strain or axisymetric
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! modify F(3,3) for axisymetric, give R/R0
         !
         F_tau(3,3) = AR/AR0
         F_t(3,3) = AR_t/AR0
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 4 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            call mdet(F_tau,detF_tau)
            call mdet(F_t,detF_t)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
            F_t = ((detFc_tau/detF_tau)**third)*F_t
         endif
         call mdet(F_tau,detF)


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point 
         !
         call integ(props,nprops,dtime,nVisco,
     +        F_tau,ER,A_t,
     +        T_tau,D_tau,A_tau,deltaA,
     +        SpUUMod,SpUPhiMod,SpPhiUMod,SpPhiPhiMod,
     +        stat)
         !
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         do n=1,nVisco
            k = (n-1)*6
            svars(1+k+jj) = A_tau(n,1,1)
            svars(2+k+jj) = A_tau(n,2,2)
            svars(3+k+jj) = A_tau(n,3,3)
            svars(4+k+jj) = A_tau(n,2,3)
            svars(5+k+jj) = A_tau(n,1,3)
            svars(6+k+jj) = A_tau(n,1,2)
         enddo
         jj = jj + nlSdv ! setup for the next intPt



         ! Time stepping algorithim based on the constitutive response.
         !  Here based on the change of the viscous stretch increment.
         !
         Limit = 0.1d0
         umeror = deltaA/Limit
         !
         if(umeror.le.half) then
            pnewdt = 1.5d0
         elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
            pnewdt = 1.25d0
         elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
            pnewdt = 0.75d0
         else
            pnewdt = half
         endif



         ! Compute/update the displacement residual vector
         !
         SmatAx(1,1) = T_tau(1,1)
         SmatAx(2,1) = T_tau(2,2)
         SmatAx(3,1) = T_tau(1,2)
         SmatAx(4,1) = T_tau(3,3)
         !
         BmatAx = zero
         do kk=1,nNode
            BmatAx(1,1+nDim*(kk-1)) = dshC(kk,1)
            BmatAx(2,2+nDim*(kk-1)) = dshC(kk,2)
            BmatAx(3,1+nDim*(kk-1)) = dshC(kk,2)
            BmatAx(3,2+nDim*(kk-1)) = dshC(kk,1)
            BmatAx(4,1+nDim*(kk-1)) = sh(kk)/(AR/(two*Pi))
         enddo
         !
         bodyForce = zero ! the body force vector may be specified here
         !
         BodyForceResAX = zero
         do kk=1,nNode
            BodyForceResAx(1+nDim*(kk-1),1) = sh(kk)*bodyForce(1)
            BodyForceResAx(2+nDim*(kk-1),1) = sh(kk)*bodyForce(2)
         enddo
         !
         Ru = Ru + detmapJC*w(intpt)*AR*
     +        (
     +        -matmul(transpose(BmatAx),SmatAx)
     +        + BodyForceResAx
     +        )



         ! Compute/update the electric potential residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         !
         bodyCharge = zero ! the free charge density may be specified here
         !
         Dvec(1,1) = D_tau(1,1)
         Dvec(2,1) = D_tau(2,1)
         !
         Rphi = Rphi + detmapJC*w(intpt)*AR*
     +        (
     +        matmul(dshC,Dvec)
     +        + transpose(Nvec)*bodyCharge
     +        )




         ! Compute/update the displacement tangent matrix
         !
         GmatAx = zero
         do kk=1,nNode
            GmatAx(1,1+nDim*(kk-1)) = dshC(kk,1)
            GmatAx(2,2+nDim*(kk-1)) = dshC(kk,1)
            GmatAx(3,1+nDim*(kk-1)) = dshC(kk,2)
            GmatAx(4,2+nDim*(kk-1)) = dshC(kk,2)
            GmatAx(5,1+nDim*(kk-1)) = sh(kk)/(AR/(two*Pi))
         enddo

         G0matAx = zero
         do kk=1,nNode
            G0matAx(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0matAx(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0matAx(3,1+nDim*(kk-1)) = dshC0(kk,2)
            G0matAx(4,2+nDim*(kk-1)) = dshC0(kk,2)
            G0matAX(5,1+nDim*(kk-1)) = sh0(kk)/ARc
         enddo

         AmatAx = zero
         AmatAx(1,1) = SpUUMod(1,1,1,1)
         AmatAx(1,2) = SpUUMod(1,1,2,1)
         AmatAx(1,3) = SpUUMod(1,1,1,2)
         AmatAx(1,4) = SpUUMod(1,1,2,2)
         AmatAx(1,5) = SpUUMod(1,1,3,3)
         AmatAx(2,1) = SpUUMod(2,1,1,1)
         AmatAx(2,2) = SpUUMod(2,1,2,1)
         AmatAx(2,3) = SpUUMod(2,1,1,2)
         AmatAx(2,4) = SpUUMod(2,1,2,2)
         AmatAx(2,5) = SpUUMod(2,1,3,3)
         AmatAx(3,1) = SpUUMod(1,2,1,1)
         AmatAx(3,2) = SpUUMod(1,2,2,1)
         AmatAx(3,3) = SpUUMod(1,2,1,2)
         AmatAx(3,4) = SpUUMod(1,2,2,2)
         AmatAx(3,5) = SpUUMod(1,2,3,3)
         AmatAx(4,1) = SpUUMod(2,2,1,1)
         AmatAx(4,2) = SpUUMod(2,2,2,1)
         AmatAx(4,3) = SpUUMod(2,2,1,2)
         AmatAx(4,4) = SpUUMod(2,2,2,2)
         AmatAx(4,5) = SpUUMod(2,2,3,3)
         AmatAx(5,1) = SpUUMod(3,3,1,1)
         AmatAx(5,2) = SpUUMod(3,3,2,1)
         AmatAx(5,3) = SpUUMod(3,3,1,2)
         AmatAx(5,4) = SpUUMod(3,3,2,2)
         AmatAx(5,5) = SpUUMod(3,3,3,3)

         QmatAx = zero
         QmatAx(1,1) = third*(AmatAx(1,1)+AmatAx(1,4)+AmatAx(1,5)) 
     +        - (two/three)*T_tau(1,1)
         QmatAx(2,1) = third*(AmatAx(2,1)+AmatAx(2,4)+AmatAx(2,5))
     +        - (two/three)*T_tau(1,2)
         QmatAx(3,1) = third*(AmatAx(3,1)+AmatAx(3,4)+AmatAx(3,5))
     +        - (two/three)*T_tau(1,2)
         QmatAx(4,1) = third*(AmatAx(4,1)+AmatAx(4,4)+AmatAx(4,5))
     +        - (two/three)*T_tau(2,2)
         QmatAx(5,1) = third*(AmatAx(5,1)+AmatAx(5,4)+AmatAx(5,5))
     +        - (two/three)*T_tau(3,3)
         QmatAx(1,4) = QmatAx(1,1)
         QmatAx(2,4) = QmatAx(2,1)
         QmatAx(3,4) = QmatAx(3,1)
         QmatAx(4,4) = QmatAx(4,1)
         QmatAx(5,4) = QmatAx(5,1)
         QmatAx(1,5) = QmatAx(1,1)
         QmatAx(2,5) = QmatAx(2,1)
         QmatAx(3,5) = QmatAx(3,1)
         QmatAx(4,5) = QmatAx(4,1)
         QmatAx(5,5) = QmatAx(5,1)
            

         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            ! This is the tangent using the F-bar method with the
            !  4 node fully integrated linear element
            !
            Kuu = Kuu + detMapJC*w(intpt)*AR*
     +           (
     +           matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +           + matmul(transpose(GmatAx),matmul(QmatAx,
     +           (G0matAx-GmatAx)))
     +           )
         else
            !
            ! This is the tangent NOT using the F-bar method with all
            !  other elements
            !
            Kuu = Kuu + detMapJC*w(intpt)*AR*
     +           (
     +           matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +           )
         endif



         ! Compute/update the electric potential tangent matrix
         !
         Amatphi(1,1) = SpPhiPhiMod(1,1)
         Amatphi(1,2) = SpPhiPhiMod(1,2)
         Amatphi(2,1) = SpPhiPhiMod(2,1)
         Amatphi(2,2) = SpPhiPhiMod(2,2)
         !
         Kphiphi = Kphiphi - detmapJC*w(intPt)*AR*
     +        (
     +        matmul(matmul(dshC,Amatphi),transpose(dshC))
     +        )


         

         ! Compute/update the electric potential - displacement tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         !
         AmatPhiU = zero
         AmatPhiU(1,1) = SpPhiUMod(1,1,1)
         AmatPhiU(1,2) = SpPhiUMod(1,2,1)
         AmatPhiU(1,3) = SpPhiUMod(1,1,2)
         AmatPhiU(1,4) = SpPhiUMod(1,2,2)
         AmatPhiU(1,5) = SpPhiUMod(1,3,3)
         AmatPhiU(2,1) = SpPhiUMod(2,1,1)
         AmatPhiU(2,2) = SpPhiUMod(2,2,1)
         AmatPhiU(2,3) = SpPhiUMod(2,1,2)
         AmatPhiU(2,4) = SpPhiUMod(2,2,2)
         AmatPhiU(2,5) = SpPhiUMod(2,3,3)
         !
         Kphiu = Kphiu - detMapJC*w(intpt)*AR*
     +        (
     +        matmul(matmul(dshC,AmatPhiU),GmatAx)
     +        )


         ! Compute/update the displacement - electric potential tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         !
         AmatUPhi = zero
         AmatUPhi(1,1) = SpUPhiMod(1,1,1)
         AmatUPhi(2,1) = SpUPhiMod(2,1,1)
         AmatUPhi(3,1) = SpUPhiMod(1,2,1)
         AmatUPhi(4,1) = SpUPhiMod(2,2,1)
         AmatUPhi(5,1) = SpUPhiMod(3,3,1)
         AmatUPhi(1,2) = SpUPhiMod(1,1,2)
         AmatUPhi(2,2) = SpUPhiMod(2,1,2)
         AmatUPhi(3,2) = SpUPhiMod(1,2,2)
         AmatUPhi(4,2) = SpUPhiMod(2,2,2)
         AmatUPhi(5,2) = SpUPhiMod(3,3,2)
         !
         Kuphi = Kuphi + detMapJC*w(intpt)*AR*
     +        (
     +        matmul(matmul(transpose(GmatAx),AmatUPhi),transpose(dshC))
     +        )


      enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface terms here
      !
      !
      if(ndload.gt.0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i=1,ndload
            !
            ! based on my convention the face which the boundary condition
            !  acts on is the ``label''
            !
            face = jdltyp(i,1)  ! label
            flux = adlmag(i,1)  ! BC magnitude
            !
            ! For surface charge, flux = omega
            ! For surface tension, flux = gamma
            !
            if((face.ge.1).and.(face.le.4)) then
               !
               ! surface charge applied
               !
               if(face.eq.1) then
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,1)-coordsC(1,2))**two) + 
     +                 ((coordsC(2,1)-coordsC(2,2))**two))
                  Rphi(1,1) = Rphi(1,1) + Pi*flux*Le*
     +                 third*(two*coordsC(1,1)+coordsC(1,2))
                  Rphi(2,1) = Rphi(2,1) + Pi*flux*Le*
     +                 third*(coordsC(1,1)+two*coordsC(1,2))
                  !
                  ! Modify the tangent matrix
                  !
                  Kphiu(1,1) = Kphiu(1,1) - (Pi*flux/Le)*
     +                 (coordsC(1,1)-coordsC(1,2))*
     +                 third*(two*coordsC(1,1)+coordsC(1,2)) -
     +                 Pi*flux*Le*third*two
                  Kphiu(1,2) = Kphiu(1,2) - (Pi*flux/Le)*
     +                 (coordsC(2,1)-coordsC(2,2))*
     +                 third*(two*coordsC(1,1)+coordsC(1,2))
                  Kphiu(1,3) = Kphiu(1,3) - (Pi*flux/Le)*
     +                 (coordsC(1,2)-coordsC(1,1))*
     +                 third*(two*coordsC(1,1)+coordsC(1,2)) -
     +                 Pi*flux*Le*third*one
                  Kphiu(1,4) = Kphiu(1,4) - (Pi*flux/Le)*
     +                 (coordsC(2,2)-coordsC(2,1))*
     +                 third*(two*coordsC(1,1)+coordsC(1,2))
                  Kphiu(2,1) = Kphiu(2,1) - (Pi*flux/Le)*
     +                 (coordsC(1,1)-coordsC(1,2))*
     +                 third*(coordsC(1,1)+two*coordsC(1,2)) -
     +                 Pi*flux*Le*third*one
                  Kphiu(2,2) = Kphiu(2,2) - (Pi*flux/Le)*
     +                 (coordsC(2,1)-coordsC(2,2))*
     +                 third*(coordsC(1,1)+two*coordsC(1,2))
                  Kphiu(2,3) = Kphiu(2,3) - (Pi*flux/Le)*
     +                 (coordsC(1,2)-coordsC(1,1))*
     +                 third*(coordsC(1,1)+two*coordsC(1,2)) -
     +                 Pi*flux*Le*third*two
                  Kphiu(2,4) = Kphiu(2,4) - (Pi*flux/Le)*
     +                 (coordsC(2,2)-coordsC(2,1))*
     +                 third*(coordsC(1,1)+two*coordsC(1,2))
                  !
               elseif(face.eq.2) then
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,2)-coordsC(1,3))**two) + 
     +                 ((coordsC(2,2)-coordsC(2,3))**two))
                  Rphi(2,1) = Rphi(2,1) + Pi*flux*Le*
     +                 third*(two*coordsC(1,2)+coordsC(1,3))
                  Rphi(3,1) = Rphi(3,1) + Pi*flux*Le*
     +                 third*(coordsC(1,2)+two*coordsC(1,3))
                  !
                  ! Modify the tangent matrix
                  !
                  Kphiu(2,3) = Kphiu(2,3) - (Pi*flux/Le)*
     +                 (coordsC(1,2)-coordsC(1,3))*
     +                 third*(two*coordsC(1,2)+coordsC(1,3)) -
     +                 Pi*flux*Le*third*two
                  Kphiu(2,4) = Kphiu(2,4) - (Pi*flux/Le)*
     +                 (coordsC(2,2)-coordsC(2,3))*
     +                 third*(two*coordsC(1,2)+coordsC(1,3))
                  Kphiu(2,5) = Kphiu(2,5) - (Pi*flux/Le)*
     +                 (coordsC(1,3)-coordsC(1,2))*
     +                 third*(two*coordsC(1,2)+coordsC(1,3)) -
     +                 Pi*flux*Le*third*one
                  Kphiu(2,6) = Kphiu(2,6) - (Pi*flux/Le)*
     +                 (coordsC(2,3)-coordsC(2,2))*
     +                 third*(two*coordsC(1,2)+coordsC(1,3))
                  Kphiu(3,3) = Kphiu(3,3) - (Pi*flux/Le)*
     +                 (coordsC(1,2)-coordsC(1,3))*
     +                 third*(coordsC(1,2)+two*coordsC(1,3)) -
     +                 Pi*flux*Le*third*one
                  Kphiu(3,4) = Kphiu(3,4) - (Pi*flux/Le)*
     +                 (coordsC(2,2)-coordsC(2,3))*
     +                 third*(coordsC(1,2)+two*coordsC(1,3))
                  Kphiu(3,5) = Kphiu(3,5) - (Pi*flux/Le)*
     +                 (coordsC(1,3)-coordsC(1,2))*
     +                 third*(coordsC(1,2)+two*coordsC(1,3)) -
     +                 Pi*flux*Le*third*two
                  Kphiu(3,6) = Kphiu(3,6) - (Pi*flux/Le)*
     +                 (coordsC(2,3)-coordsC(2,2))*
     +                 third*(coordsC(1,2)+two*coordsC(1,3))
                  !
               elseif(face.eq.3) then
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,3)-coordsC(1,4))**two) + 
     +                 ((coordsC(2,3)-coordsC(2,4))**two))
                  Rphi(3,1) = Rphi(3,1) + Pi*flux*Le*
     +                 third*(two*coordsC(1,3)+coordsC(1,4))
                  Rphi(4,1) = Rphi(4,1) + Pi*flux*Le*
     +                 third*(coordsC(1,3)+two*coordsC(1,4))
                  !
                  ! Modify the tangent matrix
                  !
                  Kphiu(3,5) = Kphiu(3,5) - (Pi*flux/Le)*
     +                 (coordsC(1,3)-coordsC(1,4))*
     +                 third*(two*coordsC(1,3)+coordsC(1,4)) -
     +                 Pi*flux*Le*third*two
                  Kphiu(3,6) = Kphiu(3,6) - (Pi*flux/Le)*
     +                 (coordsC(2,3)-coordsC(2,4))*
     +                 third*(two*coordsC(1,3)+coordsC(1,4))
                  Kphiu(3,7) = Kphiu(3,7) - (Pi*flux/Le)*
     +                 (coordsC(1,4)-coordsC(1,3))*
     +                 third*(two*coordsC(1,3)+coordsC(1,4)) -
     +                 Pi*flux*Le*third*one
                  Kphiu(3,8) = Kphiu(3,8) - (Pi*flux/Le)*
     +                 (coordsC(2,4)-coordsC(2,3))*
     +                 third*(two*coordsC(1,3)+coordsC(1,4))
                  Kphiu(4,5) = Kphiu(4,5) - (Pi*flux/Le)*
     +                 (coordsC(1,3)-coordsC(1,4))*
     +                 third*(coordsC(1,3)+two*coordsC(1,4)) -
     +                 Pi*flux*Le*third*one
                  Kphiu(4,6) = Kphiu(4,6) - (Pi*flux/Le)*
     +                 (coordsC(2,3)-coordsC(2,4))*
     +                 third*(coordsC(1,3)+two*coordsC(1,4))
                  Kphiu(4,7) = Kphiu(4,7) - (Pi*flux/Le)*
     +                 (coordsC(1,4)-coordsC(1,3))*
     +                 third*(coordsC(1,3)+two*coordsC(1,4)) -
     +                 Pi*flux*Le*third*two
                  Kphiu(4,8) = Kphiu(4,8) - (Pi*flux/Le)*
     +                 (coordsC(2,4)-coordsC(2,3))*
     +                 third*(coordsC(1,3)+two*coordsC(1,4))
                  !
               elseif(face.eq.4) then
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,4)-coordsC(1,1))**two) + 
     +                 ((coordsC(2,4)-coordsC(2,1))**two))
                  Rphi(4,1) = Rphi(4,1) + Pi*flux*Le*
     +                 third*(two*coordsC(1,4)+coordsC(1,1))
                  Rphi(1,1) = Rphi(1,1) + Pi*flux*Le*
     +                 third*(coordsC(1,4)+two*coordsC(1,1))
                  !
                  ! Modify the tangent matrix
                  !
                  Kphiu(4,7) = Kphiu(4,7) - (Pi*flux/Le)*
     +                 (coordsC(1,4)-coordsC(1,1))*
     +                 third*(two*coordsC(1,4)+coordsC(1,1)) -
     +                 Pi*flux*Le*third*two
                  Kphiu(4,8) = Kphiu(4,8) - (Pi*flux/Le)*
     +                 (coordsC(2,4)-coordsC(2,1))*
     +                 third*(two*coordsC(1,4)+coordsC(1,1))
                  Kphiu(4,1) = Kphiu(4,1) - (Pi*flux/Le)*
     +                 (coordsC(1,1)-coordsC(1,4))*
     +                 third*(two*coordsC(1,4)+coordsC(1,1)) -
     +                 Pi*flux*Le*third*one
                  Kphiu(4,2) = Kphiu(4,2) - (Pi*flux/Le)*
     +                 (coordsC(2,1)-coordsC(2,4))*
     +                 third*(two*coordsC(1,4)+coordsC(1,1))
                  Kphiu(1,7) = Kphiu(1,7) - (Pi*flux/Le)*
     +                 (coordsC(1,4)-coordsC(1,1))*
     +                 third*(coordsC(1,4)+two*coordsC(1,1)) -
     +                 Pi*flux*Le*third*one
                  Kphiu(1,8) = Kphiu(1,8) - (Pi*flux/Le)*
     +                 (coordsC(2,4)-coordsC(2,1))*
     +                 third*(coordsC(1,4)+two*coordsC(1,1))
                  Kphiu(1,1) = Kphiu(1,1) - (Pi*flux/Le)*
     +                 (coordsC(1,1)-coordsC(1,4))*
     +                 third*(coordsC(1,4)+two*coordsC(1,1)) -
     +                 Pi*flux*Le*third*two
                  Kphiu(1,2) = Kphiu(1,2) - (Pi*flux/Le)*
     +                 (coordsC(2,1)-coordsC(2,4))*
     +                 third*(coordsC(1,4)+two*coordsC(1,1))
                  !
               endif
               !
            elseif((face.ge.11).and.(face.le.14)) then
               !
               ! surface tension applied
               !
               if(face.eq.11) then
                  !
                  ! Surface tension on face 1
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,1)-coordsC(1,2))**two) + 
     +                 ((coordsC(2,1)-coordsC(2,2))**two))
                  Ru(1,1) = Ru(1,1)
     +                 - pi*flux*(coordsC(1,1)+coordsC(1,2))
     +                 *(coordsC(1,1)-coordsC(1,2))/Le -
     +                 pi*flux*Le
                  Ru(2,1) = Ru(2,1) 
     +                 - pi*flux*(coordsC(1,1)+coordsC(1,2))
     +                 *(coordsC(2,1)-coordsC(2,2))/Le
                  Ru(3,1) = Ru(3,1) 
     +                 - pi*flux*(coordsC(1,1)+coordsC(1,2))
     +                 *(coordsC(1,2)-coordsC(1,1))/Le -
     +                 pi*flux*Le
                  Ru(4,1) = Ru(4,1) 
     +                 - pi*flux*(coordsC(1,1)+coordsC(1,2))
     +                 *(coordsC(2,2)-coordsC(2,1))/Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kuu(1,1) = Kuu(1,1) + two*coordsC(1,1)*pi*flux/Le - 
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(1,1)-coordsC(1,2))*
     +                 (coordsC(1,1)-coordsC(1,2))*pi*flux/(Le**three)+
     +                 (coordsC(1,1)-coordsC(1,2))*pi*flux/Le
                  Kuu(1,2) = Kuu(1,2) - 
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(1,1)-coordsC(1,2))*
     +                 (coordsC(2,1)-coordsC(2,2))*pi*flux/(Le**three)+
     +                 (coordsC(2,1)-coordsC(2,2))*pi*flux/Le
                  Kuu(1,3) = Kuu(1,3) - two*coordsC(1,2)*pi*flux/Le - 
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(1,1)-coordsC(1,2))*
     +                 (coordsC(1,2)-coordsC(1,1))*pi*flux/(Le**three)+
     +                 (coordsC(1,2)-coordsC(1,1))*pi*flux/Le
                  Kuu(1,4) = Kuu(1,4) - 
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(1,1)-coordsC(1,2))*
     +                 (coordsC(2,2)-coordsC(2,1))*pi*flux/(Le**three)+
     +                 (coordsC(2,2)-coordsC(2,1))*pi*flux/Le
                  Kuu(2,1) = Kuu(2,1) + 
     +                 (coordsC(2,1)-coordsC(2,2))*pi*flux/Le -
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(2,1)-coordsC(2,2))*
     +                 (coordsC(1,1)-coordsC(1,2))*pi*flux/(Le**three)
                  Kuu(2,2) = Kuu(2,2) + 
     +                 (coordsC(1,1)+coordsC(1,2))*pi*flux/Le -
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(2,1)-coordsC(2,2))*
     +                 (coordsC(2,1)-coordsC(2,2))*pi*flux/(Le**three)
                  Kuu(2,3) = Kuu(2,3) + 
     +                 (coordsC(2,1)-coordsC(2,2))*pi*flux/Le -
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(2,1)-coordsC(2,2))*
     +                 (coordsC(1,2)-coordsC(1,1))*pi*flux/(Le**three)
                  Kuu(2,4) = Kuu(2,4) - 
     +                 (coordsC(1,1)+coordsC(1,2))*pi*flux/Le -
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(2,1)-coordsC(2,2))*
     +                 (coordsC(2,2)-coordsC(2,1))*pi*flux/(Le**three)
                  Kuu(3,1) = Kuu(3,1) - two*coordsC(1,1)*pi*flux/Le - 
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(1,2)-coordsC(1,1))*
     +                 (coordsC(1,1)-coordsC(1,2))*pi*flux/(Le**three)+
     +                 (coordsC(1,1)-coordsC(1,2))*pi*flux/Le
                  Kuu(3,2) = Kuu(3,2) - 
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(1,2)-coordsC(1,1))*
     +                 (coordsC(2,1)-coordsC(2,2))*pi*flux/(Le**three)+
     +                 (coordsC(2,1)-coordsC(2,2))*pi*flux/Le
                  Kuu(3,3) = Kuu(3,3) + two*coordsC(1,2)*pi*flux/Le - 
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(1,2)-coordsC(1,1))*
     +                 (coordsC(1,2)-coordsC(1,1))*pi*flux/(Le**three)+
     +                 (coordsC(1,2)-coordsC(1,1))*pi*flux/Le
                  Kuu(3,4) = Kuu(3,4) - 
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(1,2)-coordsC(1,1))*
     +                 (coordsC(2,2)-coordsC(2,1))*pi*flux/(Le**three)+
     +                 (coordsC(2,2)-coordsC(2,1))*pi*flux/Le
                  Kuu(4,1) = Kuu(4,1) - 
     +                 (coordsC(2,1)-coordsC(2,2))*pi*flux/Le -
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(2,2)-coordsC(2,1))*
     +                 (coordsC(1,1)-coordsC(1,2))*pi*flux/(Le**three)
                  Kuu(4,2) = Kuu(4,2) - 
     +                 (coordsC(1,1)+coordsC(1,2))*pi*flux/Le -
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(2,2)-coordsC(2,1))*
     +                 (coordsC(2,1)-coordsC(2,2))*pi*flux/(Le**three)
                  Kuu(4,3) = Kuu(4,3) - 
     +                 (coordsC(2,1)-coordsC(2,2))*pi*flux/Le -
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(2,2)-coordsC(2,1))*
     +                 (coordsC(1,2)-coordsC(1,1))*pi*flux/(Le**three)
                  Kuu(4,4) = Kuu(4,4) + 
     +                 (coordsC(1,1)+coordsC(1,2))*pi*flux/Le -
     +                 (coordsC(1,1)+coordsC(1,2))*
     +                 (coordsC(2,2)-coordsC(2,1))*
     +                 (coordsC(2,2)-coordsC(2,1))*pi*flux/(Le**three)
                  !
               elseif(face.eq.12) then
                  !
                  ! Surface tension on face 2
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,2)-coordsC(1,3))**two) + 
     +                 ((coordsC(2,2)-coordsC(2,3))**two))
                  Ru(3,1) = Ru(3,1) 
     +                 - pi*flux*(coordsC(1,2)+coordsC(1,3))
     +                 *(coordsC(1,2)-coordsC(1,3))/Le -
     +                 pi*flux*Le
                  Ru(4,1) = Ru(4,1) 
     +                 - pi*flux*(coordsC(1,2)+coordsC(1,3))
     +                 *(coordsC(2,2)-coordsC(2,3))/Le
                  Ru(5,1) = Ru(5,1) 
     +                 - pi*flux*(coordsC(1,2)+coordsC(1,3))
     +                 *(coordsC(1,3)-coordsC(1,2))/Le -
     +                 pi*flux*Le
                  Ru(6,1) = Ru(6,1) 
     +                 - pi*flux*(coordsC(1,2)+coordsC(1,3))
     +                 *(coordsC(2,3)-coordsC(2,2))/Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kuu(3,3) = Kuu(3,3) + two*coordsC(1,2)*pi*flux/Le - 
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(1,2)-coordsC(1,3))*
     +                 (coordsC(1,2)-coordsC(1,3))*pi*flux/(Le**three)+
     +                 (coordsC(1,2)-coordsC(1,3))*pi*flux/Le
                  Kuu(3,4) = Kuu(3,4) - 
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(1,2)-coordsC(1,3))*
     +                 (coordsC(2,2)-coordsC(2,3))*pi*flux/(Le**three)+
     +                 (coordsC(2,2)-coordsC(2,3))*pi*flux/Le
                  Kuu(3,5) = Kuu(3,5) - two*coordsC(1,3)*pi*flux/Le - 
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(1,2)-coordsC(1,3))*
     +                 (coordsC(1,3)-coordsC(1,2))*pi*flux/(Le**three)+
     +                 (coordsC(1,3)-coordsC(1,2))*pi*flux/Le
                  Kuu(3,6) = Kuu(3,6) - 
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(1,2)-coordsC(1,3))*
     +                 (coordsC(2,3)-coordsC(2,2))*pi*flux/(Le**three)+
     +                 (coordsC(2,3)-coordsC(2,2))*pi*flux/Le
                  Kuu(4,3) = Kuu(4,3) + 
     +                 (coordsC(2,2)-coordsC(2,3))*pi*flux/Le -
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(2,2)-coordsC(2,3))*
     +                 (coordsC(1,2)-coordsC(1,3))*pi*flux/(Le**three)
                  Kuu(4,4) = Kuu(4,4) + 
     +                 (coordsC(1,2)+coordsC(1,3))*pi*flux/Le -
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(2,2)-coordsC(2,3))*
     +                 (coordsC(2,2)-coordsC(2,3))*pi*flux/(Le**three)
                  Kuu(4,5) = Kuu(4,5) + 
     +                 (coordsC(2,2)-coordsC(2,3))*pi*flux/Le -
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(2,2)-coordsC(2,3))*
     +                 (coordsC(1,3)-coordsC(1,2))*pi*flux/(Le**three)
                  Kuu(4,6) = Kuu(4,6) - 
     +                 (coordsC(1,2)+coordsC(1,3))*pi*flux/Le -
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(2,2)-coordsC(2,3))*
     +                 (coordsC(2,3)-coordsC(2,2))*pi*flux/(Le**three)
                  Kuu(5,3) = Kuu(5,3) - two*coordsC(1,2)*pi*flux/Le - 
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(1,3)-coordsC(1,2))*
     +                 (coordsC(1,2)-coordsC(1,3))*pi*flux/(Le**three)+
     +                 (coordsC(1,2)-coordsC(1,3))*pi*flux/Le
                  Kuu(5,4) = Kuu(5,4) - 
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(1,3)-coordsC(1,2))*
     +                 (coordsC(2,2)-coordsC(2,3))*pi*flux/(Le**three)+
     +                 (coordsC(2,2)-coordsC(2,3))*pi*flux/Le
                  Kuu(5,5) = Kuu(5,5) + two*coordsC(1,3)*pi*flux/Le - 
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(1,3)-coordsC(1,2))*
     +                 (coordsC(1,3)-coordsC(1,2))*pi*flux/(Le**three)+
     +                 (coordsC(1,3)-coordsC(1,2))*pi*flux/Le
                  Kuu(5,6) = Kuu(5,6) - 
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(1,3)-coordsC(1,2))*
     +                 (coordsC(2,3)-coordsC(2,2))*pi*flux/(Le**three)+
     +                 (coordsC(2,3)-coordsC(2,2))*pi*flux/Le
                  Kuu(6,3) = Kuu(6,3) - 
     +                 (coordsC(2,2)-coordsC(2,3))*pi*flux/Le -
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(2,3)-coordsC(2,2))*
     +                 (coordsC(1,2)-coordsC(1,3))*pi*flux/(Le**three)
                  Kuu(6,4) = Kuu(6,4) - 
     +                 (coordsC(1,2)+coordsC(1,3))*pi*flux/Le -
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(2,3)-coordsC(2,2))*
     +                 (coordsC(2,2)-coordsC(2,3))*pi*flux/(Le**three)
                  Kuu(6,5) = Kuu(6,5) - 
     +                 (coordsC(2,2)-coordsC(2,3))*pi*flux/Le -
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(2,3)-coordsC(2,2))*
     +                 (coordsC(1,3)-coordsC(1,2))*pi*flux/(Le**three)
                  Kuu(6,6) = Kuu(6,6) + 
     +                 (coordsC(1,2)+coordsC(1,3))*pi*flux/Le -
     +                 (coordsC(1,2)+coordsC(1,3))*
     +                 (coordsC(2,3)-coordsC(2,2))*
     +                 (coordsC(2,3)-coordsC(2,2))*pi*flux/(Le**three)
               !
               elseif(face.eq.13) then
                  !
                  ! Surface tension on face 3
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,3)-coordsC(1,4))**two) + 
     +                 ((coordsC(2,3)-coordsC(2,4))**two))
                  Ru(5,1) = Ru(5,1) 
     +                 - pi*flux*(coordsC(1,3)+coordsC(1,4))
     +                 *(coordsC(1,3)-coordsC(1,4))/Le -
     +                 pi*flux*Le
                  Ru(6,1) = Ru(6,1) 
     +                 - pi*flux*(coordsC(1,3)+coordsC(1,4))
     +                 *(coordsC(2,3)-coordsC(2,4))/Le
                  Ru(7,1) = Ru(7,1) 
     +                 - pi*flux*(coordsC(1,3)+coordsC(1,4))
     +                 *(coordsC(1,4)-coordsC(1,3))/Le -
     +                 pi*flux*Le
                  Ru(8,1) = Ru(8,1) 
     +                 - pi*flux*(coordsC(1,3)+coordsC(1,4))
     +                 *(coordsC(2,4)-coordsC(2,3))/Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kuu(5,5) = Kuu(5,5) + two*coordsC(1,3)*pi*flux/Le - 
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(1,3)-coordsC(1,4))*
     +                 (coordsC(1,3)-coordsC(1,4))*pi*flux/(Le**three)+
     +                 (coordsC(1,3)-coordsC(1,4))*pi*flux/Le
                  Kuu(5,6) = Kuu(5,6) - 
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(1,3)-coordsC(1,4))*
     +                 (coordsC(2,3)-coordsC(2,4))*pi*flux/(Le**three)+
     +                 (coordsC(2,3)-coordsC(2,4))*pi*flux/Le
                  Kuu(5,7) = Kuu(5,7) - two*coordsC(1,4)*pi*flux/Le - 
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(1,3)-coordsC(1,4))*
     +                 (coordsC(1,4)-coordsC(1,3))*pi*flux/(Le**three)+
     +                 (coordsC(1,4)-coordsC(1,3))*pi*flux/Le
                  Kuu(5,8) = Kuu(5,8) - 
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(1,3)-coordsC(1,4))*
     +                 (coordsC(2,4)-coordsC(2,3))*pi*flux/(Le**three)+
     +                 (coordsC(2,4)-coordsC(2,3))*pi*flux/Le
                  Kuu(6,5) = Kuu(6,5) + 
     +                 (coordsC(2,3)-coordsC(2,4))*pi*flux/Le -
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(2,3)-coordsC(2,4))*
     +                 (coordsC(1,3)-coordsC(1,4))*pi*flux/(Le**three)
                  Kuu(6,6) = Kuu(6,6) + 
     +                 (coordsC(1,3)+coordsC(1,4))*pi*flux/Le -
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(2,3)-coordsC(2,4))*
     +                 (coordsC(2,3)-coordsC(2,4))*pi*flux/(Le**three)
                  Kuu(6,7) = Kuu(6,7) + 
     +                 (coordsC(2,3)-coordsC(2,4))*pi*flux/Le -
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(2,3)-coordsC(2,4))*
     +                 (coordsC(1,4)-coordsC(1,3))*pi*flux/(Le**three)
                  Kuu(6,8) = Kuu(6,8) - 
     +                 (coordsC(1,3)+coordsC(1,4))*pi*flux/Le -
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(2,3)-coordsC(2,4))*
     +                 (coordsC(2,4)-coordsC(2,3))*pi*flux/(Le**three)
                  Kuu(7,5) = Kuu(7,5) - two*coordsC(1,3)*pi*flux/Le - 
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(1,4)-coordsC(1,3))*
     +                 (coordsC(1,3)-coordsC(1,4))*pi*flux/(Le**three)+
     +                 (coordsC(1,3)-coordsC(1,4))*pi*flux/Le
                  Kuu(7,6) = Kuu(7,6) - 
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(1,4)-coordsC(1,3))*
     +                 (coordsC(2,3)-coordsC(2,4))*pi*flux/(Le**three)+
     +                 (coordsC(2,3)-coordsC(2,4))*pi*flux/Le
                  Kuu(7,7) = Kuu(7,7) + two*coordsC(1,4)*pi*flux/Le - 
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(1,4)-coordsC(1,3))*
     +                 (coordsC(1,4)-coordsC(1,3))*pi*flux/(Le**three)+
     +                 (coordsC(1,4)-coordsC(1,3))*pi*flux/Le
                  Kuu(7,8) = Kuu(7,8) - 
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(1,4)-coordsC(1,3))*
     +                 (coordsC(2,4)-coordsC(2,3))*pi*flux/(Le**three)+
     +                 (coordsC(2,4)-coordsC(2,3))*pi*flux/Le
                  Kuu(8,5) = Kuu(8,5) - 
     +                 (coordsC(2,3)-coordsC(2,4))*pi*flux/Le -
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(2,4)-coordsC(2,3))*
     +                 (coordsC(1,3)-coordsC(1,4))*pi*flux/(Le**three)
                  Kuu(8,6) = Kuu(8,6) - 
     +                 (coordsC(1,3)+coordsC(1,4))*pi*flux/Le -
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(2,4)-coordsC(2,3))*
     +                 (coordsC(2,3)-coordsC(2,4))*pi*flux/(Le**three)
                  Kuu(8,7) = Kuu(8,7) - 
     +                 (coordsC(2,3)-coordsC(2,4))*pi*flux/Le -
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(2,4)-coordsC(2,3))*
     +                 (coordsC(1,4)-coordsC(1,3))*pi*flux/(Le**three)
                  Kuu(8,8) = Kuu(8,8) + 
     +                 (coordsC(1,3)+coordsC(1,4))*pi*flux/Le -
     +                 (coordsC(1,3)+coordsC(1,4))*
     +                 (coordsC(2,4)-coordsC(2,3))*
     +                 (coordsC(2,4)-coordsC(2,3))*pi*flux/(Le**three)
                  !
               elseif(face.eq.14) then
                  !
                  ! Surface tension on face 4
                  !
                  ! Modify the displacement residual
                  !
                  Le = dsqrt(((coordsC(1,4)-coordsC(1,1))**two) + 
     +                 ((coordsC(2,4)-coordsC(2,1))**two))
                  Ru(7,1) = Ru(7,1) 
     +                 - pi*flux*(coordsC(1,4)+coordsC(1,1))
     +                 *(coordsC(1,4)-coordsC(1,1))/Le -
     +                 pi*flux*Le
                  Ru(8,1) = Ru(8,1) 
     +                 - pi*flux*(coordsC(1,4)+coordsC(1,1))
     +                 *(coordsC(2,4)-coordsC(2,1))/Le
                  Ru(1,1) = Ru(1,1) 
     +                 - pi*flux*(coordsC(1,4)+coordsC(1,1))
     +                 *(coordsC(1,1)-coordsC(1,4))/Le -
     +                 pi*flux*Le
                  Ru(2,1) = Ru(2,1) 
     +                 - pi*flux*(coordsC(1,4)+coordsC(1,1))
     +                 *(coordsC(2,1)-coordsC(2,4))/Le
                  !
                  ! Modify the tangent matrix
                  !
                  Kuu(7,7) = Kuu(7,7) + two*coordsC(1,4)*pi*flux/Le - 
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(1,4)-coordsC(1,1))*
     +                 (coordsC(1,4)-coordsC(1,1))*pi*flux/(Le**three)+
     +                 (coordsC(1,4)-coordsC(1,1))*pi*flux/Le
                  Kuu(7,8) = Kuu(7,8) - 
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(1,4)-coordsC(1,1))*
     +                 (coordsC(2,4)-coordsC(2,1))*pi*flux/(Le**three)+
     +                 (coordsC(2,4)-coordsC(2,1))*pi*flux/Le
                  Kuu(7,1) = Kuu(7,1) - two*coordsC(1,1)*pi*flux/Le - 
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(1,4)-coordsC(1,1))*
     +                 (coordsC(1,1)-coordsC(1,4))*pi*flux/(Le**three)+
     +                 (coordsC(1,1)-coordsC(1,4))*pi*flux/Le
                  Kuu(7,2) = Kuu(7,2) - 
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(1,4)-coordsC(1,1))*
     +                 (coordsC(2,1)-coordsC(2,4))*pi*flux/(Le**three)+
     +                 (coordsC(2,1)-coordsC(2,4))*pi*flux/Le
                  Kuu(8,7) = Kuu(8,7) + 
     +                 (coordsC(2,4)-coordsC(2,1))*pi*flux/Le -
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(2,4)-coordsC(2,1))*
     +                 (coordsC(1,4)-coordsC(1,1))*pi*flux/(Le**three)
                  Kuu(8,8) = Kuu(8,8) + 
     +                 (coordsC(1,4)+coordsC(1,1))*pi*flux/Le -
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(2,4)-coordsC(2,1))*
     +                 (coordsC(2,4)-coordsC(2,1))*pi*flux/(Le**three)
                  Kuu(8,1) = Kuu(8,1) + 
     +                 (coordsC(2,4)-coordsC(2,1))*pi*flux/Le -
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(2,4)-coordsC(2,1))*
     +                 (coordsC(1,1)-coordsC(1,4))*pi*flux/(Le**three)
                  Kuu(8,2) = Kuu(8,2) - 
     +                 (coordsC(1,4)+coordsC(1,1))*pi*flux/Le -
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(2,4)-coordsC(2,1))*
     +                 (coordsC(2,1)-coordsC(2,4))*pi*flux/(Le**three)
                  Kuu(1,7) = Kuu(1,7) - two*coordsC(1,4)*pi*flux/Le - 
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(1,1)-coordsC(1,4))*
     +                 (coordsC(1,4)-coordsC(1,1))*pi*flux/(Le**three)+
     +                 (coordsC(1,4)-coordsC(1,1))*pi*flux/Le
                  Kuu(1,8) = Kuu(1,8) - 
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(1,1)-coordsC(1,4))*
     +                 (coordsC(2,4)-coordsC(2,1))*pi*flux/(Le**three)+
     +                 (coordsC(2,4)-coordsC(2,1))*pi*flux/Le
                  Kuu(1,1) = Kuu(1,1) + two*coordsC(1,1)*pi*flux/Le - 
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(1,1)-coordsC(1,4))*
     +                 (coordsC(1,1)-coordsC(1,4))*pi*flux/(Le**three)+
     +                 (coordsC(1,1)-coordsC(1,4))*pi*flux/Le
                  Kuu(1,2) = Kuu(1,2) - 
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(1,1)-coordsC(1,4))*
     +                 (coordsC(2,1)-coordsC(2,4))*pi*flux/(Le**three)+
     +                 (coordsC(2,1)-coordsC(2,4))*pi*flux/Le
                  Kuu(2,7) = Kuu(2,7) - 
     +                 (coordsC(2,4)-coordsC(2,1))*pi*flux/Le -
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(2,1)-coordsC(2,4))*
     +                 (coordsC(1,4)-coordsC(1,1))*pi*flux/(Le**three)
                  Kuu(2,8) = Kuu(2,8) - 
     +                 (coordsC(1,4)+coordsC(1,1))*pi*flux/Le -
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(2,1)-coordsC(2,4))*
     +                 (coordsC(2,4)-coordsC(2,1))*pi*flux/(Le**three)
                  Kuu(2,1) = Kuu(2,1) - 
     +                 (coordsC(2,4)-coordsC(2,1))*pi*flux/Le -
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(2,1)-coordsC(2,4))*
     +                 (coordsC(1,1)-coordsC(1,4))*pi*flux/(Le**three)
                  Kuu(2,2) = Kuu(2,2) + 
     +                 (coordsC(1,4)+coordsC(1,1))*pi*flux/Le -
     +                 (coordsC(1,4)+coordsC(1,1))*
     +                 (coordsC(2,1)-coordsC(2,4))*
     +                 (coordsC(2,1)-coordsC(2,4))*pi*flux/(Le**three)
                  !
               endif
               !
            else
               write(*,*) 'Unknown face=',face
               write(80,*) 'Unknown face=',face
               call xit
            endif

         enddo ! loop over ndload
      endif ! ndload.gt.0 or not
      !
      ! End loop over surface terms
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Rphi,Kuu,Kuphi,Kphiu,Kphiphi,
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine UAX4

************************************************************************

      subroutine U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS,nVisco)
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 u(nNode,3),du(nNode,ndofel),phi(nNode),coordsC(mcrd,nNode)
      real*8 uNew(nNode,ndofel),uOld(nNode,ndofel),u_t(nNode,ndofel)
      real*8 dPhi(nNode)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,nVisco,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,nIntS
      integer nlSdv,kk,IntS,faceFlag

      real*8 Iden(3,3),Le,theta0,phi0,Ru(3*nNode,1),Rphi(nNode,1)
      real*8 Kuu(3*nNode,3*nNode),Kphiphi(nNode,nNode),sh0(nNode)
      real*8 dshxi(nNode,3),dsh0(nNode,3),dshC0(nNode,3),detMapJ0C
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),D_tau(3,1)
      real*8 sh(nNode),detMapJ,dsh(nNode,3),detMapJC,Limit,umeror
      real*8 dshC(nNode,3),F_tau(3,3),A_t(nVisco,3,3),A_tau(nVisco,3,3)
      real*8 F_t(3,3),detF_tau,xi(nInt,3),detF,T_tau(3,3),bodyForce(3)
      real*8 SpUUMod(3,3,3,3),bodyCharge,Gmat(9,3*nNode),Amat(9,9)
      real*8 Smat(6,1),Bmat(6,3*nNode),BodyForceRes(3*nNode,1),flux
      real*8 G0mat(9,3*nNode),Qmat(9,9),dA,Nvec(1,nNode),AmatPhiU(3,9)
      real*8 xLocal(nIntS),yLocal(nIntS),zLocal(nIntS),wS(nIntS),detF_t
      real*8 Kuphi(3*nNode,nNode),Kphiu(nNode,3*nNode),AmatUPhi(9,3)
      real*8 SpPhiUMod(3,3,3),SpPhiPhiMod(3,3),SpUPhiMod(3,3,3),detMapJ0
      real*8 ER(3,1),deltaA


      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Get element parameters
      !
      nlSdv  = jprops(2) ! number of local sdv's per integ point


      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rphi  = zero
      Kuu = zero
      Kphiphi = zero
      Kuphi = zero
      Kphiu = zero
      Energy = zero



      ! Obtain nodal displacements and electric potentials
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         k = k + 1
         phi(i) = Uall(k)
         dPhi(i) = DUall(k,1)
      enddo



      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo



      ! Impose any time-stepping changes on the increments of
      !  electric potential or displacement if you want
      !
      ! electric potential increment
      !
      do i=1,nNode
         if(dabs(dPhi(i)).gt.1.d6) then
            pnewdt = 0.5
            return
         endif
      enddo
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,7))**two) + 
     +     ((coordsC(2,1)-coordsC(2,7))**two) +
     +     ((coordsC(3,1)-coordsC(3,7))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.d0*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo



      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Here, check for hourglass stabilization and get
      !  the deformation gradient for use in the `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures, 33, 3277-3296.
      !
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.8) then
         call calcShape3DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         write(80,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      call mapShape3D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif



      ! Map shape functions from local to global current coordinate system
      !
      call mapShape3D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Calculate the deformation gradient at the element centriod
      !  at the the begining and end of the increment for use in 
      !  the `F-bar' method
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      call mdet(Fc_tau,detFc_tau)
      call mdet(Fc_t,detFc_t)
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nInt.eq.1) then
         call xint3D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
      elseif(nInt.eq.8) then
         call xint3D8pt(xi,w,nIntPt) ! 8-pt integration, nInt=8 above
      else
         write(*,*) 'Invalid number of int points, nInt=',nInt
         write(80,*) 'Invalid number of int points, nInt=',nInt
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step.
            !  Give initial conditions.
            !
            do n=1,nVisco
               A_t(n,:,:) = Iden
            enddo
            !
         else
            !
            ! this is not the first increment, read old values
            !
            do n=1,nVisco
               k = (n-1)*6
               A_t(n,1,1) = svars(1+k+jj)
               A_t(n,2,2) = svars(2+k+jj)
               A_t(n,3,3) = svars(3+k+jj)
               A_t(n,2,3) = svars(4+k+jj)
               A_t(n,3,2) = A_t(n,2,3)
               A_t(n,1,3) = svars(5+k+jj)
               A_t(n,3,1) = A_t(n,1,3)
               A_t(n,1,2) = svars(6+k+jj)
               A_t(n,2,1) = A_t(n,1,2)
            enddo
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.8) then
            call calcShape3DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.8'
            write(80,*) 'Incorrect number of nodes: nNode.ne.8'
            call xit
         endif


         ! Map shape functions from local to global reference coordinate system
         !
         call mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Map shape functions from local to global current coordinate system
         !
         call mapShape3D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Obtain the referential electric field at this integration point
         !
         ER = zero
         do k=1,nNode
            do i=1,nDim
               ER(i,1) = ER(i,1) + phi(k)*dsh(k,i)
            enddo
         enddo


         ! Obtain, and modify the deformation gradient at this integration
         !  point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.  Also, take care of plane-strain or axisymetric
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 8 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.8).and.(nInt.eq.8)) then
            call mdet(F_tau,detF_tau)
            call mdet(F_t,detF_t)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
            F_t = ((detFc_tau/detF_tau)**third)*F_t
         endif
         call mdet(F_tau,detF)


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the time integration at this integ. point to compute
         !  all the specific forms and parameters needed for the solution
         !
         call integ(props,nprops,dtime,nVisco,
     +        F_tau,ER,A_t,
     +        T_tau,D_tau,A_tau,deltaA,
     +        SpUUMod,SpUPhiMod,SpPhiUMod,SpPhiPhiMod,
     +        stat)
         !
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         do n=1,nVisco
            k = (n-1)*6
            svars(1+k+jj) = A_tau(n,1,1)
            svars(2+k+jj) = A_tau(n,2,2)
            svars(3+k+jj) = A_tau(n,3,3)
            svars(4+k+jj) = A_tau(n,2,3)
            svars(5+k+jj) = A_tau(n,1,3)
            svars(6+k+jj) = A_tau(n,1,2)
         enddo
         jj = jj + nlSdv ! setup for the next intPt


         ! Time stepping algorithim based on the constitutive response.
         !  Here based on the change of the viscous stretch increment.
         !
         Limit = 0.1d0
         umeror = deltaA/Limit
         !
         if(umeror.le.half) then
            pnewdt = 1.5d0
         elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
            pnewdt = 1.25d0
         elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
            pnewdt = 0.75d0
         else
            pnewdt = half
         endif


         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(3,3)
         Smat(4,1) = T_tau(1,2)
         Smat(5,1) = T_tau(2,3)
         Smat(6,1) = T_tau(1,3)
         !
         Bmat = zero
         do kk=1,nNode
            Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,3+nDim*(kk-1)) = dshC(kk,3)
            Bmat(4,1+nDim*(kk-1)) = dshC(kk,2)
            Bmat(4,2+nDim*(kk-1)) = dshC(kk,1)
            Bmat(5,2+nDim*(kk-1)) = dshC(kk,3)
            Bmat(5,3+nDim*(kk-1)) = dshC(kk,2)
            Bmat(6,1+nDim*(kk-1)) = dshC(kk,3)
            Bmat(6,3+nDim*(kk-1)) = dshC(kk,1)
         enddo
         !
         bodyForce = zero ! the body force vector may be specified here
         !
         BodyForceRes = zero
         do kk=1,nNode
            BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*bodyForce(1)
            BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*bodyForce(2)
            BodyForceRes(3+nDim*(kk-1),1) = sh(kk)*bodyForce(3)
         enddo
         !
         Ru = Ru + detmapJC*w(intpt)*
     +        (
     +        -matmul(transpose(Bmat),Smat)
     +        + BodyForceRes
     +        )


         ! Compute/update the electric potential residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         !
         bodyCharge = zero ! the free charge density may be specified here
         !
         Rphi = Rphi + detmapJC*w(intpt)*
     +        (
     +        matmul(dshC,D_tau)
     +        + transpose(Nvec)*bodyCharge
     +        )



         ! Compute/update the displacement tangent matrix
         !
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
            Gmat(3,3+nDim*(kk-1)) = dshC(kk,1)
            Gmat(4,1+nDim*(kk-1)) = dshC(kk,2)
            Gmat(5,2+nDim*(kk-1)) = dshC(kk,2)
            Gmat(6,3+nDim*(kk-1)) = dshC(kk,2)
            Gmat(7,1+nDim*(kk-1)) = dshC(kk,3)
            Gmat(8,2+nDim*(kk-1)) = dshC(kk,3)
            Gmat(9,3+nDim*(kk-1)) = dshC(kk,3)
         enddo

         G0mat = zero
         do kk=1,nNode
            G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(3,3+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(4,1+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(5,2+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(6,3+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(7,1+nDim*(kk-1)) = dshC0(kk,3)
            G0mat(8,2+nDim*(kk-1)) = dshC0(kk,3)
            G0mat(9,3+nDim*(kk-1)) = dshC0(kk,3)
         enddo

         Amat = zero
         Amat(1,1) = SpUUMod(1,1,1,1)
         Amat(1,2) = SpUUMod(1,1,2,1)
         Amat(1,3) = SpUUMod(1,1,3,1)
         Amat(1,4) = SpUUMod(1,1,1,2)
         Amat(1,5) = SpUUMod(1,1,2,2)
         Amat(1,6) = SpUUMod(1,1,3,2)
         Amat(1,7) = SpUUMod(1,1,1,3)
         Amat(1,8) = SpUUMod(1,1,2,3)
         Amat(1,9) = SpUUMod(1,1,3,3)
         Amat(2,1) = SpUUMod(2,1,1,1)
         Amat(2,2) = SpUUMod(2,1,2,1)
         Amat(2,3) = SpUUMod(2,1,3,1)
         Amat(2,4) = SpUUMod(2,1,1,2)
         Amat(2,5) = SpUUMod(2,1,2,2)
         Amat(2,6) = SpUUMod(2,1,3,2)
         Amat(2,7) = SpUUMod(2,1,1,3)
         Amat(2,8) = SpUUMod(2,1,2,3)
         Amat(2,9) = SpUUMod(2,1,3,3)
         Amat(3,1) = SpUUMod(3,1,1,1)
         Amat(3,2) = SpUUMod(3,1,2,1)
         Amat(3,3) = SpUUMod(3,1,3,1)
         Amat(3,4) = SpUUMod(3,1,1,2)
         Amat(3,5) = SpUUMod(3,1,2,2)
         Amat(3,6) = SpUUMod(3,1,3,2)
         Amat(3,7) = SpUUMod(3,1,1,3)
         Amat(3,8) = SpUUMod(3,1,2,3)
         Amat(3,9) = SpUUMod(3,1,3,3)
         Amat(4,1) = SpUUMod(1,2,1,1)
         Amat(4,2) = SpUUMod(1,2,2,1)
         Amat(4,3) = SpUUMod(1,2,3,1)
         Amat(4,4) = SpUUMod(1,2,1,2)
         Amat(4,5) = SpUUMod(1,2,2,2)
         Amat(4,6) = SpUUMod(1,2,3,2)
         Amat(4,7) = SpUUMod(1,2,1,3)
         Amat(4,8) = SpUUMod(1,2,2,3)
         Amat(4,9) = SpUUMod(1,2,3,3)
         Amat(5,1) = SpUUMod(2,2,1,1)
         Amat(5,2) = SpUUMod(2,2,2,1)
         Amat(5,3) = SpUUMod(2,2,3,1)
         Amat(5,4) = SpUUMod(2,2,1,2)
         Amat(5,5) = SpUUMod(2,2,2,2)
         Amat(5,6) = SpUUMod(2,2,3,2)
         Amat(5,7) = SpUUMod(2,2,1,3)
         Amat(5,8) = SpUUMod(2,2,2,3)
         Amat(5,9) = SpUUMod(2,2,3,3)
         Amat(6,1) = SpUUMod(3,2,1,1)
         Amat(6,2) = SpUUMod(3,2,2,1)
         Amat(6,3) = SpUUMod(3,2,3,1)
         Amat(6,4) = SpUUMod(3,2,1,2)
         Amat(6,5) = SpUUMod(3,2,2,2)
         Amat(6,6) = SpUUMod(3,2,3,2)
         Amat(6,7) = SpUUMod(3,2,1,3)
         Amat(6,8) = SpUUMod(3,2,2,3)
         Amat(6,9) = SpUUMod(3,2,3,3)
         Amat(7,1) = SpUUMod(1,3,1,1)
         Amat(7,2) = SpUUMod(1,3,2,1)
         Amat(7,3) = SpUUMod(1,3,3,1)
         Amat(7,4) = SpUUMod(1,3,1,2)
         Amat(7,5) = SpUUMod(1,3,2,2)
         Amat(7,6) = SpUUMod(1,3,3,2)
         Amat(7,7) = SpUUMod(1,3,1,3)
         Amat(7,8) = SpUUMod(1,3,2,3)
         Amat(7,9) = SpUUMod(1,3,3,3)
         Amat(8,1) = SpUUMod(2,3,1,1)
         Amat(8,2) = SpUUMod(2,3,2,1)
         Amat(8,3) = SpUUMod(2,3,3,1)
         Amat(8,4) = SpUUMod(2,3,1,2)
         Amat(8,5) = SpUUMod(2,3,2,2)
         Amat(8,6) = SpUUMod(2,3,3,2)
         Amat(8,7) = SpUUMod(2,3,1,3)
         Amat(8,8) = SpUUMod(2,3,2,3)
         Amat(8,9) = SpUUMod(2,3,3,3)
         Amat(9,1) = SpUUMod(3,3,1,1)
         Amat(9,2) = SpUUMod(3,3,2,1)
         Amat(9,3) = SpUUMod(3,3,3,1)
         Amat(9,4) = SpUUMod(3,3,1,2)
         Amat(9,5) = SpUUMod(3,3,2,2)
         Amat(9,6) = SpUUMod(3,3,3,2)
         Amat(9,7) = SpUUMod(3,3,1,3)
         Amat(9,8) = SpUUMod(3,3,2,3)
         Amat(9,9) = SpUUMod(3,3,3,3)


         Qmat = zero
         Qmat(1,1) = third*(Amat(1,1)+Amat(1,5)+Amat(1,9)) 
     +        - (two/three)*T_tau(1,1)
         Qmat(2,1) = third*(Amat(2,1)+Amat(2,5)+Amat(2,9))
     +        - (two/three)*T_tau(2,1)
         Qmat(3,1) = third*(Amat(3,1)+Amat(3,5)+Amat(3,9))
     +        - (two/three)*T_tau(3,1)
         Qmat(4,1) = third*(Amat(4,1)+Amat(4,5)+Amat(4,9))
     +        - (two/three)*T_tau(1,2)
         Qmat(5,1) = third*(Amat(5,1)+Amat(5,5)+Amat(5,9))
     +        - (two/three)*T_tau(2,2)
         Qmat(6,1) = third*(Amat(6,1)+Amat(6,5)+Amat(6,9))
     +        - (two/three)*T_tau(3,2)
         Qmat(7,1) = third*(Amat(7,1)+Amat(7,5)+Amat(7,9))
     +        - (two/three)*T_tau(1,3)
         Qmat(8,1) = third*(Amat(8,1)+Amat(8,5)+Amat(8,9))
     +        - (two/three)*T_tau(2,3)
         Qmat(9,1) = third*(Amat(9,1)+Amat(9,5)+Amat(9,9))
     +        - (two/three)*T_tau(3,3)
         Qmat(1,5) = Qmat(1,1)
         Qmat(2,5) = Qmat(2,1)
         Qmat(3,5) = Qmat(3,1)
         Qmat(4,5) = Qmat(4,1)
         Qmat(5,5) = Qmat(5,1)
         Qmat(6,5) = Qmat(6,1)
         Qmat(7,5) = Qmat(7,1)
         Qmat(8,5) = Qmat(8,1)
         Qmat(9,5) = Qmat(9,1)
         Qmat(1,9) = Qmat(1,1)
         Qmat(2,9) = Qmat(2,1)
         Qmat(3,9) = Qmat(3,1)
         Qmat(4,9) = Qmat(4,1)
         Qmat(5,9) = Qmat(5,1)
         Qmat(6,9) = Qmat(6,1)
         Qmat(7,9) = Qmat(7,1)
         Qmat(8,9) = Qmat(8,1)
         Qmat(9,9) = Qmat(9,1)
         

         if((nNode.eq.8).and.(nInt.eq.8)) then
            !
            ! This is the tangent using the F-bar method with the
            !  8 node fully integrated linear element
            !
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           + matmul(transpose(Gmat),matmul(Qmat,(G0mat-Gmat)))
     +           )
         else
            !
            ! This is the tangent NOT using the F-bar method with all
            !  other elements
            !
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           )
         endif



         ! Compute/update the electric potential tangent matrix
         !
         Kphiphi = Kphiphi - detmapJC*w(intPt)*
     +        (
     +        matmul(matmul(dshC,SpPhiPhiMod),transpose(dshC))
     +        )


         ! Compute/update the electric potential/displacement tangent matrix
         !   Strictly speaking, use of 'F-bar' will affect this tangent;
         !   however, the effect is expected to be small and is neglected here.
         !
         AmatPhiU = zero
         AmatPhiU(1,1) = SpPhiUMod(1,1,1)
         AmatPhiU(1,2) = SpPhiUMod(1,2,1)
         AmatPhiU(1,3) = SpPhiUMod(1,3,1)
         AmatPhiU(1,4) = SpPhiUMod(1,1,2)
         AmatPhiU(1,5) = SpPhiUMod(1,2,2)
         AmatPhiU(1,6) = SpPhiUMod(1,3,2)
         AmatPhiU(1,7) = SpPhiUMod(1,1,3)
         AmatPhiU(1,8) = SpPhiUMod(1,2,3)
         AmatPhiU(1,9) = SpPhiUMod(1,3,3)
         AmatPhiU(2,1) = SpPhiUMod(2,1,1)
         AmatPhiU(2,2) = SpPhiUMod(2,2,1)
         AmatPhiU(2,3) = SpPhiUMod(2,3,1)
         AmatPhiU(2,4) = SpPhiUMod(2,1,2)
         AmatPhiU(2,5) = SpPhiUMod(2,2,2)
         AmatPhiU(2,6) = SpPhiUMod(2,3,2)
         AmatPhiU(2,7) = SpPhiUMod(2,1,3)
         AmatPhiU(2,8) = SpPhiUMod(2,2,3)
         AmatPhiU(2,9) = SpPhiUMod(2,3,3)
         AmatPhiU(3,1) = SpPhiUMod(3,1,1)
         AmatPhiU(3,2) = SpPhiUMod(3,2,1)
         AmatPhiU(3,3) = SpPhiUMod(3,3,1)
         AmatPhiU(3,4) = SpPhiUMod(3,1,2)
         AmatPhiU(3,5) = SpPhiUMod(3,2,2)
         AmatPhiU(3,6) = SpPhiUMod(3,3,2)
         AmatPhiU(3,7) = SpPhiUMod(3,1,3)
         AmatPhiU(3,8) = SpPhiUMod(3,2,3)
         AmatPhiU(3,9) = SpPhiUMod(3,3,3)
         !
         Kphiu = Kphiu - detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(dshC,AmatPhiU),Gmat)
     +        )




         ! Compute/update the displacement/electric potential tangent matrix
         !
         AmatUPhi = zero
         AmatUPhi(1,1) = SpUPhiMod(1,1,1)
         AmatUPhi(2,1) = SpUPhiMod(2,1,1)
         AmatUPhi(3,1) = SpUPhiMod(3,1,1)
         AmatUPhi(4,1) = SpUPhiMod(1,2,1)
         AmatUPhi(5,1) = SpUPhiMod(2,2,1)
         AmatUPhi(6,1) = SpUPhiMod(3,2,1)
         AmatUPhi(7,1) = SpUPhiMod(1,3,1)
         AmatUPhi(8,1) = SpUPhiMod(2,3,1)
         AmatUPhi(9,1) = SpUPhiMod(3,3,1)
         AmatUPhi(1,2) = SpUPhiMod(1,1,2)
         AmatUPhi(2,2) = SpUPhiMod(2,1,2)
         AmatUPhi(3,2) = SpUPhiMod(3,1,2)
         AmatUPhi(4,2) = SpUPhiMod(1,2,2)
         AmatUPhi(5,2) = SpUPhiMod(2,2,2)
         AmatUPhi(6,2) = SpUPhiMod(3,2,2)
         AmatUPhi(7,2) = SpUPhiMod(1,3,2)
         AmatUPhi(8,2) = SpUPhiMod(2,3,2)
         AmatUPhi(9,2) = SpUPhiMod(3,3,2)
         AmatUPhi(1,3) = SpUPhiMod(1,1,3)
         AmatUPhi(2,3) = SpUPhiMod(2,1,3)
         AmatUPhi(3,3) = SpUPhiMod(3,1,3)
         AmatUPhi(4,3) = SpUPhiMod(1,2,3)
         AmatUPhi(5,3) = SpUPhiMod(2,2,3)
         AmatUPhi(6,3) = SpUPhiMod(3,2,3)
         AmatUPhi(7,3) = SpUPhiMod(1,3,3)
         AmatUPhi(8,3) = SpUPhiMod(2,3,3)
         AmatUPhi(9,3) = SpUPhiMod(3,3,3)
         !
         Kuphi = Kuphi - detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(transpose(Gmat),AmatUPhi),transpose(dshC))
     +        )


      enddo
      !
      ! End the loop over integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface terms
      !
      if(ndload.gt.0) then
         !
         write(*,*) '3D surface terms not programmed here'
         write(80,*) '3D surface terms not programmed here'
         call xit
         !
      endif ! ndload.gt.0 or not
      !
      ! End loop over surface  terms
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Rphi,Kuu,Kuphi,Kphiu,Kphiphi,
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------

      return
      end subroutine U3D8

************************************************************************

      subroutine integ(props,nprops,dtime,nVisco,
     +        F_tau,ER,A_t,
     +        T_tau,D_tau,A_tau,deltaA,
     +        SpUUMod,SpUPhiMod,SpPhiUMod,SpPhiPhiMod,
     +        stat)
      
      implicit none

      integer i,j,k,l,m,n,stat,nprops,nVisco

      real*8 Iden(3,3),F_tau(3,3),FT(3,3),FinvT(3,3),C_tau(3,3),detC
      real*8 Cinv(3,3),trC,I1bar,detF,Finv(3,3),ER(3,1),E(3,1),vec(3,1)
      real*8 I6,TR_tau(3,3),T_tau(3,3),D_tau(3,1),dTRdF(3,3,3,3),Gshear
      real*8 SpUUmod(3,3,3,3),A_t(nVisco,3,3),A_tau(nVisco,3,3),Kbulk,G0
      real*8 SpUPhiMod(3,3,3),SpPhiUMod(3,3,3),SpPhiPhiMod(3,3),permit
      real*8 props(nprops),Gneq(nVisco),tau(nVisco),dtime,AC(3,3),trAC
      real*8 Tneq(nVisco,3,3),Fbar(3,3),AFt(3,3),trFAFt,FAFt(3,3),deltaA
      real*8 FA(3,3),lamL,dGdF(3,3),aux,lamBar

      real*8 zero,one,two,three,third,half,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0,nine=9.d0)


      ! Identity matrix
      !
      call onem(Iden)


      ! Obtain material parameters
      !
      G0     = props(1)
      lamL   = props(2)
      Kbulk  = props(3)
      permit = props(4)
      do n=1,nVisco
         Gneq(n) = props((n-1)*2+5)
         tau(n)  = props((n-1)*2+6)
      enddo


      ! Compute the determinant, the inverse, the transpose, 
      !  and the inverse transpose of the deformation gradient
      !
      call matInv3D(F_tau,Finv,detF,stat)
      FT = transpose(F_tau)
      FinvT = transpose(Finv)
      Fbar = (detF**(-two/three))*F_tau


      ! Compute the right Cauchy-Green tensor, its inverse, 
      !  its determinant, and its trace
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      call matInv3D(C_tau,Cinv,detC,stat)
      trC = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)
 

      ! Compute the trace of the distortional right Cauchy-Green tensor
      !
      I1bar = (detF**(-two/three))*trC
      
      
      ! Compute the current electric field
      !
      E = matmul(FinvT,ER)
      
      
      ! Compute quantities related to electric field
      !
      vec = matmul(Cinv,ER)
      I6 = E(1,1)*E(1,1) + E(2,1)*E(2,1) + E(3,1)*E(3,1)


      ! Compute the shear modulus
      !
      if(lamL.le.zero) then
         Gshear = G0
      else
         lamBar = dsqrt(third*I1bar)
         aux = (three - (lamBar/lamL)**two)/(one - (lamBar/lamL)**two)
         Gshear = third*G0*aux
      endif


      ! Compute the derivative of the shear modulus with
      !  respect to the deformation gradient for later use
      !
      if(lamL.le.zero) then
         dGdF = zero
      else
         dGdF = ((4.d0*G0*lamL*lamL)/
     +        (nine*(lamBar*lamBar - lamL*lamL)**two))*F_tau
      endif
 

      ! Compute the equilibrium 1st Piola stress
      !
      TR_tau = (detF**(-two/three))*Gshear*(F_tau - third*trC*FinvT)
     +     + Kbulk*detF*(detF - one)*FinvT
     +     + permit*detF*(matmul(E,transpose(matmul(Cinv,ER))) 
     +                    - half*I6*FinvT)
 

      ! Compute the current electric displacement
      !
      D_tau = permit*E


      ! Update state variables
      !
      deltaA = zero
      do n=1,nVisco
         A_tau(n,:,:) = (A_t(n,:,:)
     +        + (dtime/tau(n))*(detF**(two/three))*Cinv)
     +        /(one + dtime/tau(n))
         deltaA = max(deltaA,dsqrt(sum((A_tau(n,:,:)-A_t(n,:,:))*
     +        (A_tau(n,:,:)-A_t(n,:,:)))))
      enddo


      ! Compute the factor A*Cbar, its trace, and the
      !  nonequilibrium stresses
      !
      do n=1,nVisco
         AC = (detF**(-two/three))*matmul(A_tau(n,:,:),C_tau)
         trAC = AC(1,1) + AC(2,2) + AC(3,3)
         Tneq(n,:,:) = (Gneq(n)/detF)*(matmul(matmul(Fbar,A_tau(n,:,:)),
     +        transpose(Fbar)) - third*trAC*Iden)
      enddo



      ! Compute the total Cauchy stress
      !
      T_tau = (one/detF)*matmul(TR_tau,transpose(F_tau))
      do n=1,nVisco
         T_tau = T_tau + Tneq(n,:,:)
      enddo



      ! Calculate the equilibrium material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + (detF**(-two/three))*Gshear*
     +                 (
     +                 (-two/three)*F_tau(i,j)*Finv(l,k)
     +                 + (two/9.d0)*trC*Finv(j,i)*Finv(l,k)
     +                 + Iden(i,k)*Iden(j,l)
     +                 + third*trC*Finv(l,i)*Finv(j,k)
     +                 - (two/three)*Finv(j,i)*F_tau(k,l)
     +                 )
     +                 + dGdF(k,l)*
     +                 (
     +                 (detF**(-two/three))*
     +                 (F_tau(i,j) - third*trC*Finv(j,i))
     +                 )
     +                 + detF*Kbulk*
     +                 (
     +                 (detF-one)*Finv(j,i)*Finv(l,k)
     +                 + detF*Finv(j,i)*Finv(l,k)
     +                 - (detF-one)*Finv(l,i)*Finv(j,k)
     +                 )
     +                 + detF*permit*
     +                 (
     +                 - Finv(l,i)*vec(j,1)*E(k,1) 
     +                 - E(i,1)*Finv(j,k)*vec(l,1)
     +                 - E(i,1)*E(k,1)*Cinv(j,l)
     +                 + Finv(j,i)*E(k,1)*vec(l,1)
     +                 + E(i,1)*vec(j,1)*Finv(l,k)
     +                 + half*I6*Finv(l,i)*Finv(j,k)
     +                 - half*I6*Finv(j,i)*Finv(l,k)
     +                 )
               enddo
            enddo
         enddo
      enddo


      ! Calculate the equilibrium spatial tangent modulus
      !
      SpUUMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpUUMod(i,j,k,l) = SpUUMod(i,j,k,l) +
     +                      (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo




      ! Modify the spatial tangent modulus to include the viscous terms
      !
      dTRdF = zero
      do n=1,nVisco
         !
         ! Compute some terms used in the tangent modulus
         !
         FAFt = matmul(matmul(F_tau,A_t(n,:,:)),transpose(F_tau))
         trFAFt = FAFt(1,1) + FAFt(2,2) + FAFt(3,3)
         FA = matmul(F_tau,A_t(n,:,:))
         AFt = matmul(A_t(n,:,:),transpose(F_tau))
         !
         ! Compute dTRdF for the viscous mechanisms
         !
         do i=1,3
            do j=1,3
               do k=1,3
                  do l=1,3
                     dTRdF(i,j,k,l) = dTRdF(i,j,k,l) +
     +                    detF**(-two/three)*
     +                    (
     +                    Gneq(n)/(one + (dtime/tau(n)))*
     +                    (
     +                    -(two/three)*Finv(l,k)*FA(i,j)
     +                    + (two/nine)*Finv(l,k)*trFAFt*Finv(j,i)
     +                    + Iden(i,k)*A_t(n,l,j)
     +                    - (two/three)*FA(k,l)*Finv(j,i)
     +                    + third*trFAFt*Finv(l,i)*Finv(j,k)
     +                    )
     +                    )
                  enddo
               enddo
            enddo
         enddo
         !
      enddo ! loop over n
      !
      ! Update the spatial tangent modulus to include the
      !  viscous mechansims
      !
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpUUMod(i,j,k,l) = SpUUMod(i,j,k,l) +
     +                       (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo



      ! Calculate the spatial stress/electric potential modulus modulus
      !
      SpUPhiMod = zero
      do i=1,3
         do j=1,3
            do l=1,3
               SpUPhiMod(i,j,l) = SpUPhiMod(i,j,l) +
     +              permit*(Iden(i,l)*E(j,1) + 
     +                       E(i,1)*Iden(j,l) - 
     +                       Iden(i,j)*E(l,1))
            end do
         end do
      end do


      ! Calculate the spatial electric displacement/electric potential modulus
      !
      SpPhiPhiMod = permit*Iden


      ! Calculate the spatial electric displacement/strain modulus
      !
      SpPhiUMod = zero
      do j=1,3
         do k=1,3
            do l=1,3
               SpPhiUMod(j,k,l) = SpPhiUMod(j,k,l) -
     +              permit*(Iden(j,k)*E(l,1)
     +              + E(k,1)*Iden(j,l)
     +              - E(j,1)*Iden(l,k))
            enddo
         enddo
      enddo



      return
      end subroutine integ

************************************************************************

      subroutine AssembleElement(nDim,nNode,ndofel,
     +     Ru,Rphi,Kuu,Kuphi,Kphiu,Kphiphi,
     +     rhs,amatrx)
      
      !
      ! Subroutine to assemble the local elements residual and tangent
      !

      implicit none

      integer i,j,k,l,m,n,A11,A12,B11,B12,nDim,nNode,nDofEl,nDofN

      real*8 Ru(nDim*nNode,1),Rphi(nNode,1),Kuu(nDim*nNode,nDim*nNode)
      real*8 Kphiphi(nNode,nNode),Kuphi(nDim*nNode,nNode),rhs(ndofel,1)
      real*8 Kphiu(nNode,nDim*nNode),amatrx(ndofel,ndofel)


      ! Total number of degrees of freedom per node
      !
      nDofN = nDofEl/nNode


      ! init
      !
      rhs(:,1) = 0.d0
      amatrx = 0.d0

      if(nDim.eq.2) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1) = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            !
            ! electric potential
            !
            rhs(A11+2,1) = Rphi(i,1)
         enddo
         !
         ! Assemble the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11) = Kuu(A12,B12)
               amatrx(A11,B11+1) = Kuu(A12,B12+1)
               amatrx(A11+1,B11) = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               !
               ! electric potential
               !
               amatrx(A11+2,B11+2) = Kphiphi(i,j)
               !
               ! displacement - electric potential
               !
               amatrx(A11,B11+2) = Kuphi(A12,j)
               amatrx(A11+1,B11+2) = Kuphi(A12+1,j)
               !
               ! electric potential - displacement
               !
               amatrx(A11+2,B11) = Kphiu(i,B12)
               amatrx(A11+2,B11+1) = Kphiu(i,B12+1)
               !
            enddo
         enddo
         !
      elseif(nDim.eq.3) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1)   = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            rhs(A11+2,1) = Ru(A12+2,1)
            !
            ! electric potential
            !
            rhs(A11+3,1) = Rphi(i,1)
            !
         enddo
         !
         ! Assembly the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11)     = Kuu(A12,B12)
               amatrx(A11,B11+1)   = Kuu(A12,B12+1)
               amatrx(A11,B11+2)   = Kuu(A12,B12+2)
               amatrx(A11+1,B11)   = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
               amatrx(A11+2,B11)   = Kuu(A12+2,B12)
               amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
               amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
               !
               ! electric potential
               !
               amatrx(A11+3,B11+3) = Kphiphi(i,j)
               !
               ! displacement - electric potential
               !
               amatrx(A11,B11+3) = Kuphi(A12,j)
               amatrx(A11+1,B11+3) = Kuphi(A12+1,j)
               amatrx(A11+2,B11+3) = Kuphi(A12+2,j)
               !
               ! electric potential - displacement
               !
               amatrx(A11+3,B11) = Kphiu(i,B12)
               amatrx(A11+3,B11+1) = Kphiu(i,B12+1)
               amatrx(A11+3,B11+2) = Kphiu(i,B12+2)
               !
            enddo
         enddo
         !
      else
         write(*,*) 'How did you get nDim=',nDim
         call xit
      endif

      return
      end subroutine AssembleElement

!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint2D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2), w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0


      return
      end subroutine xint2D1pt
      
!************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt

************************************************************************

      subroutine xint3D1pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 2 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(1,3),w(1)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 8.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0

      return
      end subroutine xint3D1pt
     
************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(8,3),w(8)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element


      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
      end subroutine calcShape2DLinear

************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 8-node linear 3D element as shown
      !
      !      8-----------7
      !     /|          /|       zeta
      !    / |         / |       
      !   5-----------6  |       |     eta
      !   |  |        |  |       |   /
      !   |  |        |  |       |  /
      !   |  4--------|--3       | /
      !   | /         | /        |/
      !   |/          |/         O--------- xi
      !   1-----------2        origin at cube center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt,i,j

      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3)
      real*8 d2shxi(8,3,3),xi,eta,zeta

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      !
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      !
      ! The second derivatives
      !
      d2shxi = zero
      d2shxi(1,1,2) = eighth*(one - zeta)
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(1,1,3) = eighth*(one - eta)
      d2shxi(1,3,1) = d2shxi(1,1,3)
      d2shxi(1,2,3) = eighth*(one - xi)
      d2shxi(1,3,2) = d2shxi(1,2,3)
      d2shxi(2,1,2) = -eighth*(one - zeta)
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(2,1,3) = -eighth*(one - eta)
      d2shxi(2,3,1) = d2shxi(2,1,3)
      d2shxi(2,2,3) = eighth*(one + xi)
      d2shxi(2,3,2) = d2shxi(2,2,3)
      d2shxi(3,1,2) = eighth*(one - zeta)
      d2shxi(3,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,3) = -eighth*(one + eta)
      d2shxi(3,3,1) = d2shxi(2,1,3)
      d2shxi(3,2,3) = -eighth*(one + xi)
      d2shxi(3,3,2) = d2shxi(2,2,3)
      d2shxi(4,1,2) = -eighth*(one - zeta)
      d2shxi(4,2,1) = d2shxi(2,1,2)
      d2shxi(4,1,3) = eighth*(one + eta)
      d2shxi(4,3,1) = d2shxi(2,1,3)
      d2shxi(4,2,3) = -eighth*(one - xi)
      d2shxi(4,3,2) = d2shxi(2,2,3)
      d2shxi(5,1,2) = eighth*(one + zeta)
      d2shxi(5,2,1) = d2shxi(2,1,2)
      d2shxi(5,1,3) = -eighth*(one - eta)
      d2shxi(5,3,1) = d2shxi(2,1,3)
      d2shxi(5,2,3) = -eighth*(one - xi)
      d2shxi(5,3,2) = d2shxi(2,2,3)
      d2shxi(6,1,2) = eighth*(one + zeta)
      d2shxi(6,2,1) = d2shxi(2,1,2)
      d2shxi(6,1,3) = eighth*(one - eta)
      d2shxi(6,3,1) = d2shxi(2,1,3)
      d2shxi(6,2,3) = -eighth*(one + xi)
      d2shxi(6,3,2) = d2shxi(2,2,3)
      d2shxi(7,1,2) = eighth*(one + zeta)
      d2shxi(7,2,1) = d2shxi(2,1,2)
      d2shxi(7,1,3) = eighth*(one + eta)
      d2shxi(7,3,1) = d2shxi(2,1,3)
      d2shxi(7,2,3) = eighth*(one + xi)
      d2shxi(7,3,2) = d2shxi(2,2,3)
      d2shxi(8,1,2) = -eighth*(one + zeta)
      d2shxi(8,2,1) = d2shxi(2,1,2)
      d2shxi(8,1,3) = -eighth*(one + eta)
      d2shxi(8,3,1) = d2shxi(2,1,3)
      d2shxi(8,2,3) = eighth*(one - xi)
      d2shxi(8,3,2) = d2shxi(2,2,3)
      
      return
      end subroutine calcShape3DLinear

************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da

************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 8-node
      !  linear and 20-node quadratic 3D elements.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode)
      real*8 mapJ(3,3),mapJ_inv(3,3),detmapJ

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      ! The second derivatives may be calculated.
      !

      return
      end subroutine mapShape3D

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

****************************************************************************
