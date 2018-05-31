*Deck Explor
      Subroutine Explor(V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Copyright (c) 1988,1990,1992,1993,1995,1998,2003,2009,2016
C                Gaussian, Inc.  All Rights Reserved.
C
C     This is part of the Gaussian(R) 16 program.  It is based on
C     the Gaussian(R) 09 system (copyright 2009, Gaussian, Inc.),
C     the Gaussian(R) 03 system (copyright 2003, Gaussian, Inc.),
C     the Gaussian(R) 98 system (copyright 1998, Gaussian, Inc.),
C     the Gaussian(R) 94 system (copyright 1995, Gaussian, Inc.),
C     the Gaussian 92(TM) system (copyright 1992, Gaussian, Inc.),
C     the Gaussian 90(TM) system (copyright 1990, Gaussian, Inc.),
C     the Gaussian 88(TM) system (copyright 1988, Gaussian, Inc.),
C     the Gaussian 86(TM) system (copyright 1986, Carnegie Mellon
C     University), and the Gaussian 82(TM) system (copyright 1983,
C     Carnegie Mellon University). Gaussian is a federally registered
C     trademark of Gaussian, Inc.
C
C     This software contains proprietary and confidential information,
C     including trade secrets, belonging to Gaussian, Inc.
C
C     This software is provided under written license and may be
C     used, copied, transmitted, or stored only in accord with that
C     written license.
C
C     The following legend is applicable only to US Government contracts
C     under DFARS:
C
C                        RESTRICTED RIGHTS LEGEND
C
C     Use, duplication or disclosure by the US Government is subject to
C     restrictions as set forth in subparagraph (c)(1)(ii) of the Rights
C     in Technical Data and Computer Software clause at DFARS
C     252.227-7013.
C
C     Gaussian, Inc., 340 Quinnipiac St., Bldg. 40, Wallingford CT 06492
C
C     The following legend is applicable only to US Government contracts
C     under FAR:
C
C                        RESTRICTED RIGHTS LEGEND
C
C     Use, reproduction and disclosure by the US Government is subject
C     to restrictions as set forth in subparagraph (c) of the Commercial
C     Computer Software - Restricted Rights clause at FAR 52.227-19.
C
C     Gaussian, Inc., 340 Quinnipiac St., Bldg. 40, Wallingford CT 06492
C
C
C                       INTRODUCTION TO THIS LINK
C1l123
C1Explor
C1IRC
C     This link performs reaction path following calculations for Gaussian.
C     This link is intended to be used as an alternative to Link 115.
C     Among the algorithms contained herein are (see references throughout
C     the link):
C                1. Euler Integration;
C                2. Local Quadradic Approximation (LQA);
C                3. Hessian-Based Predictor-Corrector (HPC);
C                4. Second-Order Gonzales-Schlegel (GS2);
C                5. Damped-Velocity Verlet (DVV); and
C                6. Euler-Based Predictor-Corrector (EPC).
C
C     The default integrator is the HPC method of Hratchian and Schlegel.
C     An overview of reaction path methods, as well as optimization
C     techniques, is available in:
C           Hratchian and Schlegel, in Theory and Applications of
C           Computational Chemistry: The First 40 Years, Ed. C. E. Dykstra,
C           G. Frenking, K. S. Kim, and G. Scuseria (Elsevier, Amsterdam,
C           2005) 195- 249.
C     Papers specifically related to schemes available in this link
C     include:
C
C
C     Link 123 also contains code to carry out coordinate driving along any
C     chosen Hessian eigenvector. At present, the step direction is
C     determined from the Hessian eigenvectors at the initial geometry
C     only. By default, coordinate driving calculations are carried out in
C     a non-mass-weighted Cartesian coordinate space.
C
C     By default reaction pathways are carried out in mass-weighted
C     Cartesian coordinates, which corresponds to the Fukui intrinsic
C     reaction coordinate (IRC). To mass-weigh the coordinate system, the
C     code uses the built-in atomic masses. The user may also specify
C     particular isotope values for each of the atoms via the typical
C     methods.
C
C     IOps are documented in L101.
C?
C
      Integer MaxIOp
      Parameter (MaxIOp=DEFMAXIOP,MxPath=2,LURVAFl=1024,IRwGen=501,
     $  IRWExp=575,IOZRed=665,IRwPES=781,IRwFit=784,IFtSt=9,IFtMax=12,
     $  NCorOp=16,LRwGen=1000)
      Common /MUnit/IUnit(20)
      Common /IOp/    IOp(MaxIOp)
      Common /IO/     In,IOut,IPunch
#include "commonmol.inc"
C
C     Local Rxn Path data for L123 is kept in common/RxnDat/ (RwF 575).
C     This common block is used in Routines Explore and BkAjst. The
C     contents of this common block include:
C           StpSiz      The current integration step size.
C           DVVv0       The damping constant for Damped-Velocity-Verlet
C                       jobs.
C           CritGd      A stopping criteria flag for gradient-based tests.
C           CritAn      A stopping criteria flag for step-angle-based tests.
C           TsPot       The potential energy of the TS.
C           TsPotQ      The QM component of the potential energy of the TS.
C                       This is used only in QM:MM runs using the
C                       QuadMac-based integration machinery.
C           Potent      This is a 2x2 array that holds the potential energy
C                       of the last two points on the two possible paths.
C           PotSav      This is the saved energy from the previous instance
C                       of this link.
C           EQM         This is a 2x2 array that holds the QM component of
C                       the potential energy of the last two points on the
C                       two possible paths. It is only used in QM:MM runs
C                       using the QuadMac-based integration machinery.
C           EPS         The requested epsilon for the modified-BS
C                       integrator.
C           IMthd       The integer integration method flag.
C           NRdStat     The number of read-in statistical parameters that
C                       need to be reported at the end of the job.
C           ICurWh      This is an integer flag that indicates what Routine
C                       Explor and RxnCyc are to do with the current point
C                       data (arrays such as CCur, GrdCur, HesCur, HVcCur,
C                       and HVlCur). In the first half of Routine RxnCyc,
C                       ICurWh is used to determine what to do with the
C                       data in these arrays after they've been populated
C                       from the RWF file. ICurWh can then be set as needed
C                       by corrector integration routines below RxnCyc. And
C                       then, at the end of Routine RxnCyc, ICurWh is used
C                       to decide how to fill CCur which Routine Explor
C                       uses as the geometry at which the next set of PES
C                       calculations are done by the rest of Gaussian.
C                       ICurWh is initialized to 0 in Routine Explor at the
C                       start of the job and each time IStep2 is re-set.
C                       ICurWh=0 indicates that the current PES data is
C                       part of the actual exploration and should be put
C                       into the standard arrays; ICurWh=1 indicates that
C                       RxnCyc should NOT copy the current point arrays
C                       into the standard data set.
C           NAtmQX      The number of QM atoms in a QM:MM runs using the
C                       QuadMac-based integration machinery.
C           IRxFwd      The integer flag indicating whether the integrated
C                       path follows the forward, reverse, or downhill
C                       directions.
C           MaxStp      The maximum number of integration steps allowed per
C                       path.
C           IStep       The current integration step number.
C           MxOpt       The maximum number of optimization cycles allowed
C                       for a given step in a GS2 reaction path
C                       caclulations, or the maximum number of corrector
C                       convergence cycles in EulerPC and HPC reaction path
C                       calculations.
C           IStep2      The current cycle step number for integration step
C                       IStep. This is the counter for the iterative
C                       solvers.
C           NPath       The total number of requested paths.
C           IPath       The current path number.
C           IHsStp      This is a counter that gets re-set to 0 before
C                       exiting Link 123 whenever an analytic Hessian is
C                       requested. As a result, it is also used at the top
C                       of the link when a Hessian-based integrator is
C                       being used to determine whether an analytic Hessian
C                       is on the RWF or an updated Hessian needs to be
C                       evaluated at the current point.
C           IPG         Integer flags related to the molecular point group.
C           ICorFl      An integer switch used by PreDWI with HPC and
C                       EulerPC integration schemes for keeping track of
C                       special versions of re-correction.
C           IURVAFl     For jobs requesting an URVA data file, this array
C                       holds the Holerith string of the URVA filename. The
C                       first alement of the array gives the number of
C                       characters in the name and the remainder of the
C                       array holds the string.
C           NTotGd      A counter that keeps track of the total number of
C                       analytic gradient evaluations.
C           NTotHs      A counter that keeps track of the total number of
C                       analytic Hessian evaluations.
C
      Common/RxnDat/StpSiz,DVVv0,CritGd,CritAn,TsPot,TsPotQ,Potent(2,2),
     $  PotSav,EQM(2,2),EPS,IMthd,NRdStat,ICurWh,NAtmQX,IRxFwd,MaxStp,
     $  IStep,MxOpt,IStep2,NPath,IPath,IHsStp,IPG(4),ICorFl,
     $  IURVAFl(LURVAFl),NTotGd,NTotHs,IDmRxD
      Logical DEBUG,First,Restrt,HvHess,DoHEVc,Update,DoHess,SkpEGH,
     $  Compct,IntHes,DwnHil,CorDon,SumRlE,SumCnv,SvPES,WrURVA,DoSvC,
     $  DoNRTst,OK
      Dimension MxPnt(MxPath),NPhase(4),ICorOp(NCorOp),Gen(LRwGen),
     $  XX(1),JJ(1),V(MDV)
      Save XX,JJ,Zero,One
      Data XX/0.d0/,JJ/0/,Zero/0.d0/,One/1.0d0/
 1000 Format(1x,17('IRC-'),('IRC'))
 1040 Format(' ')
 2000 Format(1x,'RESTARTING THE IRC JOB')
 5000 Format(1x,'Total number of points:               ',I6,
     $  /,1x,'Total number of gradient calculations:',I6,
     $  /,1x,'Total number of Hessian calculations: ',I6)
 8000 Format(1x,'IStep ',I5,' IStep2 ',I5,
     $  ' Update activated in Explor.')
 8001 Format(1x,'IStep ',I5,' IStep2 ',I5,
     $  ' Update de-activated in Explor.')
 8010 Format(1x,'IStep ',I5,' IStep2 ',I5,
     $  ' Link 123 is requesting second derivatives.')
 8011 Format(1x,'IStep ',I5,' IStep2 ',I5,
     $  ' Link 123 is turning off second derivative evalution.')
 9000 Format(1x,'NAtomQ = ',I5,2x,'NAtmQ1 = ',I5,2x,'NAtmQX = ',I5)
 9100 Format(1x,'Path direction must be the same in the original and',
     $  ' restart jobs.')
 9900 Format(1x,'Explor: JobFlg=',I3)
C
C     Initialize the link and some key flags.
C
      Call Drum(V,MDV)
      Call TStamp(1,'Explor: Top')
      IChain = 0
      LenCom = 16+InToWP(20+LURVAFl)
      First = Mod(IOp(38),10).eq.1
      IOP10 = IOp(10)
      IPrint = IOp(33)
      Call ILSW(2,21,IPFlag)
      If(IPFlag.eq.0.and.IPrint.eq.0) then
        IPrint = Max(IPrint,1)
      else if(IPFlag.eq.2.and.IPrint.eq.0) then
        IPrint = -1
        endIf
      DEBUG = .False.
      Restrt = IOp(35).eq.1
      IOP36  = IOP(36)
      IRdInp = IOp(129)
      IHesRst= IOp(143)
      DoHess = .False.
      SkpEGH = .False.
      DoSvC = IOp(161).eq.1
CC-INSERT     
C     
C     IOp(169) controls the browsing file type
C     = 0 Original browsing file(Gaussian default)
C     = 1 New browsing file
C     =11 New browsing file (Hessian stored in addition)
      IOp169 = IOp(169)
CC-INSERT 
      Call FileIO(2,IRWGen,LRwGen,Gen,0)
      If(IPrint.ge.0) Write(IOut,1000)
C
C     Turn off symmetry right away.
C
      Call ILSW(1,26,1)
C
C     For restarted calculations, restore the job state from the checkpoint
C     file. After a bit more initialization below, we will wrap-up the
C     restart work and leave Link 123.
C
      If(Restrt) then
        Call ChkPnt(1,V,MDV)
        Call RWMol(2,0)
        endIf
C
C     Fill NAt3, NAt3TT, NFroz, and NAtomQ. Note that NAtomQ is set here as
C     though the QM:MM machinery in this link is being used. After
C     interogating the IOPs, Routine ExplOp re-sets NAtomQ to 0 if QM:MM
C     codes are not to be used.
C
      NAt3   = 3*NAtoms
      NAt3TT = (NAt3*(NAt3+1))/2
      Call FrzAtm(.False.,NAtoms,MicOpt,NFroz,V)
      Call OnMask(NAtoms,.False.,.False.,MicOpt,NAtomQ,V)
C
C     Check if first run. If this is the first time (First=TRUE) initialize
C     some variables and check to see if there is any data to read from the
C     input file. If this is not the first time into Link 123
C     (First=FALSE), we checkpoint and then pick-up the common block from
C     the RWF.
C
      If(First) then
        ICurWh = 0
        IPath = 1
        IStep = 2
        IStep2 = 0
        IHsStp = 0
        ICorFl = 0
        NTotGd = 0
        NTotHs = 0
        Call DfStat(In,IOut,IPrint,IRdInp,IOZRed,NAtoms,NRdStat,V,MDV)
      else
        If(IOp36.eq.0) Call ChkPnt(0,V,MDV)
        Call TStamp(1,'Explor: IO-1')
        Call FileIO(2,-IRWExp,LenCom,StpSiz,0)
        Call TStamp(1,'Explor: IO-1 Done')
        NAtomQ = NAtmQX
        Call GetNOp(I,J)
        endIf
C
C     Set-up the link's option flags based primarily on the IOps using
C     Routine ExplOp. Then, print the job parameters. If this is a restart
C     job, then we save some key variables that appear in the common block
C     in case they have been changed for this run. We then restore these
C     variables to their old values for allocating space, and let BkAjst
C     take care of setting the final flags to the appropriate values, etc.
C
      If(Restrt) then
        IMthd1 = IMthd
        NPath1 = NPath
        IRxFw1 = IRxFwd
        MxStp1 = MaxStp
        MxOpt1 = MxOpt
      else
        IMthd1 = 0
        NPath1 = 0
        IRxFw1 = 0
        MxStp1 = 0
        MxOpt1 = 1
        endIf
      Call ExplOp(IOut,IOp,First,Restrt,StpSiz,DVVv0,NRdStat,NCorOp,
     $  NAtmQX,CritGd,CritAn,EPS,NAtoms,NAtomQ,IMthd,IEmp,IHesPrd,
     $  IHesCor,IHsUdt,ICorOp,IPFrq,IRxFwd,NPhase,IEig,NStat,MxOpt,
     $  IQMMM,NReduc,IOnDir,Compct,IntHes,NPath,MaxStp,MasFlg,DXConv,
     $  ConvF,SumRlE,SumCnv,SvPES,DoNRTst,WrURVA,IMthd1,NPath1,IRxFw1,
     $  MxStp1)
      If(First.or.Restrt) Call InfoRx(Restrt,MaxStp,StpSiz,NPath,iRxFwd,
     $  IOp10,IHesPrd,IHesCor,IMthd,NCorOp,ICorOp,MxOpt,IHsUdt,IEig)
      If(Restrt) then
        ITemp  = NPath
        NPath  = NPath1
        NPath1 = ITemp
        ITemp  = IRxFwd
        IRxFwd = IRxFw1
        IRxFw1 = ITemp
        ITemp  = MaxStp
        MaxStp = MxStp1
        MxStp1 = ITemp
        ITemp  = MxOpt
        MxOpt  = MxOpt1
        MxOpt1 = ITemp
        endIf
C
C     Take care of NAtomQ, NAtQ3, NAtQ3T, HvHess, DoHEVc, Update, and
C     DwnHil flags. Note that NAtomQ is the number of atoms that should be
C     considered for treatments for operations specifically carried out
C     over the QM sub-system.
C
      If(IQMMM.eq.1.and.Compct) then
        NAtQ3 = 3*NAtomQ
        NAtQ3T = (NAtQ3*(NAtQ3+1))/2
      else
        NAtQ3 = 0
        NAtQ3T = 0
        endIf
      DwnHil = IRxFwd.eq.3
      HvHess = (IHesPrd.ne.-1).or.(IHesCor.gt.0)
      DoHEVc = HvHess.and.(IMthd.eq.2.or.IMthd.eq.3.or.IMthd.eq.4.or.
     $  IMthd.eq.10)
      Update = HvHess.and.IHsStp.ne.0
      If(IPrint.ge.2.and.Update) then
        Write(IOut,8000) IStep,IStep2
      else if(IPrint.ge.2) then
        Write(IOut,8001) IStep,IStep2
        endIf
      IURVA = IUnit(18)
      If(First.and.WrURVA) Call InURVA(IOut,IURVA,IPrint,.True.,HvHess,
     $  LURVAFl,IURVAFl,NAtoms,0,Zero,Zero,JJ,XX,XX,XX,
CC-INSERT     
     $  XX,XX,XX,IOp169)
CC-INSERT 
C
C     Allocate space.
C
      Call Ini123(IOut,IPrint,IIRdSt,IMskON,IAtMas,IRtMas,ITsGrd,ITsGdQ,
     $  ITsCrt,ITsVec,ITsFc,ITsEVc,ITsLamb,IIDone,IAGrad,IGradQ,IHsCar,
     $  IHsEVc,IHsEVl,IStats,ICart,IGrdSv,IVel,IGrdCu,IHesCu,IHVcCu,
     $  IHVlCu,IvUpda,IwUpda,IFroze,IScr1,LenArr,IEnd,IBlck1,IBlck2,
     $  IBlck3,IBlck4,IBlck5,IBlck6,LBlck1,LBlck2,LBlck3,LBlck4,LBlck5,
     $  LBlck6,IOfBk1,IOfBk2,IOfBk3,IOfBk4,IOfBk5,IOfBk6,NAtoms,NAtomQ,
     $  LenHes,LnHEVc,LnHEVl,LenSq,MaxStp,MxOpt,NPath,NStat,NCorOp,
     $  ICorOp,.True.,HvHess,DoHEVc,DwnHil,IFtMax,IQMMM,NReduc,IntHes,
     $  NFroz,MDV)
C
C     Initialize the PES data structure. The portion of the PES object held
C     in memory are put in place after the memory allocated by Routine
C     Ini123.
C
      IPESOb = IEnd
      LenObj = 0
      If(First) then
        If(NPath.gt.MxPath) Call GauErr('Explor: NPath > MxPath.')
        Call ISet(MaxStp,NPath,MxPnt)
        If(SvPES) Call PESIRC(IOut,IPrint,0,IRwPES,NAt3,NAt3TT,HvHess,
     $    NPath,MxPnt,0,0,.False.,XX,.False.,XX(1),.False.,XX,.False.,
     $    XX,.False.,XX(1),XX,LenObj,MDV-IPESOb+1,V(IPESOb),OK)
        IEnd = IPESOb+LenObj
      else
        If(SvPES) Call PESIRC(IOut,IPrint,2,IRwPES,NAt3,NAt3TT,HvHess,
     $    NPath,MxPnt,0,0,.False.,XX,.False.,XX(1),.False.,XX,.False.,
     $    XX,.False.,XX(1),XX,LenObj,MDV-IPESOb+1,V(IPESOb),OK)
        IEnd = IPESOb+LenObj
      endIf
C
C     Pick-up various blocks of data from the Link 123 RWF. Note that we
C     try to limit the N**2 reading if possible.
C
      If(First.and.IOp10.ne.5) then
        DoHess = HvHess.or.(.not.DwnHil)
        Goto 900
      else if(.not.First) then
        Call TStamp(1,'Explor: IO-2')
        Call TStamp(1,'Explor: Block 1 I/O')
        If(LBlck1.gt.0) Call FileIO(2,-IRWExp,LBlck1,V(IBlck1),
     $    LenCom+IOfBk1)
        Call TStamp(1,'Explor: Block 2 I/O')
        If(LBlck2.gt.0) Call FileIO(2,-IRWExp,LBlck2,V(IBlck2),
     $    LenCom+IOfBk2)
        Call TStamp(1,'Explor: Block 4 I/O')
        If(LBlck4.gt.0) Call FileIO(2,-IRWExp,LBlck4,V(IBlck4),
     $    LenCom+IOfBk4)
        Call TStamp(1,'Explor: Block 5 I/O')
        If(LBlck5.gt.0.and..not.Restrt) then
          If(HvHess) then
            ITmpLn = 2*LenHes
            IOff1 = (IPath-1)*ITmpLn
            IOff2 = LenCom+IHsCar+IOff1-1
            Call FileIO(2,-IRWExp,ITmpLn,V(IHsCar+IOff1),IOff2)
            endIf
          If(DoHEVc) then
            IOff1 = (IPath-1)*LnHEVc
            IOff2 = LenCom+IHsEVc+IOff1-1
            Call FileIO(2,-IRWExp,LnHEVc,V(IHsEVc+IOff1),IOff2)
            IOff1 = (IPath-1)*LnHEVl
            IOff2 = LenCom+IHsEVl+IOff1-1
            Call FileIO(2,-IRWExp,LnHEVl,V(IHsEVl+IOff1),IOff2)
            endIf
        else if(LBlck5.gt.0) then
          Call FileIO(2,-IRWExp,LBlck5,V(IBlck5),LenCom+IOfBk5)
          endIf
        Call TStamp(1,'Explor: Block 6 I/O')
        If(LBlck6.gt.0.) Call FileIO(2,-IRWExp,LBlck6,V(IBlck6),
     $    LenCom+IOfBk6)
        Call TStamp(1,'Explor: IO-2 Done')
        endIf
C
C     Handle restarting and checkpointing. If we're restarting a failed
C     job, then check for changes in MaxStp, NPath, and StpSiz and then
C     carefully resize arrays as needed.
C
      If(Restrt) then
        Write(IOut,2000)
        If(IRxFwd.ne.IRxFw1) then
          Write(IOut,9100)
          Call GauErr('Link 123: Invalid restart parameters.')
          endIf
        NPL = IGetVL(V(IIDone),NPath)
        If(NPL.gt.MaxStp) then
          NDone = NPath
        else
          NDone = 0
          endIf
        If(NPath.eq.NPath1.and.MaxStp.eq.MxStp1.and.NDone.eq.NPath) then
          If(IPrint.ge.0) Write(IOut,1040)
          Call Summry(IOut,V(IIDone),SumRlE,SumCnv,TSPot,V(IStats),
     $      NPath,MaxStp,NStat,V(IIRdSt),V(LenArr+1),MDV-LenArr)
          If(IPrint.ge.0) Write(IOut,1040)
          IChain = 1
          Goto 900
        else if(NPath.ne.NPath1.or.MaxStp.ne.MxStp1) then
          Call BkAjst(IOut,IPrint,IRWExp,LenCom,LenArr,IRxFw1,MxStp1,
     $      NPath1,NAtoms,NAtomQ,NStat,NCorOp,ICorOp,HvHess,DoHEVc,
     $      DwnHil,IFtMax,IQMMM,NReduc,IntHes,NFroz,ICart,SkpEGH,V,MDV)
          Call OTIRC1(IOut,-1,NAt3,MxStp1,NPath1,0,XX,Zero,XX,Zero,
     $      V(IEnd),MDV-IEnd+1)
          endIf
        DoHess = DoHess.or.IHesRst.eq.1
        If(.True..or.SkpEGH) then
          Call Ini123(IOut,IPrint,IIRdSt,IMskON,IAtMas,IRtMas,ITsGrd,
     $      ITsGdQ,ITsCrt,ITsVec,ITsFc,ITsEVc,ITsLamb,IIDone,IAGrad,
     $      IGradQ,IHsCar,IHsEVc,IHsEVl,IStats,ICart,IGrdSv,IVel,IGrdCu,
     $      IHesCu,IHVcCu,IHVlCu,IvUpda,IwUpda,IFroze,IScr1,LenArr,IEnd,
     $      IBlck1,IBlck2,IBlck3,IBlck4,IBlck5,IBlck6,LBlck1,LBlck2,
     $      LBlck3,LBlck4,LBlck5,LBlck6,IOfBk1,IOfBk2,IOfBk3,IOfBk4,
     $      IOfBk5,IOfBk6,NAtoms,NAtomQ,LenHes,LnHEVc,LnHEVl,LenSq,
     $      MaxStp,MxOpt,NPath,NStat,NCorOp,ICorOp,.True.,HvHess,DoHEVc,
     $      DwnHil,IFtMax,IQMMM,NReduc,IntHes,NFroz,MDV)
        else
          IPoint = ICart+NAt3*MaxStp*(IPath-1)+NAt3*(IStep-2)
          Call AMove(NAt3,V(IPoint),C)
          If(IStep2.le.1) then
            If(IHesPrd.gt.0) DoHess = DoHess.or.
     $        Mod(IStep-2,IHesPrd).eq.0
          else
            If(IHesCor.gt.0) DoHess = DoHess.or.
     $        Mod(IStep2-1,IHesCor).eq.0
            endIf
          IChain = 0
          Goto 900
          endIf
        endIf
C
C     Take case of some key arrays that will be needed in the lower-level
C     routines:
C           1. Load atomic masses & root-masses;
C           2. Initialize IDone if appropriate;
C           3. Load ONIOM masking array;
C
      If(MasFlg.eq.2) then
        Call ASet(NAtoms,One,V(IAtMas))
      else
        Call AMove(NAtoms,AtmWgt,V(IAtMas))
        endIf
      Call ASqrt(0,NAtoms,V(IAtMas),V(IRtMas))
      If(IPath.eq.1.and.IStep.eq.2.and..not.Restrt) Call IClear(NPath,
     $  V(IIDone))
      If(IQMMM.eq.1.and.Compct) then
        Call OnMask(NAtoms,.True.,.True.,MicOpt,NAtmQ1,V(IMskON))
        If(NAtomQ.ne.NAtmQ1) then
          Write(IOut,9000) NAtomQ,NAtmQ1,NAtmQX
          Call GauErr('Explor: NAtomQ.ne.NAtmQ1')
          endIf
        endIf
C
C     Print job parameters.
C
      If(IPrint.ge.2.and..not.(First.or.Restrt)) Call InfoRx(Restrt,
     $  MaxStp,StpSiz,NPath,iRxFwd,IOP10,IHesPrd,IHesCor,IMthd,NCorOp,
     $  ICorOp,MxOpt,IHsUdt,IEig)
C
C     Everything is loaded. Go into RxnCyc and do the work.
C
 100  Call TStamp(1,'Before RxnCyc')
      Call RxnCyc(In,IOut,IPrint,IRwFit,Multip,SvPES,DoNRTst,WrURVA,
     $  LURVAFl,IURVAFl,SkpEGH,ICurWh,StpSiz,DVVv0,IStep,IStep2,IPath,
     $  IPG,ICorFl,CritGd,CritAn,TsPot,TsPotQ,Potent,PotSav,EQM,EPS,
     $  IMthd,IEmp,IHsUdt,NCorOp,ICorOp,IRxFwd,NPhase,IPFrq,IEig,NStat,
     $  MxOpt,NReduc,IQMMM,NPath,MaxStp,Gen,C,V(IIRdSt),V(IMskON),
     $  V(IAtMas),V(IRtMas),V(ITSGrd),V(ITsGdQ),V(ITSFc),V(ITSCrt),
     $  V(ITSVec),V(ITsEVc),V(ITSLamb),V(IAGrad),V(IGradQ),V(IGrdCu),
     $  V(IGrdSv),V(IHsCar),V(IHsEVc),V(IHsEVl),V(IHesCu),V(IHVcCu),
     $  V(IHVlCu),V(IvUpda),V(IwUpda),V(IStats),V(ICart),V(IVel),
     $  V(IIDone),V(IFroze),V(IScr1),NAtoms,NAt3,LenHes,LnHEVc,LnHEVl,
     $  LenSq,NAt3TT,NAtomQ,NAtQ3,NAtQ3T,IAn,Restrt,IOP10,JobFlg,HvHess,
     $  DoHEVc,Update,Compct,IntHes,DwnHil,CorDon,IFtSt,IFtMax,MicOpt,
     $  IAtTyp,AtChMM,IAtRes,NResid,IResNa,IResNu,DXConv,ConvF,LenObj,
C     $  V(IPESOb),DoSvC,V(IEnd),MDV-IEnd+1)
CC-INSERT     
     $ IOp169,V(IPESOb),DoSvC,V(IEnd),MDV-IEnd+1)
CC-INSERT 
      Call TStamp(1,'After RxnCyc')
C
C     Save the PES data structure.
C
      If(SvPES) then
        Call PESIRC(IOut,IPrint,1,IRwPES,NAt3,NAt3TT,HvHess,NPath,MxPnt,
     $    0,0,.False.,XX,.False.,XX(1),.False.,XX,.False.,XX,.False.,
     $    XX(1),XX,LenObj,MDV-IPESOb+1,V(IPESOb),OK)
        endIf
C
C     Analyze JobFlg to see what needs to be done and use ILSW(24) to
C     control the order of derivatives calculated next.
C          JobFlg = 0  Get energy and derivatives for C.
C          JobFlg = 1  Already have energy and derivatives for C => loop
C                      back to call to RxnCyc.
C          JobFlg = 2  Integration converged or ran out of cycles =>
C                      if we need to calculate the next path, we just
C                      loop back to the call to RxnCyc since we already
C                      have the structure, energy, and derivatives at
C                      the TS.
C
      SkpEGH = .False.
      If(JobFlg.eq.0.or.JobFlg.eq.1) then
        IStep = IStep + 1
        If(CorDon) then
          ICurWh = 0
          IStep2 = 1
          If(IHesPrd.gt.0) DoHess = DoHess.or.Mod(IStep-2,IHesPrd).eq.0
        else
          IStep2 = IStep2 + 1
          If(IHesCor.gt.0) DoHess = DoHess.or.Mod(IStep2-1,IHesCor).eq.0
          endIf
        If(JobFlg.eq.1) then
          Call TStamp(1,'Explor: Block 3 I/O')
          If(LBlck3.gt.0) Call FileIO(2,-IRWExp,LBlck3,V(IBlck3),
     $      LenCom+IOfBk3)
          Call TStamp(1,'Explor: After Block 3 I/O')
          Goto 100
          endIf
        IChain = 0
      else if(JobFlg.eq.2) then
        ICurWh = 0
        JPath = IPath + 1
        IStep2 = 1
        DoHess = DoHess.or.IHesRst.eq.1
        If(JPath.le.NPath) then
          NxDone = IIGet(V(IIDone),JPath)
          If(NxDone.eq.0) then
            IStep=2
            If(IPrint.ge.0) Write(IOut,1040)
            Call TStamp(1,'Explor: Block 3 I/O')
            If(LBlck3.gt.0) Call FileIO(2,-IRWExp,LBlck3,V(IBlck3),
     $        LenCom+IOfBk3)
            Call TStamp(1,'Explor: After Block 3 I/O')
            Call TStamp(1,'Explor: Block 5 I/O')
            If(LBlck5.gt.0) then
              If(HvHess) then
                ITmpLn = 2*LenHes
                IOff1 = (IPath-1)*ITmpLn
                IOff2 = LenCom+IHsCar+IOff1-1
                Call FileIO(1,-IRWExp,ITmpLn,V(IHsCar+IOff1),IOff2)
                endIf
              If(DoHEVc) then
                IOff1 = (IPath-1)*LnHEVc
                IOff2 = LenCom+IHsEVc+IOff1-1
                Call FileIO(1,-IRWExp,LnHEVc,V(IHsEVc+IOff1),IOff2)
                IOff1 = (IPath-1)*LnHEVl
                IOff2 = LenCom+IHsEVl+IOff1-1
                Call FileIO(1,-IRWExp,LnHEVl,V(IHsEVl+IOff1),IOff2)
                endIf
              Call TStamp(1,'Explor: After Block 5 I/O')
              endIf
          else
            SkpEGH = .True.
            IStep = NxDone + 1
            Call TStamp(1,'Explor: Block 3 I/O')
            If(LBlck3.gt.0) Call FileIO(2,-IRWExp,LBlck3,V(IBlck3),
     $        LenCom+IOfBk3)
            Call TStamp(1,'Explor: After Block 3 I/O')
            Call TStamp(1,'Explor: Block 5 I/O')
            If(LBlck5.gt.0) then
              If(HvHess) then
                ITmpLn = 2*LenHes
                IOff1 = (IPath-1)*ITmpLn
                IOff2 = LenCom+IHsCar+IOff1-1
                Call FileIO(1,-IRWExp,ITmpLn,V(IHsCar+IOff1),IOff2)
                IOff1 = (JPath-1)*ITmpLn
                IOff2 = LenCom+IHsCar+IOff1-1
                Call FileIO(2,-IRWExp,ITmpLn,V(IHsCar+IOff1),IOff2)
                endIf
              If(DoHEVc) then
                IOff1 = (IPath-1)*LnHEVc
                IOff2 = LenCom+IHsEVc+IOff1-1
                Call FileIO(1,-IRWExp,LnHEVc,V(IHsEVc+IOff1),IOff2)
                IOff1 = (IPath-1)*LnHEVl
                IOff2 = LenCom+IHsEVl+IOff1-1
                Call FileIO(1,-IRWExp,LnHEVl,V(IHsEVl+IOff1),IOff2)
                IOff1 = (JPath-1)*LnHEVc
                IOff2 = LenCom+IHsEVc+IOff1-1
                Call FileIO(2,-IRWExp,LnHEVc,V(IHsEVc+IOff1),IOff2)
                IOff1 = (JPath-1)*LnHEVl
                IOff2 = LenCom+IHsEVl+IOff1-1
                Call FileIO(2,-IRWExp,LnHEVl,V(IHsEVl+IOff1),IOff2)
                endIf
              Call TStamp(1,'Explor: After Block 5 I/O')
              endIf
            endIf
          IPath = JPath
          Goto 100
        else
          If(IPrint.ge.0) Write(IOut,1040)
          Call Summry(IOut,V(IIDone),SumRlE,SumCnv,TSPot,V(IStats),
     $      NPath,MaxStp,NStat,V(IIRdSt),V(LenArr+1),MDV-LenArr)
          If(IPrint.ge.0) Write(IOut,1040)
          IChain = 1
          endIf
      else
        Write(IOut,9900) JobFlg
        Call GauErr('Logic Error in Explor: Unknown JobFlg.')
        endIf
  900 If(DoHess) then
        If(DEBUG.or.IPrint.ge.2) Write(IOut,8010) IStep,IStep2
        IHsStp = 0
        Call ILSW(1,24,2)
      else
        If(DEBUG.or.IPrint.ge.2) Write(IOut,8011) IStep,IStep2
        IHsStp = 1
        Call ILSW(1,24,1)
        endIf
      If(IChain.eq.0) then
        NTotGd = NTotGd+1
        If(DoHess) NTotHs = NTotHs+1
      else
        NTotSt = IIGet(V(IIDone),1)-1
        If(NPath.eq.2) NTotSt = NTotSt+IIGet(V(IIDone),2)-1
        Write(IOut,5000) NTotSt,NTotGd,NTotHs
        endIf
C
C     Save things to the RWF.
C
      Call TStamp(1,'Explor: IO-3')
      Call FileIO(1,-IRWGen,LRwGen,Gen,0)
      If(ITQry(IRWExp).lt.(LenCom+LenArr)) Call CondDf(IRWExp,
     $  LenCom+LenArr)
      Call FileIO(1,-IRWExp,LenCom,StpSiz,0)
      Call TStamp(1,'Explor: Block 1 I/O')
      If(LBlck1.gt.0) Call FileIO(1,-IRWExp,LBlck1,V(IBlck1),
     $  LenCom+IOfBk1)
      Call TStamp(1,'Explor: Block 2 I/O')
      If(LBlck2.gt.0) Call FileIO(1,-IRWExp,LBlck2,V(IBlck2),
     $  LenCom+IOfBk2)
      Call TStamp(1,'Explor: Block 3 I/O')
      If(LBlck3.gt.0.and.
     $  (IPath.eq.1.and.IStep.eq.3.and.IStep2.eq.1.and..not.DwnHil))
     $  Call FileIO(1,-IRWExp,LBlck3,V(IBlck3),LenCom+IOfBk3)
      Call TStamp(1,'Explor: Block 4 I/O')
      If(LBlck4.gt.0) Call FileIO(1,-IRWExp,LBlck4,V(IBlck4),
     $  LenCom+IOfBk4)
      Call TStamp(1,'Explor: Block 5 I/O')
      If(LBlck5.gt.0) then
        If(HvHess) then
          ITmpLn = 2*LenHes
          IOff1 = (IPath-1)*ITmpLn
          IOff2 = LenCom+IHsCar+IOff1-1
          Call FileIO(1,-IRWExp,ITmpLn,V(IHsCar+IOff1),IOff2)
          endIf
        If(DoHEVc) then
          IOff1 = (IPath-1)*LnHEVc
          IOff2 = LenCom+IHsEVc+IOff1-1
          Call FileIO(1,-IRWExp,LnHEVc,V(IHsEVc+IOff1),IOff2)
          IOff1 = (IPath-1)*LnHEVl
          IOff2 = LenCom+IHsEVl+IOff1-1
          Call FileIO(1,-IRWExp,LnHEVl,V(IHsEVl+IOff1),IOff2)
          endIf
        endIf
      Call TStamp(1,'Explor: Block 6 I/O')
      If(LBlck6.gt.0) Call FileIO(1,-IRWExp,LBlck6,V(IBlck6),
     $  LenCom+IOfBk6)
      Call TStamp(1,'Explor: IO-3 Done')
      If(IPrint.ge.0) Write(IOut,1000)
      Call ChainX(IChain,V,MDV)
      Return
      End
*Deck BkAjst
      Subroutine BkAjst(IOut,IPrint,IRWExp,LenCom,LenArr,IRxFw1,MxStp1,
     $  NPath1,NAtoms,NAtomQ,NStat,NCorOp,ICorOp,HvHess,DoHEVc,DwnHil,
     $  IFtMax,IQMMM,NReduc,IntHes,NFroz,ICart,SkpEGH,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine adjusts the data on RWF file number IRWExp for Link 123.
C     This is needed for restarted jobs where the number of steps or paths
C     requested by the user has changed.
C
C     Before calling this routine, V should already have the permanent
C     arrays loaded.
C
C     Arrays re-sized here include:
C           (1) Stats
C           (2) Cart
C           (3) Vel
C     The storage of these vectors in memory must be in this order.
C
      Parameter (LURVAFl=1024)
      Common/RxnDat/StpSiz,DVVv0,CritGd,CritAn,TsPot,TsPotQ,Potent(2,2),
     $  PotSav,EQM(2,2),EPS,IMthd,NRdStat,ICurWh,NAtmQX,IRxFwd,MaxStp,
     $  IStep,MxOpt,IStep2,NPath,IPath,IHsStp,IPG(4),ICorFl,
     $  IURVAFl(LURVAFl),NTotGd,NTotHs,IDmRxD
      Logical HvHess,DoHEVc,DwnHil,IntHes,SkpEGH
      Dimension ICorOp(NCorOp),V(MDV)
 2000 Format(1x,'Max # of points and cannot be decreased in a restart.',
     $  /,3x,'Old Max Points=',I6,/,3x,'New Max Points=',I6)
 2010 Format(1x,'Number of paths cannot be changed in a restart.',/,
     $  3x,'Old # of Paths=',I6,/,3x,'New # of Paths=',I6)
 2020 Format(1x,'Path direction flag cannot be changed in a restart.',/,
     $  3x,'Old IRxFwd=',I6,/,3x,'New IRxFwd=',I6)
 5000 Format(1x,'IPath=',I3,/,3x,A,' from position ',I10,
     $  ' to position ',I10)
C
C     Check to see if MaxStp and NPath are staying the same as before, in
C     which case no additional work is needed here and we can just exit.
C     If there are changes, ensure that MaxStp and NPath are not being
C     decreased.  Then, initialize the flag SkpEGH, which instructs
C     level-two routines to skip steps related to loading energy, gradient,
C     and Hessians.
C
      If(MaxStp.eq.MxStp1.and.NPath.eq.NPath1) Return
      If(MaxStp.gt.MxStp1) then
        Write(IOut,2000) MaxStp,MxStp1
        Call GauErr('L123 Restart: MaxStp>MxStp1')
        endIf
      If(NPath.gt.NPath1) then
        Write(IOut,2010) NPath,NPath1
        Call GauErr('L123 Restart: NPath>NPath1')
        endIf
      If(IRxFwd.ne.IRxFw1) then
        Write(IOut,2020) IRxFwd,IRxFw1
        Call GauErr('L123 Restart: IRxFwd.ne.IRxFw1')
        endIf
      SkpEGH = .False.
C
C     Since one or more of the key job parameters have changed, it may be
C     necessary to resize and/or shift elements of some arrays.  There are
C     three basic cases allowed here:
C           Case (1) - both MaxStp and NPath change: Arrays that must be
C                      modified are IDone, Stats, Cart, and Vel.
C           Case (2) - only MaxStp changes: Arrays that must be modified
C                      are Stats, Cart, and Vel.
C           Case (3) - only NPath changes: Arrays that must be modified are
C                      IDone, Stats, Cart, and Vel.
C
C     Begin by computing array pointers for old and new MaxStp and NPath.
C
      Call Ini123(IOut,IPrint,IIRdSt,IMskON,IAtMas,IRtMas,ITsGrd,ITsGdQ,
     $  ITsCrt,ITsVec,ITsFc,ITsEVc,ITsLamb,IIDone,IAGrad,IGradQ,IHsCar,
     $  IHsEVc,IHsEVl,IStats,ICart,IGrdSv,IVel,IGrdCu,IHesCu,IHVcCu,
     $  IHVlCu,IvUpda,IwUpda,IFroze,IScr1,LenArr,IEnd,IBlck1,IBlck2,
     $  IBlck3,IBlck4,IBlck5,IBlck6,LBlck1,LBlck2,LBlck3,LBlck4,LBlck5,
     $  LBlck6,IOfBk1,IOfBk2,IOfBk3,IOfBk4,IOfBk5,IOfBk6,NAtoms,NAtomQ,
     $  LenHes,LnHEVc,LnHEVl,LenSq,MaxStp,MxOpt,NPath,NStat,NCorOp,
     $  ICorOp,.True.,HvHess,DoHEVc,DwnHil,IFtMax,IQMMM,NReduc,IntHes,
     $  NFroz,MDV)
      Call Ini123(IOut,IPrint,JIRdSt,JMskON,JAtMas,JRtMas,JTsGrd,JTsGdQ,
     $  JTsCrt,JTsVec,JTsFc,JTsEVc,JTsLamb,JIDone,JAGrad,JGradQ,JHsCar,
     $  JHsEVc,JHsEVl,JStats,JCart,JGrdSv,JVel,JGrdCu,JHesCu,JHVcCu,
     $  JHVlCu,JvUpda,JwUpda,JFroze,JScr1,LnArr1,JEnd,JBlck1,JBlck2,
     $  JBlck3,JBlck4,JBlck5,JBlck6,MBlck1,MBlck2,MBlck3,MBlck4,MBlck5,
     $  MBlck6,JOfBk1,JOfBk2,JOfBk3,JOfBk4,JOfBk5,JOfBk6,NAtoms,NAtomQ,
     $  LenHes,LnHEVc,LnHEVl,LenSq,MxStp1,MxOpt,NPath,NStat,NCorOp,
     $  ICorOp,.True.,HvHess,DoHEVc,DwnHil,IFtMax,IQMMM,NReduc,IntHes,
     $  NFroz,MDV)
      Call TstCor(LnArr1,MDV,'BkAjst')
      LIStat = NStat*MaxStp
      LICart = 3*NAtoms*MaxStp
      LIVel  = 3*NAtoms*MaxStp
      LJStat = NStat*MxStp1
      LJCart = 3*NAtoms*MxStp1
      LJVel  = 3*NAtoms*MxStp1
      If(NPath.ne.NPath1.and.MaxStp.ne.MxStp1) then
        Call GauErr('BkAjst: Case(1) NYI!')
      else if(MaxStp.ne.MxStp1) then
        Do 201 i = NPath,1,-1
          IIVel = IVel+(i-1)*LIVel-1
          JJVel = JVel+(i-1)*LJVel-1
          If(IPrint.ge.2) Write(IOut,5000) i,'Moving Vel',IIVel+1,
     $      JJVel+1
          Call AMove1(LIVel,IIVel,JJVel,V)
 201      Continue
        Do 202 i = NPath,1,-1
          IICart = ICart+(i-1)*LICart-1
          JJCart = JCart+(i-1)*LJCart-1
          If(IPrint.ge.2) Write(IOut,5000) i,'Moving Cart',IICart+1,
     $      JJCart+1
          Call AMove1(LICart,IICart,JJCart,V)
 202      Continue
        Do 203 i = NPath,1,-1
          IIStat = IStats+(i-1)*LIStat-1
          JJStat = JStats+(i-1)*LJStat-1
          If(IPrint.ge.2) Write(IOut,5000) i,'Moving Stats',IIStat+1,
     $      JJStat+1
          Call AMove1(LIStat,IIStat,JJStat,V)
 203      Continue
      else if(NPath.ne.NPath1) then
        Call GauErr('BkAjst: Case(3) NYI!')
      else
        Call GauErr('BkAjst: Unknown restart state in Link 123.')
        endIf
C
C     Now that resizing and shifting is complete, set IRxFwd, MaxStp,
C     NPath, IStep, and IPath to the new values.  Then, if pathway 1
C     completed in the previous run, we need to set IStep for the first
C     path to IDone(1) and instruct the second-level routines to skip the
C     step where energy, gradient, and Hessian are read. Note that IStep
C     must be incremented by one in all cases.
C
      MaxStp = MxStp1
      IRxFwd = IRxFw1
      IPath = 1
      NPath = NPath1
      ICart = JCart
      IDone = IIGet(V(JIDone),1)
      If(IDone.gt.0) then
        IStep = IDone
        SkpEGH = .True.
        endIf
      IStep = IStep+1
C
C     Resize the RWF file and write out the common block and saved arrays.
C
      LenArr = LnArr1
      Call Conddf(IRWExp,LenCom+LenArr)
      Call FileIO(1,-IRWExp,LenCom,StpSiz,0)
      Call FileIO(1,IRWExp,LenArr,V,0)
      Return
      End
*Deck CorDrv
      Subroutine CorDrv(IOut,NAt3,NAtoms,RtMass,Cart1,Cart2,StpSiz,
     $  TSEVec,HsEVec,HvHess,IEig,RxnFwd,AGrad,IStep)
      Implicit Real*8(A-H,O-Z)
C
C     This routine takes the cartesian coordinates of a stationary point
C     and steps in the direction of the specified eigenvector of the second
C     derivative matrix.
C
C     IEig .......... the number specifying which eigenvector to drive set
C                     by iop(1/123)
C
C     StpSiz ........ Set by iop(1/39), see Routine ExplOp.
C
C
      Dimension RtMass(NAtoms),Cart1(NAt3),Cart2(NAt3),
     $  TSeVec(NAt3,NAt3),HsEVec(NAt3,NAt3),AGrad(NAt3,2)
      Logical RxnFwd,HvHess
      Save Zero
      Data Zero/0.d0/
 1000 Format('Negating the step direction.')
C
C     First we need to mass-weigh the cartesian coordinates and get a unit
C     vector in the direction of TSeVec.  We use RxnFwd to decide whether
C     or not to change sign.
C
      Call MWeigV(.False.,NAtoms,Cart1,RtMass,Cart1)
C
C     Here we use the AGrad array for eigenvector storage rather than the
C     actual gradient.  Then, if this is the first time in (IStep=2) then
C     flip the sign on AGrad if RxnFwd is .FALSE.
C
      If(IStep.eq.2) then
        If(HvHess) then
          Call AUnit(NAt3,HsEVec(1,IEig),AGrad)
          If(.not.RxnFwd) Call ANeg(NAt3,AGrad,AGrad)
        else
          Call AUnit(NAt3,TSEVec(1,IEig),AGrad)
          endIf
      else if(HvHess) then
        IBest = 1
        Best = SProd(NAt3,HsEVec(1,1),AGrad(1,2))
        Do 100 i = 2,NAt3
          Temp = SProd(NAt3,HsEVec(1,i),AGrad(1,2))
          If(Abs(Temp).gt.Abs(Best)) then
            IBest = i
            Best = Temp
            endIf
  100     Continue
        Call AMove(NAt3,HsEVec(1,IBest),AGrad)
        If(Best.lt.Zero) then
          Write(IOut,1000)
          Call ANeg(NAt3,AGrad,AGrad)
          endIf
      else
        Call AMove(NAt3,TSEVec(1,IEig),AGrad)
        endIf
      If(.not.HvHess.and..not.RxnFwd) Call ANeg(NAt3,AGrad,AGrad)
C
C     The new coordinates are now found by:
C         x_i  =  x_i-1  +  vec/||vec||  *  deltaS
C
      Call ACASB(NAt3,Cart1,AGrad,Cart2,StpSiz)
C
C     Undo mass-weighting.
C
      Call MWeigV(.True.,NAtoms,Cart1,RtMass,Cart1)
      Call MWeigV(.True.,NAtoms,Cart2,RtMass,Cart2)
      Return
      End
*Deck DfStat
      Subroutine DfStat(In,IOut,IPrint,IFlag,IOZRed,NAtoms,NExStat,
     $  IRdStat,MaxMem)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to define and fill the array IRdStat, which
C     contains a list of statistics to record during the job and print with
C     the final job summary.
C
C     INPUT ARGUMENTS
C     IOut        The standard output file unit number.
C     IFlag       A flag indicating what statistics should be defined in
C                 IRdStat. Available options are:
C                 1   ... Read user-defined stats from the input stream.
C                1x   ... Report all Cartesian coordinates.
C               1xx   ... If redundant internals are on the RWF, then report
C                         their values along the IRC.
C              1xxx   ... If redundant internals are on the RWF, then report
C                         values for bond coordinates along the IRC.
C             1xxxx   ... If redundant internals are on the RWF, then report
C                         values for angle coordinates along the IRC.
C            1xxxxx   ... If redundant internals are on the RWF, then report
C                         values for dihedral coordinates along the IRC.
C     IOZRed      The RWF file number where, if needed, we should look for
C                 internal coordinate definitions.
C     NAtoms      The number of atoms.
C     MaxMem      The maximum size of array IRdStat in units of floating
C                 point values.
C
C     OUTPUT ARGUMENTS
C     NExStat     The number of extra statistics defined in IRdStat.
C     IRdStat     The array of definitions for additional statistics. Its
C                 structure is the same as IZMat or IB used in internal
C                 coordinate definitios elsewhere in Gaussian.
C
C        - H.P. Hratchian, December, 2006.
C                          Key updates made in February, 2009.
C
      Integer QParse
      Logical DEBUG,HvZRed,AllRed
      Parameter (NStCol=DEFMAXREDINDEX)
      Parameter (LenLin=1000,MxPTab=1000,IType=1,Index1=IType+1,
     $  Index2=Index1+1,Index3=Index2+1,Index4=Index3+1,
     $  LResult=Index4+1)
      Dimension Line(LenLin),IPTab(MxPTab),IResult(LResult),
     $  IRdStat(NStCol,*),NTRedS(4)
 1000 Format(1x,'Reading IRC summary table definitions from the input ',
     $  'file.')
 2000 Format(1x,'Number of defined statistics: =',I10)
 9000 Format(1x,'DfStat: QParse returned with an error. IStat=',I3)
 9010 Format(1x,'DfStat: Invalid IType found. IType=',I10)
 9020 Format(1x,'DfStat: Found Cartesian label. Index1=',I10,' Index2=',
     $  I10)
 9030 Format(1x,'DfStat: Found curvature element request. Element # ',
     $  I10)
C
C     Initialize things and see if there is anything to read from the input
C     file. If we need to do work here, then pull apart IFlag into its
C     constituent pieces.
C
      DEBUG = .False.
      NExStat = 0
      If(IFlag.eq.0) Return
      IFlag1 = Mod(IFlag,10)
      IFlag2 = Mod(IFlag,100)/10
      IFlag4 = Mod(IFlag,10000)/1000
      IFlag5 = Mod(IFlag,100000)/10000
      IFlag6 = Mod(IFlag,1000000)/100000
      HvZRed = ITqry(IOZRed).gt.0
      AllRed = IFlag4.eq.1.and.IFlag5.eq.1.and.IFlag6.eq.1
C
C     If requested look for internal coordinates, and read whichever types
C     of internal coordinates are called for.
C
      If(HvZRed) then
        Call FileIO(2,-IOZRed,InToWP(4),NTRedS,0)
        NTRed = NTRedS(1)
        NTBond = NTRedS(2)
        NTAng = NTRedS(3)
        NTDih = NTRedS(4)
      else
        NTRed = 0
        NTBond = 0
        NTAng = 0
        NTDih = 0
        endIf
      If(HvZRed.and.AllRed) then
        LenAdd = InToWP(NTRed*NStCol)
        LenTmp = InToWP(NExStat*NStCol)+LenAdd
        Call TstCor(LenTmp,MaxMem,'DfStat needs more memory.')
        IOff = InToWP(4)
        Call FileIO(2,-IOZRed,LenAdd,IRdStat(1,NExStat+1),IOff)
        NExStat = NExStat+NTRed
        endIf
      If(IFlag4.eq.1.and..not.AllRed) then
        LenAdd = InToWP(NTBond*NStCol)
        LenTmp = InToWP(NExStat*NStCol)+LenAdd
        Call TstCor(LenTmp,MaxMem,'DfStat needs more memory (IFlag4).')
        IOff = InToWP(4)
        Call FileIO(2,-IOZRed,LenAdd,IRdStat(1,NExStat+1),IOff)
        NExStat = NExStat+NTBond
        endIf
      If(IFlag5.eq.1.and..not.AllRed) then
        LenAdd = InToWP(NTAng*NStCol)
        LenTmp = InToWP(NExStat*NStCol)+LenAdd
        Call TstCor(LenTmp,MaxMem,'DfStat needs more memory (IFlag5).')
        IOff = InToWP(4)+InToWP(NTBond*NStCol)
        Call FileIO(2,-IOZRed,LenAdd,IRdStat(1,NExStat+1),IOff)
        NExStat = NExStat+NTAng
        endIf
      If(IFlag6.eq.1.and..not.AllRed) then
        LenAdd = InToWP(NTDih*NStCol)
        LenTmp = InToWP(NExStat*NStCol)+LenAdd
        Call TstCor(LenTmp,MaxMem,'DfStat needs more memory (IFlag6).')
        IOff = InToWP(4)+InToWP(NTBond*NStCol)+InToWP(NTAng*NStCol)
        Call FileIO(2,-IOZRed,LenAdd,IRdStat(1,NExStat+1),IOff)
        NExStat = NExStat+NTDih
        endIf
C
C     If requested, add all Cartesian coordinates to the statistics list.
C
      If(IFlag2.eq.1) then
        LenTmp = InToWP(NStCol*(NExStat+3*NAtoms))
        Call TstCor(LenTmp,MaxMem,'DfStat needs more memory (IFlag2).')
        Do 200 i = 1,NAtoms
          Do 200 j = 1,3
            NExStat = NExStat+1
            IRdStat(1,NExStat) = i
            IRdStat(2,NExStat) = -j
            IRdStat(3,NExStat) = 0
 200        IRdStat(4,NExStat) = 0
        endIf
C
C     If requested, read a set of user-defined coordinates to add to the
C     statistics list from the input stream.
C
      If(IFlag1.eq.1) then
        Write(IOut,1000)
        Call RxDfPT(IOut,MxPTab,IPTab,IType,Index1,Index2,Index3,Index4)
        IQPDbg = 0
        If(DEBUG.or.IPrint.ge.2) IQPDbg=2
 900    Call StrIn(In,Line,LenLin,IEOF)
        If(IEOF.ne.0) Call GauErr('IRC stats definition input missing!')
        Call IClear(LResult,IResult)
        Call QPInit(IPTab,0,0,0,IQPDbg)
        IStat = QParse(IResult,IPTab,Line,LenLin)
        If(IStat.ne.0) then
          Write(IOut,9000) IStat
          Call QPErro(Line,LenLin,IOut)
          Call GauErr('Parse error in Link 123.')
          endIf
        If(IResult(IType).ne.0) then
          NExStat = NExStat+1
          Call IClear(NStCol,IRdStat(1,NExStat))
          endIf
        If(IResult(IType).eq.-1) then
          If(IResult(Index2).eq.4) then
            Call IClear(NStCol,IRdStat(1,NExStat+1))
            Call IClear(NStCol,IRdStat(1,NExStat+2))
            IRdStat(1,NExStat)   = IResult(Index1)
            IRdStat(1,NExStat+1) = IResult(Index1)
            IRdStat(1,NExStat+2) = IResult(Index1)
            IRdStat(2,NExStat)   = -1
            IRdStat(2,NExStat+1) = -2
            IRdStat(2,NExStat+2) = -3
            NExStat = NExStat+2
          else if(IResult(Index2).ge.1.and.IResult(Index2).le.3) then
            IRdStat(1,NExStat)   = IResult(Index1)
            IRdStat(2,NExStat)   = -IResult(Index2)
          else
            Write(IOut,9020) IResult(Index1),IResult(Index2)
            Call GauErr('DfStat: Invalid xyz label.')
            endIf
        else if(IResult(IType).ge.1.and.IResult(IType).le.3) then
          Call IMove(NStCol,IResult(Index1),IRdStat(1,NExStat))
        else if(IResult(IType).eq.10) then
          IRdStat(1,NExStat) = -1
        else if(IResult(IType).eq.11) then
          NCurv = 3*NAtoms-7
          If(IResult(Index1).eq.0) then
            NExStat = NExStat-1
            Do 300 i = 1,NCurv
              NExStat = NExStat+1
              Call IClear(NStCol,IRdStat(1,NExStat))
              IRdStat(1,NExStat) = -2
              IRdStat(2,NExStat) = i
 300          Continue
          else if(IResult(Index1).gt.0.and.
     $      IResult(Index1).lt.NCurv) then
            IRdStat(1,NExStat) = -2
            IRdStat(2,NExStat) = IResult(Index1)
          else
            Write(IOut,9030) IResult(Index1)
            Call GauErr('DfStat: Invalid element found.')
            endIf
        else if(IResult(IType).ne.0) then
          Write(IOut,9010) IResult(IType)
          endIf
        If(IResult(IType).ne.0) Goto 900
        endIf
C
C     Print the result of this routine.
C
      If(DEBUG.or.IPrint.ge.2) then
        Write(IOut,2000) NExStat
        Call IMtOuS(IOut,'IRdStat:',0,IRdStat,NStCol,NExStat,NStCol,
     $    NExStat)
        endIf
      Return
      End
*Deck DVV
      Subroutine DVV(IOut,IPrint,NAt3,NAtoms,IStep,MaxStp,TStep,TStep0,
     $  TStep1,v0,GScale,RtMass,Cart0,Cart1,Cart2,Vel,Grad,Grad0,Accel,
     $  VPrime,TSVec,Scr1,RxnFwd,DwnHil)
      Implicit Real*8(A-H,O-Z)
C
C     This routine carries out damped velocity Verlet reaction path
C     following (H. P.  Hratchian and H. B. Schlegel, JPCA, 2002, 106,
C     165-169).  In this routine, Eq. (2) in the Hratchian and Schlegel
C     paper is solved in three parts:
C         1.  v_i-1 = v'_i-1 + .5 a_i dt**2
C         2.  v'_i  = v_i-1 + .5 a_i-1 dt
C         3.  x_i   = x_i-1 + v'_i dt
C
C        - H.P. Hratchian, Fall 2004.
C
      Dimension RtMass(NAtoms),Cart0(NAt3),Cart1(NAt3),Cart2(NAt3),
     $  Vel(NAt3,MaxStp),Grad(Nat3),Grad0(NAt3),Accel(NAt3),
     $  VPrime(NAt3),TSVec(NAt3),Scr1(NAt3,2)
      Logical TScale,RxnFwd,DwnHil
      Save Del0,Half,Pt025,One,Three
      Data Del0/.003d0/,Half/.5d0/,Pt025/.025d0/,One/1.d0/,Three/3.d0/
 1000 Format(' DVV: v0 = ',f10.7)
 1100 Format(' DVV: Current Time Step = ',f10.7)
 1110 Format(' DVV: New Time Step     = ',f10.7)
C
      TScale = .True.
      If(IPrint.ge.0) Write(IOut,1000) v0
C
C     Begin by mass-weighing everything and converting the gradient from
C     Hartrees/bohr to acceleration units (bohr amu**.5/fs**2).  Note: The
C     array Vel is kept in MW coordinates throughout l123.
C
      Call MWeigV(.False.,NAtoms,Cart0,RtMass,Cart0)
      Call MWeigV(.False.,NAtoms,Cart1,RtMass,Cart1)
      If(IStep.eq.2.and..not.DwnHil) then
        Call AClear(NAt3,Accel)
      else
        Call MWeigV(.True.,NAtoms,Grad,RtMass,Accel)
        Call AScale(NAt3,-GScale,Accel,Accel)
        endIf
C
C     If this is the first step, set Vel(1) equal to the transition
C     vector.  Otherwise, finalize the velocity using the current
C     acceleration.  Also, we need to save the current velocity in
C     VPrime for the time step variation evaluation.
C
      If(IStep.eq.2) then
        If(DwnHil) then
          Call AMove(NAt3,Accel,Vel(1,IStep-1))
        else
          Call AMove(NAt3,TSVec,Vel(1,IStep-1))
          If(RxnFwd) Call AScale(NAt3,-One,Vel(1,IStep-1),
     $      Vel(1,IStep-1))
          endIf
      else
        Call AMove(NAt3,Vel(1,IStep-1),VPrime)
        Call ACasB(NAt3,Vel(1,IStep-1),Accel,Vel(1,IStep-1),Half*TStep)
        endIf
C
C     Compute v'_i and x_i.  For storage, v'_i is kept in Vel.  At this
C     point we also damp the velocity.
C
      VelScl = v0/VLen(NAt3,Vel(1,IStep-1))
      Call AScale(NAt3,VelScl,Vel(1,IStep-1),Vel(1,IStep-1))
      Call ACasB(NAt3,Vel(1,IStep-1),Accel,Vel(1,IStep),Half*TStep)
      Call ACasB(NAt3,Cart1,Vel(1,IStep),Cart2,TStep)
C
C     For the time variation step, compute x' and store it in Scr1.
C
      If(TScale .and. IStep.gt.3) then
        tPrime = TStep+TStep0
        Call ACasB(NAt3,Cart0,Vel(1,IStep-2),Scr1,tPrime)
        Call MWeigV(.True.,NAtoms,Grad0,RtMass,Scr1(1,2))
        Call AScale(NAt3,-GScale,Scr1(1,2),Scr1(1,2))
        Call ACasB(NAt3,Scr1,Scr1(1,2),Scr1,Half*tPrime**2)
        Call ASub(NAt3,Cart2,Scr1,Scr1)
        Call MWeigV(.True.,NAtoms,Scr1,RtMass,Scr1)
        Error = VLen(NAt3,Scr1)
        AThird = One/Three
        TStep1 = TStep*(Del0/Error)**AThird
        If(TStep1.gt.Three) TStep1 = Three
        If(TStep1.lt.Pt025) TStep1 = Pt025
        If(IPrint.ge.0) then
          Write(IOut,1100) TStep
          Write(IOut,1110) TStep1
          endIf
      else
        TStep1 = TStep
        endIf
C
C     Undo the mass-wieghting the coordinates in Cart1 and Cart2.
C
      Call MWeigV(.True.,NAtoms,Cart0,RtMass,Cart0)
      Call MWeigV(.True.,NAtoms,Cart1,RtMass,Cart1)
      Call MWeigV(.True.,NAtoms,Cart2,RtMass,Cart2)
      Return
      End
*Deck EulDWI
      Subroutine EulDWI(IOut,IPrint,IRwFit,ITSStp,ISvLrg,IRwMM,IMMMod,
     $  IDWIp,DoOni,HUpFit,NAtoms,NAt3,NAtQ3,NFtPts,MaskON,IStep,RtMass,
     $  Efit,x1,dx,gfit,IPtEM,IPtGM,StDvEM,StDvGM,CarLrg,Scr1,TotStp,
     $  NSteps,ActStp,TSVec,DwnHil,DoD3Num,D3Mat,D3Vec,LShpOb,ShpObj,V,
     $  MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine carries out the corrector Euler integration on a distance
C     weighted interpolate surface (generated by utility routines).
C
C     When the energy and gradient are needed on the DWI surface,
C     DoOni=.TRUE. indicates that the ONIOM(QM:MM) codes should be used to
C     generate the fitted surface.  DoOni=.FALSE. indicates that no
C     ONIOM(QM:MM) codes should be used in the DWI fitting routines.
C
C     The input argument ITSStp controls how the first step off of the
C     transition structure point is handled. In all cases, the transition
C     vector is used to define the tangent at the transition structure.
C     Allowed values are:
C           1     Run the requested number of Euler steps in the standard
C                 way. Only the first Euler step taken uses the transition
C                 vector.
C           2     Take a large step off of the transition structure point
C                 along the transition vector. This step is taken to be
C                 half of the total requested step size given by TotStp.
C           3     This the same as ITSStp=2 in concept. The only difference
C                 is that the first step off of the transition structure is
C                 taken as one-third the total requested step size given by
C                 TotStp.
C           4     This the same as ITSStp=2 in concept. The only difference
C                 is that the first step off of the transition structure is
C                 taken as one-fourth the total requested step size given
C                 by TotStp.
C
C     Input dummy argument ISvLrg is set to 1 or 2 to indicate whether
C     argument CarLrg is loaded with the geometry corresponding to the
C     largest DWI standard deviation in the energy or gradient.
C
C        - H. P. Hratchian, Spring 2003, 2011.
C
      Logical DoOni,HUpFit,DwnHil,DoD3Num
      Dimension XX(1),MaskON(*),RtMass(NAtoms),x1(NAt3),dx(NAt3),
     $  gFit(NAt3),CarLrg(*),TSVec(NAt3),Scr1(NAt3,7),D3Mat(*),D3Vec(*),
     $  ShpObj(*),V(MDV)
      Save XX,Two,Three,Four
      Data XX/0.d0/,Two/2.d0/,Three/3.d0/,Four/4.d0/
 9000 Format(1x,'EulDWI: ITSStp=',I1)
C
C     Set-up some initial stuff...
C
      IV = 1
      IGrad2 = IV
      IV     = IGrad2 + NAt3
      MDV1 = MDV-IV+1
      Call TstCor(IV-1,MDV,'EulDWI')
      If(ISvLrg.ne.1.and.ISvLrg.ne.2)
     $  Call GauErr('EulDWI: ISvLrg ne 1 or 2.')
      IStart = 1
      Call AMove(NAt3,x1,dx)
C
C     Evaluate the energy and gradient on the fitted surface at the initial
C     point.
C
      If(DoOni) then
        Call OniDWI(IOut,IPrint-1,IRwMM,IMMMod,IRwFit,HUpFit,IDWIp,2,
     $    NAt3,NAtQ3,MaskON,RtMass,dx,Efit,gfit,Scr1(1,5),LShpOb,ShpObj,
     $    V(IV),MDV1)
      else
        Call ShepFt(IOut,IPrint,1,IRwFit,HUpFit,IDWIp,2,NFtPts,NAt3,
     $    Efit,StdDvE,dx,Gfit,XX,StdDvG,DoD3Num,D3Mat,D3Vec,LShpOb,
     $    ShpObj,V(IV),MDV1)
        endIf
C
C     Check to see if we are sitting on the TS. If so, the transition
C     vector is the current path tangent and we therefore overwrite gfit
C     with TSVec. We may also want to take one large Euler step using the
C     transition vector before beginning the normal Euler integration with
C     very small step size. If so, then ActStp accordingly re-set.
C
      TotStX = TotStp
      If(IStep.eq.3.and..not.DwnHil) then
        Call AMove(NAt3,TSVec,gfit)
        If(ITSStp.ge.2) then
          Call AUnit(NAt3,gfit,gfit)
          If(ITSStp.eq.2) then
            TSStp = TotStp/Two
          else if(ITSStp.eq.3) then
            TSStp = TotStp/Three
          else if(ITSStp.eq.4) then
            TSStp = TotStp/Four
          else
            Write(IOut,9000) ITSStp
            Call GauErr('EulDWI: Invalid ITSStp.')
            endIf
          Call ACASB(NAt3,dx,gfit,dx,-TSStp)
          TotStX = TotStp-TSStp
          If(DoOni) then
            Call OniDWI(IOut,IPrint-1,IRwMM,IMMMod,IRwFit,HUpFit,IDWIp,
     $        2,NAt3,NAtQ3,MaskON,RtMass,dx,Efit,gfit,Scr1(1,5),LShpOb,
     $        ShpObj,V(IV),MDV1)
          else
            Call ShepFt(IOut,IPrint,1,IRwFit,HUpFit,IDWIp,2,NFtPts,NAt3,
     $        Efit,StdDvE,dx,Gfit,XX,StdDvG,DoD3Num,D3Mat,D3Vec,LShpOb,
     $        ShpObj,V(IV),MDV1)
            If(StdDvE.gt.StDvEM) then
              IPtEM = 1
              StDvEM = StdDvE
              If(ISvLrg.eq.1) Call AMove(NAt3,dx,CarLrg)
              endIf
            If(StdDvG.gt.StDvGM) then
              IPtGM  = 1
              StDvGM = StdDvG
              If(ISvLrg.eq.2) Call AMove(NAt3,dx,CarLrg)
              endIf
            endIf
          IStart = 2
          endIf
        endIf
      ActStp = TotStX/GFloat(NSteps)
C
C     Carry out Euler integration.
C
      Do 10 i=IStart,NSteps
        Call AUnit(NAt3,gfit,gfit)
        Call ACASB(NAt3,dx,gfit,dx,-ActStp)
        If(DoOni) then
          Call OniDWI(IOut,IPrint-1,IRwMM,IMMMod,IRwFit,HUpFit,IDWIp,2,
     $      NAt3,NAtQ3,MaskON,RtMass,dx,Efit,gfit,Scr1(1,5),LShpOb,
     $      ShpObj,V(IV),MDV1)
        else
          Call ShepFt(IOut,IPrint,1,IRwFit,HUpFit,IDWIp,2,NFtPts,NAt3,
     $      Efit,StdDvE,dx,Gfit,XX,StdDvG,DoD3Num,D3Mat,D3Vec,LShpOb,
     $      ShpObj,V(IV),MDV1)
          If(StdDvE.gt.StDvEM) then
            IPtEM  = i
            StDvEM = StdDvE
            If(ISvLrg.eq.1) Call AMove(NAt3,dx,CarLrg)
            endIf
          If(StdDvG.gt.StDvGM) then
            IPtGM  = i
            StDvGM = StdDvG
            If(ISvLrg.eq.2) Call AMove(NAt3,dx,CarLrg)
            endIf
          endIf
        Call AMove(NAt3,dx,Scr1)
   10   Continue
      Return
      End
*Deck ExplOp
      Subroutine ExplOp(IOut,IOp,First,Restrt,StpSiz,DVVv0,NRdStat,
     $  NCorOp,NAtmQX,CritGd,CritAn,EPS,NAtoms,NAtomQ,IMthd,IEmp,
     $  IHesPrd,IHesCor,IHsUdt,ICorOp,IPFrq,IRxFwd,NPhase,IEig,NStat,
     $  MxOpt,IQMMM,NReduc,IOnDir,Compct,IntHes,NPath,MaxStp,MasFlg,
     $  DXConv,ConvF,SumRlE,SumCnv,SvPES,DoNRTst,WrURVA,IMthd1,NPath1,
     $  IRxFw1,MxStp1)
      Implicit Real*8(A-H,O-Z)
C
C     Using the IOps, this routine loads all of the parameters and options
C     used by Link 123.
C
C     NAtomQ should be set prior to calling this routine as though QM:MM
C     techniques are being used.  After setting flags from the IOps, this
C     routine will re-set NAtomQ to 0 if appropriate.
C
C     This routine has three basic blocks.  First, load internal variables
C     based on the IOps.  Second, set appropriate defaults.  Third, using
C     job flags directly determined from IOps set other remaining job
C     parameters.
C
C     PROGRAMMERS NOTE: When adding a new integrator to l123, it is
C                       essential that one ensures that the variable
C                       IDeriv is set correctly in this routine!
C
C     Dummy Arguments:
C     ----- INPUT ARGUMENTS -----
C     IOut        The standard output file unit number.
C     IOp         The IOp array.
C     First       This logical flag is sent as .True. for the first entry
C                 to Link 123 and the first call to this routine, and
C                 .False. for all other entries to Link 123.  Some IOps
C                 should only be interpreted the first time in Link 123.
C                 First may not be .True. if Restrt is also .True..
C     Restrt      This logical flag is sent as .True. for the first step in
C                 a restart job, and .False. in all other cases.  For the
C                 first step in a restart job MaxStp, NPath, and IRxFwd are
C                 NOT defaulted.  .Flase. for standard runs.  Restrt may
C                 not be .True. if First is also .True..
C     NAtoms      The number of atoms.
C     NAtomQ      The number of QM atoms in QM:MM calculations.
C     IMthd1      The default value that should be used for IMthd in
C                 restart jobs (Restrt=.True.).
C     NPath1      The default value that should be used for NPath in
C                 restart jobs (Restrt=.True.).
C     IRxFw1      The default value that should be used for IRxFwd in
C                 restart jobs (Restrt=.True.).
C     MxStp1      The default value that should be used for MaxStp in
C                 restart jobs (Restrt=.True.).
C
C     ----- OUTPUT (AND IN/OUT) ARGUMENTS -----
C     IMthd       A flag indicating the chosen integration method.  See
C                 comments in Routine RdGeom under IOp(45) for details.
C     IEmp        A flag indicating if this is a calculation on an
C                 empirical surface and translation/rotation projections
C                 should be ignored.
C     IHesPrd     A flag indicating if analytic Hessians are computed or if
C                 Hessian updating is used for predictor steps. NOTE:
C                 Evaluation of the Hessian at the TS is controlled by
C                 IOp10, and control of Hessian evaluation during
C                 corrector/optimization cycles is given by IHesCor.
C     IHesCor     A flag indicating if analytic Hessians are computed or if
C                 Hessian updating is used for corrector steps. NOTE:
C                 Evaluation of the Hessian at the TS is controlled by
C                 IOp10, and control of Hessian evaluation during predictor
C                 steps is given by IHesPrd.
C     IHsUdt      Integer flag indicating the type of updating scheme to
C                 use when updated Hessians are requested.  See comments in
C                 Routine Explor under IOp(13) for details.
C     ICorOp      Option flags for the corrector portion of Predictor-
C                 Corrector IRC integrators. Consult the comments for
C                 IOp(131) in Link 101 for details. Briefly, the elements
C                 of this array are:
C                       ICorOp(1)   Re-correction iteration control flag.
C                       ICorOp(2)   Switch for using numerical third-order
C                                   terms in DWI surface fits duing
C                                   corrector integration cycles.
C                       ICorOp(3)   Power flag used in defining DWI
C                                   weighting functions. The weights will
C                                   go as 1/x**2p, where p is set equal to
C                                   this flag.
C                       ICorOp(4)   This flag controls how the first step
C                                   off of the transition structure is
C                                   handled in Routine EulDWI.
C                       ICorOp(5)   This flag controls whether or not DWI
C                                   fitting should use Hessian update
C                                   vectors when updating has been used. A
C                                   value of 1 indicates that update
C                                   vectors should be used whenever
C                                   possible; a value of 2 supresses the
C                                   use of the update vectors at all steps.
C                                   The default is 1.
C     IPFrq       An integer flag indicating is projected frequencies
C                 should be evaluated.  At present, this option is NOT
C                 fully implemented.
C     IRxFwd      An integer flag indicating whether the IRC is to be
C                 followed in the forward direction (IRxFwd=1), reverse
C                 direction (=2), both directions (=3), or downhill from
C                 the current point (=9).
C     NPhase      A 4-element integer array with phase definitions.
C     IEig        The eigenvalue to follow during coordinate driving jobs.
C     NStat       The number of statistics to be stored in the Stats array.
C     MxOpt       Maximum number of optimization cycles in GS2 reaction
C                 path calulcations, or the maximum number of coorector
C                 integration cycles when convergence is required in
C                 EulerPC and HPC reaction path calculations.
C     IQMMM       Integer option for QM:MM IRC calculations.  IQMMM=0,2 to
C                 ignore all QM:MM codes throughout Link 123.  IQMMM=1
C                 indicates that QM:MM codes should be used when
C                 appropriate, according to the values of NReduc, IOnDir,
C                 Compct, and IntHes.
C     NReduc      The number of eigen-values/vectors to request when
C                 diagonalizing the Hessian in QM:MM IRC jobs.
C     IOnDir      A control flag for using direct ONIOM(QM:MM)
C                 Hessian-vector product codes.
C     Compct      Logical flag indicating if the Hessian is compressed,
C                 which means that only the QM block of the Hessian is
C                 stored in QM:MM calculations.
C     IntHes      Logical flag indicating whether or not the full
C                 integrated Hessian is in hand for ONIOM(QM:MM)
C                 calculations.
C     NPath       The number of paths to follow.
C     MaxStp      The maximum number of steps along each path.
C     MasFlg      This flag indicates how atomic masses are to be loaded.
C                 Note that the atomic masses are NOT loaded in this
C                 routine.  MasFlg=0/1 indicates that the masses in /Mol/
C                 are to be used.  MasFlg=2 indicates that mass-weighting
C                 is turned off and atomic masses should all be set to 1.
C     DXConv      Delta-X convergence threshold for corrector schemes that
C                 test the change in geometry for convergence.
C     ConvF       Force convergence threshold for corrector schemes that
C                 test the force for convergence.
C     SumRlE      Logical flag indicating rather or not energies reported
C                 in the final summary table from Link 123 should be given
C                 relative to the initial TS energy.
C     SumCnv      Logical flag indicating whether or not statistics over
C                 coordinate values should be converted to
C                 Angstroms/Degrees when the final summary table is
C                 printed.
C     SvPES       Logical flag indicating whether data along the IRC should
C                 be saved to the PES data structure file.
C     DoNRTst     Logical flag indicating whether a Newton-Raphson step
C                 test should be carried out when second-order methods are
C                 used.
C     WrURVA      Logical flag indicating whether a URVA input file should
C                 be written or not.
C
C        - H.P. Hratchian, Spring 2003.
C          Significant Modifications: H.P. Hratchian, 2009.
C
      Dimension IOp(*),NPhase(4),ICorOp(NCorOp)
      Logical First,Restrt,HvHess,Compct,IntHes,SumRlE,SumCnv,SvPES,
     $  DoNRTst,WrURVA
      Save Ten,DfStp1,DfStp2,DefV0,DefCGd,DefCAn,DefEPS,IDfEig,IDNRed,
     $  DefDXC
      Data Ten/10.d0/,DfStp1/0.075d0/,DfStp2/0.1d0/,DefV0/0.04d0/,
     $  DefCGd/.00015d0/,DefCAn/30.d0/,DefEPS/1.0d-5/,IDfEig/1/,
     $  IDNRed/1/,DefDXC/0.01d0/
 1000 Format(1x,55('*'),/,' WARNING! HESSIAN CALCULATION',
     $  ' REQUESTED AT EACH STEP.',/,' HESSIAN NOT NEEDED FOR SELECTED',
     $  ' INTEGRATOR!',/,1x,55('*'))
 1010 Format(1x,'WARNING: Corrector options set for an explicit ',
     $  'integrator.',/,'          These options will be ignored.')
 9000 Format(1x,'NPath  =',I3,' NPath1 =',I3,/,1x,'IRxFwd=',I3,
     $  ' IRxFw1=',I3)
C
C     Some quick checks.
C
      If(First.and.Restrt)
     $  Call GauErr('ExplOp is confused...First=T & Restrt=T.')
C
C     Load IOps directly into appropriate job parameter flags.
C
      MxOpt = IOp(6)
      IGS2Cv = IOp(7)
      IHsUdt = IOp(13)
      IOp22 = IOp(22)
      IStpSz = 0
      If(First.or.Restrt) IStpSz = IOp(39)
      IEPS = IOp(40)
      MaxStp = IOp(42) + 1
      IRCTyp = IOp(44)
      If(First.or.Restrt) IMthd = Mod(IOp(45),100)
      IEmp = Mod(IOp(45),1000)/100
      IStpOr = Mod(IOp(45),10000)/1000
      ISmRlE = Mod(IOp(45),100000)/10000
      ISmCnv = Mod(IOp(45),1000000)/100000
      IWrURVA= Mod(IOp(45),10000000)/1000000
      ISvPES = Mod(IOp(45),100000000)/10000000
      INRTst = Mod(IOp(45),1000000000)/100000000
      IHesPrd = IOp(71)
      IHesCor = IOp(162)
      If(IHesCor.lt.0) Call GauErr('Explop: IHesCor<0')
      If(IHesPrd.lt.-1) then
        IHesCor = -IHesPrd
        IHesPrd = 0
        endIf
      IPFrq = IOp(73)
      IOnDir = Mod(IOp(98),10)
      NPhase(1) = IOp(101)
      NPhase(2) = IOp(102)
      NPhase(3) = IOp(103)
      NPhase(4) = IOp(104)
      IDVVv0 = 0
      If(First.or.Restrt) IDVVv0 = IOp(106)
      ICritG = 0
      ICritA = 0
      If(First) then
        ICritG = IOp(108)
        ICritA = IOp(109)
        endIf
      Compct = Mod(IOp(115),10).eq.4
      IEig = IOp(130)
      ICrOp = IOp(131)
      ICrSgn = 1
      If(ICrOp.lt.0) then
        ICrSgn = -1
        ICrOp = -ICrOp
        endIf
      Call IClear(NCorOp,ICorOp)
      ICorOp(1) = ICrSgn*Mod(ICrOp,100)
      ICorOp(2) = Mod(ICrOp,1000)/100
      ICorOp(3) = Mod(ICrOp,10000)/1000
      ICorOp(4) = Mod(ICrOp,100000)/10000
      ICorOp(5) = Mod(ICrOp,1000000)/100000
C
C     Set defaults.  Some defaults are set or not set depending on the
C     values of First and Restrt.
C
      If(MxOpt.le.0) MxOpt = 20
      If(First.or.Restrt) then
        If(IOp22.eq.0.and.First) then
          NPath = 2
          iRxFwd = 1
        else if(IOp22.eq.0.and.Restrt) then
          NPath = NPath1
          iRxFwd = IRxFw1
        else if(IOp22.eq.1) then
          NPath = 1
          iRxFwd = 1
        else if(IOp22.eq.2) then
          NPath = 1
          iRxFwd = 2
        else if(IOp22.eq.3) then
          NPath = 2
          iRxFwd = 1
        else if(IOp22.eq.9) then
          NPath = 1
          iRxFwd = 3
        else
          Call GauErr('L123 Error: Invalid IOp(22).')
          endIf
        endIf
      If(Restrt.and.(NPath.ne.NPath1.or.IRxFwd.ne.IRxFw1)) then
        Write(IOut,9000) NPath,NPath1,IRxFwd,IRxFw1
        Call GauErr('Restart jobs may not change path direction.')
        endIf
      If(MaxStp.le.1.and.Restrt) then
        MaxStp = MxStp1
      else if(MaxStp.le.1) then
        MaxStp = 11
        endIf
      If(Restrt) MaxStp = Max(MaxStp,MxStp1)
      If(IRCTyp.eq.0) IRCTyp = 3
      MasFlg = 0
      If(IRCTyp.eq.1) then
        MasFlg = 2
      else if(IRCTyp.eq.2) then
        Call GauErr('Internal coordinate IRC NYI.')
      else if(IRCTyp.ne.3) then
        Call GauErr('Unrecognized IOp(44) in ExplOp.')
        endIf
      If(MasFlg.eq.0.and.IEmp.eq.1) MasFlg = 2
      If(IOnDir.eq.5.or.IOnDir.eq.6) then
        IQMMM = 1
        NReduc = Min(IDNRed,3*NAtoms)
      else
        IQMMM = 2
        NReduc = 0
        endIf
      If(IMthd.eq.0.and.Restrt) then
        IMthd = IMthd1
      else if(IMthd.eq.0) then
        If(IQMMM.eq.1) then
          IMthd = 6
        else
          IMthd = 3
          endIf
        endIf
      If(IMthd.eq.1.or.IMthd.eq.5.or.IMthd.eq.10) then
        IDeriv = 1
      else
        IDeriv = 2
        endIf
      If((ICorOp(1).ne.0.or.ICorOp(2).ne.0.or.ICorOp(3).ne.0.or.
     $  ICorOp(4).ne.0).and.
     $ .not.(IMthd.eq.3.or.IMthd.eq.6)) Write(IOut,1010)
      If(ICorOp(1).eq.0) then
        If(IMthd.eq.3) ICorOp(1) = 2
        If(IMthd.eq.6) ICorOp(1) = 2
        endIf
      If(ICorOp(2).eq.0) ICorOp(2) = 1
      If(ICorOp(3).eq.0) ICorOp(3) = 1
      If(ICorOp(4).eq.0) ICorOp(4) = 2
      If(ICorOp(5).eq.0) ICorOp(5) = 1
      If(First) then
        If(IStpSz.eq.0) then
          If(IMthd.eq.1.or.IMthd.eq.6) then
            StpSiz = DfStp1
          else
            StpSiz = DfStp2
            endIf
          If(IMthd.eq.5) StpSiz = StpSiz*GFloat(5)
        else
          If(IStpOr.eq.0) IStpOr = 2
          StpSiz = GFloat(IStpSz)*Ten**(-IStpOr)
          endIf
        If(IDVVv0.le.0) then
          DVVv0 = DefV0
        else
          DVVv0 = GFloat(IDVVv0)/GFloat(10000)
          endIf
        If(ICritG.eq.0) then
          CritGd = DefCGd
        else
          CritGd = GFloat(ICritG)/GFloat(10000)
          endIf
        If(ICritA.eq.0) then
          CritAn = DefCAn
        else
          CritAn = GFloat(ICritA)
          endIf
        endIf
      If(ISmRlE.eq.0) ISmRlE = 1
      SumRlE = ISmRlE.eq.1
      If(ISmCnv.eq.0) ISmCnv = 1
      SumCnv = ISmCnv.eq.1
      If(IWrURVA.eq.0) IWrURVA = 2
      WrURVA = IWrURVA.eq.1
      If(ISvPES.eq.0) ISvPES = 1
      SvPES = ISvPES.eq.1
      If(INRTst.eq.0) INRTst = 1
      DoNRTst = INRTst.eq.1
      DefCnF = GeoCon(0)
      DXConv = DefDXC
      If(IGS2Cv.eq.0) then
        ConvF = DefCnF
      else if(IMthd.eq.6) then
        ConvF = GeoCon(IGS2Cv)
        DXConv = GFloat(IGS2Cv)/GFloat(1000000)
      else
        ConvF = GeoCon(IGS2Cv)
        DXConv = DXConv*Sqrt(ConvF/DefCnF)
        endIf
      If(IHsUdt.eq.0) IHsUdt = 3
      If(IHsUdt.lt.0.or.IHsUdt.gt.4)
     $  Call GauErr('L123: IOp(13) value invalid!')
      If(IEPS.eq.0) then
        If(IGS2Cv.ne.0.and.IMthd.eq.6) then
          Eps = DefEPS*GFloat(10)
        else
          EPS = DefEPS*ConvF/DefCnF
          If(IMthd.eq.6) EPS = EPS*GFloat(10)
          endIf
      else if(IEPS.gt.0) then
        EPS = GFloat(1)/(GFloat(10)**IEPS)
      else
        Call GauErr('IOp(1/40) not valid!  IOp(1/40) must be >= 0.')
        endIf
      If(IEig.eq.0) IEig = IDfEig
      If(IEig.lt.0.or.IEig.gt.NAtoms*3)
     $  Call GauErr('L123: IOp(130) value invalid!')
      If(IHesPrd.lt.0.and.IDeriv.gt.1) IHesPrd = 0
      If(IHesPrd.eq.0.and.IDeriv.le.1) IHesPrd = -1
      HvHess = IHesPrd.ne.-1
      If(HvHess.and.IDeriv.le.1.and.IMthd.ne.10) Write(IOut,1000)
      If(IPFrq.eq.0) IPFrq = 1
      If(IPFrq.lt.0.or.IPFrq.gt.2)
     $  Call GauErr('L123: IOp(73) value invalid!')
C
C     Set remaining job flags: NStat, NAtomQ, IntHes.
C
      NStat = 3+NRdStat
      IntHes = .not.Compct
      If(IQMMM.eq.2.or.IntHes) NAtomQ = 0
      NAtmQX = NAtomQ
      Return
      End
*Deck FlNSeq
      Subroutine FlNSeq(IOut,IPrint,ISqTyp,ISt,IFtMax,NSeq)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to fill the array NSeq in ModBS from Link 123.
C     NSeq contains the sequence of the number of steps that each
C     Bulrisch-Stoer iteration should take.  Note that the maximum number
C     of BS iterations is IFtMax-1, NOT IFtMax.  Therefore, the last
C     element of NSeq is filled with 0.
C
C     Multiple sequences are available.  The choice of sequence type is
C     given by the input dummy argument ISqTyp.  Allowed values for ISqTyp
C     are:
C           1     Original Bulirsch-Stoer sequence: NSeq(i) = 2*NSeq(i-2)
C           2     Sequence suggested by Deuflhard: NSeq(i) = 2*i
C           3     Modified version of Deuflhard's sequence (used here, not
C                 in the literature): NSeq(i) = 2*(i+9)
C
C     If the Bulirsch-Stoer sequence is selected (ISqTyp=1) the parameters
C     IFirst and ISecnd are used to initialize the sequence.
C
C     NOTE: IFtMax must be larger than 3.
C
      Parameter (IFirst=2,ISecnd=4,IThird=6)
      Dimension NSeq(IFtMax)
 2000 Format(1x,'NSeq:')
C
      If(IFtMax.lt.3) Call GauErr('FlNSeq: IFtMax is < 3.')
      If(ISqTyp.lt.1.or.ISqTyp.gt.3)
     $  Call GauErr('FlNSeq: Invalid ISqTyp!')
      NSeq(IFtMax) = 0
      If(ISqTyp.eq.1) then
        NSeq(1) = IFirst
        NSeq(2) = ISecnd
        NSeq(3) = IThird
      else if(ISqTyp.eq.2) then
        NSeq(1) = 2
        NSeq(2) = 4
        NSeq(3) = 6
      else if(ISqTyp.eq.3) then
        NSeq(1) = 20
        NSeq(2) = 22
        NSeq(3) = 24
      else
        Call GauErr('ISqTyp out of range in FlNSeq.')
        endIf
      If(ISqTyp.eq.1) then
        Do 10 I = 1, (ISt-1)
          New = 2*NSeq(2)
          NSeq(1) = NSeq(2)
          NSeq(2) = NSeq(3)
   10     NSeq(3) = New
        endIf
      Do 100 i = 4, (IFtMax-1)
        If(ISqTyp.eq.1) then
          NSeq(i) = 2*NSeq(i-2)
        else if(ISqTyp.eq.2) then
          NSeq(i) = 2*i
        else if(ISqTyp.eq.3) then
          NSeq(i) = 2*(i+9)
        else
          Call GauErr('FlNSeq Error 2!')
          endIf
  100   Continue
      If(IPrint.ge.2) then
        Write(IOut,2000)
        Call IMtOut(NSeq,IFtMax-1,1,IFtMax-1,1)
        endIf
      Return
      End
*Deck GetHes
      Subroutine GetHes(IOut,IPrint,IRwF2,NAtoms,IQMMM,IMMCRS,Compct,
     $  Update,NHesTp,Cart1,Cart2,Grad1,Grad2,GradQ1,GradQ2,theta,phi,
     $  vUpdat,wUpdat,MicOpt,MaskON,Hess,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to get the Hessian according to input
C     parameters.
C
C     INPUT arguments:
C           IOut,IPrint The standard output unit number and print level
C                       integer flags.
C           IRwF2       The RWF file number for the Hessian.  If IRwF2 <=
C                       0, then this routine uses the standard RWF file
C                       number, 585.  When Hessian updating is done, this
C                       variable is ignored (see argument Update).
C           NAtoms      This is the number of atoms.
C           IQMMM       This is an integer flag indicate whether or not
C                       QM:MM codes should be used in processing the
C                       Hessian; and, if so, what mode of QM:MM Hessian
C                       building is to be used.  Available options include:
C                       0     Default to option 2.
C                       1     QM:MM processing.  If Compct is .TRUE., then
C                             load only the projected QM contribution to
C                             the Hessian potential processed by Link 120
C                             and written to the RWF in file IMMCRS.
C                       2     No special QM:MM processing.  Just pick-up
C                             the Hessian from the normal place.
C           IMMCRS      This is the RWF file number for the MMCRS file.  If
C                       IMMCRS<=0, then the default RWF file number is used
C                       (715).  If IQMMM.eq.0 or .IQMMM.eq.2 or
C                       Compct=.FALSE., then IMMCRS is ignored.
C           Compct      When QM:MM processing is in use, this logical
C                       variable should be sent as .TRUE. to indicate that
C                       the QM Hessian on RWF file number IMMCRS is stored
C                       in compact form.  If QM:MM processing is not being
C                       used, this variable is ignored.
C           Update      This logical variable should be sent as .TRUE. if
C                       Hessian updating is to be done, or .FALSE. if an
C                       analytical Hessian is to be picked up from the RWF.
C           NHesTp      This is an integer flag indicating which type of
C                       Hessian update is desired.  For definitions of
C                       NHesTp, see Routine HUpdat.  If Update is .FALSE.
C                       NHesTp is completed ignored.
C           Cart1,Cart2 These are the Cartesian coordinates for the
C                       previous point (Cart1) and the current point
C                       (Cart2).  These arrays should be allocated
C                       3*NAtoms.  If Update is .FALSE. the are completely
C                       ignored.
C           Grad1,Grad2 These are the gradients for the previous point
C                       (Grad1) and the current point (Grad2).  These
C                       arrays should be allocated 3*NAtoms.  If Update is
C                       .FALSE. or if (IQMMM.eq.1 .and. Compct), they are
C                       completely ignored.
C           GradQ1,GradQ2     These are the QM component of the gradients
C                       for the previous point (GradQ1) and the current
C                       point (GradQ2).  These arrays should be allocated
C                       NAtQ3.  If Update is .FALSE.  or if (IQMMM.ne.1
C                       .and. .not.Compct), they are completely ignored.
C           vUpdat      This is the v Hessian updating vector. It is
C                       ignored if Update is FALSE. Otherwise, it should be
C                       allocated NAt3 (or NAtQ3 for DoOni cases) length.
C           wUpdat      This is the w Hessian updating vector.  It is
C                       ignored if Update is FALSE. Otherwise, it should be
C                       allocated NAt3 (or NAtQ3 for DoOni cases) length.
C           MicOpt      This is the standard MicOpt from /Mol/.  If
C                       IQMMM.ne.1 or Compct=.FALSE. this array is
C                       completely ignored.
C           MaskON      This is the ONIOM masking array that indicates
C                       which atoms are part of the QM region.  This
C                       argument is ignored for IQMMM.ne.1 .or.
C                       .not.Compct.
C           V,MDV       The standard scratch space, and its size.
C
C     IN/OUT arguments:
C           Hess        This is returned with the Hessian for the current
C                       point.  If updating is requested (via the Update
C                       argument), then Hess should be sent as an INPUT
C                       argument with the previous point's Hessian, which
C                       will be replaced by the new Hessian.
C
C
C        - H.P. Hratchian, 2008.
C
      Parameter (KRwF2=585,KMMCRS=715)
      Dimension Cart1(*),Cart2(*),Grad1(*),Grad2(*),GradQ1(*),GradQ2(*),
     $  vUpdat(*),wUpdat(*),MicOpt(*),MaskON(*),Hess(*),V(*)
      Logical DoOni,Compct,Update
 1000 Format(1x,'Enter Routine GetHes.')
 1010 Format(1x,'GetHes: IQMMM=',I3)
 1020 Format(1x,'NVarQM = ',I7,2x,'NVarMM = ',I7,2x,'NVarT  = ',I7,2x,
     $  'NAt3   = ',I7)
 1030 Format(1x,'LenQMH = ',I12,2x,'NAtQ3T = ',I12)
 1040 Format(1x,'LRWHes = ',I12,2x,'LHess  = ',I12)
C
      If(IPrint.ge.2) Write(IOut,1000)
C
C     A couple of quick checks on the flags sent by the calling routine.
C
      If(IQMMM.lt.0.or.IQMMM.gt.2) then
        Write(IOut,1010) IQMMM
        Call GauErr('GetHes: Invalid IQMMM value!')
        endIf
C
C     Figure out some of the key values and lengths.
C
      DoOni = IQMMM.eq.1.and.Compct
      IV = 1
      MDV0 = MDV
      NAtomR = NAtoms
      NAtomN = NAtoms
      If(DoOni) then
        NAtomQ = ICount(-1,NAtomR,MicOpt)
      else
        NAtomQ = NAtoms
        endIf
      NAt3 = NAtomR*3
      NAt3TT = (NAt3*(NAt3+1))/2
      NAtQ3 = NAtomQ*3
      NAtQ3T = (NAtQ3*(NAtQ3+1))/2
      If(DoOni.and.IMMCRS.gt.0) then
        JMMCRS = IMMCRS
      else if(DoOni) then
        JMMCRS = KMMCRS
      else
        JMMCRS = 0
        endIf
C
C     Get the Hessian from the RWF or update as requested.
C
      If(Update) then
        If(.not.DoOni) then
          IScr   = IV
          IV1    = IScr   + NAt3*2
          Call TstCor(IV1,MDV0,'GetHes-1')
          Call HUpdat(NHesTp,NAt3,Cart1,Cart2,Grad1,Grad2,theta,phi,
     $      vUpdat,wUpdat,Hess,V(IScr))
        else
          ICarQ1 = IV
          ICarQ2 = ICarQ1 + NAtQ3
          IScr   = ICarQ2 + NAtQ3
          IV1    = IScr   + NAtQ3*2
          Call TstCor(IV1,MDV0,'GetHes-2')
          Call CmpVc(1,NAt3,LCarQ,NAtoms,Cart1,V(ICarQ1),MaskON,1)
          If(LCarQ.ne.NAtQ3)
     $      Call GauErr('Length of compressed Cart1 .ne. NAtQ3!')
          Call CmpVc(1,NAt3,LCarQ,NAtoms,Cart2,V(ICarQ2),MaskON,1)
          If(LCarQ.ne.NAtQ3)
     $      Call GauErr('Length of compressed Cart2 .ne. NAtQ3!')
          Call HUpdat(NHesTp,NAtQ3,V(ICarQ1),V(ICarQ2),GradQ1,GradQ2,
     $      theta,phi,vUpdat,wUpdat,Hess,V(IScr))
          endIf
      else
        If(DoOni) then
          If(JMMCRS.le.0) Call GauErr('GetHes, Logical Error 1.')
          Call InMMCRS(IOut,IPrint,JMMCRS,-1,NAtomR,NAtomN,Compct,
     $      MicOpt,IGot,NVarT,NVarMM,NVarQM,NRF,LenQMH,IEQM,IPGrad,
     $      IGrad,IQMHes)
          If(IGot.ne.0) Call GauErr('GetHes: Error from InMMCRS.')
          If(NRF.ne.0) Call GauErr('GetHes: rigid fragments NYI!')
          If(NVarT.ne.NAt3) then
            Write(IOut,1020) NVarQM,NVarMM,NVarT,NAt3
            Call GauErr('GetHes: NVarT.ne.NAt3.')
            endIf
          If(LenQMH.ne.NAtQ3T) then
            Write(IOut,1030) LenQMH,NAtQ3T
            Call GauErr('GetHes: LenQMH.ne.NAt3TT.')
            endIf
          IFilNo = JMMCRS
          IOff = IQMHes
          LHess = LenQMH
        else
          If(IRwF2.gt.0) then
            IFilNo = IRwF2
          else
            IFilNo = KRwF2
            endIf
          IOff = 0
          LHess = NAt3TT
          endIf
        LRWHes = ITqry(IFilNo)
        If(LRWHes.lt.LHess) then
          Write(IOut,1040) LRWHes,LHess
          Call GauErr('GetHes: LRWHes < LHess!')
          endIf
        Call FileIO(2,-IFilNo,LHess,Hess,IOff)
        endIf
      Return
      End
*Deck GS2Cor
      Subroutine GS2Cor(IOut,IPrint,NAtoms,NAt3,NAt3TT,IStep,IStep2,
     $  MxOpt,Done,StpSiz,x1,x2,g1,g2,geigen,H2,HsEVec,HsEVal,TSVec,dx,
     $  xpivot,p,AtMass,RtMass,Scr1,RxnFwd,DwnHil,ConvF)
      Implicit real*8(a-h,o-z)
C
C     This routine carries out the constrained optimizations needed for
C     the second order method of Gonzalez and Schlegel. The primary
C     deviation from the original implementation of GS2 is that L123 carries
C     out the reation path integration in mass-weighted Cartesians rather
C     than mass-weighted internals.
C
C     This method is discussed in detail in:
C          C. Gonzalez and H.B. Schlegel, JCP, 90, 2154 (1989).
C
      Dimension x1(NAt3),x2(NAt3),g1(NAt3),g2(NAt3),geigen(NAt3),
     $  H2(NAt3TT),HsEVec(NAt3,NAt3),HsEVal(NAt3),TSVec(NAt3),dx(NAt3),
     $  xpivot(NAt3),p(NAt3),AtMass(NAtoms),RtMass(NAtoms),IRESLT(3),
     $  Scr1(NAt3,*)
      Logical DEBUG,Done,DwnHil,RxnFwd
      Save Zero,Small,Half,ATenth,One,Two,Four,CNVFMX,CnvLam,MxLStp
      Data Zero/0.d0/,Small/1.d-4/,Half/.5d0/,ATenth/0.1d0/,One/1.d0/,
     $  Two/2.d0/,Four/4.d0/,CNVFMX/1.5d0/,CnvLam/1.d-10/,MxLStp/100/
 1000 Format('         Item               Value     Threshold  ',
     $ 'Pt',I3,' Converged?')
 1002 Format(' Maximum Force            ',F8.6,5X,F8.6,5X,3A1)
 1004 Format(' RMS     Force            ',F8.6,5X,F8.6,5X,3A1)
 1006 Format(' Maximum Displacement     ',F8.6,5X,F8.6,5X,3A1)
 1008 Format(' RMS     Displacement     ',F8.6,5X,F8.6,5X,3A1)
 1080 Format(' Optimization completed.',/,
     $       '    -- Optimized point #',I4,' Found.')
 2000 Format('   Lambda converged in ',I3,' iterations.')
 2110 Format('   Initial lambda convergence testing, step ',I3,/,
     $  '     ALam1=',F15.10,' FLam1=',F15.10)
 2111 Format('   Initial lambda convergence testing, step ',I3,/,
     $  '     ALam2=',F15.10,' FLam2=',F15.10)
 2120 Format(/,'   Solving for lambda...ICount=',I10,/,
     $  '     ALam1=',F15.12,' ALam2=',F15.12,/,
     $  '     ALam=',F15.12,'  FLam=',F15.12)
 2200 Format(' GS2Cor2 found left lambda bracket.  # of steps=',I3,/,
     $  '     ALam1=',F15.10,' FLam1=',F15.10,/)
 2201 Format(' GS2Cor2 found right lambda bracket. # of steps=',I3,/,
     $  '     ALam2=',F15.10,' FLam2=',F15.10,/)
 3000 Format(' Opt Criteria Satisfied, But Only First Step!',/,
     $  4x,'--> Take At Least One More Optimization Step.')
C
      DEBUG = .False.
C
C     Begin by mass-weighting the gradient, Hessian, and coordinates.
C
      Call MWeigV(.False.,NAtoms,x1,RtMass,x1)
      Call MWeigV(.False.,NAtoms,x2,RtMass,x2)
      Call MWeigV(.True.,NAtoms,g1,RtMass,g1)
      Call MWeigV(.True.,NAtoms,g2,RtMass,g2)
      Call MWeigh(.True.,NAtoms,AtMass,H2)
C
C     Solve for the pivot point (xpivot) and form the p vector.
C
      If(IStep.eq.3.and..not.DwnHil) then
        Call AMove(NAt3,TSVec,Scr1)
        If(RxnFwd) Call AScale(NAt3,-One,Scr1,Scr1)
      else
        GMag = VLen(NAt3,g1)
        Call AScale(NAt3,-One/GMag,g1,Scr1)
        endIf
      Call ACASB(NAt3,x1,Scr1,xpivot,Half*StpSiz)
      Call ASub(NAt3,x2,xpivot,p)
C
C     Now, move p and g to the Hessian eigenvector space.
C
      Call MatMP1(1,0,1,NAt3,NAt3,NAt3,1,HsEVec,p,Scr1)
      Call AMove(NAt3,Scr1,p)
      Call MatMP1(1,0,1,NAt3,NAt3,NAt3,1,HsEVec,g2,geigen)
C
C     Begin by setting up the initial bracketing of ALam. ALam1 is the
C     left bracket and ALam2 is the right bracket. FLam1 and FLam2 are the
C     function values are ALam1 and ALam2, respectively. Once we get
C     acceptable initial brackets set up, then we use simple bisection to
C     find the root.
C
C     Do While(FLam1.gt.Zero)
      ICount = 0
      ALam1 = HsEVal(1)
  10  ICount = ICount+1
      If(ICount.gt.MxLStp)
     $  Call GauErr('L123, GS2: Initial lambda (left) bracket failed!')
      ALam1 = ALam1-Max(Abs(HsEVal(1)-HsEVal(2)),ATenth)
      FLam1 = HLamb(NAt3,geigen,p,HsEVal,ALam1,StpSiz)
      If(DEBUG.or.IPrint.ge.3) Write(IOut,2110) ICount,ALam1,FLam1
      If(FLam1.gt.Zero) Goto 10
      If(DEBUG.or.IPrint.ge.2) Write(IOut,2200) ICount,ALam1,FLam1
C
C     Do While(FLam2.lt.Small)
      ICount = 0
      ALam2 = HsEVal(1)-Two*Small
  20  ICount = ICount+1
      If(ICount.gt.MxLStp)
     $  Call GauErr('L123, GS2: Initial lambda (right) bracket failed!')
      ALam2 = ALam2+Half*Abs(HsEVal(1)-ALam2)
      FLam2 = HLamb(NAt3,geigen,p,HsEVal,ALam2,StpSiz)
      If(DEBUG.or.IPrint.ge.3) Write(IOut,2111) ICount,ALam2,FLam2
      If(FLam2.lt.Small) Goto 20
      If(DEBUG.or.IPrint.ge.2) Write(IOut,2201) ICount,ALam2,FLam2
C
C     Now, carry out the full bisection to find Lam.
C
C
C     Do While(Abs(FLam.gt.CnvLam))
      ICount = 0
  50  ICount = ICount+1
      If(ICount.gt.MxLStp) Call GauErr(
     $  'L123, GS2: LAMBDA BRACKETING FAILED!')
      ALam = (ALam1+ALam2)/Two
      FLam = HLamb(NAt3,geigen,p,HsEVal,ALam,StpSiz)
      If(DEBUG.or.IPrint.ge.2) Write(IOut,2120) ICount,ALam1,ALam2,ALam,
     $  FLam
      If(FLam.gt.CnvLam) then
        ALam2 = ALam
        Goto 50
      else if(FLam.lt.-CnvLam) then
        ALam1 = ALam
        Goto 50
      else if(DEBUG.or.IPrint.ge.0) then
        Write(IOut,2000) ICount
        endIf
C
C     Lambda has been found, so solve for dx in the Hessian eigenvector
C     space and rotate back to the Cartesian space.
C
      Do 600 i = 1,NAt3
  600   dx(i) = -(geigen(i)-ALam*p(i))/(HsEVal(i)-ALam)
      Call MatMP2(1,0,1,NAt3,NAt3,NAt3,1,HsEVec,p,Scr1)
      Call AMove(NAt3,Scr1,p)
      Call MatMP2(1,0,1,NAt3,NAt3,NAt3,1,HsEVec,dx,Scr1)
      Call AMove(NAt3,Scr1,dx)
C
C     Test to see if optimization convergence has been met. The projected
C     gradient is placed in Scr1.
C
      temp = SProd(NAt3,g2,p)/SProd(NAt3,p,p)
      Call AScale(NAt3,temp,p,Scr1)
      Call ASub(NAt3,g2,Scr1,Scr1)
      Grms = SProd(NAt3,Scr1,Scr1)
      Grms = Sqrt(Grms/NAt3)
      Gmax = ArrMax(Scr1,NAt3,.True.)
      dxrms = SProd(NAt3,dx,dx)
      dxrms = Sqrt(dxrms/NAt3)
      dxmax = ArrMax(dx,NAt3,.True.)
      CONVX = Four*CONVF
      T1 = CONVF*CNVFMX
      T2 = CONVX*CNVFMX
      If(DEBUG.or.IPrint.ge.0) then
        Write(IOut,1000) IStep-1
        Call Convgd(Gmax,T1,IRESLT)
        Write(IOut,1002) Gmax,T1,IRESLT
        Call Convgd(Grms,CONVF,IRESLT)
        Write(IOut,1004) Grms,CONVF,IRESLT
        Call Convgd(dxmax,T2,IRESLT)
        Write(IOut,1006) dxmax,T2,IRESLT
        Call Convgd(dxrms,CONVX,IRESLT)
        Write(IOut,1008) dxrms,CONVX,IRESLT
        endIf
      Done = (Gmax.le.T1).and.(Grms.le.CONVF).and.
     $  (dxmax.le.T2).and.(dxrms.le.CONVX)
      If(IStep2.gt.MxOpt) Call GauErr(' GS2 Optimization Failure.')
      If(IStep2.eq.1) then
        If(Done) Write(IOut,3000)
        Done = .False.
        endIf
C
C     If the optimization is converged, then we're done and we can return to
C     the main routine(s) to take a predictor step. Otherwise, we need to add
C     the Cartesian displacements to x2.
C
      If(Done) then
        Write(IOut,1080) IStep-1
        Call MWeigV(.True.,NAtoms,x1,RtMass,x1)
        Call MWeigV(.True.,NAtoms,x2,RtMass,x2)
        Call MWeigV(.False.,NAtoms,g1,RtMass,g1)
        Call MWeigV(.False.,NAtoms,g2,RtMass,g2)
        Call MWeigh(.False.,NAtoms,AtMass,H2)
      else
        Call AAdd(NAt3,dx,x2,x2)
        Call MWeigV(.True.,NAtoms,x1,RtMass,x1)
        Call MWeigV(.True.,NAtoms,x2,RtMass,x2)
        Call MWeigV(.False.,NAtoms,g1,RtMass,g1)
        Call MWeigV(.False.,NAtoms,g2,RtMass,g2)
        Call MWeigh(.False.,NAtoms,AtMass,H2)
        endIf
      Return
      End
*Deck InfoRx
      Subroutine InfoRx(Restrt,MaxStp,StpSiz,NPath,iRxFwd,IOP10,IHesPrd,
     $  IHesCor,IMthd,NCorOp,ICorOp,MxOpt,IHsUdt,IEig)
      Implicit Real*8(A-H,O-Z)
C
C     Print out the basic job parameters. Basically this routine writes out
C     the user's directives the first time we come into L123.
C
C        - H.P. Hratchian, Spring 2003.
C
      Common /IO/ In,IOut,IPunch
      Logical Restrt,HvHess
      Character*36 MthTyp(10),Recorr(7),HesUp(4)
      Dimension ICorOp(NCorOp)
      Save MthTyp,Recorr,HesUp,Zero
      Data MthTyp/'Euler','LQA','HPC','GS2','DVV','Euler-PC','UNKNOWN',
     $  'UNKNOWN','UNKNOWN','Coord Drive'/,
     $  Recorr/'Never','Yes','Always','Test Only','DWI Refinement',
     $  'Yes, and Exact PES','Unknown'/,
     $  HesUp/'SR1','PSB','Bofill','SQRT(Bofill)'/,Zero/0.d0/
 1000 Format(1x,72('-'))
 1010 Format(' INPUT DATA FOR L123')
 2000 Format(' GENERAL PARAMETERS:')
 2010 Format(' Follow reaction path in both directions.')
 2012 Format(' Rxn path following direction = ',A)
 2020 Format(' Maximum points per path      =',I4)
 2030 Format(' Step size                    =',f8.3,' sqrt(amu)*bohr')
 2032 Format(' Step size                    =',f8.3,' bohr')
 2034 Format(' Step size                    =',f8.3,' fs')
 3000 Format(' Integration scheme           = ',A)
 3010 Format('    Redo corrector integration= ',A)
 3011 Format('      ==> ICorOp(1) = ',I2)
 3015 Format('    DWI Weight Power          = ',I2)
 3020 Format('    DWI will use numerical third-order information.')
 3025 Format('    DWI will use Hessian update vectors when possible.')
 3026 Format('    DWI will always use full N**2 Hessians.')
 3100 Format('    Max correction cycles     =',I4)
 3500 Format('    Following eigenvector number ',I4)
 4000 Format(' Initial Hessian              = CalcFC')
 4005 Format(' Initial Hessian              = ReadFC from chk')
 4200 Format(' Hessian evaluation           = All updating')
 4201 Format(' Hessian evaluation           = Analytic Hessian always')
 4205 Format(' Hessian evaluation           = Analytic every ',I3,
     $  ' predictor steps')
 4206 Format('                              = Analytic every ',I3,
     $  ' corrector steps')
 4207 Format(' Hessian evaluation           = Update at all predictor',
     $  ' steps')
 4208 Format('                              = Update at all corrector',
     $  ' steps')
 4209 Format(' Hessian evaluation           = Analytic at all',
     $  ' predictor steps')
 4210 Format('                              = Analytic at all',
     $  ' corrector steps')
 4250 Format(' Hessian evaluation           = ',
     $  'Gradient only (Hessian at TS only).')
 4299 Format(' Hessian updating method      = ',A)
C
C     Print out trajectory options
C
      Write(IOut,1000)
      Write(IOut,1010)
      Write(IOut,1000)
C
C     Print out the general parameters.
C
      Write(IOut,2000)
      If(NPath.eq.1) then
        If(iRxFwd.eq.1) write(IOut,2012) 'Forward'
        If(iRxFwd.eq.2) write(IOut,2012) 'Reverse'
        If(iRxFwd.eq.3) write(IOut,2012) 'Downhill'
      else
        write(IOut,2010)
        endIf
      Write(IOut,2020) MaxStp-1
      If(Restrt) then
        Write(IOut,2030) StpSiz
      else if(IMthd.eq.5) then
        Write(IOut,2034) StpSiz
      else
        If(StpSiz.lt.Zero) then
          Write(IOut,2030) -StpSiz
        else
          Write(IOut,2032) StpSiz
          endIf
        endIf
C
C     Print out integration information. If the method is GS2 (IMthd=4),
C     we print out the maximum number of optimization cycles allowed per
C     IRC step. If the method is coordinate driving (IMthd=10) we print
C     out which eigenvector we will follow.
C
      Write(IOut,3000) MthTyp(IMthd)
      If(IMthd.eq.3.or.IMthd.eq.6) then
        If(ICorOp(1).ge.1.and.ICorOp(1).le.3) then
          Write(IOut,3010) Recorr(ICorOp(1))
        else if(ICorOp(1).eq.-1) then
          Write(IOut,3010) Recorr(4)
        else if(ICorOp(1).eq.10) then
          Write(IOut,3010) Recorr(5)
        else if(ICorOp(1).eq.11) then
          Write(IOut,3010) Recorr(6)
        else
          Write(IOut,3010) Recorr(7)
          Write(IOut,3011) ICorOp(1)
          endIf
        Write(IOut,3015) 2*ICorOp(3)
        If(ICorOp(2).eq.2) Write(IOut,3020)
        If(ICorOp(5).eq.1) then
          Write(IOut,3025)
        else if(ICorOp(5).eq.2) then
          Write(IOut,3026)
          endIf
        endIf
      If(IMthd.eq.4.or.((IMthd.eq.3.or.IMthd.eq.6).and.
     $  (ICorOp(1).gt.1))) Write(IOut,3100) MxOpt
      If(IMthd.eq.10) Write(IOut,3500) IEig
C
C     Print out Hessian information.
C
      HvHess = IHesPrd.ne.-1
      If(HvHess) then
        If(IOP10.eq.0.or.IOP10.eq.4) then
          Write(IOut,4000)
        else if(IOP10.eq.5) then
          Write(IOut,4005)
          endIf
        If(IHesPrd.eq.0.and.IHesCor.eq.0) then
          Write(IOut,4200)
        else if(IHesPrd.eq.1.and.IHesCor.eq.1) then
          Write(IOut,4201)
        else
          If(IHesPrd.eq.0) then
            Write(IOut,4207)
          else if(IHesPrd.eq.1) then
            Write(IOut,4209)
          else if(IHesPrd.gt.1) then
            Write(IOut,4205) IHesPrd
          else
            Call GauErr('InfoRx: Invalid IHesPrd found.')
            endIf
          If(IHesCor.eq.0) then
            Write(IOut,4208)
          else if(IHesCor.eq.1) then
            Write(IOut,4210)
          else if(IHesCor.gt.1) then
            Write(IOut,4206) IHesCor
          else
            Call GauErr('InfoRx: Invalid IHesCor found.')
            endIf
          endIf
        If(IHesPrd.ne.1.or.IHesCor.ne.1) Write(IOut,4299) HesUp(IHsUdt)
      else
        Write(IOut,4250)
        endIf
      Write(IOut,1000)
      Return
      End
*Deck Ini123
      Subroutine Ini123(IOut,IPrint,IIRdSt,IMskON,IAtMas,IRtMas,ITsGrd,
     $  ITsGdQ,ITsCrt,ITsVec,ITsFc,ITsEVc,ITsLamb,IIDone,IAGrad,IGradQ,
     $  IHsCar,IHsEVc,IHsEVl,IStats,ICart,IGrdSv,IVel,IGrdCu,IHesCu,
     $  IHVcCu,IHVlCu,IvUpda,IwUpda,IFroze,IScr1,LenArr,IEnd,IBlck1,
     $  IBlck2,IBlck3,IBlck4,IBlck5,IBlck6,LBlck1,LBlck2,LBlck3,LBlck4,
     $  LBlck5,LBlck6,IOfBk1,IOfBk2,IOfBk3,IOfBk4,IOfBk5,IOfBk6,NAtoms,
     $  NAtomQ,LenHes,LnHEVc,LnHEVl,LenSq,MaxStp,MxOpt,NPath,NStat,
     $  NCorOp,ICorOp,DoScr,HvHess,DoHEVc,DwnHil,IFtMax,IQMMM,NReduc,
     $  IntHes,NFroz,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Allocate memory for Link 123 arrays.
C
C     LenArr is the length of the array space that includes permanent
C     arrays that need to be written/read to/from the RWF. IEnd is the
C     total length of the array space needed for the current run of Link
C     123.
C
C     To allow Link 123 to read blocks of arrays as needed (which can be
C     used to avoid N**3 I/O if unnecessary), this routine returns lengths
C     and off-sets for array blocks. Three values are returned for each
C     block. A general nomenclature is used for these values; LBlck* is the
C     length of Block *, IBlck* is the position in memory where Block *
C     should be loaded, and IOfBk* is the off-set for Block * to use when
C     calling FileIO (relative to the start of the array space on the RWF).
C     There are 6 blocks:
C           Block 1: Arrays for data read from the input file and some
C                    other basic information we want to store.
C           Block 2: O(N) arrays for data specific to the TS.
C           Block 3: O(N**2) and related arrays for data specific to the TS.
C           Block 4: O(N) arrays for data generated during the nuclear
C                    propagation.
C           Block 5: O(N**2) arrays for data generated during the nuclear
C                    propagation.
C           Block 6: Arrays that can be re-sized by routine BkAjst.
C
C     If DoScr=.True. then allocate scratch arrays as well and check if
C     there is enough available memory.
C
C     If HvHess=.True. allocate for Hessian based integration else allocate
C     for gradient only integration.
C
C     If DoHEVc=.True. allocate for Hessian eigen-vector/value arrays.
C
C     A note regarding the ordering of arrays below.  In order to
C     facilitate restart jobs where the maximum number of steps might
C     changes, we want to isolate the allocation of permanent arrays that
C     depend on MaxStp to the end of the list. See Routine BkAjst for the
C     list of the arrays that belong to this group.
C
C     PROGRAMMERS NOTE: When adding new arrays to the link, one must
C                       remember to modify the calls to this routine
C                       (Ini123).  Routines in l123 that call Ini123
C                       include Explor and BkAjst.
C
C
C----------------------------------------------------------------------
C Variable       Purpose                             Size
C----------------------------------------------------------------------
C         Save data read from the input file
C
C IRdStat    This is the user-defined stats        (4,NStat-3)
C            definition array.
C
C         Save some basic information
C
C MskOni     Mask array for ONIOM(QM:MM)           (NAtomQ)
C AtMass     Atomic Masses                         (NAtoms)
C RtMass     Square Root of AtMass                 (NAtoms)
C
C         Save information at starting point
C
C TsGrad     TS gradient                           (NAt3)
C TsGrdQ     TS gradient over only QM atoms in     (NAtQ3)
C            QM:MM calculations.
C TsFc       TS force constants                    (LenHes)
C TsCart     TS coords                             (NAt3)
C TsVec      Transition vector                     (NAt3)
C TsEVec     TS Hessian eigenvectors               (LnHEVc)
C TsLamb     TS diagonalized force constants       (LnHEVl)
C
C        Save information from previous cycles
C
C IDone      The no. of steps each path takes      (INTOWP(NPath))
C AGrad      gradient                              (NAt3,2,NPath)
C GradQ      Gradient over only QM atoms in QM:MM  (NAtQ3,2,NPath)
C            calculations.
C HesCar     Cartesian Hessian                     (LenHes,2,NPath)
C HsEVec     Hessian Eigenvector matrix            (LnHEVc,NPath)
C HsEVal     Hessian Eigenvalue vector             (LnHEVl,NPath)
C Stats      Statistics (see RxStat)               (NStat,MaxStp,NPath)
C Cart       Cartesian for each step               (NAt3,MaxStp,NPath)
C GrdSav     Gradient saved from the previous step (NAt3)
C Vel        Velocity array for DVV jobs           (NAt3,MaxStp,NPath)
C
C          Scrach Arrays (DoScr=.True. only)
C
C GrdCur     Gradient for the current point        (NAt3)
C HesCur     Hessian for the current point         (LenHes)
C vUpdat     v vector for Hessian updating         (NAt3)
C wUpdat     w vector for Hessian updating         (NAt3)
C Frozen     Logical array to flag frozen atoms    (NAtoms)
C Scr1       Scratch array for vectors             (NAt3,20)
C DWICar     Scratch array for DWI fitting         (NAt3,2) or (NAtQ3,2)
C DWIE       Scratch array for DWI fitting         (2)
C DWIGrd     Scratch array for DWI fitting         (NAt3,2) or (NAtQ3,2)
C DWIHes     Scratch array for DWI fitting         (LenHes,2)
C----------------------------------------------------------------------
C
      Parameter (MaxLay=4)
      Logical DEBUG,DoScr,HvHess,DoHEVc,DwnHil,IntHes
      Dimension ICorOp(NCorOp)
 1000 Format(' Ini123 Set-up parameters:',/,
     $       ' DoScr  =',L1,2x,' HvHess =',L1,2x,' DoHEVc =',L1,
     $       ' DwnHil =',L1,2x,/,
     $       ' NAtoms =',I12,2x,' NAtomQ =',I12,2x,' LenHes =',I12,/,
     $       ' LnHEVc =',I12,2x,' LnHEVl =',I12,2x,' LenSq  =',I12,/,
     $       ' MaxStp =',I12,2x,' MxOpt  =',I12,2x,' NPath  =',I12,/,
     $       ' NStat  =',I12,2x,' NCorOp =',I12,2x,' ICorOp1=',I12,/,
     $       ' IFtMax =',I12,2x,' NReduc =',I12,2x,' NFroz  =',I12,/)
 1100 Format(' Link 123 Allocation parameters:',/,
     $       ' IIRdSt =',I12,2x,' IMskON =',I12,2x,' IAtMas =',I12,/,
     $       ' IRtMas =',I12,2x,' ITsGrd =',I12,2x,' ITsGdQ =',I12,/,
     $       ' ITsCrt =',I12,2x,' ITsVec =',I12,2x,' ITsFc  =',I12,/,
     $       ' ITsEVc =',I12,2x,' ITsLamb=',I12,2x,' IIDone =',I12,/,
     $       ' IAGrad =',I12,2x,' IGradQ =',I12,2x,' IHsCar =',I12,/,
     $       ' IHsEVc =',I12,2x,' IHsEVl =',I12,2x,' IStats =',I12,/,
     $       ' ICart  =',I12,2x,' IGrdSv =',I12,2x,' IVel   =',I12,/,
     $       ' IGrdCu =',I12,2x,' IHesCu =',I12,2x,' IHEvCu =',I12,/,
     $       ' IHVlCu =',I12,2x,' IvUpda =',I12,2x,' IwUpda =',I12,/,
     $       ' IFroze =',I12,2x,' IScr1  =',I12,2x,' LenArr =',I12,/,
     $       ' IEnd   =',I12,/)
 1200 Format(' Link 123 RWF Array Block parameters:',/,
     $       ' IBlck1 =',I12,2x,' IBlck2 =',I12,2x,' IBlck3 =',I12,/,
     $       ' IBlck4 =',I12,2x,' IBlck5 =',I12,2x,' IBlck6 =',I12,/,
     $       ' LBlck1 =',I12,2x,' LBlck2 =',I12,2x,' LBlck3 =',I12,/,
     $       ' LBlck4 =',I12,2x,' LBlck5 =',I12,/,
     $       ' LBlck6 =',I12,2x,' IOfBk1 =',I12,2x,' IOfBk2 =',I12,/,
     $       ' IOfBk3 =',I12,2x,' IOfBk4 =',I12,2x,' IOfBk5 =',I12,/,
     $       ' IOfBk6 =',I12,/)
C
      DEBUG = .False.
      NAt3   = 3*NAtoms
      NAt3TT = (NAt3*(NAt3+1))/2
      NAtQ3  = 3*NAtomQ
      NAtQ3T = (NAtQ3*(NAtQ3+1))/2
      IGFlag = 1
      If(DwnHil) IGFlag = 0
      IHFlag = 0
      If(HvHess) IHFlag = 1
      IHEVec = 0
      If(DoHEVc) IHEVec = 1
      If(IQMMM.ne.1.or.IntHes) then
        LenHes = NAt3TT
        LnUpda = NAt3
      else
        LenHes = NAtQ3T
        LnUpda = NAtQ3
        endIf
      If(IQMMM.eq.1.and.NReduc.gt.0) then
        LnHEVc = NAt3*NReduc
        LnHEVl = NReduc
        MinCVV = Max(6,4+MaxLay)
        LenSq = Max((NAt3+1)*MinCVV,NReduc*(NAt3-3*NFroz+3))
      else
        LnHEVc = NAt3*NAt3
        LnHEVl = NAt3
        LenSq = NAt3*NAt3
        endIf
C
C     Allocate Space.
C
C     Block 1: Arrays for data read from the input file and some other
C              basic information we want to store...
      INext = 1
      IIRdSt  = INext
      IMskON  = IIRdSt+InToWP(4*(NStat-3))
      If(NAtomQ.gt.0) then
        IAtMas  = IMskON  + InToWP(NAtoms)
      else
        IAtMas  = IMskON
        endIf
      IRtMas  = IAtMas  + NAtoms
      IBlck1 = INext
      LBlck1 = (IRtMas+NAtoms)-IBlck1
      IOfBk1 = INext-1
      INext  = INext + LBlck1
C
C     Block 2: O(N) arrays for data specific to the TS...
      ITsGrd  = INext
      ITsGdQ  = ITsGrd  + NAt3*IGFlag
      ITsCrt  = ITsGdQ  + NAtQ3*IGFlag
      ITsVec  = ITsCrt  + NAt3*IGFlag
      IBlck2 = INext
      LBlck2 = (ITsVec+NAt3)-IBlck2
      IOfBk2 = INext-1
      INext  = INext + LBlck2
C
C     Block 3: O(N**2) and related arrays for data specific to the TS...
      ITsFc   = INext
      ITsEVc  = ITsFc   + LenHes*IGFlag
      ITsLamb = ITsEVc  + LnHEVc*IGFlag
      IBlck3 = INext
      LBlck3 = (ITsLamb+LnHEVl*IGFlag)-IBlck3
      IOfBk3 = INext-1
      INext  = INext + LBlck3
C
C     Block 4: O(N) arrays for data generated during the nuclear propagation...
      IIDone  = INext
      IAGrad  = IIDone  + InToWP(NPath)
      IGradQ  = IAGrad  + 2*NAt3*NPath
      IBlck4 = INext
      LBlck4 = (IGradQ+2*NAtQ3*NPath)-IBlck4
      IOfBk4 = INext-1
      INext  = INext + LBlck4
C
C     Block 5: O(N**2) arrays for data generated during the nuclear
C              propagation...
      IHsCar  = INext
      IHsEVc  = IHsCar  + 2*LenHes*NPath*IHFlag
      IHsEVl  = IHsEVc  + LnHEVc*NPath*IHEVec
      IBlck5 = INext
      LBlck5 = (IHsEVl+LnHEVl*NPath*IHEVec)-IBlck5
      IOfBk5 = INext-1
      INext  = INext + LBlck5
C
C     Block 6: Arrays that can be re-sized by routine BkAjst...
      IStats  = INext
      ICart   = IStats  + NStat*MaxStp*NPath
      IGrdSv  = ICart   + NAt3*MaxStp*NPath
      IVel    = IGrdSv  + NAt3
      IBlck6 = INext
      LBlck6 = (IVel+NAt3*MaxStp*NPath)-IBlck6
      IOfBk6 = INext-1
      INext  = INext + LBlck6
C
      IEnd    = INext
      LenArr  = IEnd-1
C
C     Arrays allocated but NOT saved to the RWF.
      If(DoScr) then
        IGrdCu  = IEnd
        IHesCu  = IGrdCu  + NAt3
        IHVcCu  = IHesCu  + LenHes*IHFlag
        IHVlCu  = IHVcCu  + LnHEVc*IHEVec
        IvUpda  = IHVlCu  + LnHEVl*IHEVec
        IwUpda  = IvUpda  + LnUpda*IHFlag
        IFroze  = IwUpda  + LnUpda*IHFlag
        IScr1   = IFroze  + InToWP(NAtoms)
        IEnd    = IScr1   + 20*NAt3
      else
        IGrdCu  = IEnd
        IHesCu  = IEnd
        IHVcCu  = IEnd
        IHVlCu  = IEnd
        IvUpda  = IEnd
        IwUpda  = IEnd
        IFroze  = IEnd
        IScr1   = IEnd
        endIf
C
C     Print out the pointers
C
      If(DEBUG.or.IPrint.ge.2) Write(IOut,1000) DoScr,HvHess,DoHEVc,
     $  DwnHil,NAtoms,NAtomQ,LenHes,LnHEVc,LnHEVl,LenSq,MaxStp,MxOpt,
     $  NPath,NStat,NCorOp,ICorOp(1),IFtMax,NReduc,NFroz
      If(DEBUG.or.IPrint.ge.3) Write(IOut,1100) IIRdSt,IMskON,IAtMas,
     $  IRtMas,ITsGrd,ITsGdQ,ITsCrt,ITsVec,ITsFc,ITsEVc,ITsLamb,IIDone,
     $  IAGrad,IGradQ,IHsCar,IHsEVc,IHsEVl,IStats,ICart,IGrdSv,IVel,
     $  IGrdCu,IHesCu,IHVcCu,IHVlCu,IvUpda,IwUpda,IFroze,IScr1,LenArr,
     $  IEnd
      If(DEBUG.or.IPrint.ge.2) Write(IOut,1200) IBlck1,IBlck2,IBlck3,
     $  IBlck4,IBlck5,IBlck6,LBlck1,LBlck2,LBlck3,LBlck4,LBlck5,LBlck6,
     $  IOfBk1,IOfBk2,IOfBk3,IOfBk4,IOfBk5,IOfBk6
      Call TstCor(IEnd,MDV,'Ini123')
      Return
      End
*Deck InURVA
      Subroutine InURVA(IOut,IURVA,IPrint,Init,HvHess,LFilNm,IFilNm,
     $  NAtoms,ITStep,RxnCor,Potent,IAn,Cart,Grad,Hess,
CC-INSERT
     $ Tang,Curv,AtMass,IOp169)
CC-INSERT
      Implicit Real*8(A-H,O-Z)
C
C     This routine writes an input file for Elfi Kraka's Unified Reaction
C     Valley Approach program.
C
C     ARGUMENTS:
C           IOut        This is the standard output unit number.
C           IURVA       This is the URVA file unit number.
C           IPrint      This is the standard print level flag.
C           IFilNm      This character input argument is the name of the
C                       URVA file.
C
C     -H. P. Hratchian, 2011.
C
      Logical Debug,Init,HvHess
      Character*4096 FilNam, ChName
CC-INSERT
      Dimension Tang(*),Curv(*),AtMass(*)
      Character*20 FulHesFName,QMHesFName
      Character*4 tmpChar,tmpChar2
      Logical IfRecord,ISExist
      Integer IRec1, IRec2,LinLen,iren
CC-INSERT

      Dimension IFilNm(LFilNm),IAn(NAtoms),Cart(*),Grad(*),Hess(*)
      Save Zero
      Data Zero/0.d0/
 1000 Format(1x,'Enter Routine InURVA.',/,4x,'ITStep=',I5)
 1010 Format(1x,'Reading URVA input file name:')
 2000 Format('BEGIN')
 2010 Format('NZ,NSubs',/,2(I4))
 2020 Format('IAnZ,IZ1,IZ2,IZ3,IZ4,LBl,LAlpha,LBeta')
 2021 Format(8(I4))
 2030 Format('BL,Alpha,Beta')
 2031 Format(3(E22.15))
 2040 Format('NVar',/,I4)
 2050 Format('Values,IntVec,FPVec')
 2051 Format(E22.15,I4,E22.15)
 2060 Format('NoSym',/,'   1')
 3000 Format('NAtom,NAt3,NAt3TT',/,3(I4))
 3010 Format('CC')
 3011 Format(3(E22.15))
 3020 Format('FX_ZMat_Orientation')
 3021 Format(3(E22.15))
 3040 Format('FFX_ZMat_Orientation')
 3041 Format(3(E22.15))
 4000 Format('IPoCou,Energy,XXIRC',/,I4,2(1x,E22.15))
 8000 Format('END')

CC-INSERT
C     
C     Formats for new browsing
C     
 2005 Format('BEGIN (no Hessian)')
 2015 Format('Natoms,NatomQ',/,2(I4))
 3044 Format('Tangent vector eta (mass-weighted)')

 3045 Format('Curvature vector kappa (mass-weighted)')
 3054 Format('Number of curvature decomposition parameters: ',I4)
 3055 Format('Number of atoms involved in the decomposition:',I4)
 3056 Format(4(I4,2x))
 3057 Format('Atomic masses needed for the',
     $       ' decomposition')
 3058 Format(3(E22.15))

CC-INSERT

C
C     Set some flags and set-up the the URVA file.
C
      Debug = .False.
      If(Debug.or.IPrint.ge.2) Write(IOut,1000) ITStep
      If(.not.HvHess)
     $  Call GauErr('InURVA: HvHess is FALSE! URVA requires Hessians.')
      NAt3 = 3*NAtoms
      NAt3TT = (NAt3*(NAt3+1))/2
      If(Init) then
        IDump = 0
        Call FFGtNm(IDump,.True.,FilNam,LName,ChName,LChNam)
        If(LChNam.gt.(LFilNm-1)*NChrPW(0))
     $    Call GauErr('Name buffer too small in INURVA.')
        IFilNm(1) = 0
        Call PutCSt(ChName(1:LChNam),IFilNm(2),IFilNm(1))
        If(Debug.or.IPrint.ge.1) then
          Write(IOut,1010)
          Call StrOut(IOut,IFilNm(2),IFilNm(1),1)
          endIf
      else
        ICur = 0
        Call HolToC(IFilNm(2),ICur,IFilNm(1),ChName)
        LChNam = IFilNm(1)
        Call OpnApd(IURVA,ChName(1:LChNam))
CC-INSERT
        IF (IOp169.EQ.0) THEN
CC-INSERT
C
C       Write the basic structural data, including the z-matrix
C       representation of the geometry. For now, we always report
C       Cartesian coordinates.
C
        Write(IURVA,2000)
        Write(IURVA,2010) NAtoms,NAt3
        Write(IURVA,2020)
        Do 100 I = 1,NAtoms
          Ixyz = 3*(I-1)+1
          Write(IURVA,2021) IAn(I),0,0,0,0,Ixyz,Ixyz+1,Ixyz+2
  100     Continue
        Write(IURVA,2030)
        Write(IURVA,2031) (Cart(I),I=1,NAt3)
        Write(IURVA,2040) NAt3
        Write(IURVA,2050)
        Do 110 I = 1,NAt3
          Write(IURVA,2051) Cart(I),0,Zero
  110     Continue
        Write(IURVA,2060)
C
C       Write the Cartesian coordinates and the energy derivatives in
C       the z-matrix frame (note that z-matrix frame is always the same
C       as the Cartesian in the current Link 123 IRC code).
        Write(IURVA,3000) NAtoms,NAt3,NAt3TT
        Write(IURVA,3010)
        Write(IURVA,3011) (Cart(I),I=1,NAt3)
        Write(IURVA,3020)
        Write(IURVA,3021) (-Grad(I),I=1,NAt3)
        Write(IURVA,3040)
        Write(IURVA,3041) (Hess(I),I=1,NAt3TT)
C
C       Write the energy, END statement, and close the file.
        Write(IURVA,4000) ITStep,Potent,Abs(RxnCor)
        Write(IURVA,8000)
CC-INSERT     
      ENDIF
CC-INSERT     
CC-INSERT     
C     
C     Write new browsing file
C     
      If(IOp169.gt.9) Itmp = MOD(IOp169,10)
      IF((IOp169.eq.1).or.(Itmp.eq.1)) THEN
      Write(IURVA,2005)
      Write(IURVA,2015) NAtoms,NAtoms
C     
C     get atomic masses
C     
      Mode=0
      Do 5001, kk=1,NAtoms
      IZZ = IAn(kk)
      AtMass(kk) = EleMas(Mode,IZZ)
 5001 Continue
      Write(IURVA,3057)
      Write(IURVA,3058) (AtMass(kk),kk=1,NAtoms)
      Write(IURVA,3010)
      Write(IURVA,3011) (Cart(I),I=1,NAt3)
      Write(IURVA,3044)
      Write(IURVA,3021) (Tang(I),I=1,NAt3)
      Write(IURVA,3045)
      Write(IURVA,3041) (Curv(I),I=1,NAt3)

      Write(IURVA,4000) ITStep,Potent,Abs(RxnCor)
      Write(IURVA,8000)

      ENDIF
CC-INSERT     
      Close(Unit=IURVA)
CC-INSERTED     
      If(IOp169.gt.9) Itmp = MOD(IOp169,10)
      IF((IOp169.gt.10).or.(Itmp.eq.1))THEN

      IURVAs = IURVA + 3
      Open(Unit=IURVAs,File="URVA_s.1",STATUS="UNKNOWN")
      Write(IURVAs,9000) Abs(RxnCor)
      Write(IURVAs,*) ITStep
      Write(IOut,*) "URVA_s.1 written."
 9000 Format(E22.15)
      Close(Unit=IURVAs)

      IURVAs = IURVAs + 3
      Open(Unit=IURVAs,File="URVA_c.0",STATUS="OLD",Err=5002)
      If(1.EQ.2) then
 5002 Call GauErr("Unable to open URVA_c.0 file.")
 5003 Call GauErr("Unabel to read URVA_c.0 file.")
      endIf
      Read(Unit=IURVAs,*,Err=5003) IRec1
      Read(Unit=IURVAs,*,Err=5003) IRec2
      If(IRec2.le.0) Call GauErr("Invalid IRec2 value.")
      If(IRec1.gt.IRec2) Call GauErr("Invalid IRec1 value.")
      Write(IOut,*) "URVA_c.0 loaded."
      Close(Unit=IURVAs)
      If((IRec1.eq.0).or.(IRec1.eq.IRec2)) then
      IfRecord=.True.
      else
      IfRecord=.false.
      endIf

      Write(tmpChar,'(I0)') ITStep
C     Write(IOut,*) "YT-"//tmpChar//"-"
      LinLen = LinEnd(tmpChar)
      If(LinLen.eq.1) tmpChar2 = "000"//tmpChar(1:1)
      If(LinLen.eq.2) tmpChar2 = "00"//tmpChar(1:1)//tmpChar(2:2)
      If(LinLen.eq.3) tmpChar2 = "0"//tmpChar(1:1)//
     $ tmpChar(2:2)//tmpChar(3:3)
      If(LinLen.eq.4) tmpChar2 = tmpChar

C     Write(IOut,*) "YT- Hessian base name is"//tmpChar2//"."
      If(IfRecord) then
      iren = rename("URVA_q.1",
     $                 "URVA_qh-"//tmpChar2//".2")
      INQUIRE(File="URVA_qh-"//tmpChar2//".2",Exist=ISExist)
      If(ISExist) then
      Write(IOut,*) "QM Hessian written to "//
     $              "URVA_qh-"//tmpChar2//".2"
      else

      Call GauErr("Failed to rename URVA_q.1")
      endIf
      endIf
      If(IfRecord) then
      iren = rename("URVA_q.2",
     $                 "URVA_fh-"//tmpChar2//".2")
      INQUIRE(File="URVA_fh-"//tmpChar2//".2",Exist=ISExist)
      If(ISExist) then
      Write(IOut,*) "Full Hessian written to "//
     $              "URVA_fh-"//tmpChar2//".2"
      else
      Call GauErr("Failed to rename URVA_q.2")
      endIf
      endIf

      If((IRec1.eq.0).or.(IRec1.eq.IRec2)) then
      IRec1 = 1
      else
      IRec1 = IRec1 + 1
      endIf

      IURVAs = IURVAs + 2
      Open(Unit=IURVAs,File="URVA_c.0",STATUS="UNKNOWN")
      Write(Unit=IURVAs,*) IRec1
      Write(Unit=IURVAs,*) IRec2
      Close(Unit=IURVAs)
      Write(IOut,*) "URVA_c.0 updated."



      ENDIF
CC-INSERTED
        endIf
      Return
      End
*Deck JobTst
      Subroutine JobTst(IOut,NAtoms,NAt3,IMthd,IStep,MaxStp,StpSiz,
     $  RtMass,Cart1,Cart2,Cart3,AGrad,Scr1,CritGd,CritAn,CorDon,Done)
      Implicit Real*8(a-h,o-z)
C
C     This routine tests to see if the job is finished yet.
C
C        - H.P. Hratchian, Spring 2003.
C
      Logical CorDon,Done
      Dimension RtMass(NAtoms),Cart1(NAt3),Cart2(NAt3),Cart3(NAt3),
     $  AGrad(NAt3,2),Scr1(NAt3,2)
      Save Zero
      Data Zero/0.d0/
 1000 Format(' PES minimum detected on this side of the pathway.')
 1010 Format('    Magnitude of the gradient = ',f14.7)
 1020 Format('    Step angle = ',f8.3)
 1030 Format('    Step is smaller than ',F14.7,' a.u.')
 2000 Format(' Maximum number of steps reached.')
C
      Done = .False.
      If(.not.CorDon) Return
C
C     Check to see if the gradient stopping criteria is met.  If this is a
C     coordinate driving job (IMthd.eq.10) we skip this test.
C
      Call MWeigV(.True.,NAtoms,AGrad,RtMass,Scr1)
      GMag = VLen(NAt3,Scr1)
      If(GMag.lt.CritGd.and.CritGd.gt.Zero.and.IMthd.ne.10) then
        Write(IOut,1000)
        Write(IOut,1010) GMag
        Done = .True.
        Return
        endIf
      If(IStep.eq.3) then
        Done = IStep.gt.MaxStp
        If(Done) Write(IOut,2000)
        Return
        endIf
C
C     Check to see if the step angle stopping criteria is met.  If this is
C     a coordinate driving job (IMthd.eq.10) we skip this test.
C
      Call ASub(NAt3,Cart1,Cart2,Scr1(1,1))
      Call ASub(NAt3,Cart3,Cart2,Scr1(1,2))
      Call MWeigV(.False.,NAtoms,Scr1,RtMass,Scr1)
      Call MWeigV(.False.,NAtoms,Scr1(1,2),RtMass,Scr1(1,2))
      StpAng = AngDeg(0,NAt3,Scr1(1,1),Scr1(1,2))
      StpLen = VLen(NAt3,Scr1(1,2))
      SmlStp = StpSiz/GFloat(10)
C
C     Check if we've essentially turned around, as though we stepped over
C     the minimum and are doubling back.
      If(StpAng.lt.CritAn.and.CritAn.gt.Zero.and.IMthd.ne.10) then
        Write(IOut,1000)
        Write(IOut,1020) StpAng
        Done = .True.
C
C     Check we've taken a very small step, which indicates that we have
C     stepped over the minimum and are turning back. NOTE: This test is
C     based on the value of StpSiz, which is in units of time (not
C     distance) for DVV calculations. Therefore, this test is skipped for
C     DVV jobs.
      else if(StpLen.lt.SmlStp.and.IMthd.ne.5) then
        Write(IOut,1000)
        Write(IOut,1030) SmlStp
        Done = .True.
C
C     Check if we've taken too many steps.
      else if(IStep.gt.MaxStp) then
        Write(IOut,2000)
        Done = .True.
        endIf
      Return
      End
*Deck LQA
      Subroutine LQA(IOut,IPrint,Cart1,Cart2,Grad,HsEVec,HsEVal,StpSz,
     $  ActStp,NAt3,NAtoms,RtMass,IStep,RxnFwd,X,GPrime,AGrad,A,Alpha,
     $  Scr1,ILQAFg,DwnHil)
      Implicit real*8(a-h,o-z)
C
C     Integrate the reaction path equation on a quadratic potential
C     surface across a distance of StpSiz.  The actual step length is
C     ActStp.  This method is the LQA method of Page and McIver and
C     can be found in:
C          M. Page and J. McIver, JCP, 88, 922 (1988).
C     For an excellent and succinct treatment of LQA see also:
C          M. Page, C. Doubleday, and J. McIver, JCP, 93, 5634 (1990)
C
C        - H.P. Hratchian, Summer 2001.
C
      Dimension Cart1(NAt3),Cart2(NAt3),HsEVal(NAt3),HsEVec(NAt3,NAt3),
     $  Grad(NAt3),RtMass(NAtoms),X(NAt3),GPrime(NAt3),AGrad(NAt3),
     $  A(NAt3,NAt3),Alpha(NAt3,NAt3),Scr1(NAt3,NAt3)
      Logical DEBUG,GradEx,RxnFwd,TisInf,DwnHil
      Save Zero,Half,Six,One,Epslon,NEuler,Two,MxIter,MxCycl,TenTho
      Data Zero/0.d0/,Half/.5d0/,Six/6.d0/,One/1.d0/,Epslon/1.d-6/,
     $  NEuler/5000/,Two/2.d0/,MxIter/50/,MxCycl/100000/,TenTho/1.d4/
 1000 Format(3x,'Gradient is small and 1 (-) eigenvalue.',/,
     $  ' Taking an Euler step.')
 1010 Format(3x,'Warning:  LQA convergence difficult!  Criteria ',
     $  'decreased to 1.d-3.')
 1020 Format(3x,'LQA: T_Est iteration FAILED!!!',/,
     $  ' LQA - ActStp = ',F15.6)
 1030 Format(3x,'LQA: T_Est iteration completed in',I5,' iterations.')
 2000 Format(3x,'LQA: Resetting T_Est step size to ',F15.6)
C
C     Preliminary bookkeeping.
C
      DEBUG = .False.
      ILQAFg = 0
      StpSiz = StpSz
      GradEx = .False.
      GMag = VLen(NAt3,Grad)
C
C     If requested, print some of the input data.
C
      If(DEBUG.or.IPrint.ge.3) then
        Call OutMtS(IOut,'LQA: Cartesian coordinates',0,0,Cart1,NAt3,1,
     $    NAt3,1)
        Call OutMtS(IOut,'LQA: Cartesian gradient',0,0,Grad,NAt3,1,
     $    NAt3,1)
        Call OutMtS(IOut,'LQA: Hessian eigenvectors',0,0,HsEVec,NAt3,
     $    NAt3,NAt3,NAt3)
        endIf
C
C     Now, decide if we need to take a simple Euler step or if
C     we can use LQA (i.e. are we at the TS or not?), and take step.
C
      If(IStep.eq.2.and..not.DwnHil) then
        If(DEBUG.or.IPrint.ge.2) Write(IOut,1000)
        Call AMove(NAt3,HsEVec,AGrad)
        GradEx = .True.
        Call MWeigV(.False.,NAtoms,Cart1,RtMass,X)
        If(RxnFwd.and.IStep.eq.2) Call AScale(NAt3,-One,AGrad,AGrad)
        GMag = VLen(NAt3,AGrad)
        Call ACASB(NAt3,X,AGrad,X,StpSiz/GMag)
        ActStp = StpSiz
      else
        Call MWeigV(.True.,NAtoms,Grad,RtMass,AGrad)
        ASixth = One/Six
        Call MWeigV(.False.,NAtoms,Cart1,RtMass,X)
        Call AClear(NAt3*NAt3,Alpha)
        Call AClear(NAt3*NAt3,A)
C
C       Use StpSiz and the AGrad (the gradient) to determine an estimate
C       for T, TEst.  Use this value to integrate dS/dT (see Page and
C       McIver references above) to calculate ActSiz.  There is a minor
C       hack to the value of StpDif to ensure that this iterative block
C       runs through at least two cycles.
C
        GDist = VLen(NAt3,AGrad)
        TEst = StpSiz/GDist
        If(.not.GradEx) then
          Call MatMP1(1,0,1,NAt3,NAt3,NAt3,1,HsEVec,AGrad,GPrime)
        else
          Call AMove(NAt3,AGrad,GPrime)
          endIf
        Do 10 J=1,MxIter
          ActStp = Zero
          TTemp = Zero
          Do 20 M = 1,MxCycl
            dSSum = Zero
            Do 30 N = 1,NAt3
              ExTerm = -Two*GFloat(M)*HsEVal(N)*TEst/GFloat(NEuler)
              Sum = (GPrime(N)**2)*Exp(ExTerm)
   30         dSSum = Sum+dSSum
            dSSum = (Sqrt(dSSum))*TEst/GFloat(NEuler)
            OldStp = ActStp
            ActStp = ActStp+dSSum
            If(ActStp.ge.StpSiz) goto 5
   20       TTemp = TTemp+(TEst/NEuler)
    5     StpDif = ActStp-StpSiz
          If(J.eq.1) StpDif = Epslon+Epslon
          If(StpDif.le.Epslon.or.TTemp.gt.TenTho) goto 15
          Scal = (StpSiz-OldStp)/(ActStp-OldStp)
          TEst = TTemp-TEst/NEuler+Scal*TEst/NEuler
          If(DEBUG.or.IPrint.ge.2) Write(IOut,2000)TEst
   10     continue
        If(StpDif.le.1.d-3) then
          If(DEBUG.or.IPrint.ge.1) Write(IOut,1010)
          Goto 15
          endIf
        If(DEBUG.or.IPrint.ge.1) Write(IOut,1020) ActStp
        ILQAFg = 1
        Return
   15   If(DEBUG.or.IPrint.ge.1) Write(IOut,1030) J
        TEst = TTemp
        TisInf = .False.
C
C       Determine the matrices (Alpha and A) needed to propagate the
C       nuclear coordinates for the next step prediction.
C
        If(TisInf) then
          Do 40 N = 1,NAt3
   40       Alpha(N,N) = -One/HsEVal(N)
        else
          do 50 N = 1,NAt3
            If(Abs(HsEVal(N)*TEst).lt.Epslon) then
              Alpha(N,N) = -TEst+Half*HsEVal(N)*TEst**2-ASixth*HsEVal(N)
     $          **2*TEst**3
            else
              Alpha(N,N) = (exp(-TEst*HsEVal(N))-One)/HsEVal(N)
              endIf
   50       continue
          endIf
        Call MatMpB(1,0,1,NAt3,NAt3,NAt3,NAt3,Alpha,HsEVec,Scr1)
        Call MatMP2(1,0,1,NAt3,NAt3,NAt3,NAt3,HsEVec,Scr1,A)
C
C       Now that we have everything we need, propagate the nuclear
C       coordinates
C
        Call MatMP2(1,1,1,NAt3,NAt3,NAt3,1,A,AGrad,X)
        endIf
C
C     Some final bookkeeping
C
      Call MWeigV(.True.,NAtoms,X,RtMass,Cart2)
      Return
      End
*Deck ModBS
      Subroutine ModBS(IOut,IPrint,IRwFit,ITSStp,ISvLrg,IRwMM,IMMMod,
     $  IDWIp,DoOni,HUpFit,NAtoms,NAt3,NAtQ3,MaskON,IStep,RtMass,Efit,
     $  x1,dx,gfit,Scr1,StpSiz,EPS,TSVec,IPtEM,IPtGM,StDvEM,StDvGM,
     $  CarLrg,YZ,dY,BSFail,DwnHil,IFtSt,IFtMax,NFtPts,DoD3Num,D3Mat,
     $  D3Vec,LShpOb,ShpObj,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine carries out modified Bulirsch-Stoer integration on a
C     fitted potential energy surface. The integration here is 'modified'
C     in that Bulirsch-Stoer typically uses a modified-midpoint
C     integration, but here we use simple Euler integration.
C
C     Input dummy argument ISvLrg is set to 1 or 2 to indicate whether
C     argument CarLrg is loaded with the geometry corresponding to the
C     largest DWI standard deviation in the energy or gradient.
C
C        - H. P. Hratchian, 2003, 2008-11.
C
      Logical DoOni,HUpFit,BSFail,DwnHil,DoD3Num
      Dimension XX(1),MaskON(*),RtMass(NAtoms),x1(NAt3),dx(NAt3),
     $  gfit(NAt3),TSVec(NAt3),CarLrg(*),YZ(NAt3),dY(NAt3),
     $  Scr1(NAt3,10),D3Mat(*),D3Vec(*),ShpObj(*),V(MDV)
      Save XX,One,NEulSt
      Data XX/0.d0/,One /1.d0/,NEulSt/0/
 1000 Format(1x,'Modified Bulirsch-Stoer Extrapolation Cycles:',/,4x
     $  ,'EPS = ',F20.15)
 1030 Format(4x,'Cycle',I6,' NS=',I6)
 2000 Format(4x,'Cycle',I6,' NS=',I6,' Truncation Error = ',F12.9)
 2010 Format(1x,'Efit after extrapolation:',F12.6)
 2020 Format(1x,'BS completed in',I5,' cycles and',I8,
     $  ' integration steps.')
 3000 Format(1x,'Maximum DWI energy   std dev = ',F12.9,' at pt ',I5/,
     $  1x,'Maximum DWI gradient std dev = ',F12.9,' at pt ',I5)
 9000 Format(1x,'WARNING: Bulirsch-Stoer Method is not Converging')
C
C     Initialize things and carry out some basic error traps.
C
      IV = 1
      MDV1 = MDV
      If(ISvLrg.ne.1.and.ISvLrg.ne.2)
     $  Call GauErr('ModBS: ISvLrg ne 1 or 2.')
C
C     Allocate space in V.
C
      IXest  = IV
      INSeq  = IXest  + IFtMax
      IPEZa  = INSeq  + InToWP(IFtMax)
      IPEZb  = IPEZa  + NAt3*(IFtMax-1)
      IV     = IPEZb  + NAt3
      Call TstCor(IV-1,MDV,'ModBS')
C
C     Try BS integration with extrapolation first. If extrapolation fails,
C     then look for convergence of the simple Euler integration with
C     decreasing step sizes.
C
      Call TStamp(1,'Top ModBS')
      If(IPrint.ge.0) Write(IOut,1000) EPS
      Call FlNSeq(IOut,IPrint,1,IFtSt,IFtMax,V(INSeq))
      Do 100 I=1,(IFtMax-1)
        Call TStamp(1,'ModBS100')
        NSeq = IIGet(V(INSeq),I)
        JXest = IXest+I-1
        Call EulDWI(IOut,IPrint,IRwFit,ITSStp,ISvLrg,IRwMM,IMMMod,IDWIp,
     $    DoOni,HUpFit,NAtoms,NAt3,NAtQ3,NFtPts,MaskON,IStep,RtMass,
     $    Efit,x1,dx,gfit,IPtEM,IPtGM,StDvEM,StDvGM,CarLrg,Scr1,StpSiz,
     $    NSEQ,V(JXest),TSVec,DwnHil,DoD3Num,D3Mat,D3Vec,LShpOb,ShpObj,
     $    V(IV),MDV1)
        NEulSt = NEulSt + NSeq
        Call PEZero(I,V(IXEST),dx,yz,dy,V(IPEZa),V(IPEZb),NAt3)
C
C       Check truncation error to see if we are converged.
C
        If(I.lt.2) then
          If(IPrint.ge.1) Write(IOut,1030) I,NSeq
        else
          ErrMx1 = ArrMax(dy,NAt3,.True.)
          If(IPrint.ge.1) Write(IOut,2000) I,NSeq,ErrMx1
          ErrMax = ErrMx1/EPS
          If(ErrMax.lt.One) then
            Call AMove(NAt3,yz,dx)
C
C           Get the new gradient and energy at this extrapolated point.
C
            If(DoOni) then
              Call OniDWI(IOut,0,IRwMM,IMMMod,IRwFit,HUpFit,IDWIp,2,
     $          NAt3,NAtQ3,MaskON,RtMass,dx,Efit,gfit,Scr1(1,5),LShpOb,
     $          ShpObj,V(IV),MDV1)
            else
              Call ShepFt(IOut,IPrint,1,IRwFit,HUpFit,IDWIp,2,NFtPts,
     $          NAt3,Efit,StdDvE,dx,Gfit,XX,StdDvG,DoD3Num,D3Mat,D3Vec,
     $          LShpOb,ShpObj,V(IV),MDV1)
              If(StdDvE.gt.StDvEM) then
                IPtEM = -1
                StDvEM = StdDvE
                If(ISvLrg.eq.1) Call AMove(NAt3,dx,CarLrg)
                endIf
              If(StdDvG.gt.StDvGM) then
                IPtGM  = -1
                StDvGM = StdDvG
                If(ISvLrg.eq.2) Call AMove(NAt3,dx,CarLrg)
                endIf
              endIf
            If(IPrint.ge.1) write(IOut,2020) I,NEulSt
            If(IPrint.ge.2) then
              Write(IOut,2010) Efit
              Call DmpFrc(IOut,'gfit:',1,NAtoms,gfit)
              endIf
            Write(IOut,3000) StDvEM,IPtEM,StDvGM,IPtGM
            Return
            endIf
          endIf
  100   Continue
      Write(IOut,3000) StDvEM,IPtEM,StDvGM,IPtGM
C
C     Maximum number of steps exceeded.
C
      Write(IOut,9000)
      BSFail = .True.
      Call TStamp(1,'Bot ModBS')
      Return
      End
*Deck NRtest
      Subroutine NRtest(IOut,IPrint,ILQAFg,NAt3,NAtoms,NFroz,CritGd,
     $  Done,Cart1,Cart2,dx,Grad,GPrime,HsEVec,HsEVal,StpSiz,StpIRC,
     $  RtMass,Scr1)
      Implicit Real*8(a-h,o-z)
C
C     This routine is used to take a Newton-Raphson step to test if the IRC
C     integrator has stepped beyond the minimum.  If the NR step is shorter
C     than the last IRC step, the NR step will replace the IRC step.  If
C     CritGd is greater than zero, then the job completion flag (Done) will
C     be returned as True.
C
C     ILQAFg is a flag expected from Routine LQA.  If it is sent here as 0,
C     then the LQA step before calling NRTest ran normally.  Therefore, if
C     NRTest predicts a smaller step than LQA did there's a numerical
C     problem and the job should be stopped.
C
C        - H.P. Hratchian, Summer 2004, 2009.
C
      Logical Done
      Dimension Cart1(NAt3),Cart2(NAt3),dx(NAt3),Grad(NAt3),
     $  GPrime(NAt3),HsEVec(NAt3,*),HsEVal(*),RtMass(NAtoms),Scr1(NAt3)
      Save zero,small
      Data zero/0.d0/,small/1.d-5/
 1000 Format(1x,'NRTest: Step size=',F8.4,' Act step=',F8.4,' NR step=',
     $  F8.4)
 5000 Format(1x,'PES minimum detected within a second-order step.')
 5100 Format(1x,'NR Step is smaller than IRC step ==> NR Step',
     $  ' Replacing IRC Step!')
C
C     Begin by setting Done to .False. and mass-weigh Cart1 and Grad.
C
      Done = .False.
      Call MWeigV(.False.,NAtoms,Cart1,RtMass,Cart1)
      Call MWeigV(.True.,NAtoms,Grad,RtMass,Grad)
C
C     Take a NR step in the eigenvector space.  To do this, begin by
C     converting Grad to the eigenvector space.
C
      Call MatMP1(1,0,1,NAt3,NAt3,NAt3,1,HsEVec,Grad,GPrime)
      Call AClear(NAt3,dx)
      If(NFroz.gt.0) then
        NVar = NAt3-3*NFroz
      else
        NVar = NAt3
        endIf
      Do 10 i = 1, NVar
   10   If(Abs(HsEVal(i)).gt.small) dx(i) = -GPrime(i)/HsEVal(i)
      Call MatMP2(1,0,1,NAt3,NAt3,NAt3,1,HsEVec,dx,Scr1)
      Call AMove(NAt3,Scr1,dx)
      xNRMag = VLen(NAt3,dx)
      If(IPrint.ge.2) write(IOut,1000) StpSiz,StpIRC,xNRMag
C
C     Undo the mass-weighting of Cart1, Grad, and dx.
C
      Call MWeigV(.True.,NAtoms,Cart1,RtMass,Cart1)
      Call MWeigV(.True.,NAtoms,dx,RtMass,dx)
      Call MWeigV(.False.,NAtoms,Grad,RtMass,Grad)
C
C     If xNRMag is less than the IRC step size (StpIRC), then re-set Cart2
C     to the NR step.
C
      If(ILQAFg.eq.0) then
        TstStp = StpIRC
      else
        TstStp = StpSiz
        endIf
      If(xNRMag.lt.TstStp) then
        If(CritGd.gt.zero) then
          Write(IOut,5000)
          Done = .True.
        else
          Call AAdd(NAt3,Cart1,dx,Cart2)
          If(IPrint.ge.1) write(IOut,5100)
          endif
      else if(ILQAFg.eq.1) then
        Call GauErr('LQA Failure.')
        endIf
      Return
      End
*Deck NumD3
      Subroutine NumD3(IOut,N,x1,x2,H1,H2,D3Mat,D3Vec)
      Implicit Real*8(A-H,O-Z)
C
C     This routine generates the matrix and vector used to describe a
C     one-dimensional numerical third derivative that is evaluated from
C     changes in coordinates and Hessians. x1 and x2 are the input
C     coordinate vectors; H1 and H2 are the input Hessian matrices. D3Mat
C     and D3Vec are the output matrix and vector representations of the
C     numerical third derivative tensor.
C
      Dimension x1(N),x2(N),H1(*),H2(*),D3Mat(*),D3Vec(N)
      Save One
      Data One/1.d0/
 1000 Format(1x,'Enter NumD3.')
C
      Write(IOut,1000)
      NTT = (N*(N+1))/2
      Call ASub(NTT,H1,H2,D3Mat)
      Call ASub(N,x1,x2,D3Vec)
      Temp = One/SProd(N,D3Vec,D3Vec)
      Call AScale(N,Temp,D3Vec,D3Vec)
      Return
      End
*Deck OniDWI
      Subroutine OniDWI(IOut,IPrint,IRwMM,IMMMod,IRwFit,HUpFit,IDWIp,
     $  NPts,NAt3,NAtQ3,MaskON,RtMass,Z,Efit,Gfit,Scr,LShpOb,ShpObj,V,
     $  MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This is a head routine for carrying out a general distance weighted
C     interpolant surface for ONIOM(QM:MM) cases where the QM portion of
C     the potential surface is fitted while the MM contribution is
C     analytically evaluated using the utility MM routines.
C
C     The final (QM:MM) energy and gradient is found at the point Z and
C     returned in Efit and Gfit.  This routine currently only works in
C     Cartesian coordinates. The data for the MH sub-system should be
C     stored in the standard Shepard interpolation data structure, which
C     can be prepared using Routine PESShp.
C
C     A listing of the arguments in the call to this subroutine is given
C     below.  Array dimensions are also shown.
C
C     IOut        The standard output unit number.
C     IPrint      The standard print option flag.
C     IRwMM       The RWF file number where the MM parameter data is kept.
C     IMMMod      The RWF file number where the ONIOM model system MM data
C                 structure is stored.
C     NPts        This is the number of points to be used in the fit.
C     NAt3        This is the dimension of the coordinate space.
C     NAtQ3       This is 3 times the number of QM atoms.
C     MaskON      This is the ONIOM masking array that indicates how
C                 active-QM-only atoms are compressed.
C     RtMass      This is the array of the atomic root-masses. (NAtoms)
C     Z           This is the point at which Efit and Gfit need to be
C                 computed. (NAt3)
C     Efit        This is the fitted energy at Z.
C     Gfit        This is the fitted gradient at Z. (NAt3)
C     Scr         This is another scratch array. (NAt3,3)
C     V,MDV       The standard scratch space, and its size.
C
C        - H. P. Hratchian, 2008, 2009.
C
      Logical DEBUG,HUpFit
#include "commonmol.inc"
      Dimension XX(1),JJ(1),MaskON(*),RtMass(NAtoms),Z(NAt3),Gfit(NAt3),
     $  Scr(NAt3,3),ShpObj(*),V(MDV)
      Save Zero,XX,JJ
      Data Zero/0.d0/,XX/0.d0/,JJ/0/
 1000 Format(1x,50('-'),/,' ENTERING SUBROUTINE OniDWI')
 9000 Format(1x,'Efit = ',f14.8)
 9999 Format(1x,' LEAVING SUBROUTINE OniDWI')
C
      DEBUG = .False.
      If(DEBUG.or.IPrint.ge.2) write(IOut,1000)
      IV = 1
      MDV1 = MDV-IV+1
C
C     Do the distance-weighted interpolant surface fitting.
C
      Call AClear(NAt3,gfit)
      Call CmpVc(1,NAt3,LZqm,NAtoms,Z,Scr,MaskON,1)
      If(LZqm.ne.NAtQ3) Call GauErr('OniDWI: LZqm.ne.NAtQ3 !')
      Call ShepFt(IOut,IPrint,1,IRwFit,HUpFit,IDWIp,2,NPts,NAtQ3,EMH,
     $  StdDvE,Scr,Gfit,XX,StdDvG,.False.,XX,XX,LShpOb,ShpObj,V(IV),
     $  MDV1)
      Call UCmpVc(1,NAtQ3,NAtoms,gfit,gfit,MaskON)
      Call MWeigV(.True.,NAtoms,Z,RtMass,Z)
C
C     Run the MM evaluations.
C
      JobMM = 1
      IHVAcc = 0
      EpsMM = Zero
      MaxLEF = -1
      LMMPrm = ITqry(IRwMM)
      LModDt = ITqry(IMMMod)
      IDip   = IV
      IGrdRL = IDip   + 3
      IGrdML = IGrdRL + NAt3
      IMMPrm = IGrdML + NAt3
      IModDt = IMMPrm + LMMPrm
      IV1    = IModDt + LModDt
      MDV1   = MDV-IV1+1
      Call TstCor(IV1-1,MDV,'OniDWI')
      Call FileIO(2,-IRwMM,LMMPrm,V(IMMPrm),0)
      Call FileIO(2,-IMMMOD,LModDt,V(IModDt),0)
      ERL = Zero
      Call AClear(NAt3,V(IGrdRL))
      Call EGHMM(IOut,IPrint,0,JobMM,IHVAcc,EpsMM,NAtoms,IAtTyp,AtChMM,
     $  Z,MaxLEF,XX,V(IMMPrm),LMMPrm,JJ,ERL,V(IGrdRL),XX,V(IDip),XX,
     $  V(IV1),MDV1)
      Call AClear(NAt3,V(IGrdML))
      Call EGHMod(IOut,IPrint,0,EpsMM,NAtoms,MicOpt,Z,EML,V(IGrdML),XX,
     $  V(IModDt),V(IModDt),V(IV1),MDV1)
C
C     Finalize things.
C
      Efit = ERL-EML+EMH
      Call MWeigV(.True.,NAtoms,V(IGrdML),RtMass,V(IGrdML))
      Call MWeigV(.True.,NAtoms,V(IGrdRL),RtMass,V(IGrdRL))
      Call ASub(NAt3,gfit,V(IGrdML),gfit)
      Call AAdd(NAt3,V(IGrdRL),gfit,gfit)
      Call MWeigV(.False.,NAtoms,Z,RtMass,Z)
      If(DEBUG.or.IPrint.ge.3) then
        Write(IOut,9000) Efit
        Write(IOut,9999)
        endIf
      Return
      End
*Deck OnMask
      Subroutine OnMask(NAtoms,MskFrz,FilMsk,MicOpt,NAct,Mask)
      Implicit real*8(a-h,o-z)
C
C     This routine is used to count the number of active QM atoms in an
C     ONIOM(QM:MM) system, and to possibly build a masking array of 0s and
C     1s where 1s indicate active atomic centers.
C
C     Input dummy arguments:
C     NAtoms      The number of (real) atoms.
C     MskFrz      A logical flag that indicates whether (.TRUE.) or not
C                 (.FALSE.) frozen atoms in the QM region should be masked
C                 as active.
C     FilMsk      Whether (.TRUE.) or not (.FALSE.) the array Mask should
C                 be filled.  If not, then only NAct is set and returned.
C     MicOpt      The standard MicOpt array from Common /Mol/.
C
C     Output dummy arguments:
C     NAct        This is the number of active atoms.
C     Mask        This is the mask array, which is filled only when
C                 FilMsk=.True.  Mask is dimensioned (NAtoms).  If
C                 FilMsk=.False., Mask is entirely ignored.
C
C
C        - H.P. Hratchian, 2008.
C
      Dimension MicOpt(NAtoms),Mask(*)
      Logical MskFrz,MicFrz,FilMsk,MicQM
C
      NAct = 0
      Do 100 i = 1, NAtoms
        If((MskFrz.and.MicQM(MicOpt(i)).and..not.MicFrz(MicOpt(i))).or.
     $    (.not.MskFrz.and.MicQM(MicOpt(i)))) then
          If(FilMsk) Mask(i) = 1
          NAct = NAct+1
        else if(FilMsk) then
          Mask(i) = 0
          endIf
  100   Continue
      Return
      End
*Deck OTIRC1
      Subroutine OTIRC1(IOut,IMode,NAt3,MaxStp,NPath,ITStp,Cart,Energy,
     $  Force,RxCord,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Store trajectory data during an IRC.
C
C     The operation of this routine is controled by INPUT argument IMode.
C     Allowed values of IMode include:
C          -1     Resize the OT bucket.
C           1     Write new data to the OT bucket.
C
C     This routine is a slightly modified version of Routine OTIRC, which
C     is part of Link 115.  This modified version is intended for use in
C     Link 123.
C
C        - H.P. Hratchian.
C
#include "commonmol.inc"
      Parameter (NStatD=2)
      Dimension V(MDV), Stats(2), XX(1), TRotT(3,3)
      Dimension Cart(NAt3),Force(NAt3)
      Common /SymInf/ NOp1,NOp2,JTrans(3,8),TRot(3,3),TrVec(3)
      Save XX
      Data XX/0.0d0/
 1100 Format(1x,'OTIRC1 is resizing the OT bucket.')
 9000 Format(1x,'OTIRC1: IMode=',I6)
C
C     Set some key flags.
C
      MaxJb1 = 1
      NStat = NStatD
      IJOff1 = 0
      IType1 = 3
      IJob1 = 1
      MaxSt1 = Max(2500,MaxStp*NPath)
      ITStp1 = 0
C
C     If IMode is -1, just call RdWrOT to resize the OT bucket and then
C     return to the calling program unit.
C
      If(IMode.eq.-1) then
        Write(IOut,1100)
        Call RdWrOT(MaxSt1,MaxJb1,NStat,NAt3,IJOff1,IType1,ITStp1,IJob1,
     $    IDone,.True.,Stats,.True.,C,.False.,XX,.True.,XX,IOut,0,-1,V,
     $    MDV)
        Return
      else if(IMode.eq.1) then
        If(ITStp.eq.1) then
          Call RdWrOT(MaxSt1,MaxJb1,NStat,NAt3,IJOff1,IType1,ITStp1,
     $      IJob1,IDone,.True.,Stats,.True.,C,.False.,XX,.True.,XX,IOut,
     $      0,0,V,MDV)
          endIf
        Stats(1) = Energy
        Stats(2) = RxCord
        ICStOr = 1
        IFX = ICStOr + NAt3
        IV = IFX + NAt3
        Call TstCor(IV,MDV,'OTIRC')
        Call ILSW(2,26,NoSym)
        If(NoSym.eq.0) then
          Call AMove(9,TRot,TRotT)
          Call Trspn2(3,TRotT)
          Call RotF1(NAtoms,TRotT,Force)
          endIf
        ITStp1 = ITStp
        IDone = ITStp1
        Call RdWrOT(MaxSt1,MaxJb1,NStat,NAt3,IJOff1,IType1,ITStp1,IJob1,
     $    IDone,.True.,Stats,.True.,Cart,.False.,V,.True.,Force,IOut,0,
     $    2,V,MDV)
      else
        Write(IOut,9000) IMode
        Call GauErr('OTIRC1: Invalid IMode.')
        endIf
      Return
      End
*Deck PreDWI
      Subroutine PreDWI(IOut,IPrint,IRwFit,ICurWh,NCorOp,ICorOp,ICorFl,
     $  IRwMM,IMMMod,IQMMM,Compct,DoHEVc,NAtoms,NAt3,NAtQ3,NAtQ3T,
     $  LenHes,LnHEVc,LnHEVl,MaskON,xCur,x1,x2,gCur,grad1,grad2,HCur,
     $  HVcCur,HVlCur,H1,H2,H2EVec,H2EVal,Update,thetaUp,phiUp,vUpdat,
     $  wUpdat,ECur,E1,E2,gradQ1,gradQ2,EQM1,EQM2,Scr,dx,gFit,IStep,
     $  IStep2,MxOpt,NAt3TT,EPS,AtMass,RtMass,dxConv,CorDon,TSVec,
     $  DwnHil,IFtSt,IFtMax,CarLrg,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine carries out the corrector integration phase of
C     algorithms using modified-Bulrisch-Stoer integration on a fitted
C     distance weighted interpolant surface.
C
      Dimension ICorOp(NCorOp),MaskON(*),xCur(NAt3),x1(NAt3),x2(NAt3),
     $  dx(NAt3),gCur(NAt3),grad1(NAt3),grad2(NAt3),gFit(NAt3),
     $  scr(NAt3,15),HCur(*),HVcCur(*),HVlCur(*),H1(*),H2(*),H2EVec(*),
     $  H2EVal(*),vUpdat(*),wUpdat(*),gradQ1(*),gradQ2(*),TSVec(NAt3),
     $  AtMass(NAtoms),RtMass(NAtoms),CarLrg(*),V(MDV)
      Logical DEBUG,DoD3Num,HUpFit,Update,CorDon,DoOni,Flip,BSFail,
     $  Compct,DoHEVc,DwnHil
      Save Zero,One
      Data Zero/0.d0/,One/1.d0/
 1000 Format(3x,'Using distance weighted interpolants to fit PES.')
 1020 Format(' SUMMARY OF CORRECTOR INTEGRATION:')
 1030 Format('   Predictor End Point Energy = ',F12.6)
 1035 Format('   Old End Point Energy       = ',F12.6)
 1040 Format('   Corrected End Point Energy = ',F12.6)
 1050 Format('   Predictor End-Start Dist.  = ',F12.6)
 1060 Format('   Old End-Start Dist.        = ',F12.6)
 1070 Format('   New End-Start Dist.        = ',F12.6)
 1080 Format('   New End-Old End Dist.      = ',F12.6)
 1090 Format('   Largest Error-Start Dist.  = ',F12.6)
 1091 Format('   Largest Error-End Dist.    = ',F12.6)
 1100 Format(1x,'PreDWI: IReCr0=',I2,' IReCr1=',I2,' IReCr2=',I2,
     $  ' IReCr3=',I2)
 1200 Format(1x,'At Start of PreDWI: ECur=',F12.6,' E1=',F12.6,' E2=',
     $  F12.6)
 8000 Format(' CORRECTOR INTEGRATION CONVERGENCE:',/,
     $  '   Recorrection delta-x convergence threshold:',F12.6)
 8010 Format('   Delta-x Convergence Met')
 8011 Format('   Delta-x Convergence NOT Met')
 8020 Format(' Recorrect forced on => Unset convergence flag.')
 9100 Format(1x,'PreDWI: IReCr1=',I2)
C
C     Tell the user that we're starting the corrector routine and take care
C     of some basic initialization. This includes interpretation of
C     ICorOp, which yields a set of control flags:
C           IReCr0      This is an integer flag indicating whether
C                       recorrection is being used. The currently allowed
C                       values and their meanings are:
C                             1 ... No recorrection.
C                             2 ... Do recorrection if convergence has not
C                                   been met.
C                             3 ... Always do at least one re-correction
C                                   cycle, even if convergence has been met
C                                   in the first macro-cycle.
C           IReCr1      This is an integer flag indicating how the final
C                       structure and PES data determine on the Shepard
C                       surface should be used by the rest of the link.
C                             1 ... The final (converged) structure and
C                                   fitted surface data should be used for
C                                   the next predictor integration step.
C                             2 ... Set the necessary flags so that Link
C                                   123 evaluates the exact PES energy and
C                                   gradient at the final (converged)
C                                   structure for the purpose of comparing
C                                   the final results on the Shepard
C                                   surface to the true PES. The final
C                                   integrated IRC should be the same as if
C                                   this flag were set to 1.
C                             3 ... Same as 2, except that the exact PES
C                                   energy and gradient are used for the
C                                   next predictor integration step.
C           IReCr2      This is an integer flag which, if recorrection is
C                       being used, indicates how convergence is
C                       determined. Note that this flag is ignored if
C                       recorrection is not being used. This flag works in
C                       a bit-fashion. Therefore, one can request
C                       convergence testing over multiple criteria by
C                       combining bits. For example, if one wants
C                       convergence to be tested using options values 1 and
C                       4, the flag to pass here is 5. The currently
C                       allowed values and their meanings are:
C                             1 ... Converge on the change in the corrector
C                                   integration end-point geometry.
C                             2 ... Converge on the largest standard
C                                   deviation of the gradient encountered
C                                   on the Shepard surface during corrector
C                                   integration in Routine EulDWI.
C           IReCr3      This is an integer flag which, if recorrection is
C                       being used, indicates if iterative refinement of
C                       the Shepard surface is being used; and if so, what
C                       iterative refinement scheme is to be used.
C                             1 ... No iterative refinement is being used.
C                             2 ... Use iterative refinement.
C           IDWIp       This is an integer that defines the exponent p used
C                       in the DWI surface fitting routine. See Routine
C                       ShepFt for details.
C           DoD3Num     This is a logical flag indicating whether numerical
C                       third derivatives are to be computed with Routine
C                       NumD3 for use in the DWI surface fitting procedure.
C           ITSStp      This is an integer flag used by Routine EulDWI to
C                       control how the first step off the TS is done on
C                       the fitted surface. See Routine EulDWI for details.
C           HUpFit      This is a logical flag indicating whether Hessian
C                       update-vector based algebra is to be used in the
C                       DWI fitting procedure.
C
      DEBUG = .False.
      If(IPrint.ge.1.or.DEBUG) Write(IOut,1000)
      IV = 1
      MDV1 = MDV-IV+1
      If(IStep2.eq.1) ICorFl = 0
      DoOni = IQMMM.eq.1.and.Compct
      If(ICorOp(6).gt.0) then
        IReCr0 = Mod(ICorOp(6),10)
        IReCr1 = Mod(ICorOp(6),100)/10
        IReCr2 = Mod(ICorOp(6),1000)/100
        IReCr3 = Mod(ICorOp(6),10000)/1000
      else
        If(ICorOp(1).eq.-1) then
          IReCr0 = 1
          IReCr1 = 2
          IReCr2 = 1
          IReCr3 = 1
        else if(ICorOp(1).eq.1) then
          IReCr0 = 1
          IReCr1 = 1
          IReCr2 = 1
          IReCr3 = 1
        else if(ICorOp(1).eq.2) then
          IReCr0 = 2
          IReCr1 = 1
          IReCr2 = 1
          IReCr3 = 1
        else if(ICorOp(1).eq.3) then
          IReCr0 = 3
          IReCr1 = 1
          IReCr2 = 1
          IReCr3 = 1
        else if(ICorOp(1).eq.10) then
          IReCr0 = 2
          IReCr1 = 1
          IReCr2 = 2
          IReCr3 = 2
        else if(ICorOp(1).eq.-11) then
          IReCr0 = 1
          IReCr1 = 3
          IReCr2 = 1
          IReCr3 = 1
        else
          Call GauErr('PreDWI: Invalid ICorOp(1) value.')
          endIf
        endIf
      If(DEBUG.or.IPrint.ge.1) Write(IOut,1100) IReCr0,IReCr1,IReCr2,
     $  IReCr3
      If(IReCr0.lt.1.or.IReCr0.gt.3.or.
     $  IReCr1.lt.1.or.IReCr1.gt.3.or.
     $  IReCr2.lt.1.or.IReCr2.gt.2.or.
     $  IReCr3.lt.1.or.IReCr3.gt.2) then
        Write(IOut,1100) IReCr0,IReCr1,IReCr2,IReCr3
        Call GauErr('PreDWI: Invlaid IReCr found.')
        endIf
      If(ICorOp(2).ne.1.and.ICorOp(2).ne.2)
     $  Call GauErr('L123: Invalid ICorOp(2).')
      DoD3Num = ICorOp(2).eq.2
      IDWIp = ICorOp(3)
      ITSStp = ICorOp(4)
      HUpFit = Update.and.ICorOp(5).eq.1
      If(DEBUG.or.IPrint.ge.2) then
        Write(IOut,1200) ECur,E1,E2
        Call OutMtS(IOut,'At Start of PreDWI: xCur',0,0,xCur,NAt3,1,
     $    NAt3,1)
        Call OutMtS(IOut,'At Start of PreDWI: x1',0,0,x1,NAt3,1,NAt3,1)
        Call OutMtS(IOut,'At Start of PreDWI: x2',0,0,x2,NAt3,1,NAt3,1)
        Call OutMtS(IOut,'At Start of PreDWI: gCur',0,0,gCur,NAt3,1,
     $    NAt3,1)
        Call OutMtS(IOut,'At Start of PreDWI: grad1',0,0,grad1,NAt3,1,
     $    NAt3,1)
        Call OutMtS(IOut,'At Start of PreDWI: grad2',0,0,grad2,NAt3,1,
     $    NAt3,1)
        endIf
C
C     Check to see if we converged the corrector integration in the
C     previous macro-cycle. If so, then do what IReCr1 requires of us, set
C     CorDon to TRUE, and exit.
C
      ICurWh = 0
      If(IStep2.gt.1.and.ICorFl.eq.1) then
        If(IReCr1.eq.1) then
          Call GauErr('PreDWI: IReCr1=1 and ICorFl=1!')
        else if(IReCr1.eq.2) then
          CorDon = .True.
        else if(IReCr1.eq.3) then
          CorDon = .True.
          Call AMove(NAt3,xCur,x2)
          Call AMove(NAt3,gCur,grad2)
          Call AMove(LenHes,HCur,H2)
          If(DoHEVc) then
            Call AMove(LnHEVc,HVcCur,H2EVec)
            Call AMove(LnHEVl,HVlCur,H2EVal)
            endIf
        else
          Write(IOut,9100) IReCr1
          Call GauErr('PreDWI: IReCr1 is invalid.')
          endIf
        Return
        endIf
C
C     Figure out the total integration step size to take on the fitted
C     surface in Routine EulDWI.
C
      Call TstCor(IV+NAt3-1,MDV,'PreDWI 1')
      Call ASub(NAt3,x1,x2,V(IV))
      Call MWeigV(.False.,NAtoms,V(IV),RtMass,V(IV))
      StpSiz = VLen(NAt3,V(IV))
C
C     Take case of mass-weighting.
C
      If(DoOni) then
        If(IReCr3.ne.1)
     $    Call GauErr('PreDWI: Iter. refinement w/QM:MM NYI.')
        IAtMas = IV
        IRtMas = IAtMas + NAtoms
        IV     = IRtMas + NAtoms
        Call TstCor(IV-1,MDV,'PreDWI 2')
        Call CmpVc(0,NAtoms,LAtMas,NAtoms,AtMass,V(IAtMas),MaskON,1)
        If(LAtMas*3.ne.NAtQ3) Call GauErr('PreDWI: LAtMas*3.ne.NAtQ3!')
        Call CmpVc(0,NAtoms,LAtMas,NAtoms,RtMass,V(IRtMas),MaskON,1)
        If(LAtMas*3.ne.NAtQ3) Call GauErr('PreDWI: LAtMas*3.ne.NAtQ3!')
        NAtomQ = NAtQ3/3
        Call MWeigV(.False.,NAtoms,xCur,RtMass,xCur)
        Call MWeigV(.False.,NAtoms,x1,RtMass,x1)
        Call MWeigV(.False.,NAtoms,x2,RtMass,x2)
        Call MWeigV(.True.,NAtoms,grad1,RtMass,grad1)
        Call MWeigV(.True.,NAtoms,grad2,RtMass,grad2)
        Call MWeigV(.True.,NAtomQ,gradQ1,V(IRtMas),gradQ1)
        Call MWeigV(.True.,NAtomQ,gradQ2,V(IRtMas),gradQ2)
        Call MWeigh(.True.,NAtomQ,V(IAtMas),H1)
        Call MWeigh(.True.,NAtomQ,V(IAtmas),H2)
        If(HUpFit) then
          Call MWeigV(.True.,NAtomQ,vUpdat,V(IRtMas),vUpdat)
          Call MWeigV(.True.,NAtomQ,wUpdat,V(IRtMas),wUpdat)
          endIf
      else
        Call MWeigV(.False.,NAtoms,xCur,RtMass,xCur)
        Call MWeigV(.False.,NAtoms,x1,RtMass,x1)
        Call MWeigV(.False.,NAtoms,x2,RtMass,x2)
        Call MWeigV(.True.,NAtoms,gCur,RtMass,gCur)
        Call MWeigV(.True.,NAtoms,grad1,RtMass,grad1)
        Call MWeigV(.True.,NAtoms,grad2,RtMass,grad2)
        Call MWeigh(.True.,NAtoms,AtMass,HCur)
        Call MWeigh(.True.,NAtoms,AtMass,H1)
        Call MWeigh(.True.,NAtoms,Atmass,H2)
        If(HUpFit) then
          Call MWeigV(.True.,NAtoms,vUpdat,RtMass,vUpdat)
          Call MWeigV(.True.,NAtoms,wUpdat,RtMass,wUpdat)
          endIf
        endIf
      Call ASUnit(NAt3,x2,x1,Scr)
      If(IStep.eq.3.and..not.DwnHil) then
        Call AUnit(NAt3,TSVec,Scr(1,2))
      else
        Call AUnit(NAt3,grad1,Scr(1,2))
        endIf
      temp = SProd(NAt3,Scr,Scr(1,2))
      If(IStep.le.3.and.temp.gt.Zero) then
        Flip = .True.
        Call AScale(NAt3,-One,grad1,grad1)
        Call AScale(NAt3,-One,TSVec,TSVec)
      else
        Flip = .False.
        endIf
C
C     Set-up the data needed by the surface fitting routines used during
C     the modified-Bulrisch-Stoer corrector integration in Routines ModBS
C     and EulDWI.
C
      Call TStamp(1,'DWI Set-up')
      If(IReCr3.eq.2) then
        IFitOp = 2
        NFtPts = IStep2+1
        NDtPts = MxOpt+1
        IFtPt  = IStep2+1
      else
        IFitOp = 1
        NFtPts = 2
        NDtPts = 2
        IFtPt  = 0
        endIf
      Call Shp123(IOut,IPrint,IFitOp,IRwFit,IStep2,NAtoms,NAt3,DoOni,
     $  NAtQ3,NAtQ3T,IFtPt,NFtPts,NDtPts,MaskON,xCur,x1,x2,ECur,E1,E2,
     $  gCur,grad1,grad2,HCur,H1,H2,HUpFit,thetaUp,phiUp,vUpdat,wUpdat,
     $  EQM1,EQM2,gradQ1,gradQ2,LShpOb,V(IV),MDV1)
      IShpOb = IV
      IV     = IShpOb + LShpOb
      MDV1 = MDV-IV+1
C
C     If requested, form the numerical third derivatives on the line from
C     x1 to x2. Note that the numerical third derivatives are evaluated
C     AFTER mass-weighting has been completed.
C
      If(DoD3Num.and.DoOni)
     $  Call GauErr('PreDWI: DoD3Num.and.DoOni NYI.')
      If(DoD3Num) then
        ID3Mat = IV
        ID3Vec = ID3Mat + NAt3TT
        IV     = ID3Vec + NAt3
        MDV1 = MDV-IV+1
        Call TstCor(IV,MDV,'PreDWI 3')
        Call NumD3(IOut,NAt3,x1,xCur,H1,HCur,V(ID3Mat),V(ID3Vec))
        If(DEBUG) then
          Call OutMtS(IOut,'PreDWI, after NumD3: D3Vec=',0,0,V(ID3Vec),
     $      NAt3,1,NAt3,1)
          Call LTOutS(IOut,'PreDWI,after NumD3: D3Mat=',0,NAt3,
     $      V(ID3Mat),0)
          endIf
      else
        ID3Mat = IV
        ID3Vec = IV
        endIf
C
C     Call the Bulrish-Stoer integrator...
C
      Call TStamp(1,'Before mBS')
      IPtEM = -999
      IPtGM = -999
      StDvEM = Zero
      StDvGM = Zero
      NFail = 0
      BSFail = .False.
      Call ModBS(IOut,IPrint,IRwFit,ITSStp,2,IRwMM,IMMMod,IDWIp,DoOni,
     $  HUpFit,NAtoms,NAt3,NAtQ3,MaskON,IStep,RtMass,Efit,x1,dx,gfit,
     $  Scr(1,1),StpSiz,EPS,TSVec,IPtEM,IPtGM,StDvEM,StDvGM,CarLrg,
     $  Scr(1,11),Scr(1,12),BSFail,DwnHil,IFtSt,IFtMax,NFtPts,DoD3Num,
     $  V(ID3Mat),V(ID3Vec),LShpOb,V(IShpOb),V(IV),MDV1)
      If(BSFail) NFail = NFail+1
      If(BSFail.and.NFail.eq.5) Call GauErr('L123 - PreDWI Error!')
C
C     Compute some metrics about the just-completed corrector integration.
C     Some of these are only to be printed to the output file; others are
C     used when determining convergence below. Here is the list of values
C     computed in this block of code:
C           PE2StX      This is the distance from the predictor's end to
C                       the predictor's start points (x2-x1).
C           Ol2StX      This is the distance from the last corrector-cycle
C                       end point to the last corrector-cycle start point,
C                       which is always the same as the predictor's start
C                       point (xCur-x1).
C           St2NwX      This is the distance from the new corrector-cycle
C                       end point to the previous corrector-cycle start
C                       point, which is always the same as the predictor's
C                       start point (dx-x1).
C           Ol2NwX      This is the distance from the new corrector-cycle
C                       end point to the previous corrector-cycle end point
C                       (dx-xCur).
C
      Call TStamp(1,'After mBS')
      Call ASub(NAt3,x2,x1,Scr)
      PE2StX = VLen(NAt3,Scr)
      Call ASub(NAt3,xCur,x1,Scr)
      Ol2StX = VLen(NAt3,Scr)
      Call ASub(NAt3,dx,x1,Scr)
      St2NwX = VLen(NAt3,Scr)
      Call ASub(NAt3,dx,xCur,Scr)
      Ol2NwX = VLen(NAt3,Scr)
      If(IReCr3.eq.2) then
        Call ASub(NAt3,CarLrg,x1,Scr)
        Er2StX = VLen(NAt3,Scr)
        Call ASub(NAt3,CarLrg,x2,Scr)
        Er2NwX = VLen(NAt3,Scr)
      else
        Er2StX = Zero
        Er2NwX = Zero
        endIf
      If(IPrint.ge.1.or.DEBUG) then
        Write(IOut,1020)
        Write(IOut,1030) E2
        Write(IOut,1035) ECur
        Write(IOut,1040) Efit
        Write(IOut,1050) PE2StX
        Write(IOut,1060) Ol2StX
        Write(IOut,1070) St2NwX
        Write(IOut,1080) Ol2NwX
        If(IReCr3.eq.2) then
          Write(IOut,1090) Er2StX
          Write(IOut,1091) Er2NwX
          endIf
        If(IPrint.ge.2.or.DEBUG) then
          Call DmpFrc(IOut,'End Point Cartesians (a.u.):',1,NAtoms,xCur)
          Call DmpFrc(IOut,'End Point Gradient (a.u.):',1,NAtoms,gCur)
          Call DmpFrc(IOut,'Corrected End Point Cartesians (a.u.):',1,
     $      NAtoms,dx)
          Call DmpFrc(IOut,'Corrected Gradient (a.u.):',1,NAtoms,gfit)
          endIf
        endIf
C
C     Set CorDon according to the jobs settings and the metrics evaluated
C     above. After the initial CorDon evaluation, CorDon may be re-set if
C     IStep2=1 and IReCr0=2. Other special conditions for re-setting CorDon
C     are handled later in this routine.
C
      CorDon = .True.
      If(IReCr0.ne.1) then
        Write(IOut,8000) dxConv
        If(IReCr2.eq.1) then
          CorDon = Ol2NwX.lt.dxConv
        else if(IReCr2.eq.2) then
          Temp = GFloat(1)/GFloat(1000)
          CorDon = StDvGM.lt.Temp
        else
          Write(IOut,1100) IReCr0,IReCr1,IReCr2,IReCr3
          Call GauErr('PreDWI: Invalid IReCr2.')
          endIf
        CorDon = CorDon.and..not.BSFail
        endIf
      If(IReCr0.ne.1) then
        If(CorDon) then
          Write(IOut,8010)
          If(IStep2.eq.1.and.IReCr0.eq.3) then
            Write(IOut,8020)
            CorDon = .False.
            endIf
        else
          Write(IOut,8011)
          endIf
        endIf
      If(.not.CorDon.and.IStep2.ge.MxOpt)
     $  Call GauErr('Maximum number of corrector steps exceded.')
C
C     Depending on CorDon move dx and gfit to to the right place for use in
C     the rest of Link 123.
C
      If(CorDon) then
        Call AMove(NAt3,dx,x2)
        Call AMove(NAt3,gfit,grad2)
        E2 = Efit
      else
        Call AMove(NAt3,dx,xCur)
        Call AMove(NAt3,gfit,gCur)
        ECur = Efit
        endIf
C
C     Take care to undo gradient and/or transition vector flipping done
C     earlier, and undo mass-weighting.
C
      If(Flip) then
        Call AScale(NAt3,-One,grad1,grad1)
        Call AScale(NAt3,-One,TSVec,TSVec)
        endIf
      If(DoOni) then
        IAtMas = 1
        IRtMas = IAtMas + NAtoms
        IV     = IRtMas + NAtoms
        Call TstCor(IV,MDV,'PreDWI 4')
        Call CmpVc(0,NAtoms,LAtMas,NAtoms,AtMass,V(IAtMas),MaskON,1)
        Call CmpVc(0,NAtoms,LAtMas,NAtoms,RtMass,V(IRtMas),MaskON,1)
        Call MWeigV(.True.,NAtoms,xCur,RtMass,xCur)
        Call MWeigV(.True.,NAtoms,x1,RtMass,x1)
        Call MWeigV(.True.,NAtoms,x2,RtMass,x2)
        Call MWeigV(.False.,NAtoms,grad1,RtMass,grad1)
        Call MWeigV(.False.,NAtoms,grad2,RtMass,grad2)
        Call MWeigV(.False.,NAtomQ,gradQ1,V(IRtMas),gradQ1)
        Call MWeigV(.False.,NAtomQ,gradQ2,V(IRtMas),gradQ2)
        Call MWeigh(.False.,NAtomQ,V(IAtMas),H1)
        Call MWeigh(.False.,NAtomQ,V(IAtmas),H2)
        If(HUpFit) then
          Call MWeigV(.False.,NAtomQ,vUpdat,V(IRtMas),vUpdat)
          Call MWeigV(.False.,NAtomQ,wUpdat,V(IRtMas),wUpdat)
          endIf
      else
        Call MWeigV(.True.,NAtoms,xCur,RtMass,xCur)
        Call MWeigV(.True.,NAtoms,x1,RtMass,x1)
        Call MWeigV(.True.,NAtoms,x2,RtMass,x2)
        Call MWeigV(.False.,NAtoms,gCur,RtMass,gCur)
        Call MWeigV(.False.,NAtoms,grad1,RtMass,grad1)
        Call MWeigV(.False.,NAtoms,grad2,RtMass,grad2)
        Call MWeigh(.False.,NAtoms,AtMass,HCur)
        Call MWeigh(.False.,NAtoms,AtMass,H1)
        Call MWeigh(.False.,NAtoms,Atmass,H2)
        If(HUpFit) then
          Call MWeigV(.False.,NAtoms,vUpdat,RtMass,vUpdat)
          Call MWeigV(.False.,NAtoms,wUpdat,RtMass,wUpdat)
          endIf
        If(IReCr3.eq.2) Call MWeigV(.True.,NAtoms,CarLrg,RtMass,CarLrg)
        endIf
C
C     If CorDon has been set to TRUE, set ICorFl to 1. For some special
C     cases, we then re-set CorDon to FALSE. This block of code also sets
C     ICurWh and ensures that xCur is filled with the correct structure as
C     needed.
C
      If(CorDon) then
        ICorFl = 1
        If(IReCr1.ne.1) then
          CorDon = .False.
          Call AMove(NAt3,x2,xCur)
          ICurWh = 1
          endIf
      else
        ICurWh = 1
        If(IReCr3.eq.2) Call AMove(NAt3,CarLrg,xCur)
        endIf
      Return
      End
*Deck RStat1
      Function RStat1(IOut,UseGen,IWhich,NRdStat,NAtoms,Gen,Couplng,
     $  IRdStat,Cart)
      Implicit Real*8(A-H,O-Z)
C
C     This function computes value to be loaded into the Stats array used
C     by Link 123 that corresponds to user-defined statistic number IWhich.
C     The array IRdStat holds the user-defined stats definitions and Cart
C     holds the Cartesian coordinates for the structure of interest (note:
C     the Cartesian coordinates should be in bohr). Gen is the standard Gen
C     array of data. Couplng is the vector of path/transverse vib mode
C     coupling coefficients. If UseGen is sent as .FASLE., it is assumed
C     that Gen has not been passed and all statistics that depend on
C     elements in Gen are set to zero.
C
C     IRdStat encodes for different user-requested statistics. If
C     IRdStat(1,i) is > 0, column i of IRdStat corresponds to a standard
C     internal coordinate definition using the same format as used in IZRed
C     in Link 101. If IAtom1 < 0, it encodes for specific property values.
C     Recognized values include:
C           -1 ... Reaction Curvature (Gen(57))
C           -2 ... Reaction path coupling value. The IAtom2 value indicates
C                  which transverse mode is of interest.
C
C     -H.P. Hratchian, 2006, 2011.
C
      Logical UseGen
      Dimension Gen(*),Couplng(*),IRdStat(4,NRdStat),Cart(3,NAtoms)
      Save Zero
      Data Zero/0.d0/
 9000 Format(1x,'RStat1: Def(1)=',I6,'  Def(2)=',I6,'  Def(3)=',I6,
     $  ' Def(4)=',I6)
C
C     Determine the function value.
C
      IAtom1 = IRdStat(1,IWhich)
      IAtom2 = IRdStat(2,IWhich)
      IAtom3 = IRdStat(3,IWhich)
      IAtom4 = IRdStat(4,IWhich)
      If(UseGen) then
        RxCurv = Gen(57)
      else
        RxCurv = Zero
        endIf
      If(IAtom1.lt.0) then
        If(IAtom1.eq.-1) then
          RStat1 = RxCurv
        else if(IAtom1.eq.-2) then
          If(IAtom2.le.0.or.IAtom2.gt.3*NAtoms)
     $      Call GauErr('RStat1: Invalid IAtom2 for IAtom1=-2.')
          RStat1 = Couplng(IAtom2)
        else
          Write(IOut,9000) IAtom1,IAtom2,IAtom3,IAtom4
          Call GauErr('RStat1: Unknown property requested!')
          endIf
      else if(IAtom3.eq.0.and.IAtom2.lt.0) then
        If(IAtom2.eq.-1) then
          RStat1 = Cart(1,IAtom1)
        else if(IAtom2.eq.-2) then
          RStat1 = Cart(2,IAtom1)
        else if(IAtom2.eq.-3) then
          RStat1 = Cart(3,IAtom1)
        else
          Call GauErr('RStat 1 is confused!')
          endIf
      else if(IAtom3.eq.0) then
        RStat1 = Distan(Cart,IAtom1,IAtom2,0)
      else if(IAtom4.eq.0) then
        Call AngTst(IAtom1,IAtom2,IAtom3,Cart,RStat1)
      else
        Call RTors(IAtom1,IAtom2,IAtom3,IAtom4,Cart,RStat1)
        endIf
      Return
      End
*Deck RxDfPT
      Subroutine RxDfPT(IOut,MaxTab,Table,IType,Index1,Index2,Index3,
     $  Index4)
      Implicit Integer(A-Z)
C
C     This routine builds a parse table for reading user-defined statistics
C     definitions tracked during IRC calculations.
C
C     The parse table built here defines a grammar that the parser will use
C     to convert each user-defined statistics definition (one line of input
C     per) into a results array comprised of a definition type (position
C     IType) and up to four index values (positions Index1, Index2, Index3,
C     and Index4). The IType values include:
C            0    The line is blank.
C           -1    Cartesian coordinate. Index1 holds the atom center number
C                 and Index2 holds 1/2/3/4 for x/y/z/all.
C            1    Bond distance coordinate. Index1 and Index2 hold the two
C                 relevant atom center numbers.
C            2    Angle coordinate. Index1, Index2, and Index3 hold the
C                 three relevant atom center numbers.
C            3    Dihedral coordinate. Index1, Index2, Index3, and Index4
C                 hold the four relevant atom center numbers.
C           10    Magnitude of the curvature vector.
C           11    An element of the curvature vector. Index1 holds the
C                 value for the element of interest. If all elements are
C                 requested, Index1 is loaded with -1.
C
C     -H. P. Hratchian, 2012.
C
      Dimension Table(MaxTab)
C
C     Set-up the parse table beginning with the recurrsion into State Oper.
C
      Call PTInit(IOut,Table,MaxTab,I)
      Call PTStat(IOut,Table,MaxTab,I,'UserStats')
      Call PTTran(IOut,Table,MaxTab,I,.False.,'@AllInts',' ',0,0)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'@CartDefs',' ',0,0)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'@OtherDefs',' ',0,0)
      Call PTSpec(IOut,Table,MaxTab,I,'Nul',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'EOL','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     AllInts State.    This state handles user definitions that conform to
C                       the old style where all integer values are provided
C                       to define an internal or Cartesian coordinate.
C
      Call PTStat(IOut,Table,MaxTab,I,'AllInts')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',Index1,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',Index2,0)
      Call PTElem(IOut,Table,MaxTab,I,3,'EOL','Exi',IType,0,0,0,0,-1,0,
     $  0,Index2,0,0,0,0,4,0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',Index3,0)
      Call PTSpec(IOut,Table,MaxTab,I,'EOL','Exi',IType,1)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',Index4,0)
      Call PTSpec(IOut,Table,MaxTab,I,'EOL','Exi',IType,2)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'EOL','Exi',IType,3)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     CartDefs State.   This state handles user definitions using Cartesian
C                       coordinate specifications.
C
      Call PTStat(IOut,Table,MaxTab,I,'CartDefs')
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Cartesian','GtCart1',
     $  IType,-1)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'XYZ','GtCart1',IType,-1)
      Call PTElem(IOut,Table,MaxTab,I,2,'X','GtCart2',IType,0,0,0,0,-1,
     $  0,0,Index2,0,0,0,0,1,0,0)
      Call PTElem(IOut,Table,MaxTab,I,2,'Y','GtCart2',IType,0,0,0,0,-1,
     $  0,0,Index2,0,0,0,0,2,0,0)
      Call PTElem(IOut,Table,MaxTab,I,2,'Z','GtCart2',IType,0,0,0,0,-1,
     $  0,0,Index2,0,0,0,0,3,0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     GtCart1 ... Used by CartDefs when both atomic center and Cartesian
C                 type are being speficied as options.
      Call PTStat(IOut,Table,MaxTab,I,'GtCart1')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'(',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',Index1,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,',',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.False.,'All',' ',Index2,4)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'XYZ',' ',Index2,4)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'X',' ',Index2,1)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Y',' ',Index2,2)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Z',' ',Index2,3)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,')','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     GtCart2 ... Used by CartDefs when the Cartesian type has been
C                 specified by the operator and only the atom center number
C                 is speficied as an option.
      Call PTStat(IOut,Table,MaxTab,I,'GtCart2')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'(',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',Index1,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,')','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     OtherDefs State.  This state handles all other allowed user
C                       definitions.
C
      Call PTStat(IOut,Table,MaxTab,I,'OtherDefs')
      Call PTTran(IOut,Table,MaxTab,I,.False.,
     $  'MagnitudeCurvatureVector','Exi',IType,10)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'MagCurvatureVector',
     $  'Exi',IType,10)
      Call PTTran(IOut,Table,MaxTab,I,.True.,'|K|','Exi',IType,10)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'CurvatureVector',
     $  'GtCurv1',IType,11)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,'GtCurv1')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'(',' ',0,0)
      Call PTSpec(IOut,Table,MaxTab,I,'EOL','Exi',Index1,-1)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',Index1,0)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'All',' ',Index1,-1)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,')','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     End the parse table.
C
      Call PTEnd(IOut,Table,MaxTab,I)
      Return
      End
*Deck RxnCyc
      Subroutine RxnCyc(In,IOut,IPrint,IRwFit,Multip,SvPES,DoNRTst,
     $  WrURVA,LURVAFl,IURVAFl,SkpEGH,ICurWh,StpSiz,DVVv0,IStep,IStep2,
     $  IPath,IPG,ICorFl,CritGd,CritAn,TsPot,TsPotQ,Potent,PotSav,EQM,
     $  EPS,IMthd,IEmp,IHsUdt,NCorOp,ICorOp,IRxFwd,NPhase,IPFrq,IEig,
     $  NStat,MxOpt,NReduc,IQMMM,NPath,MaxStp,Gen,CCur,IRdStat,MaskON,
     $  AtMass,RtMass,TsGrad,TsGrdQ,TsFC,TsCart,TsVec,TsEVec,TsLamb,
     $  AGrad,GradQ,GrdCur,GrdSav,HesCar,HsEVec,HsEVal,HesCur,HVcCur,
     $  HVlCur,vUpdat,wUpdat,Stats,Cart,Vel,IDone,Frozen,Scr1,NAtoms,
     $  NAt3,LenHes,LnHEVc,LnHEVl,LenSq,NAt3TT,NAtomQ,NAtQ3,NAtQ3T,IAn,
     $  Restrt,IOP10,JobFlg,HvHess,DoHEVc,Update,Compct,IntHes,DwnHil,
     $  CorDon,IFtSt,IFtMax,MicOpt,IAtTyp,AtChMM,IAtRes,NResid,IResNa,
     $  IResNu,DXConv,ConvF,LenObj,PESObj,DoSvC,V,MDV,
CC-INSERT     
     $ IResNu,DXConv,ConvF,LenObj,IOp169,PESObj,DoSvC,V,MDV)
CC-INSERT 
      IMPLICIT REAL*8(A-H,O-Z)
C
C     This Routine is the main driver for L123. It sets up the various
C     vectors and data needed by the various reaction path following methods
C     available in the link.
C
C     PROGRAMMERS NOTE: RxnCyc sets up the transition vector in mass
C                       weighted coordinates. However, the Cartesian
C                       coordinates (Cart), gradient (AGrad), and Hessian
C                       (HesCar and TsFc) are NOT mass weighted. For DVV
C                       jobs, velocities (Vel) are kept in mass-weighted
C                       units.
C
C        - H.P. Hratchian, Spring 2003.
C          Significant improvments to L123 made by HPH, Fall 2004.
C          Further improvements made by HPH in 2009-2012.
C
      Parameter (IRwGen=501,IRWF1=584,IRWF2=585,IRW601=601,IRwMM=687,
     $  IMMMod=714,IMMCRS=715,IRwNMS=765,IRwPES=781)
      Integer FilNum
      Logical DEBUG,SvPES,DoNRTst,WrURVA,Restrt,CorDon,Done,SkpEGH,
     $  HvHess,DoHEVc,Update,Compct,IntHes,RxnFwd,DwnHil,Frozen(NAtoms),
     $  DoEmp,Do601,GetNM,ChekNM,HavROA,DoSvC,OK
      Dimension IPG(4),Potent(2,2),EQM(2,2),ICorOp(NCorOp),NPhase(4),
     $  Gen(*),CCur(NAt3),JJ(1),IRdStat(4,*),IURVAFl(LURVAFl),MaskOn(*),
     $  AtMass(NAtoms),RtMass(NAtoms),COM(3),TsGrad(NAt3),TsGrdQ(*),
     $  TsFC(*),TsCart(NAt3),TsVec(NAt3),TsEVec(LnHEVc),TsLamb(LnHEVl),
     $  AGrad(NAt3,2,NPath),GradQ(NAtQ3,2,NPath),GrdCur(NAt3),
     $  GrdSav(NAt3),HesCar(LenHes,2,*),HsEVec(NAt3,NAt3,*),
     $  HsEVal(NAt3,*),HesCur(*),HVcCur(*),HVlCur(*),vUpdat(*),
     $  wUpdat(*),Stats(NStat,MaxStp,NPath),Cart(NAt3,MaxStp,NPath),
     $  Vel(NAt3,MaxStp,NPath),IDone(NPath),IAn(NAtoms),Scr1(NAt3,20),
     $  MicOpt(NAtoms),IAtTyp(NAtoms),AtChMM(NAtoms),IAtRes(*),
     $  IResNa(*),IResNu(*),XX(1),PESObj(*),V(MDV)
      Dimension CMCoor(3),Pmom(3),RotVec(9),RJunk(1)
      Common /PhyCon/PhyCon(30)
      Common /MUnit/IUnit(20)
      Common /SYMINF/NOp1,NOp2,JTrans(3,8),TRot(3,3),TrVec(3)
      Save Zero,One,TFac,XX,JJ
      Data Zero/0.d0/,One/1.0D0/,TFac/1.d15/,XX/0.0d0/,JJ/0/
 1000 Format(1x,'Potential energy = ',F14.7)
 1010 Format(1x,'Current Structure is TS -> form Hessian eigenvectors.')
 1011 Format(28x,'Diagonalizing Hessian.')
 1012 Format(28x,'Recovering Hessian eigenvectors from chk file.')
 1020 Format(1x,'WARNING: NO IMAGINARY FREQUENCIES AT TS!')
 1021 Format(1x,'WARNING: TS HAS MORE THAN 1 IMAGINARY FREQUENCY!')
 1060 Format(10x,8('*'),' Start new reaction path calculation ',8('*'))
 1070 Format(1x,'RCFC Option Requested - Data Read From Chk File:')
 1075 Format(3x,'Energy From Chk = ',F14.7,/)
 1090 Format(1x,'SCF convergence failure!  Path Direction Number',I6,
     $  ' Point Number',I6)
 1300 Format(1x,'Setting the phase of the TS vector:',4(2x,I5))
 1500 Format(1x,'Reaction path inflection point has been passed.')
 1510 Format(3x,'Previous lowest Hessian eigenvalue=',F16.10,/,
     $  3x,'Current lowest Hessian eigenvalue =',F16.10)
 2000 Format(1x,'Point Number',I3,' in FORWARD path direction.')
 2001 Format(1x,'Point Number',I3,' in REVERSE path direction.')
 2050 Format(1x,'Calculating another point on the path.')
 2101 Format(1x,'Using Euler Reaction Path Following.')
 2102 Format(1x,'Using LQA Reaction Path Following.')
 2105 Format(1x,'Using DVV Reaction Path Following.')
 2110 Format(1x,'Using coordinate driving method.')
 2200 Format(1x,'Using modified Bulirsch-Stoer corrector integration.')
 2500 Format(1x,'Error in corrector energy =       ',F16.10,/,
     $  1x,'Magnitude of corrector gradient = ',F16.10,/,
     $  1x,'Magnitude of analytic gradient =  ',F16.10,/,
     $  1x,'Magnitude of difference =         ',F16.10,/,
     $  1x,'Angle between gradients (degrees)=  ',F8.4)
 3000 Format(1x,'Pt',I3,' Step number',I4,' out of a maximum of',I4)
 5000 Format(1x,'Calculation of FORWARD path complete.')
 5001 Format(1x,'Calculation of REVERSE path complete.')
 5010 Format(1x,'Reaction path calculation complete.')
 5020 Format(1x,'Beginning calculation of the REVERSE path.')
 5030 Format(1x,'Copying FFX from chk file, length=',I20,'.')
 5501 Format(1x,'**** End of Projected Frequency Analysis ****')
 9000 Format(1x,'RxnCyc: IMthd = ',I2)
 9100 Format(1x,'RxnCyc: Step size is ',F12.5,'.')
 9900 Format(1x,'RxnCyc: ICurWh = ',I2)
C
C     Initialization stuff
C
C     GScale:  Converts nuclear gradient to acceleration for DVV jobs.
C              GScale is also used for computing vibrational frequencies
C              when Hessians are diagonalized. (Hartree/Bohr ->
C              amu*Bohr/fs**2)
      Call TStamp(1,'Top RxnCyc')
      IV = 1
      IURVA = IUnit(18)
      DEBUG = .False.
      GetNM = .False.
      JobFlg = -1
      Done   = .False.
      RxnFwd = iRxFwd.eq.1.or.DwnHil
      If(IPath.eq.2) RxnFwd = .not.RxnFwd
      CorDon = .True.
      DoEmp = iEmp.eq.1
      GScale = One/(TFac*TFac*(PhyCon(7)**2)*PhyCon(2)/PhyCon(8))
      IDir0 = 0
C
C     Allocate and initialize a few arrays that we carve out of the scratch
C     space.
C
      ICplng = IV
      IV     = ICplng + NAt3
      Call TstCor(IV-1,MDV,'RxnCyc-Top')
      Call AClear(IV-1,V(ICplng))
C
C     Take care of frozen atom arrays.
C
      Call FrzAtm(.True.,NAtoms,MicOpt,NFroz,Frozen)
C
C     Initialize IDone and Stats arrays.
C
      If(IStep.eq.2) then
        IDone(IPath) = 0
        If(IPath.eq.1) Call AClear(NStat*MaxStp*NPath,Stats)
        endIf
C
C     Check flag for SCF convergence failure.
C
      Call ILSW(2,5,Kill)
      If(Kill.ne.1) then
        If(IPrint.ge.0.and.IStep.eq.2.and..not.Restrt.and.IPath.eq.1)
     $    Write(IOut,1060)
      else
        Call ILSW(1,5,0)
        Write(IOut,1090) IPath,IStep
        If(IPath.eq.1.and.IStep.eq.2)
     $    Call GauErr('Unrecoverable SCF failure')
        IDone(IPath) = 0
        JobFlg = 2
        Return
        endIf
C
C     Read in and fill the structure (Cart), energy (Potent), gradient
C     (AGrad), and Hessian (HesCar) for the current structure (kept in CCur
C     and Cart). If this is the first step, also fill TSCart, TSPot,
C     TSGrad, and TSFc. If this is the first step but NOT the first path to
C     be followed, simply fill Cart, Potent, AGrad, and HesCar from TSCart,
C     TSPot, TSGrad, and TSFc. If this is an intermediate step, then shift
C     the old energies, gradients, and Hesians first. At this point in the
C     code, we also diagonalize the Hessian to form its eigenvectors
C     (HsEVec and TsEVec) and eigenvalues (HsEVal and TSLamb), as
C     necessary. After fnding the eigenvectors of the Hessian at the
C     transition structure (TsEVec) we load the transition vector (TsVec).
C
C     If Symmetry has been used by the rest of GAUSSIAN, then we must
C     transform the starting Cartesian coordinates to Standard Orientation
C     so that the computed (and read in) gradient and Hessian will
C     correspond to our structure properly. We also will turn off symmetry
C     for the rest of the calculation.
C
      HsEVlP = One
      HsEVlC = One
      If(DoHEVc) HsEVlP = HsEVal(1,IPath)
      If(.not.SkpEGH) then
        If(IStep.ge.3) then
          If(IQMMM.eq.1.and.Compct) then
            EQM(2,IPath) = EQM(1,IPath)
            Call AMove(NAtQ3,GradQ(1,1,IPath),GradQ(1,2,IPath))
            endIf
          Call GtEGrd(IOut,IPrint,1,IQMMM,IRWGen,IRwF1,IMMCRS,NAtoms,
     $      NAtomQ,Compct,MicOpt,MaskON,PotCur,EQM(1,IPath),GrdCur,
     $      GradQ(1,1,IPath),V(IV),MDV-IV+1)
          If(HvHess) then
            If(IStep2.eq.1) then
              Call AMove(LenHes,HesCar(1,1,IPath),HesCur)
              Call GetHes(IOut,IPrint,IRwF2,NAtoms,IQMMM,0,Compct,
     $          Update,IHsUdt,Cart(1,IStep-2,IPath),CCur,
     $          AGrad(1,1,IPath),GrdCur,GradQ(1,2,IPath),
     $          GradQ(1,1,IPath),thetaUp,phiUp,vUpdat,wUpdat,MicOpt,
     $          MaskON,HesCur,V(IV),MDV-IV+1)
            else
              Call AMove(LenHes,HesCar(1,2,IPath),HesCur)
              Call GetHes(IOut,IPrint,IRwF2,NAtoms,IQMMM,0,Compct,
     $          Update,IHsUdt,Cart(1,IStep-2,IPath),CCur,
     $          AGrad(1,2,IPath),GrdCur,GradQ(1,2,IPath),
     $          GradQ(1,1,IPath),thetaUp,phiUp,vUpdat,wUpdat,MicOpt,
     $          MaskON,HesCur,V(IV),MDV-IV+1)
              endIf
            endIf
        else if(DwnHil) then
          If(IStep.eq.2.and.IOp10.eq.5) then
            JRWGen = FilNum(IRWGen,IUnit(9))
            JRwF1 = FilNum(IRwF1,IUnit(9))
            JMMCRS = FilNum(IMMCRS,IUnit(9))
            JRwF2 = FilNum(IRwF2,IUnit(9))
          else
            JRWGen = IRWGen
            JRwF1 = IRwF1
            JMMCRS = IMMCRS
            JRwF2 = IRwF2
            endIf
          Call GtEGrd(IOut,IPrint,1,IQMMM,JRWGen,JRwF1,JMMCRS,NAtoms,
     $      NAtomQ,Compct,MicOpt,MaskON,PotCur,EQM(1,IPath),GrdCur,
     $      GradQ(1,1,IPath),V(IV),MDV-IV+1)
          If(HvHess) then
            Call ILSW(2,24,iHvHes)
            If(IOp10.ne.5.and.iHvHes.le.1)
     $        Call GauErr('L123: Hessian Requested but not available!')
            Call GetHes(IOut,IPrint,JRwF2,NAtoms,IQMMM,0,Compct,.False.,
     $        Junk,XX,XX,XX,XX,XX,XX,thetaUp,phiUp,vUpdat,wUpdat,MicOpt,
     $        MaskON,HesCur,V(IV),MDV-IV+1)
            endIf
          Call AMove(NAt3,CCur,Cart)
        else if(IPath.eq.1) then
          If(IOP10.eq.5) then
            If(IPrint.ge.0) then
              Write(IOut,1070)
              Call PrtNam('Chk',IOut)
              endIf
            Call GtEGrd(IOut,IPrint,1,IQMMM,FilNum(IRWGen,IUnit(9)),
     $        FilNum(IRwF1,IUnit(9)),FilNum(IMMCRS,IUnit(9)),NAtoms,
     $        NAtomQ,Compct,MicOpt,MaskON,TSPot,TSPotQ,TSGrad,TsGrdQ,
     $        V(IV),MDV-IV+1)
            IRwF2C = FilNum(IRwF2,IUnit(9))
            Call GetHes(IOut,IPrint,IRwF2C,NAtoms,IQMMM,0,Compct,
     $        .False.,0,XX,XX,XX,XX,XX,XX,thetaUp,phiUp,vUpdat,wUpdat,
     $        MicOpt,MaskON,TSFc,V(IV),MDV-IV+1)
            GetNM = .True.
            If(IQMMM.eq.1) then
              LRwF2C = ITqry(IRwF2C)
              If(LRwF2C.gt.0) then
                If(IPrint.ge.0) Write(IOut,5030) LRwF2C
                Call ConDD2(IRwF2,LRwF2C)
                Call Trsfr(LRwF2C,IRwF2C,IRwF2,V(IV),MDV-IV+1)
                endIf
              endIf
            If(IPrint.ge.0) Write(IOut,1075) TSPot
          else
            Call GtEGrd(IOut,IPrint,1,IQMMM,IRWGen,IRwF1,IMMCRS,NAtoms,
     $        NAtomQ,Compct,MicOpt,MaskON,TSPot,TsPotQ,TSGrad,TsGrdQ,
     $        V(IV),MDV-IV+1)
            Call GetHes(IOut,IPrint,IRwF2,NAtoms,IQMMM,0,Compct,.False.,
     $        0,XX,XX,XX,XX,XX,XX,thetaUp,phiUp,vUpdat,wUpdat,MicOpt,
     $        MaskON,TSFc,V(IV),MDV-IV+1)
            endIf
C
C         Make sure geometry is at center of mass coordinates
C
          If(.not.DoEmp) then
            Call CntMas(NAtoms,AtMass,CCur,TotWt,COM)
            Call Transl(-1,3,NAtoms,COM,CCur)
            endIf
          Call AMove(NAt3,CCur,TsCart)
C
C         Transform to standard orientation and set TRot to a unit matrix;
C         turn off symmetry for the remainder of the calculation.
          Call ILSW(2,26,NoSym)
          If(NoSym.ne.1) Call ILSW(1,26,1)
          If(.not.DoEmp) Call ProjG(IOut,IPrint-1,1,.False.,NAtoms,
     $      TsCart,TsGrad,Scr1,XX,Scr1(1,7))
          Call AMove(NAt3,TsCart,Cart)
          PotCur = TSPot
          EQM(1,IPath) = TsPotQ
          Call AMove(NAt3,TSGrad,GrdCur)
          Call AMove(NAtQ3,TSGrdQ,GradQ(1,1,IPath))
          If(HvHess) Call AMove(LenHes,TSFc,HesCur)
        else if(IPath.eq.2) then
          Call AMove(NAt3,TSCart,Cart(1,1,IPath))
          PotCur = TSPot
          EQM(1,IPath) = TsPotQ
          Call AMove(NAt3,TSGrad,GrdCur)
          Call AMove(NAtQ3,TSGrdQ,GradQ(1,1,IPath))
          If(HvHess) Call AMove(LenHes,TSFc,HesCur)
          If(DoHEVc) then
            Call AMove(LnHEVc,TsEVec,HVcCur)
            Call AMove(LnHEVl,TSLamb,HVlCur)
            endIf
        else
          Call GauErr('L123: RxnCyc is confused!')
          endIf
        endIf
      If(NFroz.gt.0) Call ZeroD(NAtoms,Frozen,1,1,GrdCur)
C
C     Print current energy and derivatives.
C
      If(DEBUG.or.IPrint.ge.2) then
        If(IStep.eq.2.and..not.DwnHil) then
          Write(IOut,1000) TsPot
          Call DmpFrc(IOut,'TS Forces:',1,NAtoms,TsGrad)
          Call DumpFC(IOut,'TS Force constants:',NAt3,TsFC)
        else
          Write(IOut,1000) PotCur
          Call DmpFrc(IOut,'Current forces:',1,NAtoms,GrdCur)
          If(HvHess.and.(DEBUG.or.IPrint.ge.3)) Call DumpFC(IOut,
     $      'Current force constants:',NAt3,HesCur)
          endIf
        endIf
C
C     If requested, write the current force constant matrix to RWF 601.
C     What gets written to RWF 601 is a real flag followed by the full FC
C     matrix. The real flag indicates whether this is an updated or
C     analytic Hessian. The values of this flag are:
C           0     Analytic Hessian
C           1     Updated Hessian
C
      Do601 = .False.
      If(Do601.and.HvHess.and..not.(IPath.eq.2.and.IStep.eq.2)) then
        Len601 = (MaxStp*NPath+1)*(LenHes+1)
        If(IPath.eq.1) then
          IPos = (IStep-2)*(LenHes+1)
        else
          IPos = (IStep-3)*(LenHes+1)+IDone(1)*(LenHes+1)
          endIf
        If(Update) then
          UpdFlg = GFloat(1)
        else
          UpdFlg = GFloat(0)
          endIf
        If(IPath.eq.1.and.IStep.eq.2) Call ConDDF(IRW601,Len601)
        Call FileIO(1,-IRW601,1,UpdFlg,IPos)
        Call FileIO(1,IRW601,LenHes,HesCur,0)
        endIf
C
C     Diagonalize FC matrix as needed. If this is the first step but NOT
C     the first pathway, then TsVec and TsLamb have already been solved, so
C     only solve for the Hessian eigenvectors and eigenvalues as needed.
C
      Call TStamp(1,'Bef DiagFC')
      If(.not.SkpEGH) then
        If((HvHess.or.(IStep.eq.1.and.IPath.eq.1.and..not.DwnHil)).and.
     $    .not.DoEmp) then
          Call MofI(NAtoms,CCur,AtMass,CMCoor,Pmom,RotVec)
          Call RigHnd(RotVec)
          Call AClear(NAt3*6,Scr1)
          Call TRVect(IOut,0,NAtoms,AtMass,CCur,CMCoor,RotVec,NTrRot,
     $      Scr1)
          endIf
        If(IStep.eq.2.and.IPath.eq.1.and..not.DwnHil) then
          If(IPrint.ge.1) Write(IOut,1010)
          If(.not.DoEmp) then
            If(GetNM) then
              IRwNMC = FilNum(IRwNMS,IUnit(9))
              NVib = 0
              NDim = 0
              NDim0 = 0
              Call NMSave(IOut,0,-2,0,0,XX,IRwNMC,NAtoms,NVib,0,JJ,
     $          HavROA,NumROA,NDFDPol,NDim,NDim0,LE2,LFlags,JJ,XX,XX,XX)
              GetNM = NVib.gt.0
              If(GetNM) then
                IAtMas = IV
                IFlags = IAtMas + NAtoms
                IE2 = IFlags + InToWP(NAt3*LFlags)
                IV1 = IE2 + NAt3*LE2
                Call TstCor(IV1,MDV,'RxnCyc-NMSave')
                Call NMSave(IOut,0,2,0,0,XX,IRwNMC,NAtoms,NVib,0,JJ,
     $            HavROA,NumROA,NDFDPol,NDim,NDim0,LE2,LFlags,V(IFlags),
     $            V(IAtMas),V(IE2),TsEVec)
                GetNM = ChekNM(NAtoms,NVib,AtMass,V(IAtMas))
                endIf
              endIf
            If(GetNM) then
              If(IPrint.ge.1) Write(IOut,1012)
              Call RawFrq(NVib,V(IE2),TsLamb)
              Call MWeigA(.False.,.True.,NAtoms,NVib,AtMass,TsEVec,
     $          TsEVec)
            else
              If(IPrint.ge.1) Write(IOut,1011)
              IFFX   = IV
              IOrtho = IFFX   + LenHes
              IVV    = IOrtho + LenSq
              IV1    = IVV    + LenSq
              Call TstCor(IV1,MDV,'RxnCyc-VibFq2, 1')
              Call AMove(LenHes,TsFC,V(IFFX))
              Call VibFq2(In,IOut,-1,0,NAtoms,0,IAn,
     $          Cart(1,IStep-1,IPath),NAt3,V(IFFX),AtMass,V(IOrtho),
     $          V(IVV),TsEVec,Scr1,Scr1(1,3),PhyCon,NVib,NImag,.False.,
     $          IPG,NFroz,Frozen,0,RJunk,RJunk,TRot,RJunk,.False.,
     $          NTrRot,NReduc,IAtTyp,AtChMM,MicOpt,.False.,ITemp1,
     $          ITemp2,.True.,XNorm,ITemp3,ITemp4,IntHes,IDir0,ITemp5,
     $          V(IV1),MDV-IV1+1)
              Call AMove(NVib,Scr1(1,3),TsLamb)
              endIf
            Call FixPh1(NAt3,NAt3,NVib,0,0,0,0,TsEVec,Scr1,Scr1)
          else
            NVib = NAt3
            Call DiagD(TsFc,TsEVec,TsLamb,NAt3,Scr1,Scr1(1,2),NAt3,
     $        .False.)
            endIf
          Call AMove(NAt3,TsEVec,TsVec)
          If(NFroz.gt.0) Call ZeroD(NAtoms,Frozen,1,1,TsVec)
          If(IPrint.ge.2) Call OutMtS(IOut,'TsVec:',0,0,TsVec,NAt3,1,
     $      NAt3,1)
          If(TsLamb(1).ge.Zero) write(IOut,1020)
          If(TsLamb(2).lt.Zero) write(IOut,1021)
          If(NPhase(1).ne.0) then
            Write(IOut,1300) NPhase(1),NPhase(2),NPhase(3),NPhase(4)
            Call MWeigV(.False.,NAtoms,TsVec,RtMass,TsVec)
            Call FxPhas(IOut,NPhase(1),NPhase(2),NPhase(3),NPhase(4),
     $        TsVec,CCur,NAtoms)
            Call ANeg(NAt3,TsVec,TsVec)
            Call MWeigV(.True.,NAtoms,TsVec,RtMass,TsVec)
            endIf
          Call AMove(NAt3,TsVec,TsEVec)
          If(DoHEVc) then
            Call AMove(LnHEVc,TsEVec,HVcCur)
            Call AMove(LnHEVl,TsLamb,HVlCur)
            endIf
        else if(DoHEVc.and..not.(IStep.eq.2.and.IPath.eq.2)) then
          If(.not.DoEmp) then
            IFFX   = IV
            IOrtho = IFFX   + LenHes
            IVV    = IOrtho + LenSq
            IV1    = IVV    + LenSq
            Call TstCor(IV1,MDV,'RxnCyc-VibFq2, 2')
            Call AMove(LenHes,HesCur,V(IFFX))
            Call VibFq2(In,IOut,-1,0,NAtoms,0,IAn,Cart(1,IStep-1,IPath),
     $        NAt3,V(IFFX),AtMass,V(IOrtho),V(IVV),HVcCur,Scr1,
     $        Scr1(1,3),PhyCon,NVib,NImag,.False.,IPG,NFroz,Frozen,0,
     $        RJunk,RJunk,TRot,RJunk,.False.,NTrRot,0,IAtTyp,AtChMM,
     $        MicOpt,.False.,ITemp1,ITemp2,.True.,XNorm,ITemp3,ITemp4,
     $        IntHes,IDir0,ITemp5,V(IV1),MDV-IV1+1)
            Call AClear(NAt3,HVlCur)
            Call AMove(NVib,Scr1(1,3),HVlCur)
            Call FixPh1(NAt3,NAt3,NVib,0,0,0,0,HVcCur,Scr1,Scr1)
          else
            Call DiagD(HesCur,HVcCur,HVlCur,NAt3,Scr1,Scr1(1,2),NAt3,
     $        .False.)
            If(IPrint.ge.2) then
              Call OutMtS(IOut,'Hessian eigenvectors:',0,0,HVcCur,NAt3,
     $          NAt3,NAt3,NAt3)
              Call OutMtS(IOut,'Hessian eigenvalues:',0,0,HVlCur,NAt3,1,
     $          NAt3,1)
              endIf
            NVib = NAt3
            endIf
          endIf
        endIf
C
C     If we are re-evaluating the energy and derivatives during
C     predictor-corrector runs, report on the quality of the last point's
C     fitted energy and gradient relative to the "real" values just
C     evaluated.
C
      If(ICurWh.eq.1) then
        Call ASub(NAt3,GrdSav,GrdCur,Scr1)
        XLen1 = VLen(NAt3,GrdSav)
        XLen2 = VLen(NAt3,GrdCur)
        XLen3 = VLen(NAt3,Scr1)
        Write(IOut,2500) PotSav-PotCur,XLen1,XLen2,XLen3,AngDeg(0,NAt3,
     $    GrdSav,GrdCur)
        endIf
      If(.not.SkpEGH) then
        If(IStep.ge.3.and.IStep2.eq.1) then
          Potent(2,IPath) = Potent(1,IPath)
          Call AMove(NAt3,AGrad(1,1,IPath),AGrad(1,2,IPath))
          If(HvHess) Call AMove(LenHes,HesCar(1,1,IPath),
     $      HesCar(1,2,IPath))
          endIf
        If(ICurWh.eq.0) then
          Potent(1,IPath) = PotCur
          Call AMove(NAt3,GrdCur,AGrad(1,1,IPath))
          If(HvHess) Call AMove(LenHes,HesCur,HesCar(1,1,IPath))
          If(DoHEVc) then
            Call AMove(LnHEVc,HVcCur,HsEVec(1,1,IPath))
            Call AMove(LnHEVl,HVlCur,HsEVal(1,IPath))
            endIf
          endIf
        endIf
C
C     If we have been diagonalizing the Hessian, check to see if we have
C     crossed the inflection point on the reaction path.
C
      If(DoHEVc) HsEVlC = HsEVal(1,IPath)
      If(ICurWh.eq.0.and.DoHEVc.and.(HsEVlP.lt.Zero).and.
     $  (HsEVlC.gt.Zero)) then
        Write(IOut,1500)
        Write(IOut,1510) HsEVlP,HsEVlC
        endIf
C
C     Compute the "stiff-ness" of the IRC integration. This code is turned
C     on and off by logical variable DoRStf. DoRStf is hard-wired for now.
C     In the future, we may want to connect this logical to an IOp.
C
C     DoRStf = .False.
C     If(DoRStf.and.IStep.ge.3) then
C     Call RxStif(IOut,NAt3,NAt3TT,SpecRd,AGrad,HesCar,Scr1,Scr2)
C     endIf
C
C     Given the supplied step size value (StpSiz) determine if the user
C     desires this step size on the mass-weighted PES or on the non-MW PES.
C     If the user supplied the step size for the non-MW PES, then determine
C     the appropriate step size on the MW PES for use in the rest of L123.
C     This only needs to be worked out in the first step of the job, and it
C     is carried out in Routine StpSet. This is skipped if we are running
C     DVV since in this case StpSiz is a time step value. We also skip
C     StpSet if we are doing coordinate driving.
C
      If(IPath.eq.1.and.IStep.lt.3.and.IMthd.ne.5.and.IMthd.ne.10.and.
     $  .not.DoEmp) Call StpSet(IOut,IPrint,NAt3,NAtoms,TSVec,
     $  AGrad(1,1,IPath),StpSiz,RtMass,Scr1,DwnHil)
      If(StpSiz.lt.Zero) then
        Write(IOut,9100) StpSiz
        Call GauErr('RxnCyc: Invalid step size found!')
        endIf
C
C     Now, we do the actual integration of the reaction path following
C     equations. First, check to make sure that a valid/available method
C     has been chosen.
C
      If((IMthd.le.0).or.(IMthd.gt.6.and.IMthd.lt.10).or.
     $  (IMthd.gt.10)) then
        Write(IOut,9000) IMthd
        Call GauErr('L123 ERROR: ILLEGAL VALUE OF IMTHD!')
        endIf
      If(IPath.eq.1) then
        ITStep = IStep-2
        If(NPath.eq.2) ITStep = ITStep + Max(0,IDone(2)-1)
      else if(IPath.eq.2) then
        ITStep = Max(0,IDone(1)-1)+IStep-2
      else
        Call GauErr('RxnCyc: IPath.ne.1 or 2.')
        endIf
      If(IPrint.ge.0.and.IStep.ge.3) write(IOut,3000) ITStep,IStep2,
     $  MxOpt
C
C     If we're using some sort of corrector routine in the selected method,
C     we need to call it now.
C
      Call TStamp(1,'Bef Corrector')
      If(IStep.ge.3.and.(IStep-1).gt.IDone(IPath)) then
        If(IMthd.eq.3.or.IMthd.eq.6) then
          If(IPrint.ge.1) Write(IOut,2200)
          ICarLg = IV
          IV     = ICarLg + NAt3
          Call TstCor(IV-1,MDV,'RxnCyc before PreDWI.')
          Call PreDWI(IOut,IPrint,IRwFit,ICurWh,NCorOp,ICorOp,ICorFl,
     $      IRwMM,IMMMod,IQMMM,Compct,DoHEVc,NAtoms,NAt3,NAtQ3,NAtQ3T,
     $      LenHes,LnHEVc,LnHEVl,MaskON,CCur,Cart(1,IStep-2,IPath),
     $      Cart(1,IStep-1,IPath),GrdCur,AGrad(1,2,IPath),
     $      AGrad(1,1,IPath),HesCur,HVcCur,HVlCur,HesCar(1,2,IPath),
     $      HesCar(1,1,IPath),HsEVec(1,1,IPath),HsEVal(1,IPath),Update,
     $      thetaUp,phiUp,vUpdat,wUpdat,PotCur,Potent(2,IPath),
     $      Potent(1,IPath),GradQ(1,2,IPath),GradQ(1,1,IPath),
     $      EQM(2,IPath),EQM(1,IPath),Scr1,Scr1(1,15),Scr1(1,16),IStep,
     $      IStep2,MxOpt,NAt3TT,EPS,AtMass,RtMass,dxConv,CorDon,TSVec,
     $      DwnHil,IFtSt,IFtMax,V(ICarLg),V(IV),MDV-IV+1)
        else if(IMthd.eq.4) then
          CorStp = Abs(Stats(2,IStep-1,IPath)-Stats(2,IStep-2,IPath))
          Call GS2Cor(IOut,IPrint,NAtoms,NAt3,LenHes,IStep,IStep2,MxOpt,
     $      CorDon,CorStp,Cart(1,IStep-2,IPath),Cart(1,IStep-1,IPath),
     $      AGrad(1,2,IPath),AGrad(1,1,IPath),Scr1(1,1),
     $      HesCar(1,1,IPath),HsEVec(1,1,IPath),HsEVal(1,IPath),TSVec,
     $      Scr1(1,2),Scr1(1,3),Scr1(1,4),AtMass,RtMass,Scr1(1,5),
     $      RxnFwd,DwnHil,ConvF)
          endIf
        endIf
C
C     If corrector cycles are complete (or not being used), then evaluate
C     statistics.
C
CC-INSERT     
C     
C     Allocate space for curvature and path direction vectors in new
C     browsing file
C     
      ITang = IV
      ICurv = ITang + NAt3
      IScr = ICurv + NAt3
      IV1 = IScr
      Call TstCor(IV1-1,MDV,'RxnCyc-Decomp')
CC-INSERT

      If((IStep.ge.3).and.((IStep-1).gt.IDone(IPath)).and.CorDon) then
        If(DoSvC) then
          IPtNum = IStep - 1
          If(.not.RxnFwd) IPtNum = IPtNum + MaxStp
          Call SvCalc(IOut,IPrint,IPtNum,-1,V(IV),MDV1)
          endIf
        If(IPFrq.eq.2.and.HvHess) then
          If(.not.DoEmp) then
            NVib = 0
            Call LenE2(.False.,0,LE2Fix,NumROA,NROADy,LE2,LFlags)
C            IFlags = IV
C            IE2 = IFlags + InToWP(NAt3*LFlags)
C            IV1 = IE2 + NAt3*LE2
CC-INSERT     
            IFlags = IV1
            IE2 = IFlags + InToWP(NAt3*LFlags)
            IV2 = IE2 + NAt3*LE2
CC-INSERT     
C     Call TstCor(IV1-1,MDV,'RxnCyc-IRCVib')
CC-INSERT     
            Call TstCor(IV2-1,MDV,'RxnCyc-IRCVib')
CC-INSERT     
            Call IRCVib(IOut,IPrint,.False.,0,NAtoms,NAt3,NAt3TT,.True.,
     $      Zero,Zero,Multip,NTrRot,PhyCon,Gen,IAn,IAtTyp,AtMass,
     $      AtChMM,MicOpt,IAtRes,NResid,IResNa,IResNu,NVib,NFroz,
     $      Frozen,Cart(1,IStep-1,IPath),AGrad(1,1,IPath),
     $      HesCar(1,1,IPath),V(IE2),LFlags,V(IFlags),V(ICplng),ZPE,
C     $ ETherm,HTherm,GTherm,V(IV1),MDV-IV1+1)
CC-INSERT     
     $      ETherm,HTherm,GTherm,V(IV2),MDV-IV2+1)
CC-INSERT    

          else
            Call RxCurv(IOut,NAt3,AGrad(1,1,IPath),HesCar(1,1,IPath),
     $        RCurv,V(IV),MDV-IV+1)
            Gen(57) = RCurv
            endIf
          Write(IOut,5501)
        else
          Gen(57) = Zero
          endIf
C
C       Save vital statistics to the Stats array.
C
        Call RxStat(IOut,IPrint,Gen,IRdStat,Stats(1,1,IPath),
     $    Cart(1,IStep-2,IPath),Cart(1,IStep-1,IPath),Potent(1,IPath),
     $    V(ICplng),RtMass,Scr1,NAtoms,NAt3,IStep-1,IDone,IPath,RxnFwd,
     $    2,NStat,NPath,StpSiz,IAn)
        endIf
C
C     Print out information and stats for the latest step. Also, if
C     requested, append the latest step's data to a URVA input file.
C
      If(.not.(IStep.eq.2.and.IPath.eq.2).and.CorDon) then
        Call RxPrnt(IOut,IPrint,Stats(1,1,IPath),Cart(1,IStep-1,IPath),
     $    NAtoms,NAt3,IStep,IStep2,IDone,IPath,NPath,NStat,IAn)
CC-INSERT
C     
C     Calculate eta and kappa, needed for new browsing file

C     
      Call RxCurvXX(IOut,NAt3,AGrad(1,1,IPath),HesCar(1,1,IPath),
      $ RCurv,Stats(2,IStep-1,IPath),Natoms,AtMass,RtMass,V(ITang),
      $ V(ICurv),V(IScr))
CC-INSERT

        If(WrURVA) Call InURVA(IOut,IURVA,IPrint,.False.,HvHess,LURVAFl,
     $    IURVAFl,NAtoms,ITStep,Stats(2,IStep-1,IPath),Potent(1,IPath),
     $    IAn,Cart(1,IStep-1,IPath),AGrad(1,1,IPath),HesCar(1,1,IPath),
CC-INSERT     
     $    V(ITang),V(ICurv),AtMass,IOp169)
CC-INSERT 
        endIf
C
C     Now, take a (predictor) step or wrap-up the job if we're done!
C
      Call TStamp(1,'Before Predictor')
      If(IStep.ge.3) Call JobTst(IOut,NAtoms,NAt3,IMthd,IStep,MaxStp,
     $  StpSiz,RtMass,Cart(1,IStep-3,IPath),Cart(1,IStep-2,IPath),
     $  Cart(1,IStep-1,IPath),AGrad(1,1,IPath),Scr1,CritGd,CritAn,
     $  CorDon,Done)
      If(.not.Done.and.CorDon) then
        If(IPrint.ge.0) Write(IOut,2050)
        If(IPrint.ge.1) then
          If(RxnFwd) then
            Write(IOut,2000) IStep-1
          else
            Write(IOut,2001) IStep-1
            endIf
          endIf
        If(IMthd.eq.1.or.IMthd.eq.4.or.IMthd.eq.6) then
          If(IPrint.ge.1) Write(IOut,2101)
          Call RxnEul(IOut,IPrint,NAt3,NAtoms,IStep,StpSiz,
     $      Cart(1,IStep-1,IPath),Cart(1,IStep,IPath),AGrad(1,1,IPath),
     $      TSVec,RtMass,RxnFwd,DwnHil)
          ActStp = StpSiz
        else if(IMthd.eq.2.or.IMthd.eq.3) then
          If(IPrint.ge.1) Write(IOut,2102)
          IA     = IV
          IAlpha = IA     + NAt3*NAt3
          ILQASc = IAlpha + NAt3*NAt3
          IV1    = ILQASc + NAt3*NAt3
          Call TstCor(IV1,MDV,'RxnCyc-LQA')
          Call LQA(IOut,IPrint,Cart(1,IStep-1,IPath),
     $      Cart(1,IStep,IPath),AGrad(1,1,IPath),HsEVec(1,1,IPath),
     $      HsEVal(1,IPath),StpSiz,ActStp,NAt3,NAtoms,RtMass,IStep,
     $      RxnFwd,Scr1(1,1),Scr1(1,2),Scr1(1,3),V(IA),V(IAlpha),
     $      V(ILQASc),ILQAFg,DwnHil)
          If(DoNRTst.and.(HsEVal(1,IPath).gt.Zero)) Call NRtest(IOut,
     $      IPrint,ILQAFg,NAt3,NAtoms,NFroz,CritGd,Done,
     $      Cart(1,IStep-1,IPath),Cart(1,IStep,IPath),Scr1,
     $      AGrad(1,1,IPath),Scr1(1,2),HsEVec(1,1,IPath),
     $      HsEVal(1,IPath),StpSiz,ActStp,RtMass,Scr1(1,3))
        else if(IMthd.eq.5) then
          If(IPrint.ge.1) write(IOut,2105)
          If(IStep.gt.3) then
            TStep0 = Stats(3,IStep-1,IPath)-Stats(3,IStep-2,IPath)
            Call AMove(NAt3,Cart(1,IStep-2,IPath),Scr1)
          else
            TStep0 = Zero
            endIf
          If(IStep.eq.2.and.IPath.eq.2) StpSiz = Stats(3,2,1)
          Call DVV(IOut,IPrint,NAt3,NAtoms,IStep,MaxStp,StpSiz,TStep0,
     $      TStep1,DVVv0,GScale,RtMass,Scr1,Cart(1,IStep-1,IPath),
     $      Cart(1,IStep,IPath),Vel(1,1,IPath),AGrad(1,1,IPath),
     $      AGrad(1,2,IPath),Scr1(1,2),Scr1(1,3),TSVec,Scr1(1,4),RxnFwd,
     $      DwnHil)
        else if(IMthd.eq.10) then
          If(IPrint.ge.1) write(IOut,2110)
          Call CorDrv(IOut,NAt3,NAtoms,RtMass,Cart(1,IStep-1,IPath),
     $      Cart(1,IStep,IPath),StpSiz,TSEVec,HsEVec(1,1,IPath),HvHess,
     $      IEig,RxnFwd,AGrad(1,1,IPath),IStep)
        else
          Write(IOut,9000) IMthd
          Call GauErr('L123 ERROR: ILLEGAL VALUE OF IMTHD!')
          endIf
        Call RxStat(IOut,IPrint,Gen,IRdStat,Stats(1,1,IPath),
     $    Cart(1,IStep-1,IPath),Cart(1,IStep,IPath),Potent(1,IPath),
     $    V(ICplng),RtMass,Scr1,NAtoms,NAt3,IStep,IDone,IPath,RxnFwd,0,
     $    NStat,NPath,StpSiz,IAn)
        If(IMthd.eq.5) StpSiz = TStep1
        endIf
      Call TStamp(1,'After Predictor')
C
C     If the integration along this path is complete (Done=.True.), then
C     take care of the final set of stats, set JobFlg to 2, and return to
C     Explor.
C
      If(Done) then
        If(IPrint.ge.0) then
          If(RxnFwd) then
            Write(IOut,5000)
          else
            Write(IOut,5001)
            endIf
          If(IPath.ge.NPath) then
            Write(IOut,5010)
          else
            Write(IOut,5020)
            endIf
          endIf
        IDone(IPath) = IStep - 1
        ITemp = NAt3*25
        RxnCoo = Stats(2,IStep-1,IPath)
        Call TStamp(1,'Bef OTIRC1 A')
        Call OTIRC1(IOut,1,NAt3,MaxStp,NPath,ITStep+1,
     $    Cart(1,IStep-1,IPath),Potent(1,IPath),AGrad(1,1,IPath),RxnCoo,
     $    Scr1,ITemp)
        If(SvPES) Call PESIRC(IOut,IPrint,3,IRwPES,NAt3,NAt3TT,HvHess,
     $    NPath,JJ,IPath,IStep-1,.True.,Cart(1,IStep-1,IPath),.True.,
     $    Potent(1,IPath),.True.,AGrad(1,1,IPath),HvHess,
     $    HesCar(1,1,IPath),.True.,RxnCoo,XX,LenObj,0,PESObj,OK)
        Call TStamp(1,'Aft OTIRC1 A')
        JobFlg = 2
        Return
        endIf
C
C     If we're carrying out a multi-step corrector step, then we may need
C     to index IStep to IStep-1 for book keeping. This is also where we
C     write the IRC data to the checkpoint file.
C
      If(.not.CorDon) then
        IStep = IStep-1
      else if((IPath.eq.1.or.IStep.ne.2)
     $  .and.(IStep-1).gt.IDone(IPath)) then
        ITemp = NAt3*25
        RxnCoo = Stats(2,IStep-1,IPath)
        Call TStamp(1,'Bef OTIRC1 B')
        Call OTIRC1(IOut,1,NAt3,MaxStp,NPath,ITStep+1,
     $    Cart(1,IStep-1,IPath),Potent(1,IPath),AGrad(1,1,IPath),RxnCoo,
     $    Scr1,ITemp)
        If(SvPES) Call PESIRC(IOut,IPrint,3,IRwPES,NAt3,NAt3TT,HvHess,
     $    NPath,JJ,IPath,IStep-1,.True.,Cart(1,IStep-1,IPath),.True.,
     $    Potent(1,IPath),.True.,AGrad(1,1,IPath),HvHess,
     $    HesCar(1,1,IPath),.True.,RxnCoo,XX,LenObj,0,PESObj,OK)
        Call TStamp(1,'Aft OTIRC1 B')
        endIf
      If(ICurWh.eq.0) then
        Call AMove(NAt3,Cart(1,IStep,IPath),CCur)
      else if(ICurWh.ne.1) then
        Write(IOut,9900) ICurWh
        Call GauErr('RxnCyc: Invalid ICurWh.')
        endIf
C
C     Move PotCur to PotSav and GrdCur to GrdSav so that it is available in
C     the next instance of the link.
C
      PotSav = PotCur
      Call AMove(NAt3,GrdCur,GrdSav)
C
C     Signal IRC not done.
C
      JobFlg = 0
      Return
      End
*Deck RxnEul
      Subroutine RxnEul(IOut,IPrint,NAt3,NAtoms,IStep,StpSiz,Cart1,
     $  Cart2,Grad,TSVec,RtMass,RxnFwd,DwnHil)
      Implicit Real*8(a-h,o-z)
C
C     This routine carries out Euler integration of the Rxn. Path Following
C     Equation.
C
C        - H.P. Hratchian, Spring 2003.
C
      Dimension Cart1(NAt3),Cart2(NAt3),Grad(NAt3),TSVec(NAt3),
     $  RtMass(NAt3)
      Logical RxnFwd,DwnHil
      Save One
      Data One/1.d0/
 1000 Format(' Gradient is small - USING TRANSITION VECTOR.')
C
C     Mass Weigh the Cartesians and gradient.
C
      If(IPrint.ge.2) then
        Call OutMtS(IOut,'RxnEul: Cart1 before MW''ing=',0,0,Cart1,NAt3,
     $    1,NAt3,1)
        Call OutMtS(IOut,'RxnEul: Grad before MW''ing=',0,0,Grad,NAt3,1,
     $    NAt3,1)
        endIf
      Call MWeigV(.False.,NAtoms,Cart1,RtMass,Cart1)
      Call MWeigV(.True.,NAtoms,Grad,RtMass,Grad)
      If(IPrint.ge.2) then
        Call OutMtS(IOut,'RxnEul: Cart1 after MW''ing=',0,0,Cart1,NAt3,
     $    1,NAt3,1)
        Call OutMtS(IOut,'RxnEul: Grad after MW''ing=',0,0,Grad,NAt3,1,
     $    NAt3,1)
        endIf
C
C     Now, do the integration.
C
      GMag = VLen(NAt3,Grad)
      If(.not.DwnHil.and.IStep.eq.2) then
        If(IPrint.ge.2) Call OutMtS(IOut,'RxnEul: TSVec=',0,0,TSVec,
     $    NAt3,1,NAt3,1)
        If(IPrint.ge.0) Write(IOut,1000)
        If(.not.RxnFwd) Call AScale(NAt3,-One,TSVec,TSVec)
        GMag = VLen(NAt3,TSVec)
        Call ACASB(NAt3,Cart1,TSVec,Cart2,-StpSiz/GMag)
        If(.not.RxnFwd) Call AScale(NAt3,-One,TSVec,TSVec)
      else
        Call ACASB(NAt3,Cart1,Grad,Cart2,-StpSiz/GMag)
        endIf
      If(IPrint.ge.2) Call OutMtS(IOut,
     $  'RxnEul: Cart2 before de-MW''ing=',0,0,Cart2,NAt3,1,NAt3,1)
C
C     Un-Mass Weigh the cartesians and gradient before leaving.
C
      Call MWeigV(.True.,NAtoms,Cart1,RtMass,Cart1)
      Call MWeigV(.True.,NAtoms,Cart2,RtMass,Cart2)
      Call MWeigV(.False.,NAtoms,Grad,RtMass,Grad)
      If(IPrint.ge.2) Call OutMtS(IOut,
     $  'RxnEul: Cart2 after de-MW''ing=',0,0,Cart2,NAt3,1,NAt3,1)
      Return
      End
*Deck RxPrnt
      Subroutine RxPrnt(IOut,IPrint,Stats,Cart,NAtoms,NAt3,IStep,IStep2,
     $  IDone,IPath,NPath,NStat,IAn)
      Implicit Real*8(A-H,O-Z)
C
C     This routine prints out summary information after each step is
C     completed.
C
C        - H.P. Hratchian, Spring 2003.
C
C     Stats contains information from each step.
C       Stats(1,i) = EPot   = Potential Energy
C       Stats(2,i) = TRxnC  = Reaction coordinate (amu**(1/2)*bohr)
C
      Parameter (MxAtPr=50)
      Dimension Stats(NStat,*),Cart(NAt3),IDone(NPath),IAn(NAtoms),JJ(1)
      Common /PhyCon/  PhyCon(30)
      Save JJ
      Data JJ/0/
 1000 Format(1X,'Point Number:',I4,10X,'Path Number:',I4)
 1020 Format(25X,'INITIAL STRUCTURE')
 1050 Format(19X,' CURRENT STRUCTURE')
 2000 Format('   CHANGE IN THE REACTION COORDINATE = ',F10.5,/,
     $       '   NET REACTION COORDINATE UP TO THIS POINT = ',F10.5,/,
     $       '  # OF POINTS ALONG THE PATH = ',I3,/,
     $       '  # OF STEPS = ',I3,/)
C
      If(IPrint.ge.0) Write(IOut,1000) IStep-2,IPath
C
C     If this is the initial step, just print out the inital coordinates
C     and then exit the routine.
C
      If(IStep.eq.2) then
        If(IPrint.ge.2) then
          Write(IOut,1020)
          Call CorPrS(IOut,'Cartesian Coordinates (Ang):',0,NAtoms,IAn,
     $      .False.,JJ,Cart,PhyCon(1))
          endIf
        Return
        endIf
C
C     Print out coordinates of the newest structure.
C
      If(IPrint.ge.2.or.(IPrint.ge.1.and.NAtoms.le.MxAtPr)) then
        Write(IOut,1050)
        Call CorPrS(IOut,'Cartesian Coordinates (Ang):',0,NAtoms,IAn,
     $    .False.,JJ,Cart,PhyCon(1))
        endIf
C
C     Print out some of the stats...
C
      If(NStat.lt.2) Call GauErr('RxPrnt: NStat too small !')
      ITStp = IStep-2
      If(IPath.gt.1) then
        Do 100 I = 1, (IPath-1)
  100     ITStp = IDone(I)+ITStp-1
        endIf
      If(IPrint.ge.0) Write(IOut,2000)
     $  Abs(Stats(2,IStep-1)-Stats(2,IStep-2)),Abs(Stats(2,IStep-1)),
     $  ITStp,IStep2
      Return
      End
*Deck RxStat
      Subroutine RxStat(IOut,IPrint,Gen,IRdStat,Stats,Cart0,Cart1,
     $  Potent,Couplng,RtMass,Scr,NAtoms,NAt3,IStep,IDone,IPath,RxnFwd,
     $  IFlag,NStat,NPath,StpSiz,IAnZ)
      Implicit Real*8(A-H,O-Z)
C
C     This routine computes and prints some useful statistics at each step
C     for L123.
C
C        - H.P. Hratchian, Spring 2003.
C
C
C     IFlag = 0  Normal Step (Standard print out - needed for GaussView animation.)
C     IFlag = 1  Predictor Step.
C     IFlag = 2  Completed Corrector Step.
C     IFlag = 3  Job Done - only wrap up stats for last point.
C
C
C     Stats saves certain information from each step
C       Stats(1,i)       = EPot   = Potential Energy
C       Stats(2,i)       = TRxnC  = Reaction coordinate (amu**(1/2)*bohr)
C       Stats(3,i)       = TotTim = Total time step for DVV (fs)
C       Stats(4:NStat,i) = User-defined values.
C
      Dimension Gen(*),IRdStat(4,*),Stats(NStat,*),Cart0(NAt3),
     $  Cart1(NAt3),Couplng(*),RtMass(NAtoms),Scr(NAt3,4),IDone(NPath),
     $  IAnZ(NAtoms),JJ(1)
      Common /PhyCon/  PhyCon(30)
      Logical RxnFwd
      Character*36 StpTyp(0:2)
      Save JJ,StpTyp,Zero
      Data JJ/0/,StpTyp/'*****PARTIAL STEP SUMMARY*****',
     $  '*****PREDICTOR STEP SUMMARY*****',
     $  '*****CORRECTOR STEP SUMMARY*****'/,Zero/0.d0/
 1000 Format(1x,'INFO FROM RxStat')
 1005 Format(19X,a)
 1050 Format(19X,' NEW STRUCTURE')
 1051 Format(19x,' PREDICTED STRUCTURE')
 1052 Format(19x,' CORRECTED STRUCTURE')
 9999 Format(' ')
C
C     Begin by printing out what type of step we're reporting statistics on
C     if the user has asked for extra printing (i.e. IPrint>1).
C
C     Note: We don't print out anything either way if the job is done since
C           all we're doing is filling in the some final data.
C
      If(IPrint.ge.2.and.IFlag.ne.3) then
        Write(IOut,9999)
        Write(IOut,1000)
        Write(IOut,1005) StpTyp(IFLag)
        endIf
C
C     Save the potential energy for the last position (IFlag=3).
C
      If(IFlag.eq.3) then
        Stats(1,IStep-1) = Potent
        Return
        endIf
C
C     Print out coordinates of the newest structure.
C
      If(IPrint.ge.2) then
        If(IFlag.eq.1) then
          Write(IOut,1051)
        else if(IFlag.eq.2) then
          write(IOut,1052)
        else
          write(IOut,1050)
          endIf
        Call CorPrS(IOut,'Cartesian Coordinates (Ang):',0,NAtoms,
     $    IAnZ,.False.,JJ,Cart1,PhyCon(1))
        Write(IOut,9999)
        endIf
C
C     Calculate the standard stats.
C
      Call MWeigV(.False.,NAtoms,Cart1,RtMass,Scr)
      Call MWeigV(.False.,NAtoms,Cart0,RtMass,Scr(1,2))
      Call ASub(NAt3,Scr,Scr(1,2),Scr)
      RxnCoo = Sqrt(sprod(NAt3,Scr,Scr))
      If(.not.RxnFwd) RxnCoo = -RxnCoo
      If(IStep.eq.2) then
        TRxnC = RxnCoo
        TotTim = StpSiz
      else
        TRxnC = Stats(2,IStep-1)+RxnCoo
        TotTim = Stats(3,IStep-1)+StpSiz
        endIf
C
C     Now, fill the array Stats...Remember that we're filling in the energy
C     for the last step and that we don't fill the energy if this is a
C     completed corrector step (IFlag = 2) - fill corrector energies when
C     the next predictor step is being reported.
C
      If(NStat.lt.3) Call GauErr('RxStat: NStat too small !')
      If(IFlag.eq.2) then
        Stats(1,IStep) = Potent
      else
        Stats(1,IStep-1) = Potent
        endIf
      Stats(2,IStep) = TRxnC
      Stats(3,IStep) = TotTim
C
C     Fill elements in the array Stats corresponding to user-defined stats.
C
      Do 100 i = 4,NStat
        Stats(i,IStep) = RStat1(IOut,.True.,i-3,NStat-3,NAtoms,Gen,
     $    Couplng,IRdStat,Cart1)
  100   Continue
C
C     When IStep=2, we take care of Stats(:,2) using all of the standard
C     code here, but we also take care of Stats(:,1). The special care
C     required for filling Stats(:,1) is done in this block. Note that the
C     TS statistics are only taken care of if IFlag.ne.2 since some of Gen
C     and other data arrays may well refer to the current point being
C     considered by corrector routines, etc.
C
      If(IStep.eq.2.and.IFlag.ne.2) then
        Stats(2,1) = Zero
        Stats(3,1) = Zero
        Do 200 i = 4,NStat
          Stats(i,1) = RStat1(IOut,.False.,i-3,NStat-3,NAtoms,Gen,
     $      Couplng,IRdStat,Cart0)
  200     Continue
        endIf
C
C     Print out statistical information.
C
      ITStp = IStep-1
      If(IPath.gt.1) then
        Do 900 I = 1, (IPath-1)
          ITStp = IDone(I) + ITStp
          If(I.ge.2) ITStp = ITStp-1
  900     Continue
        endIf
      Return
      End
*Deck RxStif
      Subroutine RxStif(IOut,NAt3,NAt3TT,SpecRd,Grad,Hess,Scr1,Scr2)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to evaluate the spectral radius of the
C     pre-factor matrix in the linearized ODE form of the (steepest
C     descent) reaction path equation [see Burger and Yang, JCP, 124,
C     224102 (2006), Eq. 6].
C
C     The dummy argument SpecRd is returned from this routine as the
C     spectral radius of the matrix:
C                                    __                           __
C                                    |                             |
C                                    |    H      g Transpose[Hg]   |
C           SpecRd = Spectral Radius |  ----- - -----------------  |
C                                    |   |g|          |g|^3        |
C                                    |_                           _|
C
C     Grad and Hess are the gradient and Hessian, respectively.
C
C     Scr1 is a scratch array that is NAt3*4 long.
C
C     Scr2 is a scratch array that is NAt3*NAt3 long.
C
C
C        - H.P. Hratchian, January 2007.
C
C
      Dimension Grad(NAt3),Hess(NAt3TT),Scr1(NAt3,4),Scr2(NAt3,NAt3,*)
 1000 Format(1x,'Entering Routine RxStif.')
C
      Write(IOut,1000)
C
C     Build the matrix we will diagonalize.
C
      Temp = GFloat(1)/VLen(NAt3,Grad)
      Call MultLV(.True.,1,0,NAt3,NAt3,Hess,Grad,Scr1)
      Call MatMpB(1,0,1,NAt3,1,NAt3,NAt3,Grad,Scr1,Scr2)
      Call AScale(NAt3*NAt3,Temp*Temp*Temp,Scr2,Scr2)
      Call MatOut(Scr2,NAt3,NAt3,NAt3,NAt3)
      Call Square(Hess,Scr2(1,1,2),NAt3,NAt3,0)
      Call AScale(NAt3*NAt3,Temp,Scr2(1,1,2),Scr2(1,1,2))
      Call ASub(NAt3*NAt3,Scr2(1,1,2),Scr2,Scr2)
      Call MatOut(Scr2,NAt3,NAt3,NAt3,NAt3)
C
C     Get the eigenvalues and spectral radius of the matrix in Scr2.
C
      Call RG(NAt3,NAt3,Scr2,Scr1,Scr1(1,2),0,Scr2(1,1,2),Scr1(1,3),
     $  Scr1(1,4),IErr)
      SpecRd = ArrMax(Scr1,NAt3,.True.)
      Return
      End
*Deck Shp123
      Subroutine Shp123(IOut,IPrint,IFitOp,IRwFit,IStep2,NAtoms,NAt3,
     $  DoOni,NAtQ3,NAtQ3T,IFtPt,NFtPts,NDtPts,MaskON,CCur,Cart1,Cart2,
     $  ECur,E1,E2,GCur,Grad1,Grad2,HCur,Hess1,Hess2,HUpFit,thetaUp,
     $  phiUp,vUpdat,wUpdat,EQM1,EQM2,GradQ1,GradQ2,LShpOb,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine prepares the input data for a distance weighted
C     interpolant fit for Link 123.  Basically, this routine is responsible
C     for transfering the data form used throughout l123 into the form
C     necessary for the DWI routines.
C
C     The work in this routine typically involves a series of AMove calls.
C     In the case where ONIOM(QM:MM) IRC code is indicated, the ONIOM mask
C     array is used to load the DWI arrays.
C
C     Input Dummy Arguments:
C           IOut              The standard output file unit number.
C           IPrint            The standard print integer flag.
C           IFitOp            This input option flag indicates whether the
C                             current set of PES data is to be appended to
C                             the fitting data set (IFitOp=2) or if the
C                             fitting set should be initialized and only
C                             the first and second point information passed
C                             in is put into the fitting set (IFitOp=1).
C           IRwFit            This is the RWF file number where Shepard
C                             fitting data should be stored via Routine
C                             PESShp.
C           IStep2            This is the current corrector iteration step
C                             number. Note that for all IFitOp values, when
C                             IStep2=1, the fitting data set is initialized
C                             and only the first and second input data
C                             points are loaded into the data set.
C           NAtoms            The number of "real" atoms.
C           NAt3              3*NAtoms.
C           DoOni             This is a logical flag that indicates whether
C                             or not QM:MM codes should be used in
C                             processing the DWI fitting.
C           NAtQ3             3*NAtoms(QM centers).  This is ignored if
C                             DoOni=.FASLE..
C           NAtQ3T            NAtQ3*(NAtQ3+1)/2.  This is ignored if
C                             DoOni=.FALSE..
C           IFtPt             This input integer gives the data point
C                             number where the current input data should be
C                             added in the PES fitting data set. NOTE that
C                             IFtPt is ignored if IFitOp=1.
C           NFtPts            This input integer gives the number of points
C                             that are to be used in the later calls to
C                             fitting routines.
C           NDtPts            This input integer gives the maximum number
C                             of points in the data set of points that can
C                             be used in PES fitting.
C           MaskON            This is the ONIOM masking integer array.
C                             This is ignored if DoOni=.FALSE..
C           Cart1,Cart2       Cartesian coordinates of the two points that
C                             will be used in the DWI fitting.  The length
C                             of each set of coordinates should be over the
C                             full set of coordinates, NAt3.
C           Grad1,Grad2       Gradients of the two points that will be used
C                             in the DWI fitting.  The length of each set
C                             of gradients should be over the full set of
C                             coordinates, NAt3.
C           Hess1,Hess2       Hessian matrices of the two points that will
C                             be used in the DWI fitting.  The length of
C                             each set of force constants should either be
C                             NAt3TT or NAtQ3T.  NAt3TT is used if the
C                             current DWI surface is being fitted over all
C                             atoms; NAtQ3T is used if the the DWI surface
C                             is being fit over only QM atoms (i.e., in
C                             QM:MM IRC runs).
C           HUpFit            This is an input logical argument that is
C                             sent as TRUE if the Hessian update vectors
C                             and scalars are being used for point 2
C                             Hessian input. If HUpFit is FALSE, thetaUp,
C                             phiUp, vUpdat, and wUpdat are ignored.
C           thetaUp,phiUp     These are the Hessian update scalars. The
C                             scalar thetaUp is the prefactor for vector
C                             vUpdat; scalar phiUp is the prefactor for
C                             vector wUpdat.  Note that these arguments are
C                             ignored for HUpFit=FALSE.
C           vUpdat,wUpdat     These are the Hessian update vectors v and w.
C                             The length of these vectors is NAt3 or NAtQ3
C                             for DoOni=TRUE.  Note that these vectors are
C                             ignored for HUpFit=FALSE.
C           EQM1,EQM2         These are the QM contributions to the QM:MM
C                             energies for the two points being fitted.
C                             These are ignored if DoOni=.FALSE..
C           GradQ1,GradQ2     Gradients of the two points that will be used
C                             in the DWI fitting when QM:MM codes are being
C                             employed.  These arrays are ignored if
C                             DoOni=.FALSE..
C           V,MDV             V is a scratch array of length MDV.
C
C        - H. P. Hratchian, 2006.
C        - H. P. Hratchian, 2008 - HUpFits for QM:MM IRC.
C
      Logical DEBUG,PESOk,HUpFit,DoOni
      Character*(1) CJunk
      Dimension JJ(1),XX(1),MaskON(*),CCur(*),Cart1(*),Cart2(*),GCur(*),
     $  Grad1(*),Grad2(*),HCur(*),Hess1(*),Hess2(*),vUpdat(*),wUpdat(*),
     $  GradQ1(*),GradQ2(*),V(MDV)
      Save Zero,JJ,XX,CJunk
      Data Zero/0.d0/,JJ/0/,XX/0.d0/,CJunk/' '/
 2000 Format(1x,'Shp123: E1=',F25.8,' E2=',F25.8)
 2100 Format(1x,'Shp123: ECur=',F25.8)
 8000 Format(1x,'At the end of Routine Shp123, the Shepard fitting ',
     $  'data object definitions:')
 9000 Format(1x,'Shp123: IFitOp=',I2)
 9010 Format(1x,'Shp123: MDtPts=',I10,' NDtPts=',I10)
C
C     Initialize things, and do some checks.
C
      DEBUG = .False.
      LShpOb = 0
      If(DoOni.and.NAtQ3T.ne.(NAtQ3*(NAtQ3+1)/2))
     $  Call GauErr('Shp123: Problem with NAtQ3 and NAtQ3T!')
      If(IFitOp.ne.1.and.IFitOp.ne.2) then
        Write(IOut,9000) IFitOp
        Call GauErr('Shp123: Invalid IFitOp.')
        endIf
      If(NDtPts.lt.2.and.IFitOp.eq.1)
     $  Call GauErr('Shp123 is confused: NDtPts < 2 but IFitOp=1.')
      If(IFtPt.gt.NDtPts) Call GauErr('Shp123: IFtPt > NDtPts.')
C
C     Begin by setting LenC.
C
      If(DoOni) then
        LenC = NAtQ3
      else
        LenC = NAt3
        endIf
C
C     If requested, print out the input data.
C
      If(DEBUG.or.IPrint.ge.2) then
        Write(IOut,2000) E1,E2
        Call OutMtS(IOut,'Shp123 - Cart1:',0,0,Cart1,Lenc,1,lenc,1)
        Call OutMtS(IOut,'Shp123 - Cart2:',0,0,Cart2,Lenc,1,lenc,1)
        Call OutMtS(IOut,'Shp123 - Grad1:',0,0,Grad1,Lenc,1,lenc,1)
        Call OutMtS(IOut,'Shp123 - Grad2:',0,0,Grad2,Lenc,1,lenc,1)
        Call LtOutS(IOut,'Shp123 - Hess1:',0,LenC,Hess1,0)
        Call LtOutS(IOut,'Shp123 - Hess2:',0,LenC,Hess2,0)
        Write(IOut,2100) ECur
        Call OutMtS(IOut,'Shp123 - CCur:',0,0,CCur,Lenc,1,lenc,1)
        Call OutMtS(IOut,'Shp123 - GCur:',0,0,GCur,Lenc,1,lenc,1)
        Call LtOutS(IOut,'Shp123 - HCur:',0,LenC,HCur,0)
        endIf
C
C     Based on the value of IFitOp, initialize or load the Shepard
C     interpolation data structure.
C
      IShpOb = 1
      IHUpRef = 0
      If(IStep2.eq.1.or.IFitOp.eq.1) then
        Call PESShp(IOut,IPrint,0,IRwFit,LenC,2,NDtPts,0,.False.,XX,
     $    .False.,Zero,.False.,XX,.False.,XX,.False.,Zero,Zero,XX,XX,
     $    .False.,IHUpRef,XX,JJ,LShpOb,MDV-IShpOb+1,V(IShpOb),PESOk)
      else if(IFitOp.eq.2) then
        MDtPts = 0
        Call PESShp(IOut,IPrint,2,IRwFit,LenC,2,MDtPts,0,.False.,XX,
     $    .False.,Zero,.False.,XX,.False.,XX,.False.,Zero,Zero,XX,XX,
     $    .False.,IHUpRef,XX,JJ,LShpOb,MDV,V(IShpOb),PESOk)
        If(MDtPts.ne.NDtPts) then
          Write(IOut,9010) MDtPts,NDtPts
          Call GauErr('Shp123: MDtPts.ne.NDtPts.')
          endIf
      else
        Write(IOut,9000) IFitOp
        Call GauErr('Shp123: Invalid IFitOp.')
        endIf
      IV     = IShpOb + LShpOb
      Call TstCor(IV-1,MDV,'Shp123-1')
C
C     Add new data to the data structure.
C
      If(IStep2.eq.1.or.IFitOp.eq.1) then
        If(DoOni) then
          ICart = IV
          IV1   = ICart + LenC
          Call TstCor(IV1-1,MDV,'Shp123-2')
          Call CmpVc(1,NAtoms,LDWIC1,NAtoms,Cart1,V(ICart),MaskON,1)
          If(LDWIC1.ne.LenC) Call GauErr('DWI123: Masking error 1.')
          IHUpRef = 1
          If(HUpFit) then
            Call TstCor(IV1+LenC-1,MDV,'Shp123-3')
            Call MultLV(.True.,1,0,LenC,LenC,Hess1,V(ICart),V(IV1))
            endIf
          Call PESShp(IOut,IPrint,3,IRwFit,LenC,2,NDtPts,1,.True.,
     $      V(ICart),.True.,EQM1,.True.,GradQ1,.True.,Hess1,.False.,
     $      Zero,Zero,XX,XX,HUpFit,IHUpRef,V(IV1),JJ,LShpOb,MDV,
     $      V(IShpOb),PESOk)
          Call CmpVc(1,NAtoms,LDWIC1,NAtoms,Cart2,V(ICart),MaskON,1)
          If(LDWIC1.ne.LenC) Call GauErr('DWI123: Masking error 2.')
          If(HUpFit) Call MultLV(.True.,1,0,LenC,LenC,Hess1,V(ICart),
     $      V(IV1))
          Call PESShp(IOut,IPrint,3,IRwFit,LenC,2,NDtPts,2,.True.,
     $      V(ICart),.True.,EQM2,.True.,GradQ2,.True.,Hess2,HUpFit,
     $      thetaUp,phiUp,vUpdat,wUpdat,HUpFit,IHUpRef,V(IV1),JJ,LShpOb,
     $      MDV,V(IShpOb),PESOk)
        else
          IHUpRef = 1
          If(HUpFit) then
            Call TstCor(IV+LenC-1,MDV,'Shp123-3')
            Call MultLV(.True.,1,0,LenC,LenC,Hess1,Cart1,V(IV))
            endIf
          Call PESShp(IOut,IPrint,3,IRwFit,LenC,2,NDtPts,1,.True.,Cart1,
     $      .True.,E1,.True.,Grad1,.True.,Hess1,.False.,Zero,Zero,XX,XX,
     $      HUpFit,IHUpRef,V(IV),JJ,LShpOb,MDV,V(IShpOb),PESOk)
          If(HUpFit) Call MultLV(.True.,1,0,LenC,LenC,Hess1,CCur,V(IV))
          Call PESShp(IOut,IPrint,3,IRwFit,LenC,2,NDtPts,2,.True.,CCur,
     $      .True.,ECur,.True.,GCur,.True.,HCur,HUpFit,thetaUp,phiUp,
     $      vUpdat,wUpdat,HUpFit,IHUpRef,V(IV),JJ,LShpOb,MDV,V(IShpOb),
     $      PESOk)
          endIf
      else if(IFitOp.eq.2) then
        If(HUpFit) Call MultLV(.True.,1,0,LenC,LenC,Hess1,CCur,V(IV))
        Call PESShp(IOut,IPrint,3,IRwFit,LenC,2,NDtPts,NFtPts,.True.,
     $    CCur,.True.,ECur,.True.,GCur,.True.,HCur,HUpFit,thetaUp,phiUp,
     $    vUpdat,wUpdat,HUpFit,IHUpRef,V(IV),JJ,LShpOb,MDV,V(IShpOb),
     $    PESOk)
      else
        Write(IOut,9000) IFitOp
        Call GauErr('Shp123: Invalid IFitOp.')
        endIf
C
C     If needed (IFitOp=2), we save the PES data structure to the file.
C
      Call PESShp(IOut,IPrint,1,IRwFit,LenC,2,MDtPts,0,.False.,XX,
     $  .False.,Zero,.False.,XX,.False.,XX,.False.,Zero,Zero,XX,XX,
     $  .False.,IHUpRef,XX,JJ,LShpOb,MDV,V(IShpOb),PESOk)
C
C     If the print level is high enough, print the definition table for
C     V(IShpOb).
C
      If(DEBUG.or.IPrint.ge.2) then
        IDef = 0
        Write(IOut,8000)
        Call PESDef(IOut,IPrint,99,IDef,' ',0,0,0,V(IShpOb),JJ,CJunk)
        endIf
      Return
      End
*Deck StpSet
      Subroutine StpSet(IOut,IPrint,NAt3,NAtoms,TSVec,AGrad,StpSiz,
     $  RtMass,Scr1,DwnHil)
      Implicit real*8(a-h,o-z)
C
C     This routine is used to set StpSiz such that the step taken on a
C     mass-weighted PES is equal to the user's desired step size on a non
C     mass-weighted PES.
C
C     NOTE: The user's desired step size is given in the input file by
C           IOP(1/8).  If the value given is <0, the value of StpSiz set
C           in this routine is equal to -IOP(1/8).  I.e., a negative input
C           for IOP(1/8) indicates that the user has given the desired
C           step size for a mass-weighted PES.
C
C        - H.P. Hratchian, Summer 2004.
C
      Dimension TSVec(NAt3),AGrad(NAt3),RtMass(NAtoms),Scr1(NAt3)
      Logical DwnHil
      Save Zero
      Data Zero/0.d0/
 1000 Format(1x,'Entering Routine StpSet.')
 1010 Format(1x,'StpSet: ALenMW=',F8.4,'  ALenCr=',F8.4)
 2000 Format(1x,'Supplied step size of ',f8.4,' sqrt(amu)*bohr used ',
     $  'for integration.')
 2010 Format(1x,'Supplied step size of ',f8.4,' bohr.',/,
     $  4x,'Integration on MW PES will use step size of ',f8.4,
     $  ' sqrt(amu)*bohr.')
C
      If(IPrint.ge.3) then
        Write(IOut,1000)
        Call OutMtS(IOut,'StpSet: RtMass=',0,0,RtMass,NAtoms,1,NAtoms,1)
        endIf
C
C     Save the input step value (StpSiz) so that we can print it out later.
C
      OldStp = StpSiz
C
C     Test to see if StpSiz is negative.  If so, the user has provided the
C     desired step size in sqrt(amu)*bohr and we just set StpSiz = -StpSiz.
C
      If(StpSiz.lt.Zero) then
        StpSiz = -StpSiz
        Write(IOut,2000) StpSiz
C
C     Since the value of StpSiz is positive, the user-provided step size
C     value (StpSiz) is in bohr.  Using the transition vector (TSVec), we
C     determine the step size on the mass weighted PES in sqrt(amu)*bohr.
C     Note that if this is a DOWNHILL calculation we use the gradient
C    (AGrad) and not the transition vector.
C
      else
        If(DwnHil) then
          ALenMW = VLen(NAt3,AGrad)
          If(IPrint.ge.3) Call OutMtS(IOut,'AGrad(MW) in StpSet:',0,0,
     $      AGrad,NAt3,1,NAt3,1)
          Call MWeigV(.True.,NAtoms,AGrad,RtMass,Scr1)
        else
          ALenMW = VLen(NAt3,TSVec)
          If(IPrint.ge.3) Call OutMtS(IOut,'TSVec(MW) in StpSet:',0,0,
     $      TSVec,NAt3,1,NAt3,1)
          Call MWeigV(.True.,NAtoms,TSVec,RtMass,Scr1)
          endIf
        ALenCr = VLen(NAt3,Scr1)
        StpSiz = StpSiz*ALenMW/ALenCr
        If(IPrint.ge.3) then
          Call OutMtS(IOut,'TSVec/AGrad(non-MW) in StpSet:',0,0,Scr1,
     $      NAt3,1,NAt3,1)
          Write(IOut,1010) ALenMW,ALenCr
          endIf
        Write(IOut,2010) OldStp,StpSiz
        endIf
      Return
      End
*Deck Summry
      Subroutine Summry(IOut,IDone,DoRelE,DoCnv,TSPot,Stats,NPath,
     $  MaxStp,NStat,IRdStat,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine writes out the summary of information for reaction path
C     calculations.
C
C        - H.P. Hratchian, Spring 2003.
C
      Parameter (NCorCl=DEFMAXREDINDEX)
      Logical Revers,DoRelE,DoCnv
      Dimension IDone(NPath),Stats(NStat,MaxStp,NPath),
     $  IRdStat(NCorCl,*),V(MDV),JJ(1)
      Save Zero,JJ
      Data Zero/0.d0/,JJ/0/
 1000 Format(1x,74('-'),/,
     $       4x,'Summary of reaction path following',/,
     $       1x,74('-'))
 1020 Format(1x,74('-'))
 2000 Format(1x,'Energies reported relative to the TS energy of ',F18.6)
 2010 Format(1x,'NOTE: Coordinate values are given in Bohr/Rad.')
C
C     Set up some flags and allocate space in V for the version of Stats
C     that will be printed.
C
      IEPnt1 = 0
      IEPnt2 = 0
      If(NPath.eq.1) then
        ICycle = IAbs(IDone(1))
        Revers = Stats(2,ICycle,1).lt.0
        If(Revers) then
          IEPnt2 = IDone(1)
          IOpTab = 2
        else
          IEPnt1 = IDone(1)
          IOpTab = 3
          endIf
        NRows = IDone(1)
        NStart = 1
      else
        Revers = .False.
        IEPnt1 = IDone(1)
        IEPnt2 = IDone(2)
        IOpTab = 1
        NRows = IDone(1)+IDone(2)-1
        NStart = 2
        endIf
C
C     If there are any coordinate values to print in the summary table
C     below, then print out a table with the coordinate definitions here.
C     The space allocated for the final stats list, IStat1, is used here to
C     temporarily store the TS coordinate values.
C
      ITS = 1
      IR  = ITS+NStat-3
      IP  = IR+NStat-3
      IStat1 = IP+NStat-3
      LabCol = IStat1 + NStat*MaxStp*NPath
      IV     = LabCol
      Call TstCor(IV-1,MDV,'Summry')
      If(NStat.ge.4) then
        ICur = 0
        Do 200 i = 4,NStat
          V(ITS+i-4) = Stats(i,1,1)
          If(NPath.eq.2) then
            V(IP+i-4) = Stats(i,IEPnt1,1)
            V(IR+i-4) = Stats(i,IEPnt2,2)
          else if(Revers) then
            V(IR+i-4) = Stats(i,IEPnt1,2)
          else
            V(IP+i-4) = Stats(i,IEPnt1,1)
            endIf
 200      Continue
        Call CorTab(IOut,IOpTab,4,NStat-3,IRdStat,V(ITS),V(IR),V(IP),
     $    V(IStat1+NStat))
        endIf
C
C     Fill the summary table row labesl (stored in RowCol). Then, fill the
C     summary table column labels (stored in LabCol).  We also print out
C     the definitions of any user-defined stats.  This allows for
C     convenient reading and interpretation of the summary file by the
C     user.
C
      MaxLab = MDV-LabCol+1
      ICur = 0
      Call PutCSt('Energy',V(LabCol),ICur)
      Call PutDel(2,V(LabCol),ICur)
      Call PutCSt('RxCoord',V(LabCol),ICur)
      Call PutDel(2,V(LabCol),ICur)
      Call CorLab(-1,.True.,JJ,JJ,ICur,2,MaxLab)
      Do 400 i = 4,NStat
        Call CorLab(0,.False.,IRdStat(1,i-3),V(LabCol),ICur,2,MaxLab)
  400   Continue
C
C     Fill the summary table data matrix V(IStat1).
C
      If(DoRelE) then
        RelE = TSPot
      else
        RelE = Zero
        endIf
      Do 500 ICol = 1,NStat
        IRow = 1
        Do 510 IPath = NPath,1,-1
          ICycle = IAbs(IDone(IPath))
          If(Stats(2,ICycle,IPath).lt.0) then
            Do 520 IStep = ICycle,1,-1
              If(ICol.eq.1) then
                V(IStat1-1+(ICol-1)*NRows+IRow) =
     $            Stats(ICol,IStep,IPath)-RelE
              else if(ICol.ge.4.and.DoCnv) then
                V(IStat1-1+(ICol-1)*NRows+IRow) =
     $            Stats(ICol,IStep,IPath)*CorCnv(IRdStat(1,ICol-3))
              else
                V(IStat1-1+(ICol-1)*NRows+IRow) =
     $            Stats(ICol,IStep,IPath)
                endIf
              IRow = IRow+1
  520         Continue
          else
            Do 530 IStep = NStart,ICycle
              If(ICol.eq.1) then
                V(IStat1-1+(ICol-1)*NRows+IRow) =
     $            Stats(ICol,IStep,IPath)-RelE
              else if(ICol.ge.4.and.DoCnv) then
                V(IStat1-1+(ICol-1)*NRows+IRow) =
     $            Stats(ICol,IStep,IPath)*CorCnv(IRdStat(1,ICol-3))
              else
                V(IStat1-1+(ICol-1)*NRows+IRow) =
     $            Stats(ICol,IStep,IPath)
                endIf
              IRow = IRow+1
  530         Continue
            endIf
  510     Continue
  500   Continue
C
C     For now, we never print the third column of the Stats array.  For
C     this reason we need to zap it from ScrMat.
C
      Do 600 i = 4,NStat
        II = (i-2)*NRows+IStat1
        IJ = (i-1)*NRows+IStat1
        Call AMove(NRows,V(IJ),V(II))
  600   Continue
      NCols = NStat-1
C
C     Print the summary table.
C
      If(DoRelE) Write(IOut,2000) TSPot
      If(.not.DoCnv) Write(IOut,2010)
      Write(IOut,1000)
      Call MatPt2(IOut,1,V(IStat1),1,1,1,NRows,NCols,NRows,NCols,0,1,JJ,
     $  V(LabCol),0,JJ,0,0,0,JJ,JJ)
      Write(IOut,1020)
      Return
      End
*Deck UCmpVc
      Subroutine UCmpVc(IOp,LenA,LenMsk,A,B,Mask)
      Implicit Real*8(A-H,O-Z)
C
C     This routine uncompresses a vector using a mask array as its guide.
C     Mask is an integer array with 0/1 for .TRUE./.FALSE.  A and B may be
C     the same.
C
C     The dummy arguments are:
C           IOp         This argument is a control flag.  Its allowed values are:
C                         0  The elements of Mask are in one-to-one
C                            correspondence with the elements of B.
C                         1  The elements of Mask flag for three elements
C                            in B at a time.  For example, Mask is over
C                            atomic centers while B is over Cartesian
C                            coordinates.
C           LenA        This is the length of the compressed vector, A.
C           LenMsk      This is the length of the masking logical array,
C                       Mask.
C           A           This is the compressed vector, dimensioned LenA.
C           B           This is the uncompressed vector.  Its length is
C                       determined from LenMsk and IOp.
C           Mask        This is the masking vector.  Its length is LenMsk.
C
C     - H. P. Hratchian, 2007.
C
      Dimension A(LenA),B(*),Mask(LenMsk)
C
C     Figure out the length of B (LenB), and then uncompress A.
C
      If(IOp.eq.0) then
        LenB = LenMsk
        IStpAB = 1
      else if(IOp.eq.1) then
        LenB = LenMsk*3
        IStpAB = 3
      else
        Call GauErr('Invalid IOp in UCmpVc.')
        endIf
C
      IPoinA = LenA-IStpAB+1
      IPoinB = LenB-IStpAB+1
      Do 100 I = LenMsk,1,-1
      If(Mask(I).eq.1) then
        Call AMove(IStpAB,A(IPoinA),B(IPoinB))
        IPoinA = IPoinA-IStpAB
        else
          Call AClear(IStpAB,B(IPoinB))
          endIf
        IPoinB = IPoinB-IStpAB
  100   Continue
      Return
      End

*Deck RxCurXX
      Subroutine RxCurvXX(IOut,NDim,Grad,Hess,RCurv,RxnCor,
     $ Natoms,AtMass,RtMass,VTang,VCurv,VScr)
      Implicit Real*8(A-H,O-Z)
C     
C     This routine computes the reaction path curvature vecture.
C     
C     -H.P. Hratchian, 2012. (original routine does mass-weighting
C     outside)
C     
C     Modified by EK, June 2015
C     
      Dimension Grad(*),Hess(*)
      Dimension AtMass(*), RtMass(*)
      Dimension VTang(*),VCurv(*),VScr(*)
      Save One
      Data One/1.d0/
 2000 Format(/,1x,'S, Scaled RP curvature(S) = ',F14.8,3x,F14.8,/)
C     
C     Do some basic initializations and mass-weight the input gradient
C     and Hessian
C     
      idump = 0
      Call MWeigV(.True.,NAtoms,Grad,RtMass,Grad)
      Call MWeigh(.True.,NAtoms,AtMass,Hess)
C     
C     Set up the tangent vector, originally V(ITang),
C     and the gradient magnitude, GMag
C     
      GMag = VLen(NDim,Grad)
      S1 = One/GMag
      Call AScale(NDim,-S1,Grad,VTang)
C     
C     Solve for the curvature vector,originally V(ICurv),
C     and its magnitudea CurvMag,
C     

      Call MultLV(.True.,1,0,NDim,NDim,Hess,VTang,VScr)
      S2 = SProd(NDim,VTang,VScr)
      Call ACASB(NDim,VScr,VTang,VCurv,-S2)
      Call AScale(NDim,-S1,VCurv,VCurv)
C     
C     Scaling of curvature vector
c
C     Call AScale(NDim,Scale,V(ICurv),V(ICurv))
C     Call OutMtS(IOut,'Curvature Vector',0,0,V(ICurv),NDim,1,NDim,1)
C     
      RCurv = VLen(NDim,VCurv)
      if (RxnCor.lt.1.d6)
     $ Write(IOut,2000) RxnCor, RCurv
C     
C     output of tangent and curvature
C     
      if (RxnCor.lt.1.d6.and.idump.eq.1) then
C     
      Write(IOut,6010)
 6010 Format(' Mass-weighted gradient ')
      Write(IOut,6000)
     $ (VTang(jj),jj=1,NDim)
      Write(IOut,6020)
 6020 Format(' Mass-weighted curvature vector ')
      Write(IOut,6000)
     $ (VCurv(kk),kk=1,NDim)
 6000 Format(3(F20.8))
      endif
C     
C     Undo mass-weighting of gradient and Hessian
C     
      Call MWeigV(.False.,NAtoms,Grad,RtMass,Grad)
      Call MWeigh(.False.,NAtoms,AtMass,Hess)
C     
      Return
      End
      
*Deck EleMas
      Function EleMas(Mode,IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     -----------------------------------------------------------------------
c     --- Mode = 0 : mass of most abundant isotope
c     --- .ne. 0 : averaged mass of isotopes
c     --- negative mass if IZ is out of range.
c     -----------------------------------------------------------------------
      parameter (maxza=120)
      dimension amas1(maxza),amas2(maxza)
c     most abundant isotopic mass
      data (amas1(i),i=1,60) /
      1 1.007825037d0, 4.002603250d0, 7.016004500d0,
      2 9.012182500d0, 11.009305300d0, 12.000000000d0,
      3 14.003074008d0, 15.994914640d0, 18.998403250d0,

      4 19.992439100d0, 22.989769700d0, 23.985045000d0,
      5 26.981541300d0, 27.976928400d0, 30.973763400d0,
      6 31.972071800d0, 34.968852729d0, 39.962383100d0,
      7 38.963707900d0, 39.962590700d0, 44.955913600d0,
      8 47.947946700d0, 50.943962500d0, 51.940509700d0,
      9 54.938046300d0, 55.934939300d0, 58.933197800d0,
      a 57.935347100d0, 62.929599200d0, 63.929145400d0,
      1 68.925580900d0, 73.921178800d0, 74.921595500d0,
      2 79.916520500d0, 78.918336100d0, 83.911506400d0,
      3 84.911700000d0, 87.905600000d0, 88.905400000d0,
      4 89.904300000d0, 92.906000000d0, 97.905500000d0,
      5 98.906300000d0, 101.903700000d0, 102.904800000d0,
      6 105.903200000d0, 106.905090000d0, 113.903600000d0,
      7 114.904100000d0, 117.901800000d0, 120.903800000d0,
      8 129.906700000d0, 126.900400000d0, 131.904200000d0,
      9 132.905429000d0, 137.905000000d0, 138.906100000d0,
      a 139.905300000d0, 140.907400000d0, 141.907500000d0/
       data (amas1(i),i=61,maxza) /
      1 144.912700000d0, 151.919500000d0, 152.920900000d0,
      2 157.924100000d0, 158.925000000d0, 163.928800000d0,
      3 164.930300000d0, 165.930400000d0, 168.934400000d0,
      4 173.939000000d0, 174.940900000d0, 179.946800000d0,
      5 180.948000000d0, 183.951000000d0, 186.956000000d0,
      6 189.958600000d0, 192.963300000d0, 194.964800000d0,
      7 196.966600000d0, 201.970600000d0, 204.974500000d0,
      8 207.976600000d0, 208.980400000d0, 208.982500000d0,
      9 210.987500000d0, 222.017500000d0, 223.019800000d0,
      a 226.025400000d0, 227.027800000d0, 232.038200000d0,
      1 231.035900000d0, 238.050800000d0, 237.048000000d0,
      2 242.058700000d0, 243.061400000d0, 246.067400000d0,
      3 247.070200000d0, 249.074800000d0, 252.082900000d0,
      4 252.082700000d0, 255.090600000d0, 259.101000000d0,
      5 262.109700000d0, 261.108700000d0, 262.114100000d0,
      6 266.121900000d0, 264.124700000d0, 265.000000000d0,
      7 268.138800000d0, 269.000000000d0, 272.000000000d0,
      8 277.000000000d0, 280.000000000d0, 280.000000000d0,
      9 280.000000000d0, 280.000000000d0, 280.000000000d0,
      a 280.000000000d0, 280.000000000d0, 280.000000000d0/
c      averaged isotopic mass
       data (amas2(i),i=1,60) /
      1 1.007940000d0, 4.002602000d0, 6.941000000d0,
      2 9.012183100d0, 10.811000000d0, 12.010700000d0,
      3 14.006700000d0, 15.999400000d0, 18.998403163d0,
      4 20.179700000d0, 22.989769280d0, 24.305000000d0,
      5 26.981538500d0, 28.085500000d0, 30.973761998d0,
      6 32.065000000d0, 35.453000000d0, 39.948000000d0,
      7 39.098300000d0, 40.078000000d0, 44.955908000d0,
      8 47.867000000d0, 50.941500000d0, 51.996100000d0,
      9 54.938044000d0, 55.845000000d0, 58.933194000d0,
      a 58.693400000d0, 63.546000000d0, 65.380000000d0,
      1 69.723000000d0, 72.640000000d0, 74.921595000d0,
      2 78.971000000d0, 79.904000000d0, 83.798000000d0,
      3 85.467800000d0, 87.620000000d0, 88.905840000d0,
      4 91.224000000d0, 92.906370000d0, 95.950000000d0,

      5 98.907200000d0, 101.070000000d0, 102.905500000d0,
      6 106.420000000d0, 107.868200000d0, 112.414000000d0,
      7 114.818000000d0, 118.710000000d0, 121.760000000d0,
      8 127.600000000d0, 126.904470000d0, 131.293000000d0,
      9 132.905451960d0, 137.327000000d0, 138.905470000d0,
      a 140.116000000d0, 140.907660000d0, 144.242000000d0/
       data (amas2(i),i=61,maxza) /
      1 144.900000000d0, 150.360000000d0, 151.964000000d0,
      2 157.250000000d0, 158.925350000d0, 162.500000000d0,
      3 164.930330000d0, 167.259000000d0, 168.934220000d0,
      4 173.054000000d0, 174.966800000d0, 178.490000000d0,
      5 180.947880000d0, 183.840000000d0, 186.207000000d0,
      6 190.230000000d0, 192.217000000d0, 195.084000000d0,
      7 196.966569000d0, 200.590000000d0, 204.383300000d0,
      8 207.200000000d0, 208.980400000d0, 208.982400000d0,
      9 209.987100000d0, 222.017600000d0, 223.019700000d0,
      a 226.024500000d0, 227.027700000d0, 232.037700000d0,
      1 231.035880000d0, 238.028910000d0, 237.048200000d0,
      2 239.064200000d0, 243.061400000d0, 247.070400000d0,
      3 247.070300000d0, 251.079600000d0, 252.083000000d0,
      4 257.059100000d0, 258.098400000d0, 259.101000000d0,
      5 262.109700000d0, 261.108800000d0, 262.114100000d0,
      6 266.121900000d0, 264.120100000d0, 265.000000000d0,
      7 268.138800000d0, 269.000000000d0, 272.000000000d0,
      8 277.000000000d0, 280.000000000d0, 280.000000000d0,
      9 280.000000000d0, 280.000000000d0, 280.000000000d0,
      a 280.000000000d0, 280.000000000d0, 280.000000000d0/
      save amas1,amas2

      if(IZ .lt. 1 .or. IZ .gt. maxza) then

      EleMas = -1.d0

      else if (Mode .eq. 0) then

      EleMas = amas1(IZ)

      else

      EleMas = amas2(IZ)

      end if

      return
      End
