PROGRAM BINARY
! *****************************************************************
!   MOLECULAR DYNAMICS SIMULATIONS FOR BINARY SYSTEMS
!   PROGRAM DATE   : 10th March 2001  ( THERMAL DIFFUSE METHOD)
!     COMPONENT 1 AND 2 = PERMEATING SPECIES (NONINTERACTING IN THIS
!     PROGRAM (I.E. IDEAL GASES))
!     COMPONENT 3 = SILICA (OXYGENS)
!     LENGTH SCALES IN ANGSTROMS
!     TIME SCALES IN PICOSECONDS
!     RMASS = G/GMOL (MOLECULAR WEIGHT)
!     THE GAS CONSTANT R = 0.8314 (G/GMOL).(A**2)/(PS**2)/K
!     ENERGY IS IN UNITS PER MOLE
!     INPUT 'CHEM'POT.' PROVIDE IDEAL GAS DENSITIES AS n(I)(i.g. SPECIES)
!              = EXP(CMUE(I)) (PARTICLES)/(ANGSTROMS**3)
!     DENSITY PROFILES ARE IN (NUMBER OF PARTICLES)/(ANGSTROMS**3)
!     FLUXES ARE IN (NUMBER OF PARTICLES)/(PICOSECONDS)
! *****************************************************************
!
      IMPLICIT REAL*8 (A - H,O - Z)
      PARAMETER(IZSI = 7500)
      PARAMETER(IZO = 15000)
      PARAMETER(IZ1 = 2000)
      PARAMETER(IZF = 4)
      PARAMETER(IZS = 2)
      PARAMETER(IZT = 6)
      DIMENSION VV(100)
      DIMENSION UOLD(IZF,IZ1),UNEW(IZF,IZ1),GDEN(IZF,500),TDEN(IZF,500)
      DIMENSION TMDEN(IZF,500),C(2,200,4)
      DIMENSION ISIL(1000),NO(IZS),NNTE(IZF)
      DIMENSION NMOV(IZF),NFAIL(IZF),NRMV(IZF),NADD(IZF),NMPP1(IZF)
      DIMENSION N(IZF),NM1(IZF)
      DIMENSION SIG(IZT),EPS(IZT),RMASS(IZF),RMASSI(IZF),RKCON(4)
      DIMENSION FACTT(IZF),SIGG(IZF,IZT),RCUT(IZF,IZT),COLLT(1000)
      DIMENSION EPSS(IZF,IZT),REPS(IZF,IZT),SQCUT(IZF,IZT)
      DIMENSION RLST(IZF,IZT),SQLST(IZF,IZT),CCUT1(IZF,IZT)
      DIMENSION CCUT2(IZF,IZT),SIGG3(IZF,IZT),FORFT(IZF,IZT)
      DIMENSION SQSIG(IZT,IZT),SQSIGI(IZT,IZT),SQMIN(IZF,IZT)
      DIMENSION AA(IZF),  AV(IZF), VVFACT(IZF)
      DIMENSION BINI(IZF),BINII(IZF)
      DIMENSION X0M(IZF,IZ1), Y0M(IZF,IZ1), Z0M(IZF,IZ1)
      DIMENSION X0(IZF,IZ1), Y0(IZF,IZ1), Z0(IZF,IZ1)
      DIMENSION X0O(IZF,IZ1),  Y0O(IZF,IZ1),  Z0O(IZF,IZ1)
      DIMENSION XNEW(IZF,IZ1),  YNEW(IZF,IZ1),  ZNEW(IZF,IZ1)
      DIMENSION BUFFXO(IZF,IZ1),  BUFFYO(IZF,IZ1),  BUFFZO(IZF,IZ1)
      DIMENSION BUFFXN(IZF,IZ1),  BUFFYN(IZF,IZ1),  BUFFZN(IZF,IZ1)
      DIMENSION VX(IZF,IZ1),  VY(IZF,IZ1),  VZ(IZF,IZ1), VP(IZF,IZ1)
      DIMENSION VNEW(10000),TFUNC(10000),RDFSUM(IZF)
      DIMENSION FX(IZF,IZ1),  FY(IZF,IZ1),  FZ(IZF,IZ1)
      DIMENSION AX(IZF,IZ1),  AY(IZF,IZ1),  AZ(IZF,IZ1)
      DIMENSION AXNEW(IZF,IZ1),  AYNEW(IZF,IZ1),  AZNEW(IZF,IZ1)
      DIMENSION DEN(IZF,101,500),TEMPR(IZF,101,500)
      DIMENSION NVI(IZF,IZ1),NVII(IZF,IZ1)
      DIMENSION NVOLI(IZF),NVOLII(IZF),NCV(IZF)
      DIMENSION RWMIN(IZF,IZT),RWMAX(IZF,IZT)
      DIMENSION LIST(500000),     NABOR(10000)
      DIMENSION RLISTCUT(IZF,IZT),RLIST(IZF,IZT)
      DIMENSION X1(IZSI),Y1(IZSI),Z1(IZSI)
      DIMENSION X2(IZO),Y2(IZO),Z2(IZO)
      DIMENSION NEIGH(IZO,10),NEIGHO(IZO,10)
      DIMENSION DIST(IZO,10),DISTO(IZO,10)
      DIMENSION NEIOO(IZO,10),NEIGOO(IZO,10)
      DIMENSION DISO(IZO,10),DISOO(IZO,10)
      DIMENSION NEO(IZSI),NES(IZO),NNM(IZF)
      DIMENSION BUFFX(IZF,IZ1),BUFFY(IZF,IZ1),BUFFZ(IZF,IZ1)
      DIMENSION BUFFXX(IZF,IZ1),BUFFYY(IZF,IZ1),BUFFZZ(IZF,IZ1)
      DIMENSION BUFFXT(IZF,IZ1),BUFFYT(IZF,IZ1),BUFFZT(IZF,IZ1)
      DIMENSION COXO(IZF,IZ1),COYO(IZF,IZ1),COZO(IZF,IZ1)
      DIMENSION SMTAU(IZF,500),DDELVV(500)
      DIMENSION SM1XY(IZF,201),SM1Z(IZF,201)
      DIMENSION TCX0(IZF,201),TCY0(IZF,201),TCZ0(IZF,201)
      DIMENSION TX0(IZF,500,201),TY0(IZF,500,201),TZ0(IZF,500,201)
      DIMENSION DRDF(IZF,IZS,200)
      COMMON /SILIC/ XOR1(IZSI),YOR1(IZSI),ZOR1(IZSI)
      COMMON /OXYG/ XOR2(IZO),YOR2(IZO),ZOR2(IZO)
      COMMON /SCO/ XS(IZS,10000),  YS(IZS,10000),  ZS(IZS,10000)
      COMMON /KPORE/ IDOUB,ISINGO,NCON,MCUT,N1,N2
      COMMON /PPORE/ CLT,SIO,RISF,PORR,PORIN,FREM,ROUTSQ,SQO
      COMMON /IRAND1/ ISXX,ISYY,ISZZ
      COMMON /NEIG/ NEIOOR(IZO,10)
! ****   READ 'input.dat'
!
      OPEN(UNIT = 5,FILE = 'siocoorddia.dat',STATUS = 'unknown', &
           FORM = 'FORMATTED',ACCESS = 'SEQUENTIAL')
      OPEN(UNIT = 10,FILE = 'sioinfktrian1.dat',STATUS = 'unknown', &
           FORM = 'FORMATTED',ACCESS = 'SEQUENTIAL')
      OPEN(UNIT = 14,FILE = 'sioinfktrian1.out',STATUS = 'unknown', &
           FORM = 'FORMATTED',ACCESS = 'SEQUENTIAL')
      READ(10,*) IDIFF,NF,NS,NCON
      READ(10,*) NSAVE0,NSM
      NFNS = NF + NS
      READ(10,*) TOL,NITER
      FIT = 1.0/FLOAT(NITER)
      READ(10,*) MXLONG, MXSHRT, IEQ, NPRINT
      do L = 1,NF
         READ(10,*) RMASS(L)
         WRITE(14,*) RMASS(L)
      end do
      do L = 1,NFNS
         READ(10,*) SIG(L),  EPS(L)
         WRITE(14,*) SIG(L),EPS(L)
      end do
      READ(10,*) RCUTT,  ROUT,  ZL
      READ(10,*) DLTT, TEMP
      READ(10,*) ISX,  ISY,  ISZ
      READ(10,*) ISXX,ISYY,ISZZ
      READ(10,*) ISV, ISW
      READ(10,*) IAX1,  IAY1,  IAZ1
      READ(10,*) IAX2,  IAY2,  IAZ2
      READ(10,*) PORR,PORIN,FREM
      READ(10,*) B0,RK0,rkang
      B0CO = 0.5*B0
      WRITE(14,*) B0,RK0
      ROUTSQ = ROUT*ROUT/(1.62*1.62)
! *******************************
! ** INPUT DATA FOR SILICA CUBE
! ****  READ RUN PARAMETERS AND CHARACTERISTICS OF SYSTEM
      READ(5,*) CLT
      READ(5,*) SIO
      READ(5,*) RISF
      PORRSQ = PORR*PORR
      PORINS = PORIN*PORIN
      SQO = SIO*SIO
! *********************************
      READ(5,*) N1
      do I = 1,N1
         READ(5,*) X1(I),Y1(I),Z1(I)
      end do
      READ(5,*) N2
      do I = 1,N2
         READ(5,*) X2(I),Y2(I),Z2(I)
      end do
      NSI = N1
      NOXX = N2
      CLX = CLT
      CLY = CLT
      CLZ = CLT
      CLMX = CLX/2.0
      CLNX = -CLX/2.0
      CLMY = CLY/2.0
      CLNY = -CLY/2.0
      CLMZ = CLZ/2.0
      CLNZ = -CLZ/2.0
!   OXYGEN NEIGHBORS TO A GIVEN SILICON ATOM
      do 999 J = 1,NSI
      NEO(J) = 0
      do 99 I = 1,NOXX
         DELX = X1(J) - X2(I)
         if (DELX > CLMX) DELX = DELX - CLX
         if (DELX <= CLNX) DELX = DELX + CLX
         if (DABS(DELX) > SIO) GOTO 99
         DELY = Y1(J) - Y2(I)
         if (DELY > CLMY) DELY = DELY - CLY
         if (DELY <= CLNY) DELY = DELY + CLY
         if (DABS(DELY) > SIO) GOTO 99
         SQXY = DELY*DELY + DELX*DELX
         if (SQXY > SQO) GOTO 99
         DELZ = Z1(J) - Z2(I)
         if (DELZ > CLMZ) DELZ = DELZ - CLZ
         if (DELZ <= CLNZ) DELZ = DELZ + CLZ
         if (DABS(DELZ) > SIO) GOTO 99
         SQXYZ = DELZ*DELZ + SQXY
         if (SQXYZ > SQO) GOTO 99
         NEO(J) = NEO(J) + 1
         NEIGH(J,NEO(J)) = I
         DIST(J,NEO(J)) = SQXYZ
99    CONTINUE
! ** ORDER THE DISTANCES
      NN = NEO(J)
      do 77 IP = 1,1000
         DD = DIST(J,1)
         IMIN = 1
         do 100 M = 2,NN
            if (DIST(J,M) < DD) then
               DD = DIST(J,M)
               IMIN = M
            end if
100      CONTINUE
         NEIGHO(J,IP) = NEIGH(J,IMIN)
         NEIOOR(J,IP) = NEIGH(J,IMIN)
         DISTO(J,IP) = DIST(J,IMIN)
         NEIGH(J,IMIN) = NEIGH(J,NN)
         DIST(J,IMIN) = DIST(J,NN)
         NN = NN - 1
         if (NN == 0) GOTO 771
77    CONTINUE
771   CONTINUE
999   CONTINUE
!   SILICON NEIGHBORS TO A GIVEN OXYGEN ATOM
      do 9995 I = 1,NOXX
      NES(I) = 0
      do 995 J = 1,NSI
         DELX = X1(J) - X2(I)
         if (DELX > CLMX) DELX = DELX - CLX
         if (DELX <= CLNX) DELX = DELX + CLX
         if (DABS(DELX) > SIO) GOTO 995
         DELY = Y1(J) - Y2(I)
         if (DELY > CLMY) DELY = DELY - CLY
         if (DELY <= CLNY) DELY = DELY + CLY
         if (DABS(DELY) > SIO) GOTO 995
         SQXY = DELY*DELY + DELX*DELX
         if (SQXY > SQO) GOTO 995
         DELZ = Z1(J) - Z2(I)
         if (DELZ > CLMZ) DELZ = DELZ - CLZ
         if (DELZ <= CLNZ) DELZ = DELZ + CLZ
         if (DABS(DELZ) > SIO) GOTO 995
         SQXYZ = DELZ*DELZ + SQXY
         if (SQXYZ > SQO) GOTO 995
         NES(I) = NES(I) + 1
         NEIOO(I,NES(I)) = J
         DISO(I,NES(I)) = SQXYZ
995   CONTINUE
! ** ORDER THE DISTANCES
      NN = NES(I)
      do 775 IP = 1,1000
         DD = DISO(I,1)
         IMIN = 1
         do 1010 M = 2,NN
            if (DISO(I,M) < DD) then
               DD = DISO(I,M)
               IMIN = M
            end if
1010     CONTINUE
         NEIGOO(I,IP) = NEIOO(I,IMIN)
         DISOO(I,IP) = DISO(I,IMIN)
         NEIOO(I,IMIN) = NEIOO(I,NN)
         DISO(I,IMIN) = DISO(I,NN)
         NN = NN - 1
         if (NN == 0) GOTO 776
775   CONTINUE
776   CONTINUE
9995  CONTINUE
      do 888 J = 1,NSI
         XOR1(J) = X1(J)
         YOR1(J) = Y1(J)
         ZOR1(J) = Z1(J)
888   CONTINUE
      do 881 J = 1,NOXX
         XOR2(J) = X2(J)
         YOR2(J) = Y2(J)
         ZOR2(J) = Z2(J)
881   CONTINUE
!
!***************************************
!
! OUTPUT FILE
!
      WRITE(14,*) MXLONG
!
! ****   SET UP ENERGY CUTOFF LIMIT
!
      DLTS2 = 0.5*DLTT*DLTT
      DLTSQ = DLTT*DLTT
      DLT01 = 0.1*DLTT
      DLT2I = 1.0/(2.0*DLTT)
      PI = 4.0D0*atan(1.D0)
      PI05 = 0.5*PI
      PI15 = 1.5*PI
      PI2 = 2.0*PI
      RGAS = 0.8314
      RTI = 1.0/(RGAS*TEMP)
      IPRINT = 0
      LTIM = MXLONG/10
      RCMAX = RCUTT*4.0
      CL = ZL
      CLP = CL/2.0
      DELZC = CLP/101.0
      ZLINV = 1.0/ZL
      CLINV = 1.0/CL
      CLN = -CLP
      CL2M = CLN
      VCELL = CL**3
      WRITE(14,*)'CL = ',CL,'ROUT = ',ROUT
!
      C13 = 1.0/3.0
      RK0 = RK0*RGAS
      rkang = rkang*RGAS
      do 1241 L = 1,NF
      do 1241 K = 1,NFNS
      SIGG(L,K) = 0.5*(SIG(L) + SIG(K))
      SIGG3(L,K) = (2.0**C13)*SIGG(L,K)**2
      EPSS(L,K) = sqrt(EPS(L)*EPS(K))
      REPS(L,K) = 4.0*EPSS(L,K)/TEMP
      EPSS(L,K) = RGAS*EPSS(L,K)
1241  CONTINUE
      do L = 1,NF
      do M = 1,NS
         RWMIN(L,M) = SIGG3(L,NF + M)
         RWMAX(L,M) = 2.0*RWMIN(L,M)
      end do
      end do
      do 22491 L = 1,NF
      RMASSI(L) = 1.0/RMASS(L)
      FACTT(L) = sqrt(RGAS*TEMP/RMASS(L))
      VVFACT(L) = RMASS(L)/(3.0*RGAS)
      do 2241 K = 1,NFNS
         SQSIG(L,K) = SIGG(L,K)**2
         SQSIGI(L,K) = 1.0/SQSIG(L,K)
         RCUT(L,K) = RCUTT*SIGG(L,K)
         RLIST(L,K) = 1.3*RCUT(L,K)
         SQCUT(L,K) = RCUT(L,K)**2
         RLISTCUT(L,K) = RLIST(L,K)**2
         SQMIN(L,K) = 0.25*SQSIG(L,K)
         RLST(L,K)   = 1.3* RCUT(L,K)
         SQLST(L,K) = RLST(L,K)**2
         SQRC2  = SQSIG(L,K)  / SQCUT(L,K)
         SQRC6  = SQRC2 * SQRC2 * SQRC2
         SQRC7 = SQRC6/RCUT(L,K)
         SQRC13 = SQRC6*SQRC6/RCUT(L,K)
         CCUT1(L,K) = EPSS(L,K)*48.0*(SQRC13 - 0.5*SQRC7)
         CCUT2(L,K) = EPSS(L,K)*SQRC6*(28.0 - (52.0*SQRC6))
         FORFT(L,K) =  48. *EPSS(L,K)*SQSIGI(L,K)
2241  CONTINUE
22491 CONTINUE
      RLCMIN = 0.25*((0.3*RCUT(1,3))**2)
      RDFMAX = 10000.0
      do L = 1,NF
      RDFSUM(L) = 0.0
      do K = 1,NFNS
      if (RDFMAX > SQCUT(L,K)) then
      RDFMAX = SQCUT(L,K)
      end if
      end do
      end do
      DELRDF = sqrt(RDFMAX)/100.0
      DLRDFI = 1.0/DELRDF
      AFRR = (4.0*PI/3.0)*(DELRDF**3)
      do JP = 1,100
      DDELVV(JP) = AFRR*(3.0*FLOAT(JP)*(FLOAT(JP) - 1.0) + 1.0)
      do L = 1,NF
      do K = 1,NS
      DRDF(L,K,JP) = 0.0
      end do
      end do
      end do
! ***********************************
      DVN = 0.001
      do LV = 1,5000
      VNEW(LV) = FLOAT(LV)*DVN
      VNSQ = VNEW(LV)**2
      VNEW(LV) = VNEW(LV)*sqrt(2.0D00)
      TFUNC(LV) = (VNSQ + 1.0)*exp(-VNSQ)
      end do
! **********************************
      ARGMAX = 50.0
      ARGMIN = -50.0
      RCMAX = RCUTT*4.0
      RNCST = 4.46216
      SQSIG(NFNS,NFNS) = SIG(NFNS)*SIG(NFNS)
!
      RMUL2 = 16807.D0
      RMO1 = 2147483648.D0
      RIMO1 = 1.D0/RMO1
      RMO2 = 2147483647.D0
      RMOA   = RMO1 / (2. * PI)
      RIMOA = 1.0D0/RMOA
!
! ****   SET UP DIMENSIONS
!
      NZ = 100
      DEZ = CL/(2.0*FLOAT(NZ - 1))
      DEZI = 1.0/DEZ
      DEZIZ = DEZI/2.0
!
! ****   SET REAL TIME STEP AND ITS MULTIPLE
!
      DLT2   = DLTT / 2.0
      DTSQ   = DLTT * DLTT
      DTSQ2  = DTSQ / 2.0
      DTSQ4  = DTSQ / 4.0
       do 5566 L = 1,NF
      do 5566 KZ = 1,NZ
      GDEN(L,KZ) = 0.0
5566  CONTINUE
      SMINT = DLTT*NSAVE0
       NSM14 = NSM/4
      NSM12 = NSM/2
      NSM34 = 3*NSM14
       do L = 1,NF
      do 50 I = 1,NSM
        SMTAU(L,I) = (I - 1.)*SMINT
50    CONTINUE
      do 120 I = 1,NSM
        SM1XY(L,I) = 0.
        SM1Z(L,I) = 0.
120   CONTINUE
       end do
       do L = 1,2
          do K = 1,101
            do M = 1,4
               C(L,K,M) = 0.0
            end do
           end do
       end do
          do LPP = 1,4
          RKCON(LPP) = 0.0
          end do
          TRY = 0.0
        do 88888 KPO = 1,NCON
      NPV0 = 0
      ISAVE0 = 0
      JSAVE0 = 0
! *************************************
! THE SILICA ATOM COORDINATES
       CALL SIOPOR
       NO(1) = IDOUB
       NO(2) = ISINGO
       WRITE(14,*) IDOUB,ISINGO
       WRITE(14,*) ISXX,ISYY,ISZZ
! *******************************
! **   INITIAL FLUID PARTICLE POSITIONS
       do L = 1,2
       N(L) = 0
       N(L + 2) = 0
       I2 = 0
       do I = 1,10000,2
       I2 = I2 + 1
       do LL = 1,100000000
       TRY = TRY + 1.0
       ISX = mod(RMUL2*ISX,RMO2)
       ISX = mod(RMUL2*ISX,RMO2)
       RN1 = ISX*RIMO1
       X0MM = CL*(0.5 - RN1)
       ISY = mod(RMUL2*ISY,RMO2)
       ISY = mod(RMUL2*ISY,RMO2)
       RN2 = ISY*RIMO1
       Y0MM = CL*(0.5 - RN2)
       ISW = mod(RMUL2*ISW,RMO2)
       ISW = mod(RMUL2*ISW,RMO2)
       RN3 = ISW*RIMO1
       Z0MM = CL*(0.5 - RN3)
       ISX = mod(RMUL2*ISX,RMO2)
       ISX = mod(RMUL2*ISX,RMO2)
       RN4 = ISX*RIMO1
       ISY = mod(RMUL2*ISY,RMO2)
       ISY = mod(RMUL2*ISY,RMO2)
       RN5 = ISY*RIMO1
       CA = RN4
       SA = sqrt(1.D0 - RN4*RN4)
       CB = cos(PI2*RN5)
       SB = sin(PI2*RN5)
       BL1 = 0.5*B0*SA*CB
       BL2 = 0.5*B0*SA*SB
       BL3 = 0.5*B0*CA
       X0(L,I) = X0MM + BL1
       X0(L,I + 1) = X0MM - BL1
       Y0(L,I) = Y0MM + BL2
       Y0(L,I + 1) = Y0MM - BL2
       Z0(L,I) = Z0MM + BL3
       Z0(L,I + 1) = Z0MM - BL3
       X0(L + 2,I2) = X0MM
       Y0(L + 2,I2) = Y0MM
       Z0(L + 2,I2) = Z0MM
       do M = 1,NS
       do LS = 1,NO(M)
       do K = 1,2
       RZ = Z0(L,I + K - 1) - ZS(M,LS)
       RZ = RZ - anint(RZ*CLINV)*CL
       RX = X0(L,I + K - 1) - XS(M,LS)
       RY = Y0(L,I + K - 1) - YS(M,LS)
       RX = RX - anint(RX*CLINV)*CL
       RY = RY - anint(RY*CLINV)*CL
       SQ = RX*RX + RY*RY + RZ*RZ
       if (SQ <= SQSIG(1,4 + M)) GOTO 50121
       end do
       RZ = Z0(L + 2,I2) - ZS(M,LS)
       RZ = RZ - anint(RZ*CLINV)*CL
       RX = X0(L + 2,I2) - XS(M,LS)
       RY = Y0(L + 2,I2) - YS(M,LS)
       RX = RX - anint(RX*CLINV)*CL
       RY = RY - anint(RY*CLINV)*CL
       SQ = RX*RX + RY*RY + RZ*RZ
       if (SQ <= SQSIG(3,4 + M)) GOTO 50121
       end do
       end do
       N(L) = N(L) + 2
       N(L + 2) = N(L + 2) + 1
       if (N(L) >= 40) GOTO 50123
       GOTO 50122
50121  CONTINUE
       end do
50122  CONTINUE
       end do
50123  CONTINUE
       end do
! ****************************************
! ****  INITIAL VELOCITY FOR MOVING PARTICLE
! ****
        VSUM = 0.0
        VSUMSQ = 0.0
      do 44220 L = 1,NF
         NM1(L) = N(L) - 1
       if (N(L) == 0) GOTO 44220
      do 44221 J = 1,N(L)
!
! ** INITIAL VELOCITY FROM MAXWELLIAN DISTRIBUTION *****
!
      IAX1 = mod(RMUL2*IAX1,RMO2)
      IAX1 = mod(RMUL2*IAX1,RMO2)
      IAX2 = mod(RMUL2*IAX2,RMO2)
      IAX2 = mod(RMUL2*IAX2,RMO2)
      IAY1 = mod(RMUL2*IAY1,RMO2)
      IAY1 = mod(RMUL2*IAY1,RMO2)
      IAY2 = mod(RMUL2*IAY2,RMO2)
      IAY2 = mod(RMUL2*IAY2,RMO2)
      IAZ1 = mod(RMUL2*IAZ1,RMO2)
      IAZ1 = mod(RMUL2*IAZ1,RMO2)
      IAZ2 = mod(RMUL2*IAZ2,RMO2)
      IAZ2 = mod(RMUL2*IAZ2,RMO2)
!
      RNX1 = IAX1*RIMO1
      RNX2 = IAX2*RIMO1
      RNY1 = IAY1*RIMO1
      RNY2 = IAY2*RIMO1
      RNZ1 = IAZ1*RIMO1
      RNZ2 = IAZ2*RIMO1
!
      VX(L,J) = FACTT(L)*sqrt(-2.D0*log(RNX1))*cos(PI2*RNX2)
      VY(L,J) = FACTT(L)*sqrt(-2.D0*log(RNY1))*cos(PI2*RNY2)
      VZ(L,J) = FACTT(L)*sqrt(-2.D0*log(RNZ1))*cos(PI2*RNZ2)
!
      VP(L,J) = sqrt(VX(L,J)**2 + VY(L,J)**2 + VZ(L,J)**2)
      VSUM = VSUM + VP(L,J)
      VSUMSQ = VSUMSQ + VVFACT(L)*VP(L,J)*VP(L,J)
44221 CONTINUE
      do 37010 J = 1,N(L)
      FX(L,J) = 0.0
      FY(L,J) = 0.0
      FZ(L,J) = 0.0
37010 CONTINUE
44220  CONTINUE
       WRITE(14,*) 'VSUMSQ = ',VSUMSQ
!
            SUMU = 0.0
! ****  EVALUATE FORCES DUE TO OXYGEN ATOMS
      do 59182 L = 1,NF
      do 59183 J = 1,N(L)
      do 59184 M = 1,NS
      do 34780 II = 1,NO(M)
      RZ = Z0(L,J) - ZS(M,II)
      RZ = RZ - anint(RZ*CLINV)*CL
      RX = X0(L,J) - XS(M,II)
      RY = Y0(L,J) - YS(M,II)
      RX = RX - anint(RX*CLINV)*CL
      RY = RY - anint(RY*CLINV)*CL
      SQ = RX*RX + RY*RY + RZ*RZ
! *********
!
                  if (SQ > SQCUT(L,NF + M)) GOTO 34780
                     RR1 = sqrt(SQ)
                     RR2 = SQSIG(L,NF + M) / SQ
                     RR6 = RR2 * RR2 * RR2
           POTEN = 4.  *EPSS(L,NF + M)* RR6*(RR6 - 1.)
           POTEN = POTEN+ (CCUT1(L,NF + M)*RR1) + CCUT2(L,NF + M)
           FORCE = FORFT(L,NF + M)*  RR6 * (RR6 - 0.5) * RR2
           FORCE = FORCE - (CCUT1(L,NF + M)/RR1)
                     SUMU  = SUMU + POTEN
                  FIJX  = FORCE * RX
                  FIJY  = FORCE * RY
                  FIJZ  = FORCE * RZ
                     FX(L,J) = FX(L,J) + FIJX
                     FY(L,J) = FY(L,J) + FIJY
                     FZ(L,J) = FZ(L,J) + FIJZ
34780 CONTINUE
59184 CONTINUE
59183 CONTINUE
59182 CONTINUE
! **  EVALUATE INTRAMOLECULAR C=O STRETCHING FORCES
      do L = 1,2
      I2 = 0
      do J = 1,N(L),2
      I2 = I2 + 1
! ** RELATIVE SEPARATION C and O
! ***
      XTEMP3 = X0(L + 2,I2)
      DX1 = X0(L,J) - XTEMP3
      if (DX1 > CLP) then
      XTEMP3 = XTEMP3 + CL
      DX1 = X0(L,J) - XTEMP3
      end if
      if (DX1 < CLN) then
      XTEMP3 = XTEMP3 - CL
      DX1 = X0(L,J) - XTEMP3
      end if
      XTEMP3 = X0(L + 2,I2)
      DX2 = X0(L,J + 1) - XTEMP3
      if (DX2 > CLP) then
      XTEMP3 = XTEMP3 + CL
      DX2 = X0(L,J + 1) - XTEMP3
      end if
      if (DX2 < CLN) then
      XTEMP3 = XTEMP3 - CL
      DX2 = X0(L,J + 1) - XTEMP3
      end if
! ***
      YTEMP3 = Y0(L + 2,I2)
      DY1 = Y0(L,J) - YTEMP3
      if (DY1 > CLP) then
      YTEMP3 = YTEMP3 + CL
      DY1 = Y0(L,J) - YTEMP3
      end if
      if (DY1 < CLN) then
      YTEMP3 = YTEMP3 - CL
      DY1 = Y0(L,J) - YTEMP3
      end if
      YTEMP3 = Y0(L + 2,I2)
      DY2 = Y0(L,J + 1) - YTEMP3
      if (DY2 > CLP) then
      YTEMP3 = YTEMP3 + CL
      DY2 = Y0(L,J + 1) - YTEMP3
      end if
      if (DY2 < CLN) then
      YTEMP3 = YTEMP3 - CL
      DY2 = Y0(L,J + 1) - YTEMP3
      end if
! ***
      ZTEMP3 = Z0(L + 2,I2)
      DZ1 = Z0(L,J) - ZTEMP3
      if (DZ1 > CLP) then
      ZTEMP3 = ZTEMP3 + CL
      DZ1 = Z0(L,J) - ZTEMP3
      end if
      if (DZ1 < CLN) then
      ZTEMP3 = ZTEMP3 - CL
      DZ1 = Z0(L,J) - ZTEMP3
      end if
      ZTEMP3 = Z0(L + 2,I2)
      DZ2 = Z0(L,J + 1) - ZTEMP3
      if (DZ2 > CLP) then
      ZTEMP3 = ZTEMP3 + CL
      DZ2 = Z0(L,J + 1) - ZTEMP3
      end if
      if (DZ2 < CLN) then
      ZTEMP3 = ZTEMP3 - CL
      DZ2 = Z0(L,J + 1) - ZTEMP3
      end if
! ***
      R21 = sqrt(DX1*DX1 + DY1*DY1 + DZ1*DZ1)
      FORCS = RK0*(1.D0 - (B0CO/R21))
      FIJX1 = -DX1*FORCS
      FIJY1 = -DY1*FORCS
      FIJZ1 = -DZ1*FORCS
            FX(L,J) = FX(L,J) + FIJX1
            FY(L,J) = FY(L,J) + FIJY1
            FZ(L,J) = FZ(L,J) + FIJZ1
            FX(L + 2,I2) = FX(L + 2,I2) - FIJX1
            FY(L + 2,I2) = FY(L + 2,I2) - FIJY1
            FZ(L + 2,I2) = FZ(L + 2,I2) - FIJZ1
      R21 = sqrt(DX2*DX2 + DY2*DY2 + DZ2*DZ2)
      FORCS = RK0*(1.D0 - (B0CO/R21))
      FIJX2 = -DX2*FORCS
      FIJY2 = -DY2*FORCS
      FIJZ2 = -DZ2*FORCS
            FX(L,J + 1) = FX(L,J + 1) + FIJX2
            FY(L,J + 1) = FY(L,J + 1) + FIJY2
            FZ(L,J + 1) = FZ(L,J + 1) + FIJZ2
            FX(L + 2,I2) = FX(L + 2,I2) - FIJX2
            FY(L + 2,I2) = FY(L + 2,I2) - FIJY2
            FZ(L + 2,I2) = FZ(L + 2,I2) - FIJZ2
       end do
       end do
!
! ******************************
! **** EVALUATE O=C=O BOND ANGLE BENDING
      do L = 1,2
      I2 = 0
      do J = 1,N(L),2
      I2 = I2 + 1
      rx1 = X0(L + 2,I2)
      rx2 = X0(L,J)
      rx3 = X0(L,J + 1)
      ry1 = Y0(L + 2,I2)
      ry2 = Y0(L,J)
      ry3 = Y0(L,J + 1)
      rz1 = Z0(L + 2,I2)
      rz2 = Z0(L,J)
      rz3 = Z0(L,J + 1)
      rxl1 = rx1 - rx2
      ryl1 = ry1 - ry2
      rzl1 = rz1 - rz2
      rxl2 = rx3 - rx1
      ryl2 = ry3 - ry1
      rzl2 = rz3 - rz1
            rxl1 = rxl1 - anint(rxl1*CLINV)*CL
            ryl1 = ryl1 - anint(ryl1*CLINV)*CL
            rzl1 = rzl1 - anint(rzl1*CLINV)*CL
            rxl2 = rxl2 - anint(rxl2*CLINV)*CL
            ryl2 = ryl2 - anint(ryl2*CLINV)*CL
            rzl2 = rzl2 - anint(rzl2*CLINV)*CL
      rl1sq = sqrt(rxl1*rxl1 + ryl1*ryl1 + rzl1*rzl1)
      rl2sq = sqrt(rxl2*rxl2 + ryl2*ryl2 + rzl2*rzl2)
        rl1in = 1.0/rl1sq
        rl2in = 1.0/rl2sq
      rl1l2 = (rxl1*rxl2 + ryl1*ryl2 + rzl1*rzl2)*(rl1in*rl2in)
         pp1x = -rl2in*((rxl2*rl2in*rl1l2) - (rxl1*rl1in))
         pp1y = -rl2in*((ryl2*rl2in*rl1l2) - (ryl1*rl1in))
         pp1z = -rl2in*((rzl2*rl2in*rl1l2) - (rzl1*rl1in))
         pp2x = -rl1in*((rxl2*rl2in) - (rxl1*rl1in*rl1l2))
         pp2y = -rl1in*((ryl2*rl2in) - (ryl1*rl1in*rl1l2))
         pp2z = -rl1in*((rzl2*rl2in) - (rzl1*rl1in*rl1l2))
         pp3x = -pp1x - pp2x
         pp3y = -pp1y - pp2y
         pp3z = -pp1z - pp2z
      xktheta = rkang
       FX(L + 2,I2) = FX(L + 2,I2) +xktheta*pp3x
       FY(L + 2,I2) = FY(L + 2,I2) +xktheta*pp3y
       FZ(L + 2,I2) = FZ(L + 2,I2) +xktheta*pp3z
       FX(L,J) = FX(L,J) +xktheta*pp2x
       FY(L,J) = FY(L,J) +xktheta*pp2y
       FZ(L,J) = FZ(L,J) +xktheta*pp2z
       FX(L,J + 1) = FX(L,J + 1) +xktheta*pp1x
       FY(L,J + 1) = FY(L,J + 1) +xktheta*pp1y
       FZ(L,J + 1) = FZ(L,J + 1) +xktheta*pp1z
      end do
        end do
! ****         NEW ACCELERATIONS  &  NEW VELOCITY
!
           WRITE(14,*) 'SUMMMU = ',SUMU
            do 40081 L = 1,NF
           if (N(L) == 0) GOTO 40081
            do 40080 I = 1, N(L)
               AX(L,I) = FX(L,I)*RMASSI(L)
               AY(L,I) = FY(L,I)*RMASSI(L)
               AZ(L,I) = FZ(L,I)*RMASSI(L)
               BUFFXN(L,I) = X0(L,I)
               BUFFYN(L,I) = Y0(L,I)
               BUFFZN(L,I) = Z0(L,I)
               X0O(L,I) = X0(L,I)
               Y0O(L,I) = Y0(L,I)
               Z0O(L,I) = Z0(L,I)
               BUFFXO(L,I) = BUFFXN(L,I) - DLTT*VX(L,I)
               BUFFYO(L,I) = BUFFYN(L,I) - DLTT*VY(L,I)
               BUFFZO(L,I) = BUFFZN(L,I) - DLTT*VZ(L,I)
40080       CONTINUE
40081     CONTINUE
!
! ****   MOLECULAR DYNAMICS SIMULATIONS FOR
!              L-J  (12:6) POTENTIAL
!
! ********************************
! ****   START MD
      KT = 1
      KTIM = 0
        VSUMSQ = 0.0
       do 55066 L = 1,NF
      do 55066 KT = 1,101
      do 55066 KZ = 1,NZ
      DEN(L,KT,KZ) = 0.0
      TEMPR(L,KT,KZ) = 0.0
55066  CONTINUE
         ITAVG = 0
       KCOUNT = 0
       KKPPI = 0
      do 77777 ILONG = 1, MXLONG
      do 77776 ITSTEP = 1,MXSHRT
       KCOUNT = KCOUNT + 1
! ******************************************
! **   CHECK ON POSSIBLE LIST UPDATE
       NEWLIST = 0
       DMAX = 0.0
       do L = 1,NF
       do I = 1,N(L)
       DDX = X0(L,I) - X0O(L,I)
       DDY = Y0(L,I) - Y0O(L,I)
       DDZ = Z0(L,I) - Z0O(L,I)
               if (DDX > CLP) DDX = DDX - CL
               if (DDX <= CLN) DDX = DDX + CL
               if (DDY > CLP) DDY = DDY - CL
               if (DDY <= CLN) DDY = DDY + CL
               if (DDZ > CLP) DDZ = DDZ - CL
               if (DDZ <= CLN) DDZ = DDZ + CL
       DDSQ = DDX*DDX + DDY*DDY + DDZ*DDZ
       DMAX = MAX(DDSQ,DMAX)
7708   CONTINUE
       end do
       end do
       if (DMAX > RLCMIN) then
       KKPPI = KKPPI + 1
       NEWLIST = 1
       end if
! ***********************************
       if ((ILONG == 1).AND.(ITSTEP == 1)) GOTO 20111
       if ((ITSTEP == MXSHRT).OR.(NEWLIST == 1)) GOTO 20111
       GOTO 20222
20111  CONTINUE
          NT = 0
          do L = 1,NF
          NT = NT + N(L)
          end do
          NTT = 0
          do M = 1,NS
          NTT = NTT + NO(M)
          end do
          NTT = NT + NTT
          NLIST = 0
         do 15000 I = 1,NT
             NABOR(I) = NLIST + 1
             NTE = 0
             do L = 1,NF
             NTEI = NTE
             NTE = NTE + N(L)
             if (I <= NTE) then
                I1 = I - NTEI
                 X = X0(L,I1)
                 Y = Y0(L,I1)
                 Z = Z0(L,I1)
                 GOTO 88111
             end if
             end do
88111 CONTINUE
          do 14999 J = I + 1 , NTT
            if (J <= NT) then
               NTE = 0
               do L = 1,NF
             NTEI = NTE
             NTE = NTE + N(L)
             if (J <= NTE) then
                J1 = J - NTEI
                 XX = X0(L,J1)
                 YY = Y0(L,J1)
                 ZZ = Z0(L,J1)
               GOTO 88222
             end if
             end do
             else
              J1 = J - NT
               NTE = 0
               do M = 1,NS
             NTEI = NTE
             NTE = NTE + NO(M)
             if (J1 <= NTE) then
                JS1 = J1 - NTEI
                 XX = XS(M,JS1)
                 YY = YS(M,JS1)
                 ZZ = ZS(M,JS1)
               GOTO 88222
             end if
             end do
             end if
88222 CONTINUE
!
            if (J <= NT) then
                  XIJ = X - XX
                  YIJ = Y - YY
                  ZIJ = Z - ZZ
               if (XIJ > CLP) XIJ = XIJ - CL
               if (XIJ <= CLN) XIJ = XIJ + CL
               if (YIJ > CLP) YIJ = YIJ - CL
               if (YIJ <= CLN) YIJ = YIJ + CL
               if (ZIJ > CLP) ZIJ = ZIJ - CL
               if (ZIJ <= CLN) ZIJ = ZIJ + CL
               SQR = XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ
            else
                  RZ = Z - ZZ
                  RZ = RZ - anint(RZ*CLINV)*CL
                  RX = X - XX
                  RY = Y - YY
                  RX = RX - anint(RX*CLINV)*CL
                  RY = RY - anint(RY*CLINV)*CL
                  SQR = RX*RX + RY*RY + RZ*RZ
          end if
!
           NTE1 = 0
           do L = 1,NF
           NTE1 = NTE1 + N(L)
           NTE2 = 0
           do K = 1,NFNS
           if (K <= NF) NTE2 = NTE2 + N(K)
           if (K > NF) NTE2 = NTE2 + NO(K - NF)
           if (I <= NTE1.AND.J <= NTE2) then
           RLISTC = RLISTCUT(L,K)
           GOTO 88999
            end if
         end do
      end do
88999 CONTINUE
!
                if (SQR > RLISTC) GOTO 14999
                      NLIST = NLIST + 1
                      LIST(NLIST) = J
!
14999 CONTINUE
15000 CONTINUE
                NABOR(NT + 1) = NLIST + 1
       do L = 1,NF
       do I = 1,N(L)
       X0O(L,I) = X0(L,I)
       Y0O(L,I) = Y0(L,I)
       Z0O(L,I) = Z0(L,I)
       end do
       end do
20222  CONTINUE
!
! *******************************************************
      SUMU = 0.0
      KTIM = KTIM + 1
      if (KTIM == LTIM) then
      KTIM = 0
      KT = KT + 1
      end if
! ****         CALCULATE NEW VELOCITIES
          ITAVG = ITAVG +1
          do 12000 L = 1,NF
           do 12050 I = 1, N(L)
!
! ****            NEW POSITIONS AT ( t + Delta (t) )
!
        COXO(L,I) = BUFFXN(L,I)
        COYO(L,I) = BUFFYN(L,I)
        COZO(L,I) = BUFFZN(L,I)
       XNEW(L,I) = 2.0*BUFFXN(L,I) - BUFFXO(L,I) + (DLTSQ * AX(L,I))
       YNEW(L,I) = 2.0*BUFFYN(L,I) - BUFFYO(L,I) + (DLTSQ * AY(L,I))
       ZNEW(L,I) = 2.0*BUFFZN(L,I) - BUFFZO(L,I) + (DLTSQ * AZ(L,I))
           VX(L,I) = DLT2I*(XNEW(L,I) - BUFFXO(L,I))
           VY(L,I) = DLT2I*(YNEW(L,I) - BUFFYO(L,I))
           VZ(L,I) = DLT2I*(ZNEW(L,I) - BUFFZO(L,I))
          if (ILONG == IEQ) then
               BUFFX(L,I) = 0.0
               BUFFY(L,I) = 0.0
               BUFFZ(L,I) = 0.0
           end if
           if (ILONG >= IEQ) then
               BUFFXT(L,I) = BUFFX(L,I) + XNEW(L,I) - COXO(L,I)
               BUFFYT(L,I) = BUFFY(L,I) + YNEW(L,I) - COYO(L,I)
               BUFFZT(L,I) = BUFFZ(L,I) + ZNEW(L,I) - COZO(L,I)
             end if
               BUFFXO(L,I) = BUFFXN(L,I)
               BUFFYO(L,I) = BUFFYN(L,I)
               BUFFZO(L,I) = BUFFZN(L,I)
               BUFFXN(L,I) = XNEW(L,I)
               BUFFYN(L,I) = YNEW(L,I)
               BUFFZN(L,I) = ZNEW(L,I)
               X0(L,I) = XNEW(L,I)
               Y0(L,I) = YNEW(L,I)
               Z0(L,I) = ZNEW(L,I)
          if (X0(L,I) > CLP) then
          KXO = (X0(L,I) + CLP)*CLINV
          X0(L,I) = X0(L,I) - FLOAT(KXO)*CL
          end if
          if (X0(L,I) <= CLN) then
          KXO = (X0(L,I) + CLN)*CLINV
          X0(L,I) = X0(L,I) - FLOAT(KXO)*CL
          end if
          if (Y0(L,I) > CLP) then
          KYO = (Y0(L,I) + CLP)*CLINV
          Y0(L,I) = Y0(L,I) - FLOAT(KYO)*CL
          end if
          if (Y0(L,I) <= CLN) then
          KYO = (Y0(L,I) + CLN)*CLINV
          Y0(L,I) = Y0(L,I) - FLOAT(KYO)*CL
          end if
          if (Z0(L,I) > CLP) then
          KZO = (Z0(L,I) + CLP)*CLINV
          Z0(L,I) = Z0(L,I) - FLOAT(KZO)*CL
          end if
          if (Z0(L,I) <= CLN) then
          KZO = (Z0(L,I) + CLN)*CLINV
          Z0(L,I) = Z0(L,I) - FLOAT(KZO)*CL
          end if
! ******
12050       CONTINUE
12000       CONTINUE
! *******************************************************************
! **  THERMAL DIFFUSE SCATTERING FROM THE OXYGEN ATOMS
             if (IDIFF > 0) then
       NNTE(1) = 0
       do JP = 1,NF - 1
       NNTE(JP + 1) = NNTE(JP) + N(JP)
       end do
       NTE = 0
      do 88192 L = 1,NF
        ILOW = 1 + NTE
        IHIGH = N(L) + NTE
        NTE = NTE + N(L)
       if (N(L) == 0) GOTO 88192
      do 88193 J = ILOW,IHIGH
                JBEG = NABOR(J)
                JEND = NABOR(J + 1) - 1
      if (JBEG > JEND) GOTO 88193
       J1 = J - NNTE(L)
          XXNEW = BUFFXN(L,J1)
          YYNEW = BUFFYN(L,J1)
          ZZNEW = BUFFZN(L,J1)
          XI0 = BUFFXO(L,J1)
          YI0 = BUFFYO(L,J1)
          ZI0 = BUFFZO(L,J1)
! **
        TCOL = 0.0
        FFDLT = DLTT
        IICOLD = 1000000
      do 84910 KCOLL = 1,1000
      NSIL = 0
        FFDLT = FFDLT - TCOL
      do 84010 JNAB = JBEG,JEND
        FDELT0 = 0.0
                 II = LIST(JNAB)
       if (II <= NT) GOTO 84010
                 II1 = II - NT
                 JJ = II1
                 do M = 1,NS
                 JJ = JJ - NO(M)
                 if (JJ <= 0) then
                 KK1 = JJ + NO(M)
                 ZSNEW = ZS(M,KK1)
                 XSNEW = XS(M,KK1)
                 YSNEW = YS(M,KK1)
                 NSCOMP = M
                 GOTO 84911
                 end if
                 end do
84911  CONTINUE
       if (II1 == IICOLD) GOTO 84010
      RZ = ZZNEW - ZSNEW
      RZ = RZ - anint(RZ*ZLINV)*ZL
      RX = XXNEW - XSNEW
      RY = YYNEW - YSNEW
      RX = RX - anint(RX*CLINV)*CL
      RY = RY - anint(RY*CLINV)*CL
      SQ = RX*RX + RY*RY + RZ*RZ
                RO = SQ
       if (RO > RWMAX(L,NSCOMP)) GOTO 84010
              RZN = ZI0 - ZSNEW
              RZN = RZN - anint(RZN*CLINV)*CL
              RXN = XI0 - XSNEW
              RYN = YI0 - YSNEW
              RXN = RXN - anint(RXN*CLINV)*CL
              RYN = RYN - anint(RYN*CLINV)*CL
              SQN = RXN*RXN + RYN*RYN + RZN*RZN
                RNN = SQN
           RWAA = RWMIN(L,NSCOMP)
        if (RO < RWAA.AND.RNN > RWAA) GOTO 48765
        if (RO > RWAA.AND.RNN < RWAA) GOTO 48765
                GOTO 84010
48765  CONTINUE
! ****** FIND EXACT COLLISION TIME
        NSIL = NSIL + 1
        ISIL(NSIL) = II1
        RWM = RWMIN(L,NSCOMP)
! **** INITIAL ESTIMATE OF FDELT
        FMM = 1.0
        if (RNN > RWM) FMM = -1.0
        VDOT = VX(L,J1)*VX(L,J1) + VY(L,J1)*VY(L,J1) + VZ(L,J1)*VZ(L,J1)
        RVDOT = RXN*VX(L,J1) + RYN*VY(L,J1) + RZN*VZ(L,J1)
        RELR = RNN - RWM
!        FDELT0 = (-RVDOT+FMM*sqrt((RVDOT*RVDOT)-VDOT*RELR))/VDOT
        FDELT = FDELT0
! **** NEWTON METHOD
        ADOT = (AX(L,J1)*AX(L,J1) + AY(L,J1)*AY(L,J1) + AZ(L,J1)*AZ(L,J1))
        ADOT4 = 0.25*ADOT
        RVDOT2 = 2.0*RVDOT
        VADOT = VX(L,J1)*AX(L,J1) + VY(L,J1)*AY(L,J1) + VZ(L,J1)*AZ(L,J1)
        VADOT3 = VADOT*3.0
        RADOT = RXN*AX(L,J1) + RYN*AY(L,J1) + RZN*AZ(L,J1) + VDOT
        RADOT2 = RADOT*2.0
        do 78399 LP = 1,NITER
        FF1 = FDELT
        FF2 = FF1*FF1
        FF3 = FF1*FF2
        FF4 = FF2*FF2
        FNEWT = ADOT4*FF4 + VADOT*FF3 + RADOT*FF2 + RVDOT2*FF1 + RELR
        DFNEWT = ADOT*FF3 + VADOT3*FF2 + RADOT2*FF1 + RVDOT2
        FDELT = FDELT - (FNEWT/DFNEWT)
        FMULT = DFNEWT*RELR
        if (FDELT > FFDLT.OR.FDELT < 0.0.OR.FMULT > 0.0) then
        FDELT0 = FDELT0 + FIT*FFDLT
        FDELT = FDELT0
        GOTO 78399
        end if
        ERR = DABS(FDELT - FF1)
        if (ERR < TOL) then
        FDELT = FDELT + ERR
        if (FDELT > FFDLT) FDELT = FFDLT
               GOTO 78398
        end if
78399  CONTINUE
         FDELT = FFDLT
78398  CONTINUE
        if (FDELT < 0.0.OR.FDELT > FFDLT) then
        WRITE(14,*) FDELT,FFDLT,FMULT,TCOL,KCOLL
         WRITE(14,*) L,J1,II1
         WRITE(14,*) XI0,YI0,ZI0
         WRITE(14,*) XXNEW,YYNEW,ZZNEW
         WRITE(14,*) XSNEW,YSNEW,ZSNEW
         WRITE(14,*) VX(L,J1),VY(L,J1),VZ(L,J1)
         WRITE(14,*) RWMIN(L,NSCOMP),RO,RNN
        WRITE(14,*) 'ROOT FAILURE'
        STOP
        end if
84710  CONTINUE
        COLLT(NSIL) = FDELT
84010  CONTINUE
         if (NSIL == 0) GOTO 88193
! ************** FIND SHORTEST COLLISION TIME
         TCOL = 1000.0
         do 79911 JT = 1,NSIL
         if (COLLT(JT) < TCOL) then
         TCOL = COLLT(JT)
         ICOL = JT
         end if
79911  CONTINUE
! ************************************
        FDLTS2 = 0.5*TCOL*TCOL
        XI0 = XI0 + TCOL*VX(L,J1) + FDLTS2*AX(L,J1)
        YI0 = YI0 + TCOL*VY(L,J1) + FDLTS2*AY(L,J1)
        ZI0 = ZI0 + TCOL*VZ(L,J1) + FDLTS2*AZ(L,J1)
        II1 = ISIL(ICOL)
        IICOLD = II1
                 JJ = II1
                 do M = 1,NS
                 JJ = JJ - NO(M)
                 if (JJ <= 0) then
                 KK1 = JJ + NO(M)
                 ZSNEW = ZS(M,KK1)
                 XSNEW = XS(M,KK1)
                 YSNEW = YS(M,KK1)
                 NSCOMP = M
                 GOTO 84921
                 end if
                 end do
84921  CONTINUE
        XIS = XI0 - XSNEW
        YIS = YI0 - YSNEW
        ZIS = ZI0 - ZSNEW
!
              XIS = XIS - anint(XIS*CLINV)*CL
              YIS = YIS - anint(YIS*CLINV)*CL
              ZIS = ZIS - anint(ZIS*CLINV)*CL
! ****************************
            ISV = mod(RMUL2 * ISV,RMO2)
               RNV = ISV * RIMO1
            do LV = 1,5000
            if (RNV >= TFUNC(LV)) then
            LVV = LV
            GOTO 88898
            end if
            end do
88898  CONTINUE
            LVV1 = LVV - 1
            DDE = TFUNC(LVV) - TFUNC(LVV1)
            DNE = VNEW(LVV) - VNEW(LVV1)
            DMUL = (RNV - TFUNC(LVV1))
            VEL = VNEW(LVV1) + (DNE/DDE)*DMUL
! ******************************
          VVE = VEL*FACTT(L)
!          WRITE(14,*) VEL,VVE
            ISX = mod(RMUL2 * ISX,RMO2)
               RNX = ISX * RIMO1
         FFAR1 = sqrt(1.0D0 - RNX)
         FFAR2 = sqrt(RNX)
          VVVE = VVE*FFAR1
            ISZ = mod(RMUL2 * ISZ,RMO2)
               RA  = ISZ * RIMOA
         VTHET = VVVE*sin(RA)
            ISY = mod(RMUL2 * ISY,RMO2)
               RNY = ISY * RIMO1
         FVV = 1.0
         if (RNY < 0.5) FVV = -1.0
         VR = FVV*VVE*FFAR2
         VPHI = VVVE*cos(RA)
        RRIS = XIS*XIS + YIS*YIS
        RIS = sqrt(ZIS*ZIS + RRIS)
         RRIS = sqrt(RRIS)
         COST = ZIS/RIS
         THET = acos(COST)
         COSP = DABS(YIS/RRIS)
         PHI = acos(COSP)
        if (YIS > 0.0.AND.XIS > 0.0) PHI = PHI
        if (YIS <= 0.0.AND.XIS > 0.0) PHI = PI - PHI
        if (YIS <= 0.0.AND.XIS <= 0.0) PHI = PHI + PI
        if (YIS > 0.0.AND.XIS <= 0.0) PHI = PI2 - PHI
        COSP = cos(PHI)
        SINT = sin(THET)
        SINP = sin(PHI)
        VZ(L,J1) = VR*COST - VTHET*SINT
        AVEL = VR*SINT + VTHET*COST
        VX(L,J1) = AVEL*SINP + VPHI*COSP
        VY(L,J1) = AVEL*COSP - VPHI*SINP
! ***
        FDDT = (FFDLT - TCOL)
        FDDTSQ = 0.5*FDDT*FDDT
        XXNEW = XI0 + FDDT*VX(L,J1) + FDDTSQ*AX(L,J1)
        YYNEW = YI0 + FDDT*VY(L,J1) + FDDTSQ*AY(L,J1)
        ZZNEW = ZI0 + FDDT*VZ(L,J1) + FDDTSQ*AZ(L,J1)
        XNEW(L,J1) = XXNEW
        YNEW(L,J1) = YYNEW
        ZNEW(L,J1) = ZZNEW
        DMMT = DLTT - FDDT
        DMMTSQ = 0.5*DMMT*DMMT
        BUFFXO(L,J1) = XI0 - DMMT*VX(L,J1) + DMMTSQ*AX(L,J1)
        BUFFYO(L,J1) = YI0 - DMMT*VY(L,J1) + DMMTSQ*AY(L,J1)
        BUFFZO(L,J1) = ZI0 - DMMT*VZ(L,J1) + DMMTSQ*AZ(L,J1)
               BUFFXN(L,J1) = XNEW(L,J1)
               BUFFYN(L,J1) = YNEW(L,J1)
               BUFFZN(L,J1) = ZNEW(L,J1)
               X0(L,J1) = XNEW(L,J1)
               Y0(L,J1) = YNEW(L,J1)
               Z0(L,J1) = ZNEW(L,J1)
           if (ILONG >= IEQ) then
               BUFFXT(L,J1) = BUFFX(L,J1) + XNEW(L,J1) - COXO(L,J1)
               BUFFYT(L,J1) = BUFFY(L,J1) + YNEW(L,J1) - COYO(L,J1)
               BUFFZT(L,J1) = BUFFZ(L,J1) + ZNEW(L,J1) - COZO(L,J1)
          end if
          if (X0(L,J1) > CLP) then
          KXO = (X0(L,J1) + CLP)*CLINV
          X0(L,J1) = X0(L,J1) - FLOAT(KXO)*CL
          end if
          if (X0(L,J1) <= CLN) then
          KXO = (X0(L,J1) + CLN)*CLINV
          X0(L,J1) = X0(L,J1) - FLOAT(KXO)*CL
          end if
          if (Y0(L,J1) > CLP) then
          KYO = (Y0(L,J1) + CLP)*CLINV
          Y0(L,J1) = Y0(L,J1) - FLOAT(KYO)*CL
          end if
          if (Y0(L,J1) <= CLN) then
          KYO = (Y0(L,J1) + CLN)*CLINV
          Y0(L,J1) = Y0(L,J1) - FLOAT(KYO)*CL
          end if
          if (Z0(L,J1) > CLP) then
          KZO = (Z0(L,J1) + CLP)*CLINV
          Z0(L,J1) = Z0(L,J1) - FLOAT(KZO)*CL
          end if
          if (Z0(L,J1) <= CLN) then
          KZO = (Z0(L,J1) + CLN)*CLINV
          Z0(L,J1) = Z0(L,J1) - FLOAT(KZO)*CL
          end if
! *********
84910 CONTINUE
88193 CONTINUE
88192 CONTINUE
! ************************************************
         end if
! *******************************************************
            if (ITAVG >= NPRINT) then
                ITAVG = 0
            end if
78992  CONTINUE
         do L = 1,NF
           do I = 1,N(L)
!
               FX(L,I) = 0.0
               FY(L,I) = 0.0
               FZ(L,I) = 0.0
               UOLD(L,I) = 0.0
             if (ILONG >= IEQ) then
             BUFFX(L,I) = BUFFXT(L,I)
             BUFFY(L,I) = BUFFYT(L,I)
             BUFFZ(L,I) = BUFFZT(L,I)
             end if
!
          end do
        end do
! ****  EVALUATE FORCES DUE TO OXYGEN ATOMS
! **  EVALUATE INTRAMOLECULAR STRETCHING FORCES
      do L = 1,2
      I2 = 0
      do J = 1,N(L),2
      I2 = I2 + 1
! ** RELATIVE SEPARATION WITHIN TRIATOMIC STRUCTURE
      XTEMP3 = X0(L + 2,I2)
      DX1 = X0(L,J) - XTEMP3
      if (DX1 > CLP) then
      XTEMP3 = XTEMP3 + CL
      DX1 = X0(L,J) - XTEMP3
      end if
      if (DX1 < CLN) then
      XTEMP3 = XTEMP3 - CL
      DX1 = X0(L,J) - XTEMP3
      end if
      XTEMP3 = X0(L + 2,I2)
      DX2 = X0(L,J + 1) - XTEMP3
      if (DX2 > CLP) then
      XTEMP3 = XTEMP3 + CL
      DX2 = X0(L,J + 1) - XTEMP3
      end if
      if (DX2 < CLN) then
      XTEMP3 = XTEMP3 - CL
      DX2 = X0(L,J + 1) - XTEMP3
      end if
! ***
      YTEMP3 = Y0(L + 2,I2)
      DY1 = Y0(L,J) - YTEMP3
      if (DY1 > CLP) then
      YTEMP3 = YTEMP3 + CL
      DY1 = Y0(L,J) - YTEMP3
      end if
      if (DY1 < CLN) then
      YTEMP3 = YTEMP3 - CL
      DY1 = Y0(L,J) - YTEMP3
      end if
      YTEMP3 = Y0(L + 2,I2)
      DY2 = Y0(L,J + 1) - YTEMP3
      if (DY2 > CLP) then
      YTEMP3 = YTEMP3 + CL
      DY2 = Y0(L,J + 1) - YTEMP3
      end if
      if (DY2 < CLN) then
      YTEMP3 = YTEMP3 - CL
      DY2 = Y0(L,J + 1) - YTEMP3
      end if
! ***
      ZTEMP3 = Z0(L + 2,I2)
      DZ1 = Z0(L,J) - ZTEMP3
      if (DZ1 > CLP) then
      ZTEMP3 = ZTEMP3 + CL
      DZ1 = Z0(L,J) - ZTEMP3
      end if
      if (DZ1 < CLN) then
      ZTEMP3 = ZTEMP3 - CL
      DZ1 = Z0(L,J) - ZTEMP3
      end if
      ZTEMP3 = Z0(L + 2,I2)
      DZ2 = Z0(L,J + 1) - ZTEMP3
      if (DZ2 > CLP) then
      ZTEMP3 = ZTEMP3 + CL
      DZ2 = Z0(L,J + 1) - ZTEMP3
      end if
      if (DZ2 < CLN) then
      ZTEMP3 = ZTEMP3 - CL
      DZ2 = Z0(L,J + 1) - ZTEMP3
      end if
! ***
      R21 = sqrt(DX1*DX1 + DY1*DY1 + DZ1*DZ1)
      FORCS = RK0*(1.D0 - (B0CO/R21))
      FIJX1 = -DX1*FORCS
      FIJY1 = -DY1*FORCS
      FIJZ1 = -DZ1*FORCS
            FX(L,J) = FX(L,J) + FIJX1
            FY(L,J) = FY(L,J) + FIJY1
            FZ(L,J) = FZ(L,J) + FIJZ1
            FX(L + 2,I2) = FX(L + 2,I2) - FIJX1
            FY(L + 2,I2) = FY(L + 2,I2) - FIJY1
            FZ(L + 2,I2) = FZ(L + 2,I2) - FIJZ1
      R21 = sqrt(DX2*DX2 + DY2*DY2 + DZ2*DZ2)
      FORCS = RK0*(1.D0 - (B0CO/R21))
      FIJX2 = -DX2*FORCS
      FIJY2 = -DY2*FORCS
      FIJZ2 = -DZ2*FORCS
            FX(L,J + 1) = FX(L,J + 1) + FIJX2
            FY(L,J + 1) = FY(L,J + 1) + FIJY2
            FZ(L,J + 1) = FZ(L,J + 1) + FIJZ2
            FX(L + 2,I2) = FX(L + 2,I2) - FIJX2
            FY(L + 2,I2) = FY(L + 2,I2) - FIJY2
            FZ(L + 2,I2) = FZ(L + 2,I2) - FIJZ2
       end do
       end do
! ******************************
! **** EVALUATE O=C=O BOND ANGLE BENDING
      do L = 1,2
      I2 = 0
      do J = 1,N(L),2
      I2 = I2 + 1
      rx1 = X0(L + 2,I2)
      rx2 = X0(L,J)
      rx3 = X0(L,J + 1)
      ry1 = Y0(L + 2,I2)
      ry2 = Y0(L,J)
      ry3 = Y0(L,J + 1)
      rz1 = Z0(L + 2,I2)
      rz2 = Z0(L,J)
      rz3 = Z0(L,J + 1)
      rxl1 = rx1 - rx2
      ryl1 = ry1 - ry2
      rzl1 = rz1 - rz2
      rxl2 = rx3 - rx1
      ryl2 = ry3 - ry1
      rzl2 = rz3 - rz1
            rxl1 = rxl1 - anint(rxl1*CLINV)*CL
            ryl1 = ryl1 - anint(ryl1*CLINV)*CL
            rzl1 = rzl1 - anint(rzl1*CLINV)*CL
            rxl2 = rxl2 - anint(rxl2*CLINV)*CL
            ryl2 = ryl2 - anint(ryl2*CLINV)*CL
            rzl2 = rzl2 - anint(rzl2*CLINV)*CL
      rl1sq = sqrt(rxl1*rxl1 + ryl1*ryl1 + rzl1*rzl1)
      rl2sq = sqrt(rxl2*rxl2 + ryl2*ryl2 + rzl2*rzl2)
        rl1in = 1.0/rl1sq
        rl2in = 1.0/rl2sq
      rl1l2 = (rxl1*rxl2 + ryl1*ryl2 + rzl1*rzl2)*(rl1in*rl2in)
         pp1x = -rl2in*((rxl2*rl2in*rl1l2) - (rxl1*rl1in))
         pp1y = -rl2in*((ryl2*rl2in*rl1l2) - (ryl1*rl1in))
         pp1z = -rl2in*((rzl2*rl2in*rl1l2) - (rzl1*rl1in))
         pp2x = -rl1in*((rxl2*rl2in) - (rxl1*rl1in*rl1l2))
         pp2y = -rl1in*((ryl2*rl2in) - (ryl1*rl1in*rl1l2))
         pp2z = -rl1in*((rzl2*rl2in) - (rzl1*rl1in*rl1l2))
         pp3x = -pp1x - pp2x
         pp3y = -pp1y - pp2y
         pp3z = -pp1z - pp2z
      xktheta = rkang
       FX(L + 2,I2) = FX(L + 2,I2) +xktheta*pp3x
       FY(L + 2,I2) = FY(L + 2,I2) +xktheta*pp3y
       FZ(L + 2,I2) = FZ(L + 2,I2) +xktheta*pp3z
       FX(L,J) = FX(L,J) +xktheta*pp2x
       FY(L,J) = FY(L,J) +xktheta*pp2y
       FZ(L,J) = FZ(L,J) +xktheta*pp2z
       FX(L,J + 1) = FX(L,J + 1) +xktheta*pp1x
       FY(L,J + 1) = FY(L,J + 1) +xktheta*pp1y
       FZ(L,J + 1) = FZ(L,J + 1) +xktheta*pp1z
      end do
        end do
! ********************
       NNTE(1) = 0
       do JP = 1,NF - 1
       NNTE(JP + 1) = NNTE(JP) + N(JP)
       end do
        NTE = 0
       RDFSUM(1) = RDFSUM(1) + FLOAT(IDOUB)
       RDFSUM(2) = RDFSUM(2) + FLOAT(ISINGO)
      do 59192 L = 1,NF
        ILOW = 1 + NTE
        IHIGH = N(L) + NTE
        NTE = NTE + N(L)
       if (N(L) == 0) GOTO 59192
      do 59193 J = ILOW,IHIGH
                JBEG = NABOR(J)
                JEND = NABOR(J + 1) - 1
      if (JBEG > JEND) GOTO 59193
         J1 = J - NNTE(L)
      do 34010 JNAB = JBEG,JEND
                 II = LIST(JNAB)
       if (II <= NT) GOTO 34010
                 II1 = II - NT
                 JJ = II1
                 do M = 1,NS
                 JJ = JJ - NO(M)
                 if (JJ <= 0) then
                 KK1 = JJ + NO(M)
                 ZSNEW = ZS(M,KK1)
                 XSNEW = XS(M,KK1)
                 YSNEW = YS(M,KK1)
                 NM = M
                 GOTO 84931
                 end if
                 end do
84931  CONTINUE
      RZ = Z0(L,J1) - ZSNEW
      RZ = RZ - anint(RZ*CLINV)*CL
      RX = X0(L,J1) - XSNEW
      RY = Y0(L,J1) - YSNEW
      RX = RX - anint(RX*CLINV)*CL
      RY = RY - anint(RY*CLINV)*CL
      SQ = RX*RX + RY*RY + RZ*RZ
! *********
!
                 NMM = NM + NF
                     if (SQ > SQCUT(L,NMM)) GOTO 34010
                        RR1 = sqrt(SQ)
                        RR2 = SQSIG(L,NMM) / SQ
                        RR6 = RR2 * RR2 * RR2
              POTEN = 4.  *EPSS(L,NMM)* RR6*(RR6 - 1.)
              POTEN = RTI*(POTEN+ (CCUT1(L,NMM)*RR1) + CCUT2(L,NMM))
              UOLD(L,J1) = UOLD(L,J1) + POTEN
              FORCE = FORFT(L,NMM)*  RR6 * (RR6 - 0.5) * RR2
              FORCE = FORCE - (CCUT1(L,NMM)/RR1)
                     SUMU  = SUMU + POTEN
                  FIJX  = FORCE * RX
                  FIJY  = FORCE * RY
                  FIJZ  = FORCE * RZ
                     FX(L,J1) = FX(L,J1) + FIJX
                     FY(L,J1) = FY(L,J1) + FIJY
                     FZ(L,J1) = FZ(L,J1) + FIJZ
              if (SQ > RDFMAX) GOTO 34010
              INTRDF = (RR1*DLRDFI) + 1.0
              DRDF(L,M,INTRDF) = DRDF(L,M,INTRDF) + 1.0
34010 CONTINUE
59193 CONTINUE
59192 CONTINUE
! ****         NEW ACCELERATIONS  &  NEW VELOCITY
!
        TOTKE = 0.0
            do 40001 L = 1,NF
            if (N(L) == 0) GOTO 40001
            do 40000 I = 1, N(L)
               AX(L,I) = FX(L,I)*RMASSI(L)
               AY(L,I) = FY(L,I)*RMASSI(L)
               AZ(L,I) = FZ(L,I)*RMASSI(L)
       TOTKE = TOTKE + VVFACT(L)*(VX(L,I)*VX(L,I))
       TOTKE = TOTKE + VVFACT(L)*(VY(L,I)*VY(L,I))
       TOTKE = TOTKE + VVFACT(L)*(VZ(L,I)*VZ(L,I))
40000       CONTINUE
40001     CONTINUE
          SUMN = 0.0
          do L = 1,NF
          SUMN = SUMN + FLOAT(N(L))
          end do
         VSUMSQ = VSUMSQ + (TOTKE/SUMN)
               IPRINT = IPRINT + 1
               if (IPRINT >= NPRINT) then
               RCC = 1.0/(FLOAT(KCOUNT))
               IPRINT = 0
               TTEM = VSUMSQ*RCC
               WRITE(14,*) SUMU,TTEM
! ********************************************
!          AVCANG = 0.0
!             RAVA1 = 0.0
!             RAVBB1 = 0.0
!             RAVBB2 = 0.0
!         do L=1,2
!           I2 = 0
!          do I=1,N(L),2
!            I2 = I2+1
!             DDXA = DABS(X0(L,I)-X0(L,I+1))
!             if (DDXA > CLP) DDXA = DDXA-CL
!             DDYA = DABS(Y0(L,I)-Y0(L,I+1))
!             if (DDYA > CLP) DDYA = DDYA-CL
!             DDZA = DABS(Z0(L,I)-Z0(L,I+1))
!             if (DDZA > CLP) DDZA = DDZA-CL
!              DDA = DDXA*DDXA+DDYA*DDYA+DDZA*DDZA
!              RRA = sqrt(DDA)
!             DDXB1 = DABS(X0(L,I)-X0(L+2,I2))
!             if (DDXB1 > CLP) DDXB1 = DDXB1-CL
!             DDYB1 = DABS(Y0(L,I)-Y0(L+2,I2))
!             if (DDYB1 > CLP) DDYB1 = DDYB1-CL
!             DDZB1 = DABS(Z0(L,I)-Z0(L+2,I2))
!             if (DDZB1 > CLP) DDZB1 = DDZB1-CL
!              DDB1 = DDXB1*DDXB1+DDYB1*DDYB1+DDZB1*DDZB1
!              RRB1 = sqrt(DDB1)
!             DDXB2 = DABS(X0(L,I+1)-X0(L+2,I2))
!             if (DDXB2 > CLP) DDXB2 = DDXB2-CL
!             DDYB2 = DABS(Y0(L,I+1)-Y0(L+2,I2))
!             if (DDYB2 > CLP) DDYB2 = DDYB2-CL
!             DDZB2 = DABS(Z0(L,I+1)-Z0(L+2,I2))
!             if (DDZB2 > CLP) DDZB2 = DDZB2-CL
!              DDB2 = DDXB2*DDXB2+DDYB2*DDYB2+DDZB2*DDZB2
!              RRB2 = sqrt(DDB2)
!              COSST = acos(-0.5*(DDA-DDB1-DDB2)/(sqrt(DDB1*DDB2)))
!             RAVA1 = RAVA1+RRA
!             RAVBB1 = RAVBB1+RRB1
!             RAVBB2 = RAVBB2+RRB2
!              AVCANG = AVCANG+COSST
!            end do
!          end do
!          AVCANG = AVCANG/FLOAT(N(1))
!          RAVA1 = RAVA1/FLOAT(N(1))
!          RAVBB1 = RAVBB1/FLOAT(N(1))
!          RAVBB2 = RAVBB2/FLOAT(N(1))
!          WRITE(14,*) AVCANG,RAVA1,RAVBB1,RAVBB2
! ***********************************
               end if
!
       if (ILONG >= IEQ) then
         RMAS1 = 16.0/44.0
         RMAS2 = 12.0/44.0
         do L = 1,2
          I2 = 0
           do I = 1,N(L),2
             I2 = I2 + 1
!
       BUFFXX(L,I) = RMAS1*(BUFFX(L,I) + BUFFX(L,I + 1)) + RMAS2*BUFFX(L + 2,I2)
       BUFFYY(L,I) = RMAS1*(BUFFY(L,I) + BUFFY(L,I + 1)) + RMAS2*BUFFY(L + 2,I2)
       BUFFZZ(L,I) = RMAS1*(BUFFZ(L,I) + BUFFZ(L,I + 1)) + RMAS2*BUFFZ(L + 2,I2)
!
          end do
        end do


! ***   EVALUATE MSD
! **  Store each position at interval NSAVE0
          ISAVE0 = ISAVE0 + 1
          if (ISAVE0 < NSAVE0) GOTO 55001
          ISAVE0 = 0
          JSAVE0 = JSAVE0 + 1
          NPV0 = NPV0 + 1
          if (NPV0 > NSM) NPV0 = NPV0 - NSM
          do L = 1,2
          K = 0
          do 51000 I = 1,N(L),2
          K = K + 1
            TX0(L,K,NPV0) = BUFFXX(L,I)
              TY0(L,K,NPV0) = BUFFYY(L,I)
               TZ0(L,K,NPV0) = BUFFZZ(L,I)
51000     CONTINUE
          NNM(L) = K
          end do
! **  Instantaneous mean square displacement
          if (JSAVE0 <= NSM) GOTO 55001
          do 53000 J = 2,NSM
            NN = NPV0 - J + 1
            if (NN <= 0) NN = NN + NSM
            do L = 1,2
            do 52000 I = 1,NNM(L)
              RXS = TX0(L,I,NPV0) - TX0(L,I,NN)
              RYS = TY0(L,I,NPV0) - TY0(L,I,NN)
              RZS = TZ0(L,I,NPV0) - TZ0(L,I,NN)
              XXYY = RXS*RXS + RYS*RYS
              ZZZ = RZS*RZS
              SM1XY(L,J) = SM1XY(L,J) + XXYY
              SM1Z(L,J) = SM1Z(L,J) + ZZZ
              if (J == NSM14) then
              RKCON(1) = RKCON(1) + 1.0
              NZBIN = (sqrt(XXYY + ZZZ)/DELZC) + 1.0
              if (NZBIN > 100) GOTO 67110
              C(L,NZBIN,1) = C(L,NZBIN,1) + 1.0
67110         CONTINUE
              end if
              if (J == NSM12) then
              RKCON(2) = RKCON(2) + 1.0
              NZBIN = (sqrt(XXYY + ZZZ)/DELZC) + 1.0
              if (NZBIN > 100) GOTO 67220
              C(L,NZBIN,2) = C(L,NZBIN,2) + 1.0
67220         CONTINUE
              end if
              if (J == NSM34) then
              RKCON(3) = RKCON(3) + 1.0
              NZBIN = (sqrt(XXYY + ZZZ)/DELZC) + 1.0
              if (NZBIN > 100) GOTO 67330
              C(L,NZBIN,3) = C(L,NZBIN,3) + 1.0
67330         CONTINUE
              end if
              if (J == NSM) then
              RKCON(4) = RKCON(4) + 1.0
              NZBIN = (sqrt(XXYY + ZZZ)/DELZC) + 1.0
              if (NZBIN > 100) GOTO 67440
              C(L,NZBIN,4) = C(L,NZBIN,4) + 1.0
67440         CONTINUE
              end if
52000       CONTINUE
          end do
53000     CONTINUE
55001   CONTINUE
       end if
77776   CONTINUE
77777   CONTINUE
88888   CONTINUE
! ****************************************
       DELR = 0.28203328406
       DELV = (4.D0/3.D0)*PI*(DELR**3)
       do L = 1,100
       RL = FLOAT(L)
       VV(L) = DELV*(3.0*RL*RL - 3.0*RL + 1.0)
       end do
        do L = 1,NF
        do K = 1,NS
        RFD = -0.5*DELRDF
        do JP = 1,100
        RFD = RFD + DELRDF
        RDF = VCELL*DRDF(L,K,JP)/(RDFSUM(K)*DDELVV(JP))
        RDF = RDF/FLOAT(N(L))
        WRITE(14,*) RFD,RDF
        end do
        end do
        end do
7000    FORMAT(I5,3(3X,D15.6))
         WRITE(14,*) KKPPI,MXLONG
       do L = 1,2
       DENLL = 0.5*FLOAT(JSAVE0 - NSM)*FLOAT(N(L))*FLOAT(NCON)
       do J = 1,NSM
       RMSXY = SM1XY(L,J)/DENLL
       RMSZ = SM1Z(L,J)/DENLL
       WRITE(14,*) SMTAU(L,J),RMSXY,RMSZ
       end do
       end do
       WRITE(14,*) '***** CONC PROFILES *****'
       WRITE(14,*) 'DELZC = ',DELZC
       do L = 1,NF
       WRITE(14,*) 'SPECIES ',L,'   NSM14',NSM14
       ZZC = -0.5*DELZC
       do K = 1,100
       CCV = C(L,K,1)/(0.5*VV(K)*RKCON(1))
       ZZC = ZZC + DELZC
       WRITE(14,*) ZZC,CCV
       end do
       WRITE(14,*) 'SPECIES ',L,'   NSM12',NSM12
       ZZC = -0.5*DELZC
       do K = 1,100
       CCV = C(L,K,2)/(0.5*VV(K)*RKCON(2))
       ZZC = ZZC + DELZC
       WRITE(14,*) ZZC,CCV
       end do
       WRITE(14,*) 'SPECIES ',L,'   NSM34',NSM34
       ZZC = -0.5*DELZC
       do K = 1,100
       CCV = C(L,K,3)/(0.5*VV(K)*RKCON(3))
       ZZC = ZZC + DELZC
       WRITE(14,*) ZZC,CCV
       end do
       WRITE(14,*) 'SPECIES ',L,'   NSM',NSM
       ZZC = -0.5*DELZC
       do K = 1,100
       CCV = C(L,K,4)/(0.5*VV(K)*RKCON(4))
       ZZC = ZZC + DELZC
       WRITE(14,*) ZZC,CCV
       end do
       end do
      AVVOL = FLOAT(NCON)*40.0/TRY
      WRITE(14,*) B0,'ACCVOL = ',AVVOL
      WRITE(14,*) ISXX,ISYY,ISZZ
!
      CLOSE(UNIT = 10,STATUS = 'KEEP')
      CLOSE(UNIT = 14,STATUS = 'KEEP')
      CLOSE(UNIT = 5,STATUS = 'KEEP')
!
      STOP
   END PROGRAM


   SUBROUTINE SIOPOR
! *********************************************************************
!      QUENCHED SILICA PORE
!      PROGRAM NOTE :
!        1. TWO COMPONENTS: SI(4+),O(2-)
!        2. CUTS A SILICA PORE AND DETERMINES PROPERTIES
!        3. DISCARDS SINGLY BONDED SURFACE OXYGEN TRIPLETS
!        4. LENGTH SCALES IN L(SIO) = 0.162 nm UNITS
! **********************************************************************
      IMPLICIT REAL*8 (A - H,O - Z)
      PARAMETER(IZS = 2)
      PARAMETER(IZ1 = 7500)
      PARAMETER(IZ2 = 15000)
      DIMENSION X1(IZ1),Y1(IZ1),Z1(IZ1)
      DIMENSION X2(IZ2),Y2(IZ2),Z2(IZ2)
      DIMENSION XISO(IZ2),YISO(IZ2),ZISO(IZ2)
      DIMENSION XIDO(IZ2),YIDO(IZ2),ZIDO(IZ2)
      DIMENSION X01(IZ1),Y01(IZ1),Z01(IZ1)
      DIMENSION X02(IZ2),Y02(IZ2),Z02(IZ2)
      DIMENSION NEIGH(IZ2,10),NEIGHO(IZ2,10)
      DIMENSION DIST(IZ2,10),DISTO(IZ2,10)
      DIMENSION NEIOO(IZ2,10),NEIGOO(IZ2,10)
      DIMENSION DISO(IZ2,10),DISOO(IZ2,10)
      DIMENSION NEO(IZ1),NES(IZ2),NUM(IZ1)
      DIMENSION KOX(25000),IMI(25000)
      DIMENSION ISIL(25000),NTRIP(25000),NSILP(IZ1)
      DIMENSION RNSI(101),RNO(101),DELV(101)
      COMMON /SILIC/ XOR1(IZ1),YOR1(IZ1),ZOR1(IZ1)
      COMMON /OXYG/ XOR2(IZ2),YOR2(IZ2),ZOR2(IZ2)
      COMMON /SCO/ XS(IZS,10000),  YS(IZS,10000),  ZS(IZS,10000)
      COMMON /KPORE/ IDOUB,ISINGO,NCON,MCUT,N1,N2
      COMMON /PPORE/ CLT,SIO,RISF,PORR,PORIN,FREM,ROUTSQ,SQO
      COMMON /IRAND1/ ISXX,ISYY,ISZZ
      COMMON /NEIG/ NEIOOR(IZ2,10)
! SET VALUES OF CONSTANTS
      PORRSQ = PORR*PORR
      PORINS = PORIN*PORIN
      PI = 4.*atan(1.D0)
      PI2 = 2.0*PI
      RMUL2 = 16807.D0
      RMO1 = 2147483648.D0
      RIMO1 = 1./RMO1
      RMO2 = 2147483647.D0
      NSI = N1
      NOXX = N2
      CLX = CLT
      CLY = CLT
      CLZ = CLT
      CLMX = CLX/2.0
      CLNX = -CLX/2.0
      CLMY = CLY/2.0
      CLNY = -CLY/2.0
      CLMZ = CLZ/2.0
      CLNZ = -CLZ/2.0
!
! MONTE CARLO CUT
!
!     RNSI(1:101) = 0.0
!     RNO(1:101) = 0.0
      X1(1:NSI) = XOR1(1:NSI)
      Y1(1:NSI) = YOR1(1:NSI)
      Z1(1:NSI) = ZOR1(1:NSI)
      X2(1:NOXX) = XOR2(1:NOXX)
      Y2(1:NOXX) = YOR2(1:NOXX)
      Z2(1:NOXX) = ZOR2(1:NOXX)
!
      ISXX = mod(RMUL2*ISXX,RMO2)
      ISXX = mod(RMUL2*ISXX,RMO2)
      RN1 = ISXX*RIMO1
      ISYY = mod(RMUL2*ISYY,RMO2)
      ISYY = mod(RMUL2*ISYY,RMO2)
      RN2 = ISYY*RIMO1
!
      do J = 1,NSI
         X1(J) = X1(J) + (2.0*RN1 - 1.0)*CLX
         Y1(J) = Y1(J) + (2.0*RN2 - 1.0)*CLY
      end do
      do J = 1,NOXX
         X2(J) = X2(J) + (2.0*RN1 - 1.0)*CLX
         Y2(J) = Y2(J) + (2.0*RN2 - 1.0)*CLY
      end do
      do J = 1,NSI
         if (X1(J) >  CLMX) X1(J) = X1(J) - CLX
         if (X1(J) <= CLNX) X1(J) = X1(J) + CLX
         if (Y1(J) >  CLMY) Y1(J) = Y1(J) - CLY
         if (Y1(J) <= CLNY) Y1(J) = Y1(J) + CLY
      end do
      do J = 1,NOXX
         if (X2(J) >  CLMX) X2(J) = X2(J) - CLX
         if (X2(J) <= CLNX) X2(J) = X2(J) + CLX
         if (Y2(J) >  CLMY) Y2(J) = Y2(J) - CLY
         if (Y2(J) <= CLNY) Y2(J) = Y2(J) + CLY
      end do
      NS = 0
! ***********************************
      NSILON = 0
      do I = 1,NSI
         NSILON = NSILON + 1
         NSILP(NSILON) = I
      end do
      NREM = FREM*FLOAT(NSILON) + 1.0
      if (NREM == 0) GOTO 3337
      do NINT = 1,10000
         ISZZ = mod(RMUL2*ISZZ,RMO2)
         ISZZ = mod(RMUL2*ISZZ,RMO2)
         RN = ISZZ*RIMO1
         MOLS = RN*FLOAT(NSILON) + 1.0
         NS = NS + 1
         NUM(NS) = NSILP(MOLS)
         NSILP(MOLS) = NSILP(NSILON)
         NSILON = NSILON - 1
         NREM = NREM - 1
         if (NREM <= 0) EXIT
      end do
3337  CONTINUE
!
      do 556 I = 1,NS
         X01(I) = X1(NUM(I))
         Y01(I) = Y1(NUM(I))
         Z01(I) = Z1(NUM(I))
556   CONTINUE
      M = 0
      do 558 I = 1,NS
      do 557 K = 1,4
         M = M + 1
         KOX(M) = NEIOOR(NUM(I),K)
557   CONTINUE
558   CONTINUE
      NN = M
      do 559 IP = 1,20000
         IID = KOX(IP)
         IC = 0
         do J = IP + 1,NN
            if (KOX(J) == IID) then
               IC = IC + 1
               IMI(IC) = J
            end if
         end do
         if (IC == 0) GOTO 559
         do IL = 1,IC
            KOX(IMI(IL)) = KOX(NN)
            NN = NN - 1
         end do
         if (NN == IP) EXIT
559   CONTINUE
      NOX = NN
      do 444 I = 1,NOX
         X02(I) = X2(KOX(I))
         Y02(I) = Y2(KOX(I))
         Z02(I) = Z2(KOX(I))
444   CONTINUE
!
! OXYGEN NEIGHBORS TO A GIVEN SILICON ATOM
      do 9909 J = 1,NS
      NEO(J) = 0
      do 909 I = 1,NOX
         DELX = X01(J) - X02(I)
         if (DELX > CLMX) DELX = DELX - CLX
         if (DELX <= CLNX) DELX = DELX + CLX
         if (DABS(DELX) > SIO) GOTO 909
         DELY = Y01(J) - Y02(I)
         if (DELY > CLMY) DELY = DELY - CLY
         if (DELY <= CLNY) DELY = DELY + CLY
         if (DABS(DELY) > SIO) GOTO 909
         SQXY = DELY*DELY + DELX*DELX
         if (SQXY > SQO) GOTO 909
         DELZ = Z01(J) - Z02(I)
         if (DELZ > CLMZ) DELZ = DELZ - CLZ
         if (DELZ <= CLNZ) DELZ = DELZ + CLZ
         if (DABS(DELZ) > SIO) GOTO 909
         SQXYZ = DELZ*DELZ + SQXY
         if (SQXYZ > SQO) GOTO 909
         NEO(J) = NEO(J) + 1
         NEIGH(J,NEO(J)) = I
         DIST(J,NEO(J)) = SQXYZ
909   CONTINUE
! ORDER THE DISTANCES
      NN = NEO(J)
      do IP = 1,1000
         DD = DIST(J,1)
         IMIN = 1
         do M = 2,NN
            if (DIST(J,M) < DD) then
               DD = DIST(J,M)
               IMIN = M
            end if
         end do
         NEIGHO(J,IP) = NEIGH(J,IMIN)
         DISTO(J,IP) = DIST(J,IMIN)
         NEIGH(J,IMIN) = NEIGH(J,NN)
         DIST(J,IMIN) = DIST(J,NN)
         NN = NN - 1
         if (NN == 0) EXIT
      end do
9909  CONTINUE
!
! SILICON NEIGHBORS TO A GIVEN OXYGEN ATOM
      ISINGO = 0
      IDOUB = 0
      do 99905 I = 1,NOX
      NES(I) = 0
      do 9905 J = 1,NS
         DELX = X01(J) - X02(I)
         if (DELX > CLMX) DELX = DELX - CLX
         if (DELX <= CLNX) DELX = DELX + CLX
         if (DABS(DELX) > SIO) GOTO 9905
         DELY = Y01(J) - Y02(I)
         if (DELY > CLMY) DELY = DELY - CLY
         if (DELY <= CLNY) DELY = DELY + CLY
         if (DABS(DELY) > SIO) GOTO 9905
         SQXY = DELY*DELY + DELX*DELX
         if (SQXY > SQO) GOTO 9905
         DELZ = Z01(J) - Z02(I)
         if (DELZ > CLMZ) DELZ = DELZ - CLZ
         if (DELZ <= CLNZ) DELZ = DELZ + CLZ
         if (DABS(DELZ) > SIO) GOTO 9905
         SQXYZ = DELZ*DELZ + SQXY
         if (SQXYZ > SQO) GOTO 9905
         NES(I) = NES(I) + 1
         NEIOO(I,NES(I)) = J
         DISO(I,NES(I)) = SQXYZ
9905  CONTINUE
! ORDER THE DISTANCES
      NN = NES(I)
      do IP = 1,1000
         DD = DISO(I,1)
         IMIN = 1
         do M = 2,NN
            if (DISO(I,M) < DD) then
               DD = DISO(I,M)
               IMIN = M
            end if
         end do
         NEIGOO(I,IP) = NEIOO(I,IMIN)
         DISOO(I,IP) = DISO(I,IMIN)
         NEIOO(I,IMIN) = NEIOO(I,NN)
         DISO(I,IMIN) = DISO(I,NN)
         NN = NN - 1
         if (NN == 0) EXIT
      end do
      if (NES(I) <= 1) then
         ISINGO = ISINGO + 1
         XISO(ISINGO) = X02(I)
         YISO(ISINGO) = Y02(I)
         ZISO(ISINGO) = Z02(I)
      else
         IDOUB = IDOUB + 1
         XIDO(IDOUB) = X02(I)
         YIDO(IDOUB) = Y02(I)
         ZIDO(IDOUB) = Z02(I)
      end if
99905 CONTINUE
!
! CHECK FOR TRIPLET SILANOLS
      ISIL(1:NS) = 0
      do J = 1,NS
         do K = 1,4
            if (NES(NEIGHO(J,K)) == 1) then
               ISIL(J) = ISIL(J) + 1
            end if
         end do
      end do
      IT = 0
      do J = 1,NS
         if (ISIL(J) == 3) then
            IT = IT + 1
            NTRIP(IT) = J
         end if
      end do
!
      XIDO(1:IDOUB) = XIDO(1:IDOUB)*1.62
      YIDO(1:IDOUB) = YIDO(1:IDOUB)*1.62
      ZIDO(1:IDOUB) = ZIDO(1:IDOUB)*1.62
      XISO(1:ISINGO) = XISO(1:ISINGO)*1.62
      YISO(1:ISINGO) = YISO(1:ISINGO)*1.62
      ZISO(1:ISINGO) = ZISO(1:ISINGO)*1.62
!
      do L = 1,IDOUB
         XS(1,L) = XIDO(L)
         YS(1,L) = YIDO(L)
         ZS(1,L) = ZIDO(L)
      end do
      do L = 1,ISINGO
         XS(2,L) = XISO(L)
         YS(2,L) = YISO(L)
         ZS(2,L) = ZISO(L)
      end do
      RETURN
   END SUBROUTINE SIOPOR

