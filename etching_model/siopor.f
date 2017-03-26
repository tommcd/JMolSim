   SUBROUTINE SIOPOR
! *********************************************************************
!      QUENCHED SILICA PORE
!      PROGRAM NOTE :
!        1. TWO COMPONENTS: SI(4+),O(2-)
!        2. CUTS A SILICA PORE AND DETERMINES PROPERTIES
!        3. DISCARDS SINGLY BONDED SURFACE OXYGEN TRIPLETS
!        4. LENGTH SCALES IN L(SIO) = 0.162 nm UNITS
! **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(IZS=2)
      PARAMETER(IZ1=7500)
      PARAMETER(IZ2=15000)
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
      PI=4.*DATAN(1.D0)
      PI2 = 2.0*PI
      RMUL2=16807.D0
      RMO1=2147483648.D0
      RIMO1=1./RMO1
      RMO2=2147483647.D0
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
      ISXX=DMOD(RMUL2*ISXX,RMO2)
      ISXX=DMOD(RMUL2*ISXX,RMO2)
      RN1=ISXX*RIMO1
      ISYY=DMOD(RMUL2*ISYY,RMO2)
      ISYY=DMOD(RMUL2*ISYY,RMO2)
      RN2=ISYY*RIMO1
!
      DO J=1,NSI
         X1(J) = X1(J)+(2.0*RN1-1.0)*CLX
         Y1(J) = Y1(J)+(2.0*RN2-1.0)*CLY
      END DO
      DO J=1,NOXX
         X2(J) = X2(J)+(2.0*RN1-1.0)*CLX
         Y2(J) = Y2(J)+(2.0*RN2-1.0)*CLY
      END DO
      DO J=1,NSI
         IF(X1(J) >  CLMX) X1(J)=X1(J)-CLX
         IF(X1(J) <= CLNX) X1(J)=X1(J)+CLX
         IF(Y1(J) >  CLMY) Y1(J)=Y1(J)-CLY
         IF(Y1(J) <= CLNY) Y1(J)=Y1(J)+CLY
      END DO
      DO J=1,NOXX
         IF(X2(J) >  CLMX) X2(J)=X2(J)-CLX
         IF(X2(J) <= CLNX) X2(J)=X2(J)+CLX
         IF(Y2(J) >  CLMY) Y2(J)=Y2(J)-CLY
         IF(Y2(J) <= CLNY) Y2(J)=Y2(J)+CLY
      END DO
      NS = 0
! ***********************************
      NSILON = 0
      DO I=1,NSI
         NSILON = NSILON+1
         NSILP(NSILON) = I
      END DO
      NREM = FREM*FLOAT(NSILON)+1.0
      IF(NREM == 0) GOTO 3337
      DO NINT = 1,10000
         ISZZ = DMOD(RMUL2*ISZZ,RMO2)
         ISZZ = DMOD(RMUL2*ISZZ,RMO2)
         RN = ISZZ*RIMO1
         MOLS = RN*FLOAT(NSILON)+1.0
         NS = NS+1
         NUM(NS) = NSILP(MOLS)
         NSILP(MOLS) = NSILP(NSILON)
         NSILON = NSILON-1
         NREM = NREM-1
         IF(NREM <= 0) GOTO 3337
      END DO
3337  CONTINUE
!
      DO 556 I=1,NS
         X01(I) = X1(NUM(I))
         Y01(I) = Y1(NUM(I))
         Z01(I) = Z1(NUM(I))
556   CONTINUE
      M=0
      DO 558 I=1,NS
      DO 557 K=1,4
         M = M+1
         KOX(M) = NEIOOR(NUM(I),K)
557   CONTINUE
558   CONTINUE
      NN = M
      DO 559 IP=1,20000
         IID = KOX(IP)
         IC = 0
         DO 10010 J=IP+1,NN
            IF(KOX(J) == IID) THEN
            IC = IC+1
            IMI(IC) = J
            END IF
10010    CONTINUE
         IF(IC == 0) GOTO 559
         DO 10110 IL=1,IC
            KOX(IMI(IL)) = KOX(NN)
            NN = NN-1
10110    CONTINUE
         IF(NN == IP) GOTO 7776
559   CONTINUE
7776  CONTINUE
      NOX = NN
      DO 444 I=1,NOX
         X02(I) = X2(KOX(I))
         Y02(I) = Y2(KOX(I))
         Z02(I) = Z2(KOX(I))
444   CONTINUE
!
! OXYGEN NEIGHBORS TO A GIVEN SILICON ATOM
      DO 9909 J=1,NS
      NEO(J) = 0
      DO 909 I=1,NOX
         DELX=X01(J)-X02(I)
         IF(DELX > CLMX) DELX=DELX-CLX
         IF(DELX <= CLNX) DELX=DELX+CLX
         IF(DABS(DELX) > SIO) GOTO 909
         DELY=Y01(J)-Y02(I)
         IF(DELY > CLMY) DELY=DELY-CLY
         IF(DELY <= CLNY) DELY=DELY+CLY
         IF(DABS(DELY) > SIO) GOTO 909
         SQXY=DELY*DELY+DELX*DELX
         IF(SQXY > SQO) GOTO 909
         DELZ=Z01(J)-Z02(I)
         IF(DELZ > CLMZ) DELZ=DELZ-CLZ
         IF(DELZ <= CLNZ) DELZ=DELZ+CLZ
         IF(DABS(DELZ) > SIO) GOTO 909
         SQXYZ=DELZ*DELZ+SQXY
         IF(SQXYZ > SQO) GOTO 909
         NEO(J) = NEO(J)+1
         NEIGH(J,NEO(J)) = I
         DIST(J,NEO(J)) = SQXYZ
909   CONTINUE
! ORDER THE DISTANCES
      NN = NEO(J)
      DO 707 IP=1,1000
         DD = DIST(J,1)
         IMIN = 1
         DO 1020 M=2,NN
            IF(DIST(J,M) < DD) THEN
               DD = DIST(J,M)
               IMIN = M
            END IF
1020     CONTINUE
         NEIGHO(J,IP) = NEIGH(J,IMIN)
         DISTO(J,IP) = DIST(J,IMIN)
         NEIGH(J,IMIN) = NEIGH(J,NN)
         DIST(J,IMIN) = DIST(J,NN)
         NN = NN-1
         IF(NN == 0) GOTO 7071
707   CONTINUE
7071  CONTINUE
9909  CONTINUE
!
! SILICON NEIGHBORS TO A GIVEN OXYGEN ATOM
      ISINGO = 0
      IDOUB = 0
      DO 99905 I=1,NOX
      NES(I) = 0
      DO 9905 J=1,NS
         DELX=X01(J)-X02(I)
         IF(DELX > CLMX) DELX=DELX-CLX
         IF(DELX <= CLNX) DELX=DELX+CLX
         IF(DABS(DELX) > SIO) GOTO 9905
         DELY=Y01(J)-Y02(I)
         IF(DELY > CLMY) DELY=DELY-CLY
         IF(DELY <= CLNY) DELY=DELY+CLY
         IF(DABS(DELY) > SIO) GOTO 9905
         SQXY=DELY*DELY+DELX*DELX
         IF(SQXY > SQO) GOTO 9905
         DELZ=Z01(J)-Z02(I)
         IF(DELZ > CLMZ) DELZ=DELZ-CLZ
         IF(DELZ <= CLNZ) DELZ=DELZ+CLZ
         IF(DABS(DELZ) > SIO) GOTO 9905
         SQXYZ=DELZ*DELZ+SQXY
         IF(SQXYZ > SQO) GOTO 9905
         NES(I) = NES(I)+1
         NEIOO(I,NES(I)) = J
         DISO(I,NES(I)) = SQXYZ
9905  CONTINUE
! ORDER THE DISTANCES
      NN = NES(I)
      DO 7705 IP=1,1000
         DD = DISO(I,1)
         IMIN = 1
         DO 1040 M=2,NN
            IF(DISO(I,M) < DD) THEN
               DD = DISO(I,M)
               IMIN = M
            END IF
1040     CONTINUE
         NEIGOO(I,IP) = NEIOO(I,IMIN)
         DISOO(I,IP) = DISO(I,IMIN)
         NEIOO(I,IMIN) = NEIOO(I,NN)
         DISO(I,IMIN) = DISO(I,NN)
         NN = NN-1
         IF(NN == 0) GOTO 7706
7705  CONTINUE
7706  CONTINUE
      IF(NES(I) <= 1) THEN
         ISINGO = ISINGO+1
         XISO(ISINGO) = X02(I)
         YISO(ISINGO) = Y02(I)
         ZISO(ISINGO) = Z02(I)
      ELSE
         IDOUB = IDOUB+1
         XIDO(IDOUB) = X02(I)
         YIDO(IDOUB) = Y02(I)
         ZIDO(IDOUB) = Z02(I)
      END IF
99905 CONTINUE
!
! CHECK FOR TRIPLET SILANOLS
      ISIL(1:NS) = 0
      DO J=1,NS
         DO K=1,4
            IF(NES(NEIGHO(J,K)) == 1) THEN
               ISIL(J) = ISIL(J)+1
            END IF
         END DO
      END DO
      IT = 0
      DO J=1,NS
         IF(ISIL(J) == 3) THEN
            IT = IT+1
            NTRIP(IT) = J
         END IF
      END DO
      DO L=1,IDOUB
         XIDO(L) = XIDO(L)*1.62
         YIDO(L) = YIDO(L)*1.62
         ZIDO(L) = ZIDO(L)*1.62
      END DO
      DO L=1,ISINGO
         XISO(L) = XISO(L)*1.62
         YISO(L) = YISO(L)*1.62
         ZISO(L) = ZISO(L)*1.62
      END DO
      DO L=1,IDOUB
         XS(1,L) = XIDO(L)
         YS(1,L) = YIDO(L)
         ZS(1,L) = ZIDO(L)
      END DO
      DO L=1,ISINGO
         XS(2,L) = XISO(L)
         YS(2,L) = YISO(L)
         ZS(2,L) = ZISO(L)
      END DO
      RETURN
   END SUBROUTINE SIOPOR

