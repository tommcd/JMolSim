         PROGRAM SILICA
C
C
C       This program simulates the silica substance
C      using (BMH) pair potential with adding three-body
C      interections by the Stillinger and Weber potential  
C
C 

C
C*****   READ 'input.dat'  *******************
        COMMON / BLOCK1 / VX, VY, VZ
        COMMON / BLOCK2 / POSITION,FR,CL,VEL,pot_energ,number5                      
    
                         
      DIMENSION POSITION(1536,3),POSITIONNEW(1536,3),VEL(1536,3)
	DIMENSION FR(1536,3),VELNEW(1536,3),ij(3),VZ(1536) 
	DIMENSION FROLD(1536,3),L(3),RDF(100000),VX(1536),VY(1536)
      double precision POSITION,POSITIONNEW,FR,FROLD,Ek,in(3)
	double precision CL,VEL,DT,mOx,mSi,time,L,DL,lenth,t1,jn(3)
      REAL*8 VX,VY,VZ,OX(1536),OY(1536),OZ(1536),T,ji(3),nj(3),ni(3)
      REAL*8 pot_energ,AIJ,AIJ1,AIJ2,pi,e,ij,t5,erfc1,ro,e1
	real*8 bettaSiSi,bettaSiO,bettaOO,ss,angleSi,angleOx,const1Ox
	real*8 const2Ox,rOxc,const1Si,const2Si,rSic,distance
	PARAMETER(pi=3.14159265358979323846d0)
	PARAMETER(AIJ=1.88*(1.0E-16))
	PARAMETER(AIJ1=0.720*(1.0E-16))
	PARAMETER(AIJ2=2.96*(1.0E-16))
	PARAMETER(ro=0.29d0*(1.0E-10))
	PARAMETER(bettaSiSi=2.53d0*(1.0E-10))
	PARAMETER(bettaSiO=2.6d0*(1.0E-10)) 
	PARAMETER(bettaOO=2.55d0*(1.0E-10))
      PARAMETER(e1=1.602176462d0*(1.0E-19))
	PARAMETER(angleSi=3.1415926535897932384d0*2.0d0*109.471d0/360.0d0)
      PARAMETER(angleOx=3.1415926535897932384d0*2.0d0*109.471d0/360.0d0)                                                                           
      PARAMETER(const1Ox=0.3E-18)
      PARAMETER(const2Ox=2.0*(1.0E-10))
      PARAMETER(rOxc=2.6*(1.0E-10))
      PARAMETER(const1Si=18.0E-18)
	PARAMETER(const2Si=2.6*(1.0E-10))
	PARAMETER(rSic=3.0*(1.0E-10))
	PARAMETER(distance=2.0*(1.0E-10))

      integer I,j,n,RDF,nn,number5
      e=(8987551787.36818d0)*e1**2
      ss=0.5d0

 
      OPEN(unit=13,FILE='betasio1.dat',STATUS='OLD')
	OPEN(unit=14,FILE='RDF.out',STATUS='NEW')
      DO 12 I=1,1536
      READ(13,*) POSITION(I,1),POSITION(I,2),POSITION(I,3)
12    continue     
      CLOSE(UNIT=13,STATUS='KEEP')
      
C------------------------------------------------------------
      CL=4.0d0*(1.62d0*2.0d0/dsqrt(3.0d0))*(1.0E-10)  
      DO 70 I=1,1536
	DO 71 J=1,3
	POSITION(i,j)=POSITION(i,j)*(1.62d0*2.0d0/dsqrt(3.0d0))*(1.0E-10)
 71    continue
 70    continue
       T = 300.0d0
  
C     Determination of the proximity for every atom
	 DO i=1,1536
	        s=0
	  DO 65 j=1,1536
		                DO n=1,3
	ij(n)=POSITION(j,n)-POSITION(i,n)  
	ij(n)=ij(n)-(4.0d0*CL)*DNINT(ij(n)/(4.0d0*CL))
	                    END DO
			IF (lenth(ij).LE.(1.627d0)) THEN
	s=s+1
	proximity(i,s)=j
	END IF

	  END DO

	 END DO
C********** Doing the melting process ********************************

C	          DO 112 ppp=1,50
C	             DO 111 count=1,1536
C	DO n=1,3
C       copnew_atom(count,n,1)=new_atom(count,n,1)
C	END DO
C       copatom(count)=atom(count)
C
C            call energy4
C            call repul_energy1(count)
C		  energytotl1=energy      

C  67	new_atom(count,1,1)=new_atom(count,1,1)+(2.0d0*rand(isx)-1.0d0)*de
C	new_atom(count,2,1)=new_atom(count,2,1)+(2.0d0*rand(isx)-1.0d0)*de
C	new_atom(count,3,1)=new_atom(count,3,1)+(2.0d0*rand(isx)-1.0d0)*de
C
C      new_atom(count,1,1)=new_atom(count,1,1) - 
C     & boxl*DNINT(new_atom(count,1,1)*boxli)
C	new_atom(count,2,1)=new_atom(count,2,1) - 
C     & boxl*DNINT(new_atom(count,2,1)*boxli)
C      new_atom(count,3,1)=new_atom(count,1,1) - 
C     & boxl*DNINT(new_atom(count,3,1)*boxli)

      
            
C	      call energy4
C            call repul_energy1(count)

C		  energytotl=energy
C            energytotl2=energytotl-energytotl1  

C	      IF (count.EQ.i) THEN
C	EEE = energytotl2/E0
C		If(EEE.lt.0.0) go to 10002
C		if(EEE.gt.50.0) go to 10001
C           rrr=rand(isx)
C	           IF (rrr.lt.dexp(-EEE)) GO TO 10002
C10001  CONTINUE
C	DO n=1,3
C       new_atom(count,n,1)=copnew_atom(count,n,1)
C	END DO
C        atom(count)=copatom(count)
C10002  CONTINUE
C	      END IF

C	      IF (count.NE.i) THEN
C		EEE = energytotl2/Etot
C			  rrr=rand(isx)          
C		If(EEE.lt.0.0) go to 10004
C		if(EEE.gt.50.0) go to 10003
C	           IF (rrr.lt.dexp(-EEE)) GO TO 10004
C10003     CONTINUE
C				DO n=1,3
C				 new_atom(count,n,1)=copnew_atom(count,n,1)
C				END DO
C				atom(count)=copatom(count)
C10004    CONTINUE
C	      END IF

C 111   CONTINUE
C 112    CONTINUE
       
C      DO 44 I=1,1536
C	DO 45 J=1,3
C	FROLD(i,j)=FR(i,j)
C 45    continue
C 44    continue
C	 DT=0.001d0*(1.0E-12)
C	 DL=0.01d0*(1.0E-10)
C	mOx=16.00d0/(6.02214199d0*(1.0e+26))
C	mSi=28.06d0/(6.02214199d0*(1.0e+26))

       time=0.0d0
C                              TIME OF CALCULATION
 50              IF ((time).LE.(600.0d0*(1.0E-12))) THEN
	  write(*,*) time
	number5=number5+1
C                     USING THE VERLET-POSITION ALGORITHM FOR OXIGEN ATOM

                        DO 36 I=1,1024
                             DO 37 J=1,3
      POSITIONNEW(I,J)=POSITION(I,J)+DT*vel(i,j)+(ss)*DT**2*FR(I,J)/mOx
                IF (POSITIONNEW(i,j).GT.t) THEN 
	POSITIONNEW(i,j)=POSITIONNEW(i,j)-t*DINT(POSITIONNEW(i,j)/t)
	          END IF
	IF (POSITIONNEW(i,j).LT.0.0d0) THEN 
	t1=-1.0d0+DINT(POSITIONNEW(i,j)/t)
	POSITIONNEW(i,j)=POSITIONNEW(i,j)-t*t1
	          END IF
        POSITION(I,J)=POSITIONNEW(I,J) 
  37                            continue                               
  36                       continue
                        
C                     USING THE VERLET-POSITION ALGORITHM FOR SILICON ATOM

                        DO 38 I=1025,1536
                             DO 39 J=1,3
      POSITIONNEW(I,J)=POSITION(I,J)+DT*vel(i,j)+(ss)*DT**2*FR(I,J)/mSi
                IF (POSITIONNEW(i,j).GT.t) THEN 
	POSITIONNEW(i,j)=POSITIONNEW(i,j)-t*DINT(POSITIONNEW(i,j)/t)
	          END IF
	IF (POSITIONNEW(i,j).LT.0.0d0) THEN 
	t1=-1.0d0+DINT(POSITIONNEW(i,j)/t)
	POSITIONNEW(i,j)=POSITIONNEW(i,j)-t*t1
	          END IF
        POSITION(I,J)=POSITIONNEW(I,J)

  39                            continue                                
  38                       continue


        call FORCE 

C                     USING THE VERLET-VELOCITY ALGORITHM FOR OXIGEN ATOM

                        DO 40 I=1,1024
                             DO 41 J=1,3
      VELNEW(I,J)=VEL(I,J)+(0.5d0)*DT*(FROLD(i,j)/mOx+FR(i,j)/mOx)
        VEL(I,J)=VELNEW(I,J) 
  41                            continue                                
  40                       continue

C                     USING THE VERLET-VELOCITY ALGORITHM FOR SILICON ATOM

                        DO 46 I=1025,1536
                             DO 47 J=1,3
      VELNEW(I,J)=VEL(I,J)+(0.5d0)*DT*(FROLD(i,j)/mSi+FR(i,j)/mSi)
        VEL(I,J)=VELNEW(I,J) 
  47                            continue                                
  46                       continue

      DO 42 I=1,1536
	DO 43 J=1,3
	FROLD(i,j)=FR(i,j)
 43    continue
 42    continue
       IF (number5.EQ.500) THEN
	Ek=0.0d0
	DO i=1,1024
	Ek=Ek+mOx*(VEL(i,1)**2+VEL(i,2)**2+VEL(i,3)**2)/2.0d0
	END DO
	DO i=1025,1536
	Ek=Ek+mSi*(VEL(i,1)**2+VEL(i,2)**2+VEL(i,3)**2)/2.0d0
	END DO
	pot_energ=0.0d0
	 DO i=1,1536
	  DO 64 j=i+1,1536
		                DO n=1,3
	ij(n)=POSITION(j,n)-POSITION(i,n)
	ij(n)=ij(n)-(4.0d0*CL)*DNINT(ij(n)/(4.0d0*CL))
	                   END DO
	  IF((lenth(ij).EQ.0.0d0).OR.(lenth(ij).GT.9.0d0*(1.0E-10))) THEN
	GO TO 64
	  END IF

	                
	    IF (i.GE.1025) THEN
      t5=16.0d0*(e)*erfc1(lenth(ij)/bettaSiSi)/(lenth(ij))
	pot_energ=pot_energ+AIJ*dexp(-lenth(ij)/ro)+t5
	    END IF
	IF ((i.LE.1024).AND.(j.LE.1024)) THEN
      t5=4.0d0*(e)*erfc1(lenth(ij)/bettaOO)/(lenth(ij))
	pot_energ=pot_energ+AIJ1*dexp(-lenth(ij)/ro)+t5
	END IF
	IF ((i.LE.1024).AND.(j.GE.1025)) THEN
      t5=-8.0d0*(e)*erfc1(lenth(ij)/bettaSiO)/(lenth(ij))
	pot_energ=pot_energ+AIJ2*dexp(-lenth(ij)/ro)+t5
	END IF
 64	  continue
	 END DO

	write(14,*) time,Ek,Ek+pot_energ
	number5=0
	END IF
	       time=time+DT
       GO TO 50
	END IF

	      DO 54 i=1,1536
	DO 55 j=i+1,1536
	  DO n=1,3
      L(n)=POSITION(j,n)-POSITION(i,n)
	L(n)=L(n)-(4.0d0*CL)*DNINT(L(n)/(4.0d0*CL))
	  END DO
       RDF(DINT(lenth(l)/DL))=RDF(DINT(lenth(l)/DL))+1

 55    continue
 54    continue
      DO n=1,DNINT(8*CL/DL)
	write(14,*) DL*n,RDF(n)
	END DO
	      CLOSE(UNIT=14,STATUS='KEEP')
 
	 STOP
         end
C------------------------------------------------------------ 

********************************************************************************
** FICHE F.24.  INITIAL VELOCITY DISTRIBUTION                                 **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** CENTRE OF MASS AND ANGULAR VELOCITIES FOR LINEAR MOLECULES    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   THE NUMBER OF MOLECULES           **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    VX(N),VY(N),VZ(N)   VELOCITIES                        **
C    ** REAL    EX(N),EY(N),EZ(N)   ORIENTATIONS                      **
C    ** REAL    OX(N),OY(N),OZ(N)   SPACE-FIXED ANGULAR VELOCITIES    **
C    ** REAL    TEMP                REDUCED TEMPERATURE               **
C    ** REAL    INERT               REDUCED MOMENT OF INERTIA         **
C    **                                                               **
C    ** SUPPLIED ROUTINES:                                            **
C    **                                                               **
C    ** SUBROUTINE COMVEL ( TEMP )                                    **
C    **    SETS THE CENTRE OF MASS VELOCITIES FOR A CONFIGURATION OF  **
C    **    LINEAR MOLECULES AT A GIVEN TEMPERATURE.                   **
C    ** SUBROUTINE ANGVEL ( TEMP, INERT )                             **
C    **    SETS THE ANGULAR VELOCITIES FOR A CONFIGURATION OF LINEAR  **
C    **    MOLECULES AT A GIVEN TEMPERATURE.                          **
C    ** REAL FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE  **
C    ** REAL FUNCTION GAUSS ( DUMMY )                                 **
C    **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
C    **    DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             **
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** WE ASSUME UNIT MOLECULAR MASS AND EMPLOY LENNARD-JONES UNITS  **
C    **       PROPERTY                      UNITS                     **
C    **       RX, RY, RZ           (EPSILON/M)**(1.0/2.0)             **
C    **       OX, OY, OZ           (EPSILON/M*SIGMA**2)**(1.0/2.0)    **
C    **       INERT                 M*SIGMA**2                        **
C    *******************************************************************

        SUBROUTINE COMVEL ( TEMP )

        COMMON / BLOCK1 / VX, VY, VZ
C    *******************************************************************
C    ** TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  **
C    **                                                               **
C    ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (UNIT) MASS.**
C    ** THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    **
C    ** MOLECULES, AND NON-LINEAR MOLECULES.                          **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           ** 
C    **                                                               **
C    ** REAL FUNCTION GAUSS ( DUMMY )                                 **
C    **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
C    **    DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             **
C    *******************************************************************
C*****************************************
C*****  INITIAL VELOCITY FOR MOVING PARTICLE        
C***** 
      real*8 PI,PI2,RMUL2,IAX1,IAX2,RMO2,IAY1,IAY2,IAZ1,IAZ2
	real*8 RIMO1,RNX1,RNX2,RNY1,RNY2,RNZ1,RNZ2
	real*8 VX(1536),VY(1536),VZ(1536),TEMP
      PI = 4.0D0*DATAN(1.0D0)
	PI2 = 2.0*PI
	IAX1=987321.0d0
	IAY1=199971.0d0
	IAZ1=673923.0d0
	IAX2=786391.0d0
	IAY2=383837.0d0
	IAZ2=921217.0d0
	RMUL2=16807.0d0
	RMO2=2147483647.0d0
	RIMO1=1.0d0/(2147483648.0d0) 
      DO 44220 I=1,1536   
C                                                                 
C*** INITIAL VELOCITY FROM MAXWELLIAN DISTRIBUTION *****     
C                                                                    
      IAX1 = DMOD(RMUL2*IAX1,RMO2)                      
      IAX1 = DMOD(RMUL2*IAX1,RMO2)                
      IAX2 = DMOD(RMUL2*IAX2,RMO2)                        
      IAX2 = DMOD(RMUL2*IAX2,RMO2)            
      IAY1 = DMOD(RMUL2*IAY1,RMO2)            
      IAY1 = DMOD(RMUL2*IAY1,RMO2)         
      IAY2 = DMOD(RMUL2*IAY2,RMO2)                        
      IAY2 = DMOD(RMUL2*IAY2,RMO2)                   
      IAZ1 = DMOD(RMUL2*IAZ1,RMO2)                    
      IAZ1 = DMOD(RMUL2*IAZ1,RMO2)          
      IAZ2 = DMOD(RMUL2*IAZ2,RMO2)                         
      IAZ2 = DMOD(RMUL2*IAZ2,RMO2)                    
C                                                                     
      RNX1 = IAX1*RIMO1                                       
      RNX2 = IAX2*RIMO1                                           
      RNY1 = IAY1*RIMO1                                             
      RNY2 = IAY2*RIMO1                              
      RNZ1 = IAZ1*RIMO1                                   
      RNZ2 = IAZ2*RIMO1                      
C
      VX(I)=dsqrt(TEMP)*DSQRT(-2.D0*DLOG(RNX1))*DCOS(PI2*RNX2)                
      VY(I)=dsqrt(TEMP)*DSQRT(-2.D0*DLOG(RNY1))*DCOS(PI2*RNY2)                  
      VZ(I)=dsqrt(TEMP)*DSQRT(-2.D0*DLOG(RNZ1))*DCOS(PI2*RNZ2)                  
44220  CONTINUE
       return
	end

        SUBROUTINE FORCE 

        COMMON / BLOCK2 / POSITION,FR,CL,VEL,pot_energ,number5

      DIMENSION POSITION(1536,3)
	DIMENSION  k(3),km(3),m(3),number(1024),lj(3),number1(1024,1024)
	DIMENSION  u(3),v(3),w(3),uv(3),uw(3),vw(3),wv(3),wu(3),vu(3)
	DIMENSION FR(1536,3),number3(1024),number4(1024),ij(3) 
      double precision POSITION,FR,k,km,m,DUMMY1,distance,t4
	double precision u,v,w,uv,uw,vw,wv,wu,CL,lenth,lj,vu,erfc
	real*8 angleOx,angleSi,const1Ox,const2Ox,rOxc,const1Si,pi,lenth3
	real*8 const2Si,rSic,multpl,BMHSW3,BMHSW2,BMHSW1,ij,t,t1,t2,t3
	real*8 AIJ,AIJ1,AIJ2,e,erfc1,VEL(1536,3),t5,pot_energ,ij1(3),ro
	real*8 bettaSiSi,bettaSiO,bettaOO,e1,energy1,lenth1,lenth2
      integer i,n,p,number,s,j,g,f,number1,number3,number4,l,number5                                   
      PARAMETER(angleSi=3.1415926535897932384d0*2.0d0*109.471d0/360.0d0)
      PARAMETER(angleOx=3.1415926535897932384d0*2.0d0*109.471d0/360.0d0)                                                                           
      PARAMETER(const1Ox=0.3E-18)
      PARAMETER(const2Ox=2.0*(1.0E-10))
      PARAMETER(rOxc=2.6*(1.0E-10))
      PARAMETER(const1Si=18.0E-18)
	PARAMETER(const2Si=2.6*(1.0E-10))
	PARAMETER(rSic=3.0*(1.0E-10))
	PARAMETER(distance=2.0*(1.0E-10))
	PARAMETER(pi=3.14159265358979323846d0)
	PARAMETER(AIJ=1.88*(1.0E-16))
	PARAMETER(AIJ1=0.720*(1.0E-16))
	PARAMETER(AIJ2=2.96*(1.0E-16))
	PARAMETER(ro=0.29d0*(1.0E-10))
	PARAMETER(bettaSiSi=2.53d0*(1.0E-10))
	PARAMETER(bettaSiO=2.6d0*(1.0E-10)) 
	PARAMETER(bettaOO=2.55d0*(1.0E-10))
      PARAMETER(e1=1.602176462d0*(1.0E-19))
	e=(8987551787.36818d0)*(e1)**2
         DO 52 i=1,1536
           DO 53 j=1,3
	FR(i,j)=0.0d0
  53       continue
  52     continue
        DO 37 i=1,1536
                DO 38 j=1,1536
	                DO n=1,3
	ij(n)=POSITION(j,n)-POSITION(i,n)
	ij(n)=ij(n)-(4.0d0*CL)*DNINT(ij(n)/(4.0d0*CL))
	                END DO
	 lenth1=lenth(ij)
	IF((lenth1.EQ.0.0d0).OR.(lenth1.GT.9.0d0*(1.0E-10))) THEN
	GO TO 38
	END IF
		IF (j.EQ.i) THEN
	GO TO 38
	END IF

             IF ((i.GE.1025).AND.(j.GE.1025)) THEN
      		  DO f=1,3
	t=AIJ*dexp(-lenth1/ro)/(ro*lenth1)
	t4=erfc1(lenth1/bettaSiSi)
	t1=(16.0d0*e/(lenth1**2*lenth1))*t4
	t3=dexp(-lenth1**2/(bettaSiSi)**2)
	t2=t3*16.0d0*e*2.0d0/((lenth1)**2*bettaSiSi*dsqrt(pi))
	FR(i,f)=-(ij(f))*(t+t1+t2)+FR(i,f)
	    	  END DO 
	      END IF
         IF ((i.GE.1025).AND.(j.LE.1024)) THEN
      		  DO f=1,3
	t=AIJ2*dexp(-lenth1/ro)/(ro*lenth1)
	t4=erfc1(lenth1/bettaSiO)
	t1=(-8.0d0*e/(lenth1**2*lenth1))*t4
	t3=dexp(-lenth1**2/(bettaSiO)**2)
	t2=-8.0d0*(e)*t3*2.0d0/((lenth1)**2*bettaSiO*dsqrt(pi))
	FR(i,f)=-(ij(f))*(t+t1+t2)+FR(i,f) 
	    	  END DO 
	   END IF

	IF ((i.LE.1024).AND.(j.LE.1024)) THEN
      		  DO f=1,3
	t=AIJ1*dexp(-lenth1/ro)/(ro*lenth1)
	t4=erfc1(lenth1/bettaOO)
	t1=(4.0d0*e/(lenth1**2*lenth1))*t4
	t3=dexp(-lenth1**2/(bettaOO)**2)
	t2=4.0d0*(e)*2.0d0/((lenth1)**2*bettaOO*dsqrt(pi))*t3
	FR(i,f)=-(ij(f))*(t+t1+t2)+FR(i,f) 
	    	  END DO
       END IF
      IF((i.LE.1024).AND.(j.GE.1025)) THEN
      		  DO f=1,3
	t=AIJ2*dexp(-lenth1/ro)/(ro*lenth1)
	t4=erfc1(lenth1/bettaSiO)
	t1=(-8.0d0*e/(lenth1**2*lenth1))*t4
	t3=dexp(-lenth1**2/(bettaSiO)**2)
	t2=-8.0d0*(e)*t3*2.0d0/((lenth1)**2*bettaSiO*dsqrt(pi))
	FR(i,f)=-(ij(f))*(t+t1+t2)+FR(i,f) 
	    	  END DO 
	      END IF
  38            continue
  37    continue
C------------------------------------------------------------
        RETURN
        END
C------------------------------------------------------------
C        Function subprogram for the force function
C
       real*8  Function  BMHSW1(j,i,k,anglei,const1i,const2i,ric,l)
         DIMENSION i(3),J(3),k(3),jk(3),ji(3),ki(3),kj(3),ij(3),ik(3)
           double precision k,j,i,anglei,const1i,const2i,jk,ji,ki,kj,ij
	     double precision ik,t2,t3,t4,t5,ric,lenth,multpl 
	     real*8 s1,s2,s3
           integer l,n
           DO 11 n=1,3
             jk(n)=k(n)-j(n)
             ji(n)=i(n)-j(n)
             ki(n)=i(n)-k(n)
             kj(n)=j(n)-k(n)
             ij(n)=j(n)-i(n)
             ik(n)=k(n)-i(n)
 11           continue  
      s1=lenth(ij)
	s2=lenth(ik)
      t1=dexp(const2i/(s1-ric)+const2i/(s2-ric))
	s3=multpl(ij,ik)
      t1=t1*(s3/(s1*s2)-dcos(anglei))
	t1=t1*const1i
      t2=((j(l)-i(l))/((s1-ric)**2)*s1)
	t2=t2+(k(l)-i(l))/((s2-ric)**2*s2)
      t3=(s3/(s1*s2)-dcos(anglei))*const2i
	t3=t3*t2
      t4=-2*(ij(l)+ik(l))/(s1*s2)
      t5=2*s3*(ij(l)/((s1**2)*s1*s2)+ik(l)/((s2**2)*s2*s1))
      BMHSW1=(-1.0d0)*t1*(t4+t5+t3)
         return
            end
C-------------------------------------------------------------
        real*8  Function  BMHSW2(k,j,i,anglej,const1j,const2j,rjc,l)
         DIMENSION i(3),J(3),k(3),jk(3),ji(3),ki(3),kj(3),ij(3),ik(3)
           double precision k,j,i,jk,ji,ki,kj,ij,ik,anglej,const1j
		 double precision t,const2j,t1,t3,t4,rjc,lenth,multpl
	     real*8 s1,s2,s3
            integer l,n
           DO 11 n=1,3
             jk(n)=k(n)-j(n)
             ji(n)=i(n)-j(n)
             ki(n)=i(n)-k(n)
             kj(n)=j(n)-k(n)
             ij(n)=j(n)-i(n)
             ik(n)=k(n)-i(n)
 11           continue 
      s1=lenth(ji)
	s2=lenth(jk)
	s3=multpl(jk,ji)
      t1=(s3/(s1*s2)-dcos(anglej))*const1j
	t=t1
	t1=t*dexp(const2j/(s2-rjc)+const2j/(s1-rjc))
      t2=-(i(l)-j(l))/((s1-rjc)**2*s1)
	t=(s3/(s1*s2)-dcos(anglej))*const2j
	t2=t2*t
	t=2.0d0*jk(l)/(s2*s1)
      t3=-ji(l)*2.0d0*s3/(s2*(s1**2)*s1)
      BMHSW2=(-1.0d0)*t1*(t2+t+t3)
              return
                 end


C-------------------------------------------------------------                 


        real*8 Function  BMHSW3(i,k,j,anglek,const1k,const2k,rkc,l)
         DIMENSION i(3),J(3),k(3),jk(3),ji(3),ki(3),kj(3),ij(3),ik(3)
           double precision k,j,i,jk,ji,ki,kj,ij,ik,anglek,const1k
		 real*8    t,const2k,t1,t2,t4,rkc,lenth,multpl
	     real*8 s1,s2,s3
           integer l,n
           DO 11 n=1,3
             jk(n)=k(n)-j(n)
             ji(n)=i(n)-j(n)
             ki(n)=i(n)-k(n)
             kj(n)=j(n)-k(n)
             ij(n)=j(n)-i(n)
             ik(n)=k(n)-i(n)
 11           continue 
      s1=lenth(ki)
	s2=lenth(kj)
	s3=multpl(ki,kj)
      t=dexp(const2k/(s1-rkc)+const2k/(s2-rkc))*const1k  
      t1=(s3/(s1*s2)-dcos(anglek))*t
	t=(s3/(s1*s2)-dcos(anglek))*const2k
      t2=-(i(l)-k(l))*t/((s1-rkc)**2*s1)
	t3=2.0d0*(j(l)-k(l))/(s1*s2)
      t4=-2.0d0*s3*(i(l)-k(l))/((s1**2)*s1*s2)
      BMHSW3=(-1.0d0)*t1*(t2+t3+t4)
	     return
             end
         
C------------------------------------------------------------        
C        Function subprogram for the vector dot-product
C  
          Double precision Function multpl(vector1,vector2)
          DIMENSION vector1(3),vector2(3)
          double precision vector1,vector2
	       integer i
             multpl=0.0d0
             DO 10 i=1,3
             multpl=multpl+(vector1(i)*vector2(i))
 10          continue
            return
              end  
              
C------------------------------------------------------------                    
C         Function subprogram for the absolut value of the vector
C
          Double precision Function lenth(vector1)
          DIMENSION vector1(3)
           double precision vector1
	     integer y
           lenth=0.0d0
           DO 10 y=1,3
             lenth=lenth+vector1(y)**2
 10           continue
             lenth=dsqrt(lenth)   
               return
                 end
        REAL*8 FUNCTION ERFC1 ( X )

C    *******************************************************************
C    ** APPROXIMATION TO THE COMPLEMENTARY ERROR FUNCTION             **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,    **
C    **    NATIONAL BUREAU OF STANDARDS, FORMULA 7.1.26               **
C    *******************************************************************

        REAL*8        A1, A2, A3, A4, A5, P

        PARAMETER ( A1 = 0.254829592d0, A2 = -0.284496736d0 )
        PARAMETER ( A3 = 1.421413741d0, A4 = -1.453152027d0 )
        PARAMETER ( A5 = 1.061405429d0, P  =  0.3275911d0   )

        REAL*8        T, X, XSQ, TP

C    *******************************************************************

        T  = 1.0 / ( 1.0 + P * X )
        XSQ = X * X

        TP = T * ( A1 + T * ( A2 + T * ( A3 + T * ( A4 + T * A5 ) ) ) )

        ERFC1 = TP*DEXP ( -XSQ )

        RETURN
        END

       REAL*8 FUNCTION energy1(ki,kj,anglek,const1k,const2k,rkc)
         DIMENSION ki(3),kj(3)
           double precision ki,kj,anglek,const1k
		 real*8    t,const2k,t1,t2,t4,rkc,lenth,multpl
	     real*8 s1,s2,s3
           integer l,n
      s1=lenth(ki)
	s2=lenth(kj)
	s3=multpl(ki,kj)
      t=dexp(const2k/(s1-rkc)+const2k/(s2-rkc))*const1k  
      energy1=(s3/(s1*s2)-dcos(anglek))**2*t
	return
	end