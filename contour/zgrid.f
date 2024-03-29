      SUBROUTINE ZGRID(ZPIJ,KNXT,IMNEW,Z,NX,NY,X1N,XLN,Y1,DXN,DYN,
     % XP,YP,ZP,N,CAY,NRNG)
C     SETS UP SQUARE GRID FOR CONTOURING , GIVEN ARBITRARILY PLACED     
C     DATA POINTS. LAPLACE INTERPOLATION IS USED.                       
C     THE METHOD USED HERE WAS LIFTED DIRECTLY FROM NOTES LEFT BY       
C     MR IAN CRAIN FORMERLY WITH THE COMP.SCIENCE DIV.                  
C     INFO ON RELAXATION SOLN OF LAPLACE EQN SUPPLIED BY DR T MURTY.    
C     FORTRAN II   OCEANOGRAPHY/EMR   DEC/68   JDT                      
C                                                                       
C     Z = 2-D ARRAY OF HGTS TO BE SET UP. POINTS OUTSIDE REGION TO BE   
C     CONTOURED SHOULD BE INITIALIZED TO 10**35 . THE REST SHOULD BE 0.0
C     NX,NY = MAX SUBSCRIPTS OF Z IN X AND Y DIRECTIONS .               
C     X1,Y1 = COORDINATES OF Z(1,1)                                     
C     DX,DY = X AND Y INCREMENTS .                                      
C     XP,YP,ZP = ARRAYS GIVING POSITION AND HGT OF EACH DATA POINT.     
C     N = SIZE OF ARRAYS XP,YP AND ZP .                                 
C                                                                       
C     MODIFICATION FEB/69   TO GET SMOOTHER RESULTS A PORTION OF THE    
C     BEAM EQN  WAS ADDED TO THE LAPLACE EQN GIVING                     
C     DELTA2X(Z)+DELTA2Y(Z) - K(DELTA4X(Z)+DELTA4Y(Z)) = 0 .            
C     K=0 GIVES PURE LAPLACE SOLUTION. K=INF. GIVES PURE SPLINE SOLUTION
C     CAY = K = AMOUNT OF SPLINE EQN (BETWEEN 0 AND INF.)               
C     NRNG...GRID POINTS MORE THAN NRNG GRID SPACES FROM THE NEAREST    
C            DATA POINT ARE SET TO UNDEFINED.                           
C                                                                       
C     MODIFICATION DEC23/69   DATA PTS NO LONGER MOVED TO GRID PTS.     
C                                                                       
C***********************************************************************
C                                                                       
C  MODIFIED JUN30/75 TO REMOVE ALL PRINTED OUTPUT
C
C
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*4 X1N,XLN,Y1,DXN,DYN,CAY
      REAL*4 Z(NX,NY)
      REAL*4 XP(1),YP(1),ZP(1)                                       
      REAL*4 ZPIJ(N)
      INTEGER KNXT(N),IMNEW(NY)
      DATA EPS, BIG/ 0.002D0, 0.9E35/
      DATA ITMAX/ 100/
C
C  MODIFICATION OF 27 SEP 77 : X COORDINATES ARE NORMALIZED.
C
      DO 101   IX=1,NX
      DO 101   IY=1,NY
  101 Z(IX,IY)=0.0
      NTP1=N+1
      FACT1=XLN-X1N
      FACT2=1.0/FACT1
      X1=0.0
      DX=FACT2*DXN
      DO 112   IJK=1,NTP1
  112 XP(IJK)=FACT2*(XP(IJK)-X1N)
      DY=DYN
C                                                                       
C     GET ZBASE WHICH WILL MAKE ALL ZP VALUES POSITIVE BY 20*(ZMAX-ZMIN)
C***********************************************************************
C                                                                       
      ZMIN=ZP(1)
      ZMAX=ZP(1)
      DO 20 K=2,N
      IF(ZP(K)-ZMAX)14,14,12
12    ZMAX=ZP(K)
14    IF(ZP(K)-ZMIN)16,20,20
16    ZMIN=ZP(K)
20    CONTINUE
      ZRANGE=ZMAX-ZMIN
      ZBASE=ZRANGE*20.-ZMIN
      HRANGE=MIN(DX*(NX-1) , DY*(NY-1))
      DERZM=2.*ZRANGE/HRANGE
C
C     SET POINTER ARRAY KNXT
C***********************************************************************
C
      DO 60 KK=1,N
      K=1+N-KK
      KNXT(K)=0
      I= (XP(K)-X1)/DX + 1.5
      IF(I*(NX+1-I))60,60,35
35    J= (YP(K)-Y1)/DY + 1.5
      IF(J*(NY+1-J))60,60,40
40    IF(Z(I,J)-BIG)45,60,60
45    KNXT(K)=N+1
      IF(Z(I,J))55,55,50
50    KNXT(K)= Z(I,J)+.5
55    Z(I,J)=K
60    CONTINUE
C
C     AFFIX EACH DATA POINT ZP TO ITS NEARBY GRID POINT.  TAKE AVG ZP IF
C     MORE THAN ONE ZP NEARBY THE GRID POINT. ADD ZBASE AND COMPLEMENT. 
C***********************************************************************
C
      DO 80 K=1,N
      IF(KNXT(K))80,80,65
65    NPT=0
      ZSUM=0.
      I= (XP(K)-X1)/DX + 1.5
      J =(YP(K)-Y1)/DY + 1.5
      KK=K
70    NPT=NPT+1
      ZSUM=ZSUM+ ZP(KK)
      KNXT(KK)=-KNXT(KK)
      KK = -KNXT(KK)
      IF(KK-N)70,70,75
75    Z(I,J) = -ZSUM/NPT-ZBASE
80    CONTINUE
C
C     INITIALLY SET EACH UNSET GRID POINT TO VALUE OF NEAREST KNOWN PT. 
C***********************************************************************
C
      DO 110 I=1,NX
      DO 110 J=1,NY
C
C   CHANGE MADE 27 AUG 81 (ERROR FOUND BY R.WINTERBURN)
C     IF(Z(I,J))110,100,110
C100  Z(I,J) = -1.E35
C
      IF(ABS(Z(I,J)).LT.1.0E-6)   Z(I,J)= -1.E35
C
 110  CONTINUE
      DO 199 ITER=1,NRNG
      NNEW=0
      DO 197 I=1,NX
      DO 197 J=1,NY
      IF(Z(I,J)+BIG)152,192,192
152   IF(J-1)162,162,153
153   IF(JMNEW)154,154,162
154   ZIJN=ABS(Z(I,J-1))
      IF(ZIJN-BIG)195,162,162
162   IF(I-1)172,172,163
163   IF(IMNEW(J))164,164,172
164   ZIJN=ABS(Z(I-1,J))
      IF(ZIJN-BIG)195,172,172
172   IF(J-NY)173,182,182
173   ZIJN=ABS(Z(I,J+1))
      IF(ZIJN-BIG)195,182,182
182   IF(I-NX)183,192,192
183   ZIJN=ABS(Z(I+1,J))
      IF(ZIJN-BIG)195,192,192
192   IMNEW(J)=0
      JMNEW=0
      GO TO 197
195   IMNEW(J)=1
      JMNEW=1
      Z(I,J)=ZIJN
      NNEW=NNEW+1
197   CONTINUE
      IF(NNEW)200,200,199
199   CONTINUE
200   CONTINUE
      DO 202 I=1,NX
      DO 202 J=1,NY
      ABZ=ABS(Z(I,J))
      IF(ABZ-BIG)202,201,201
201   Z(I,J)=ABZ
202   CONTINUE
C
C     IMPROVE THE NON-DATA POINTS BY APPLYING POINT OVER-RELAXATION
C     USING THE LAPLACE-SPLINE EQUATION  (CARRES METHOD IS USED)
C***********************************************************************
C
      DZRMSP=ZRANGE
      RELAX=1.0
      DO 4000 ITER=1,ITMAX
      DZRMS=0.
      DZMAX=0.
      NPG =0
      DO 2000 I=1,NX
      DO 2000 J=1,NY
      Z00=Z(I,J)
      IF(Z00-BIG)205,2000,2000
205   IF(Z00)2000,208,208
208   WGT=0.
      ZSUM=0
C
      IM=0
      IF(I-1)570,570,510
510   ZIM=ABS(Z(I-1,J))
      IF(ZIM-BIG)530,570,570
530   IM=1
      WGT=WGT+1.
      ZSUM=ZSUM+ZIM
      IF(I-2)570,570,540
540   ZIMM=ABS(Z(I-2,J))
      IF(ZIMM-BIG)560,570,570
560   WGT=WGT+CAY
      ZSUM=ZSUM-CAY*(ZIMM-2.*ZIM)
570   IF(NX-I)700,700,580
580   ZIP=ABS(Z(I+1,J))
      IF(ZIP-BIG)600,700,700
600   WGT=WGT+1.
      ZSUM=ZSUM+ZIP
      IF(IM)620,620,610
610   WGT=WGT+4.*CAY
      ZSUM=ZSUM+2.*CAY*(ZIM+ZIP)
620   IF(NX-1-I)700,700,630
630   ZIPP=ABS(Z(I+2,J))
      IF(ZIPP-BIG)650,700,700
650   WGT=WGT+CAY
      ZSUM=ZSUM-CAY*(ZIPP-2.*ZIP)
700   CONTINUE
C
      JM=0
      IF(J-1)1570,1570,1510
1510  ZJM=ABS(Z(I,J-1))
      IF(ZJM-BIG)1530,1570,1570
1530  JM=1
      WGT=WGT+1.
      ZSUM=ZSUM+ZJM
      IF(J-2)1570,1570,1540
1540  ZJMM=ABS(Z(I,J-2))
      IF(ZJMM-BIG)1560,1570,1570
1560  WGT=WGT+CAY
      ZSUM=ZSUM-CAY*(ZJMM-2.*ZJM)
1570  IF(NY-J)1700,1700,1580
1580  ZJP=ABS(Z(I,J+1))
      IF(ZJP-BIG)1600,1700,1700
1600  WGT=WGT+1.
      ZSUM=ZSUM+ZJP
      IF(JM)1620,1620,1610
1610  WGT=WGT+4.*CAY
      ZSUM=ZSUM+2.*CAY*(ZJM+ZJP)
1620  IF(NY-1-J)1700,1700,1630
1630  ZJPP=ABS(Z(I,J+2))
      IF(ZJPP-BIG)1650,1700,1700
1650  WGT=WGT+CAY
      ZSUM=ZSUM-CAY*(ZJPP-2.*ZJP)
1700  CONTINUE
C
      DZ=ZSUM/WGT-Z00
      NPG=NPG+1
      DZRMS=DZRMS+DZ*DZ
      DZMAX=MAX(ABS(DZ),DZMAX)
      Z(I,J)=Z00+DZ*RELAX
2000  CONTINUE
C
C     SHIFT DATA POINTS ZP PROGRESSIVELY BACK TO THEIR PROPER PLACES AS 
C     THE SHAPE OF SURFACE Z BECOMES EVIDENT.
C***********************************************************************
C
      IF(ITER-(ITER/10)*10) 3600,3020,3600
3020  DO 3400 K=1,N
      KNXT(K) =IABS(KNXT(K))
      IF(KNXT(K))3400,3400,3030
3030  X=(XP(K)-X1)/DX
      I=X+1.5
      X= X+1.-I
      Y=(YP(K)-Y1)/DY
      J=Y+1.5
      Y=Y+1.-J
      ZPXY = ZP(K)+ZBASE
      Z00 = ABS(Z(I,J))
C
      ZW=1.E35
      IF(I-1)3120,3120,3110
3110  ZW = ABS(Z(I-1,J))
3120  ZE=1.E35
      IF(I-NX)3130,3140,3140
3130  ZE = ABS(Z(I+1,J))
3140  IF(ZE-BIG)3160,3150,3150
3150  IF(ZW-BIG)3180,3170,3170
3160  IF(ZW-BIG)3200,3190,3190
3170  ZE=Z00
      ZW=Z00
      GO TO 3200
3180  ZE=2.*Z00-ZW
      GO TO 3200
3190  ZW = 2.*Z00-ZE
C
3200  ZS=1.E35
      IF(J-1)3220,3220,3210
3210  ZS = ABS(Z(I,J-1))
3220  ZN= 1.E35
      IF(J-NY)3230,3240,3240
3230  ZN = ABS(Z(I,J+1))
3240  IF(ZN-BIG)3260,3250,3250
3250  IF(ZS-BIG)3280,3270,3270
3260  IF(ZS-BIG)3300,3290,3290
3270  ZN= Z00
      ZS= Z00
      GO TO 3300
3280  ZN = 2.*Z00-ZS
      GO TO 3300
3290  ZS = 2.*Z00-ZN
C
3300  A=(ZE-ZW)*.5
      B=(ZN-ZS)*.5
      C=(ZE+ZW)*.5-Z00
      D=(ZN+ZS)*.5-Z00
      ZXY=Z00+A*X+B*Y+C*X*X+D*Y*Y
      DELZ=Z00-ZXY
      DELZM=DERZM*(ABS(X)*DX+ABS(Y)*DY)*.80
      IF(DELZ-DELZM)3355,3355,3350
3350  DELZ=DELZM
3355  IF(DELZ+DELZM)3360,3365,3365
3360  DELZ=-DELZM
3365  ZPIJ(K)=ZPXY+DELZ
3400  CONTINUE
C
      DO 3500 K=1,N
      IF(KNXT(K))3500,3500,3410
3410  NPT=0
      ZSUM = 0.
      I= (XP(K)-X1)/DX + 1.5
      J= (YP(K)-Y1)/DY + 1.5
      KK = K
3420  NPT = NPT+1
      ZSUM = ZSUM + ZPIJ(KK)
      KNXT(KK)= -KNXT(KK)
      KK = -KNXT(KK)
      IF(KK-N)3420,3420,3430
3430  Z(I,J) =  -ZSUM/NPT
3500  CONTINUE
3600  CONTINUE
C
C     TEST FOR CONVERGENCE
C***********************************************************************
C
      IF (NPG) 3605,3605,3610
 3605 DZRMS=1.0D35
      GO TO 3620
3610  DZRMS=SQRT(DZRMS/NPG)
 3620 CONTINUE
      ROOT= DZRMS/DZRMSP
      DZRMSP=DZRMS
      DZMAXF=DZMAX/ZRANGE
      IF(ITER-(ITER/10)*10-2)3715,3710,3715
3710  DZRMS8 = DZRMS
3715  IF(ITER-(ITER/10)*10)4000,3720,4000
3720  ROOT = SQRT(SQRT(SQRT(DZRMS/DZRMS8)))
      IF(ROOT-.9999)3730,4000,4000
3730  IF(DZMAXF/(1.-ROOT)-EPS)4010,4010,3740
C
C     IMPROVE THE RELAXATION FACTOR.
C***********************************************************************
C
3740  IF((ITER-20)*(ITER-40)*(ITER-60))4000,3750,4000
3750  IF(RELAX-1.-ROOT)3760,4000,4000
3760  TPY =(ROOT+RELAX-1.)/RELAX
      ROOTGS = TPY*TPY/ROOT
      RELAXN= 2./(1.+SQRT(1.-ROOTGS))
      IF(ITER-60)3780,3785,3780
3780  RELAXN= RELAXN-.25*(2.-RELAXN)
3785  RELAX = MAX(RELAX,RELAXN)
4000  CONTINUE
4010  CONTINUE
C
C     REMOVE ZBASE FROM ARRAY Z AND RETURN.
C***********************************************************************
C
      DO 4500 I=1,NX
      DO 4500 J=1,NY
      IF(Z(I,J)-BIG)4400,4500,4500
4400  Z(I,J)=ABS(Z(I,J))-ZBASE
4500  CONTINUE
C  MODIFICATION OF 27 SEP 77 : X COORDINATES ARE RESTORED TO
C  THEIR ORIGINAL VALUES.
      DO 113   IJK=1,NTP1
  113 XP(IJK)=FACT1*XP(IJK)+X1N
      RETURN
      END
