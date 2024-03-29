      SUBROUTINE CONSEG(Z,NX,NY,X1,Y1,XL,YL,ZLEV,NDECL,LWGTL,NLEV,
     & HGT,NDIV,NARC,X,Y,S)

C     PLOTS CONTOURS OF GRID Z(I,J) IN SQUARE SEGMENTS OF LENGTH LSEG
C     GRID UNITS PER SIDE.  EACH GRID SQUARE IS SUBDIVIDED INTO NDIV**2
C     SUBSQUARES USING CUBIC POLYNOMIAL INTERPOLATION. NDIV MUST BE A
C     POWER OF 2.
C     OCEANOGRAPHY EMR  DEC/1969
C

      INTEGER LWGTL(NLEV), NDECL(NLEV)
      REAL X(1), Y(1), S(1)
      REAL Z(NX,NY)
      REAL ZZ(57,57)
      REAL ZLEV(NLEV)
      COMMON /XYS/ IND10(57,6)

      IZZ= 57
      LSEG= (IZZ-1)/NDIV-2
      DX =(XL-X1)/(NX-1)
      DY =(YL-Y1)/(NY-1)
C
      IA=1
40    IB = MIN0(IA+LSEG,NX)
      LREM = NX-IB
      IF( LREM*(LSEG-LREM))60,60,50
50    IB = (1+IA+NX)/2
60    IAM = MAX0(IA-1,1)
      IBM = MIN0(IB+1,NX)
      LXM = IBM-IAM
      XX1 = X1+(IA-1)*DX
      XXL = X1+(IB-1)*DX
      IIM = 1+(IA-IAM)*NDIV
      NXX =(IB-IA)*NDIV+1
C
C
      JA=1
80    JB = MIN0(JA+LSEG,NY)
      LREM = NY-JB
      IF(LREM*(LSEG-LREM))100,100,90
90    JB =(1+JA+NY)/2
100   JAM = MAX0(JA-1,1)
      JBM = MIN0(JB+1,NY)
      LYM = JBM-JAM
      YY1 = Y1 + (JA-1)*DY
      YYL = Y1 + (JB-1)*DY
      JJM = 1+(JA-JAM)*NDIV
      NYY = (JB-JA)*NDIV+1
C
C
      DO 120 I = IAM,IBM
      II =  I-IAM+1
      DO 120 J = JAM,JBM
      JJ =  J-JAM+1
120   ZZ(II,JJ) = Z(I,J)
C
C
      DO 140 K=1,10
      NDIVK = 2**(K-1)
      IF(NDIVK-NDIV)135,150,150
135   CALL DOUBLE(ZZ,LXM*NDIVK+1,LYM*NDIVK+1)
140   CONTINUE
C
C
150   CONTINUE
      CALL CONTUR(ZZ(IIM,JJM), NXX, NYY, XX1, YY1, XXL, YYL,
     & ZLEV,NDECL, LWGTL, NLEV, HGT, NARC, X, Y, S)
C
C
      JA=JB
      IF(JA-NY)80,160,160
160   IA=IB
      IF(IA-NX)40,170,170
170   RETURN
      END
