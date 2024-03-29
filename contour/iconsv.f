      SUBROUTINE ICONSV(XP,YP,ZP,NZX,NZY,NT)
      PARAMETER (NLEV1=51)
      CHARACTER*3 XBTYPE, YBTYPE
      CHARACTER*80 TITLE,TITLEX,TITLEY
      CHARACTER*80 FILENM
      CHARACTER*40 WORD
      DIMENSION SECTOR(28),XP(1),YP(1),ZP(1)
      COMMON /PARA/LABPT,NSM,NDIV,CAY,NARC,NRNG,HGTPT,HGT,
     &             LABC(51),LWGT(51)
      COMMON /PARAC/ TITLE,  SDPLOT, RDPLOT
      COMMON /XAX/X1,XL,XLEFT,XRIGHT,XSCALE,XINC,DX,
     &            X1PL,XLPL,NX,X1GRID,XLGRID,DIVX,XVAL(100),NXVAL
      COMMON /XAXC/TITLEX, XBTYPE
      COMMON /YAX/Y1,YL,YUP,YDOWN,YSCALE,YINC,DY,
     &            Y1PL,YLPL,NY,Y1GRID,YLGRID,DIVY,YVAL(100),NYVAL
      COMMON /YAXC/TITLEY, YBTYPE
      COMMON /ZAX/ZMIN,ZMAX,ZINC,NLEV,ZLEV(NLEV1)
C
  200 FORMAT(A80)
  360 FORMAT(A40)
  400 FORMAT(F15.4)
  500 FORMAT(1X,' WARNING : ACCEPTABLE VALUES FOR NDIV ARE 1,2 AND 4 ',/,
     & '  ACTUAL VALUE IS ',I3,'. PLOTTING IS DONE WITH NDIV = 1 ')
  600 FORMAT(1X,'NX,NY TOO LARGE',/,1X,'EXECUTION',
     & 1X,'TERMINATED BECAUSE ARRAY SIZE LIMITATIONS',/)
  700 FORMAT(1X,'NLEV = ',' THIS PROGRAM ALLOWS ONLY 51 LEVELS FOR',
     & ' CONTOURING',/,'EXECUTION TERMINATED')
C
      XBTYPE='LIN'
      YBTYPE='LIN'
C
      READ(55,200)TITLE
      READ(55,200)FILENM
       CALL FILETYPE('   ',FILENM,55,17)
      READ(55,200)TITLEX
      READ(55,400)X1
      READ(55,400)XL
      READ(55,400)XLEFT
      READ(55,400)XRIGHT
      READ(55,360) WORD
      CALL AXLEN(XBTYPE,XLEFT,XRIGHT,WORD,XSCALE,XLEN,55)
C      READ(55,400)XSCALE
      READ(55,400)XINC
      READ(55,200)TITLEY
      READ(55,400)YUP
      READ(55,400)YDOWN
      READ(55,360) WORD
      CALL AXLEN(YBTYPE,YUP,YDOWN,WORD,YSCALE,YLEN,55)
C      READ(55,400)YSCALE
      READ(55,400)YINC
      Y1=YDOWN
      YL=YUP
      X1GRID=AMAX1(X1,XLEFT)
      XLGRID=AMIN1(XL,XRIGHT)
      Y1GRID=Y1
      YLGRID=YL
      READ(55,400)DUMMY
      NSVTOT=IFIX(DUMMY)
      READ(55,400)DUMMY
      READ(55,400)DIVX
      READ(55,400)DIVY
      READ(55,400)FLAGRC
      READ(55,400)DUMMY
      READ(55,400)DUMMY
      READ(55,400)DUMMY
      READ(55,400)DUMMY
      NX=IFIX(DUMMY)
      READ(55,400)DUMMY
      NY=IFIX(DUMMY)
      IF(NX*NY.LE.NZX*NZY)   GO TO 1000
      WRITE(6,600)
      STOP
 1000 CONTINUE
      READ(55,400)DUMMY
      I1OLD=0
      NX10=NX/10
      NXIND=MIN0(57,NX)
      NXKAB=MAX0(7,NX10+2)
      READ(55,400)DUMMY
      READ(55,400)CAY
      READ(55,400)DUMMY
      NRNG=IFIX(DUMMY)
C SECTION TO DETERMINE THE CONTOURLEVELS (ARRAY ZLEV).
      READ(55,400)ZMIN
      READ(55,400)ZMAX
      READ(55,400)ZINC
C TYPE 4 INFORMATION
      READ(55,400)X1PL
      READ(55,400)DUMMY
      READ(55,400)Y1PL
      READ(55,400)DUMMY
      NSM=IFIX(DUMMY)
      READ(55,400)HGTPT
      READ(55,400)HGT
      READ(55,400)DUMMY
      LABPT=IFIX(DUMMY)
      READ(55,400)DUMMY
      NDIV=IFIX(DUMMY)
      READ(55,400)DUMMY
      NARC=IFIX(DUMMY)
      READ(55,400)DUMMY
      LABC(1)=IFIX(DUMMY)
      READ(55,400)DUMMY
      LWGT1=IFIX(DUMMY)
      IF(ABS(ZINC).GT.0.0)   GO TO 1200
      NLEV=1
      ZLEV(1)=ZMIN
      GO TO 1800
 1200 CONTINUE
      NLEV=ABS((ZMAX-ZMIN)/ZINC)+1
      IF(NLEV.LE.NLEV1)   GO TO 1400
      WRITE(6,700)NLEV
      STOP
 1400 CONTINUE
      DO 1600 I=1,NLEV
      LABC(I)=LABC(1)
      ZLEV(I)=(I-1)*ZINC+ZMIN
      LWGT(I)=LWGT1
      IF(MOD(IFIX(ZLEV(I)+0.5),10).EQ.0 )LWGT(I)=LWGT1+1
 1600 CONTINUE
 1800 CONTINUE
      IF( (NDIV.NE.1) .AND.
     &    (NDIV.NE.2) .AND.
     &    (NDIV.NE.4) )     THEN
       WRITE(6,500) NDIV
       NDIV=1
      END IF
      READ(17,444)SECTOR
  444 FORMAT(5E15.4)
      NSVT=IFIX(SECTOR(1))
      IF(NSVT.NE.NSVTOT) THEN
        write(6,*) ' WRONG NUMBER OF PROFILES '
        STOP
      END IF
C
      NT=0
 2000 CONTINUE
      READ(17,*,END=2100)RNG,NP
      READ(17,*)(YP(J),ZP(J),J=NT+1,NT+NP)
      DO 2200 J=NT+1,NT+NP
      XP(J)=RNG
 2200 CONTINUE
      NT=NT+NP
      GO TO 2000
 2100 CONTINUE
C
      RETURN
      END
