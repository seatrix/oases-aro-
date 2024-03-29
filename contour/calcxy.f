      SUBROUTINE CALCXY(Z,ZP,XP,YP,NPX,NPY)
C
      COMMON /XAX/X1,XL,XLEFT,XRIGHT,XSCALE,XINC,DX,
     % X1PL,XLPL,NX,X1GRID,XLGRID,DIVX,XVAL(100),NXVAL
      COMMON /YAX/Y1,YL,YUP,YDOWN,YSCALE,YINC,DY,
     % Y1PL,YLPL,NY,Y1GRID,YLGRID,DIVY,YVAL(100),NYVAL
      COMMON /PARA/LABPT,NSM,NDIV,CAY,NARC,NRNG,HGTPT,HGT,
     % LABC(51),LWGT(51)
      DIMENSION Z(NPX,NPY),ZP(NPX,NPY),XP(NPX,NPY),YP(NPX,NPY)
C
      XSTEP=(XL-X1)/(NPX-1)
      YSTEP=(Y1-YL)/(NPY-1)
      DO 2000 IY=1,NPY
      ANG=(IY-1)*YSTEP+YL
      CC=COS(ANG)
      CS=SIN(ANG)
      DO 1000 IX=1,NPX
      RR=(IX-1)*XSTEP+X1
      XP(IX,IY)=RR*CC
      YP(IX,IY)=RR*CS
      ZP(IX,IY)=Z(IX,IY)
 1000 CONTINUE
 2000 CONTINUE
C
      NRNGX=ABS(RR*YSTEP/(2.0*DX))
      NRNGY=ABS(RR*YSTEP/(2.0*DY))
      NRNG=MAX0(NRNGX,NRNGY,3)
      RETURN
      END
