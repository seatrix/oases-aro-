      SUBROUTINE XAXLOG(XAXIS,ISKIP,HGT)

      CHARACTER*3 XBTYPE
      CHARACTER*4  ICHAR, WORD
      CHARACTER*80 TITLEX


      COMMON /XAX/ X1,XL,XLEFT,XRIGHT,XSCALE,XINC,DX,
     % XDUM3,XDUM4,IXDUM1,XDUM5,XDUM6,DIVX,XVAL(100),NXVAL
      COMMON /XAXC/ TITLEX, XBTYPE
C
      AL2=ALOG(2.0)
      ALA=ALOG(XLEFT)/AL2
      ALB=ALOG(XRIGHT)/AL2
      NXVAL=(ALB-ALA)*XINC + 1

C
      
      DO 1000   I = 1, NXVAL
      XVAL(I)= XLEFT*2.0**((I-1)/XINC)
 1000 CONTINUE

      DO 1200   I=1,NXVAL   
      IF( MOD(XVAL(I),1.0) .NE. 0.0 )   THEN
       NDEC= 1
       GO TO 1400
      END IF
 1200 CONTINUE
      NDEC= -1
 1400 CONTINUE

C
      DX=XAXIS/(ALB-ALA)
      CALL PLOT(0.,0.,3)

      DO   2000   I=1,NXVAL
      X= (ALOG(XVAL(I))/AL2 - ALA) * DX
      CALL PLOT(X,0.,2)
      CALL PLOT(X,-0.05,2)
      IF(MOD(I-1,ISKIP).NE.0)   GO TO 2000
      XN=ABS(XVAL(I))
      IF(XN.LT.10.0)   XN=1.
      R=ALOG10(XN)
C     NDEC=1
      N=ABS(R)+2+NDEC
      IF(XVAL(I).LT.0.0)   N=N+1
      X1=X-N*HGT/2.0
      Y=-0.3
      CALL NUMBER(X1,Y,HGT,XVAL(I),0.,NDEC)
 2000 CALL PLOT(X,0.,3)
      CALL PLOT(XAXIS,0.,2)
C
C  SEARCH FOR LAST NON BLANK CHARACTER OF TITLEX
C
      DO   3000   I=80,1,-1
      N=I
 3000 IF(TITLEX(I:I).NE.' ')   GO TO 4000
 4000 CONTINUE
      X=(XAXIS-N*1.5*HGT)/2.
      Y=-0.70
      CALL SYMBOL(X,Y,HGT*1.5,TITLEX,0.,N)
      RETURN
      END
