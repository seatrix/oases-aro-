      SUBROUTINE ZIN(Z,NX,NPYRD,LU)
      DIMENSION Z(NX,NPYRD)
      DO 1000    J=1,NPYRD
      READ(LU) (Z(I,J),I=1,NX)
 1000 CONTINUE
      RETURN
      END
