C
C     ARRAYS FOR DEPTH-DEPENDENT GREEN'S FUNCTIONS
C     AND WORKING ARRAYS
C
      COMPLEX CFF(NP,NPAR)
      COMPLEX CFFS(NP)
      COMMON /STRDIS/ CFF,CFFS
      COMPLEX CBUF(NP),CFILE(ISIZE)
      REAL ARG(NP),FAC(NP)
      COMMON /BUFF1/ ARG,FAC
      COMMON /BUFF2/ CBUF
      COMMON /BUFF3/ CFILE
