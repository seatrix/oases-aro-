      SUBROUTINE CHKFILE(ISTAT)
      INCLUDE 'compar.f'
      INCLUDE 'compul.f'
*     Check input FILENAME
	
      CHARACTER*80 DUMMY
      CHARACTER*45 FILENM
      CHARACTER*45 FILE

      FILE=' '
      FILENM=' '
      CALL GETTXT(1, FILENM)
      FILENAME=FILENM
*     Check if info-file exists
      CALL READHEAD(ISTAT)
      CLOSE(15,STATUS='KEEP')
      IF (ISTAT.EQ.1) GO TO 10
      IF (ISTAT.EQ.2) GO TO 15

      IF (FILENM.NE.NEWNAME) THEN
        LX=LXTRF
        MX=MXTRF
        FREQS=FCTRF
        DLFRQP=1.0/REAL(DT*NX)
        FMIN=(LX-1)*DLFRQP
        FMINTRF=FMIN
        FMAX=(MX-1)*DLFRQP
        FMAXTRF=FMAX
        NFTRF=MX-LX+1

        CALL PUTFLT(4,FMIN)
        CALL PUTFLT(5,FMAX)
        CALL PUTFLT(6,FREQS)

        NEWNAME=FILENM

      END IF

C     WRITE(6,*)' CHECKFILE:FREQS,TMP:',FREQS,TMP

      ISTAT=0
      RETURN

10    WRITE(6,*) '>>> FILE DOES NOT EXIST <<<'
      PAUSE
      ISTAT=1
      RETURN

15    WRITE(6,*) '>>> FILE NOT IN CORRECT FORMAT <<<'
      PAUSE
      ISTAT=2
      RETURN
      END
      SUBROUTINE CHKFMIN(ISTAT)
      INCLUDE 'compar.f'
      INCLUDE 'compul.f'
      CALL GETFLT(4,FMIN)
      FMIN=AMAX1(FMINTRF,NINT(FMIN/DLFRQP)*DLFRQP)
      CALL PUTFLT(4,FMIN)
      ISTAT=0
      RETURN
      END
      SUBROUTINE CHKFMAX(ISTAT)
      INCLUDE 'compar.f'
      INCLUDE 'compul.f'
      CALL GETFLT(5,FMAX)
      FMAX=AMIN1(FMAXTRF,NINT(FMAX/DLFRQP)*DLFRQP)
      CALL PUTFLT(5,FMAX)
      ISTAT=0
      RETURN
      END
