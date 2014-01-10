      PROGRAM RDOASTL
c ********************************************************
c *                       OASES                          *
c *  Ocean Acoustic and Seismic Exploration Synthetics   *
c *                   Copyright (C)                      *
c *                  Henrik Schmidt                      *
c *       Massachusetts Institute of Technology          *
c *               Cambridge, MA 02139                    *
c ********************************************************
c
C     UNIX-FORTRAN - RANGE DEPENDENT VERSION
C     VERSION 2.2, UPDATE 6-July-1999.     
C          
      parameter (maxsect = 1000)
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comnrd.f'
      INCLUDE 'comfip.f'
C
      COMPLEX SLOW
      CHARACTER*50 FILENM
      DIMENSION CONMAX(npar)
      DIMENSION X(NP2,NPAR),FF(2,NP3),PX(MODULO)              
      DIMENSION RDUP(npar),RDDOWN(npar),CYAXIS(npar),RDINC(npar)
      DIMENSION YUP(npar),YDOWN(npar),YAXIS(npar),YINC(npar)
      DIMENSION ZMIN(npar),ZMAX(npar),ZSTEP(npar)
c >>> arrays for saving source parameters
      real vold(nla,6), sdcsav(nrd)
      integer laytold(nla),laysold(nrd)
c >>> Source strength array for backscattering
      complex uh(nrd),uv(nrd),uhb(nrd),uvb(nrd)
c >>> Sector variables
      real sectl,sleft,sright,slf(maxsect),dlr(maxsect)
      integer lfs(maxsect),indrng(maxsect)
      integer irtl(nrd)
      COMPLEX CTRF(npar)
      DIMENSION FFS(2,NP),XS(NP2),AKM(100),AKM1(100),AKM2(100)     
      CHARACTER*16 TITLEY(4)
      CHARACTER*6  OPTION(2)
      CHARACTER*80 TITLE
      character*80 atitle,ctitle,stitle
      character*40 botbuf
      character*4 pident(npar)
      LOGICAL ICONTU,AUSAMP,cylgeo,backscat,fftint,fftsec,
     &        v_field,con_int
      common /refnum/ csref,isref
      EQUIVALENCE (NREC,ISPACE),(LF,NUMFR)
      EQUIVALENCE (FF(1,1),CFF(1,1)),(FFS(1,1),CFFS(1))
      EQUIVALENCE (X(1,1),CFF(1,1)),(XS(1),CFFS(1))        
      equivalence (atitle,title)
      DATA CONMAX /npar*-1E30/
      DATA OPTION /'RDOAS ','      '/
      data pident /' Szz',' w  ',' u  ',' v  ',
     &             ' Sxx',' P  ',' S  '/                     
C          
C ********************FORMATS*******************         
C          
C          
 200  FORMAT(20A4)                
 209  FORMAT(1H ,'**********************************',
     &      /1H ,'*       OASES Version ',f2.2,'       *',
     &      /1H ,'**********************************',
     &     //1H ,'Range-Dependent TL Module RDOAST' )
 210  FORMAT(//1H ,A)         
350   FORMAT(//1H ,'    DEPTH        ALPHA       BETA      ATTENA       AT        
     1TENB         RHO       ROUGHNESS'//(1H ,3F12.5,2F12.8,2F12.5))  
500   FORMAT(//1H ,'Minimum Frequency:    ',F10.2,' HZ ',
     &        /1H ,'Maximum Frequency:    ',F10.2,' HZ ',
     &        /1H ,'Number of Frequencies:',I7 )       
550   FORMAT(//1H ,3X,'NW',5X,'  IC1',5X,'  IC2',/(1X,I5,2I10))                   
  600 FORMAT(//,'  CMIN = ',G15.6,' M/S ',/,'  CMAX = ',G15.6,' M/S ')          
C          
C          
C **********************************************         
C          
c >>> default real axis contour
      con_int=.false.
      DEBUG=.FALSE. 
      DECOMP=.FALSE.
      AUSAMP=.FALSE.
      LUGRN=30
c
c Direct access file record length
c
      call dasrec(ldaun)
C
C     CHECK THAT NP IS INTEGER POWER OF 2
C
      IPOW=0
 992  IPOW=IPOW+1
      II=2**IPOW
      IF (II.LT.NP) GO TO 992
      IF (II.NE.NP) THEN
        WRITE(6,*)
        WRITE(6,*) 
     &   '>>>> FATAL: PARAMETER NP MUST BE INTEGER POWER OF 2 <<<<'
        WRITE(6,*)
     &   '>>>>        CHANGE IN FILE compar.f AND RELINK.     <<<<'      
      END IF
C          
      rdoas=.true.
      MODU=MODULO
      ITXX=NP2
      PI=4E0*ATAN(1E0)           
      AI=CMPLX(0.,1.)             
      CNUL=CMPLX(0.,0.)
      IR=1 
      LS=1
      DELTA=1.
      THETA=0.
      FOCDEP=0.
      LTYP=1
      LINA=0
      NFLAG=.FALSE.
      ICNTIN=0
c >>> force plane geometry
      icdr=1
C
C     DEFAULT: ISOTROPIC LAYERS
C
      NTISOL=0
C          
      CALL OPFILR(1,IOER)
      IF (IOER.NE.0) STOP '>>>> ERROR: .dat FILE NOT FOUND <<<<' 
      READ(1,'(A)')TITLE            
      write(6,209) version
      WRITE(6,210)TITLE           
      CALL GETOPT(IPROF,ICNTIN,ICONTU,cylgeo,backscat,fftint,
     &            v_field,con_int)
      insrct=srctyp
c >>> Multible frequencies introduced 920910
      if (doppler) then
        READ(1,*) FREQ1,FREQ2,NFREQ,offdb,vrec
        vsou=vrec
      else IF (ICNTIN.GT.0) THEN
        READ(1,*) FREQ1,FREQ2,NFREQ,offdb
      ELSE
        OFFDB=0E0
        READ(1,*) FREQ1,FREQ2,NFREQ
      END IF
      IF (FREQ1*FREQ2.EQ.0.0) THEN
        STOP '*** FREQUENCIES MUST BE NON-ZERO, ABORTING ***'
      ELSE IF (FRCONT) THEN
        IF (NFREQ.LE.1) THEN
          STOP '*** CONTOURS REQUIRE NRFR>1 YOU STUPID FOOL ***'
        END IF
        F1LOG=LOG(FREQ1)
        F2LOG=LOG(FREQ2)
        DFLOG=(F2LOG-F1LOG)/(NFREQ-1)
      ELSE
      END IF
      IF (NFREQ.GT.1) THEN
       DLFREQ=(FREQ2-FREQ1)/(NFREQ-1)
      Else
       DLFREQ=1.
      END IF
      WRITE(6,500)FREQ1,freq2,nfreq            
c >>> create file for scattering rhs
      if (sctout) then
       call opfilb(45,ioer)
       write(45) nfreq,freq1,freq2,1e0/(nfreq*dlfreq)
      end if
c >>> file with virtual source strengths
      if (v_field) then
       call opfilw(16,ioer)
      end if

c >>> open file for sector environmental parameters
      call opnenv()
C           
C     READ IN ENVIRONMENTAL DATA
C
      if (iprof.ne.0.or.plkern.or.anspec.or.icontu) then
       read(1,'(2I10)',err=771) nsect,incsec
       go to 772
 771   write(6,*) 'Input file error'
       write(6,*) 'You probably tried to run an old data file'
       write(6,*) 'The # of sectors should now be accompagnied'
       write(6,*) 'by the output increment for options Z, I, a, c.'
       stop
 772   continue
       if (incsec.lt.1) go to 771
      else
       read(1,*) nsect
       incsec=1
      end if
c
c >>> Make bottom file for drcont
c
      if (drcont) then
       call opfilw(22,ioer)
       write(22,'(a)') 'BOTTOM         1'
      end if

      rlft=0.0
      do isect=1,nsect
       read(1,*) numl,sectl
       backspace(1)
       CALL INENVI
       call savenv()
       if (drcont) then
        rght=rlft+sectl*1000.0
        dbottom=fndbot()
        write(22,'(2F12.3)') rlft,dbottom
        write(22,'(2F12.3)') rght,dbottom
        rlft=rght
       end if
      end do      
c >>> initialize first sector env. variables
       call rwdenv()
       call rclenv()
       call rwdenv()
C          
C     SOURCE AND RECEIVER DATA
C
      CALL INSRC(SD)
      isref=lays((ls-1)/2+1)
      csref=v(isref,2)
      inls=ls
      do is=1,inls
       sdcsav(is)=sdc(is)
      end do
      CALL INREC(RD,RDLOW,INTF)
c
c >>> input receiver depths for TL vs range plots
      read(1,*) rdtlup,rdtllo,nrtl
      do ii=1,nrtl
       if (nrtl.gt.1) then
        rdt=rdtlup+(ii-1)* (rdtllo-rdtlup)/(nrtl-1)
       else
        rdt=rdtlup
       end if
       diff=1e20
       do jj=1,ir
        dfr=abs(rdt-rdc(jj))
        if (dfr.lt.diff) then
         irtl(ii)=jj
         diff=dfr
        end if
       end do
      end do
       
      NUMI=NUML-1                 
C          
C          
      READ(1,*)CMININ,CMAXIN
      IF (CMININ.EQ.0) STOP '*** CMIN MUST BE NON-ZERO ***'
      IF (((CMAXIN.GT.0).AND.(CMININ.GT.0)).OR.
     1    ((CMAXIN.LT.0).AND.(CMININ.LT.0))) THEN
       NFLAG=.FALSE.
      ELSE
       IF (CMININ.LE.0) STOP '*** CMIN/CMAX CONFLICT ***'
       CMAXIN=1E8
       NFLAG=.TRUE.
      END IF

c >>> Force full spectrum in RDOAST
c       CMAXIN=1E8
c       NFLAG=.TRUE.

      READ(1,*)NWVNOin,ICW1,ICW2       
C
C     RANGE DATA
C
      READ(1,*) XLEFT,XRIGHT,XAXIS,XINC
c 
c >>> plot data
c
      IF ((DEPTAV.and.(.not.frcont)).OR.PLTL.OR.PLKERN
     &    .OR.ANSPEC.OR.TLDEP) THEN
        DO 980 I=1,NPAR
        IF (IOUT(I).EQ.0) GO TO 980
        READ(1,*) YUP(I),YDOWN(I),YAXIS(I),YINC(I)
        YIAXIS=YAXIS(I)
 980    CONTINUE
      END IF
      IF (DRCONT.OR.TLDEP.OR.ICONTU) THEN
C
C     DEPTH AXES
C
        READ(1,*) RDUP(1),RDDOWN(1),CYAXIS(1),RDINC(1)
      END IF
      
      if (drcont) then
       write(22,'(2F12.3)') xright*1000.,dbottom
       write(22,'(2F12.3)') xright*1000.,rddown(1)
       write(22,'(2F12.3)') xleft*1000.,rddown(1)
       endfile(22)
      end if
C
C     DETERMINE MAXIMUM NUMBER OF TLDEP PLOTS
C
      NTLDEP=INT(ABS((XRIGHT-XLEFT)/XINC))+1
C        
      IF (DRCONT.or.frcont) THEN
        XSCALE=ABS(XRIGHT-XLEFT)/XAXIS
        if (drcont) YSCALE=ABS(RDUP(1)-RDDOWN(1))/CYAXIS(1)
        DO 990 I=1,npar
        IF (IOUT(I).EQ.0) GO TO 990
        READ(1,*) ZMIN(I),ZMAX(I),ZSTEP(I)
 990    CONTINUE
      END IF
c
c >>> open scratch file for frequency-range contours
c
      if (frcont) then
       open(unit=72,FORM='UNFORMATTED',status='scratch')
       open(unit=73,FORM='UNFORMATTED',status='scratch')
      end if
C          
C     OPEN PLOT FILES
C
      CALL OPFILW(19,IOER)
      CALL OPFILW(20,IOER)
      WRITE(19,6010) MODU,'MODU'
C
C     OPEN CHECK FILE
C
      IF (DEBUG) CALL OPFILW(21,IOER)
C
C     PLOT OF VELOCITY PROFILE IF OPTION 'Z' WAS CHOSEN
C
      IF (IPROF.GT.0) THEN
        READ(1,*) VPLEFT,VPRIGHT,VPLEN,VPINC
        READ(1,*) DPUP,DPLO,DPLEN,DPINC
        rewind(81)
        ranst=0e0
        do isct=1,nsect
         call rclenv()
         if (mod(isct-1,incsec).eq.0.or.isct.eq.nsect) then
          ltit=lenstr(atitle)
          write(stitle,'(a,a,f9.3)') atitle(1:ltit),' - R =',ranst
          CALL PLPROF(STITLE,VPLEN,DPLEN,VPLEFT,VPRIGHT,VPINC,
     2              DPUP,DPLO,DPINC)
         end if
         ranst=ranst+sectl
        end do
      END IF
 6010 FORMAT(1H ,I8,10X,A40)
c >>> initialize first sector
       call rwdenv()
       call rclenv()
       call rwdenv()
c
c >>> open contour files
c
      IF (DRCONT.or.FRCONT.OR.ICONTU) THEN
        CALL OPFILW(28,IOER)
        CALL OPFILW(29,IOER)
      END IF
c
c >>> Initialize
c
      CALL PINIT1
      if (DEBUG) CALL PREQV(NUML,NUMI)
C
C     SLOWNESS DIAGRAM FOR TRANS. ISOTR. MEDIA
C
      IF (NTISOL.GT.0.and.IPROF.GT.0) THEN
       CALL SNSDGM(title)
      END IF

      IF (TABLERC) THEN
C >>>  Tabulated Top  Reflection Coefficient
       call settrc()
      ELSE IF (BOTTOMRC) THEN
C >>>  Tabulated Bottom  Reflection Coefficient
       call setbrc()
      END IF

c >>> open kernel buffer file for coupling
      open(18,status='scratch',form='unformatted')

      if (backscat) then
       ndir=2
c >>> open direct access file for backscattering sources
       open(9,status='scratch',access='direct',form='unformatted',
     &       recl=2*ldaun*ir)
      else
       ndir=1
      end if
c >>>
c >>> Frequency loop
c >>>
      DO 15 JJ=1,NFREQ
       srctyp=insrct
       ftime=0e0
       CALL CLTIME
       IF (FRCONT) THEN
        FREQ=EXP(F1LOG+(JJ-1)*DFLOG)
       ELSE
        FREQ=FREQ1+(JJ-1)*DLFREQ
       END IF
       write(6,309) freq
 309   FORMAT(//1H ,'Frequency: ',F10.2,' Hz')

       DSQ=2E0*PI*FREQ             
       CSQ=DSQ*DSQ                 

       IF (TABLERC) THEN
c >>>   Read Top reflection coefficient table
        call gettrc()
       ELSE IF (BOTTOMRC) THEN
c >>>   Read Bottom reflection coefficient table
        call getbrc()
       ENDIF
c
c >>> forward/backward loop
c
      do idir=1,ndir
C     
C     OPEN SCRATCH FILES FOR unassembled transfer functions
C
       do i=1,npar
        IF (IOUT(I).NE.0) THEN
          LOGNUM=30+I
          open(unit=lognum,status='scratch',form='unformatted')
        END IF
       end do

       if (idir.eq.1) then
        ncsect=0
        sleft=0
        isfst=1
        islst=nsect
        istp=1
       else
        isfst=ncsect-1
        islst=1
        istp=-1
       end if
c 
c >>> start of sector loop
c
      rewind(81)
c >>> default only symmetric sources
      nsym=1
      do isect=isfst,islst,istp
       if (idir.eq.1) then
c >>> if not first sector, save environmental param
        if (isect.ne.isfst) then
         srctyp=98
         do ila=1,numl
          do j=1,6
           vold(ila,j)=v(ila,j)
          end do
          laytold(ila)=laytyp(ila)
         end do
         do ird=1,ir
          laysold(ird)=lay(ird)
         end do
         rangs=1e3*sectl
        end if

        if (sleft.lt.xright) then
         ncsect=ncsect+1
c >>> read in environmental parameters for sector
         call rclenv()
         if (isect.eq.nsect) then
          sectl=max(sectl,xright-sleft+rstep)
         end if
         sright=min(sleft+sectl,xright+rstep)
c >>> Initialize receiver array
         call recarr()
c >>> If not first sector, set source array
         if (isect.ne.isfst) then
          call souarr()
          deltz=abs(rdc(2)-rdc(1))
c >>> insert here computation of source strengths
          rewind(18)
          call csarr(nsym,deltz,rangs,vold,laytold,laysold,
     &               uh,uv,uhb,uvb)
c >>>>  NOTE change sign of uv. No explanation yet. 960305
c           call vneg(uv,1,uv,1,2*ir)
c           call vneg(uvb,1,uvb,1,2*ir)
c <<<<<
          if (backscat) then
           write(9,rec=isect-1) (uhb(jr),jr=1,ir)
           write(9,rec=isect-1+(nsect-1)*ir) (uvb(jr),jr=1,ir)
          end if
          if (v_field) then
           write(16,*)
           write(16,'(i10,a)') idir,  '   # Prop. dir.: 1 forw; 2 backw'
           write(16,'(i10,a)') isect,  '   # Sector #'
           write(16,'(f10.3,a)') sleft,'   # Range in km '
           write(16,'(a10,4a14)') ' Depth (m)',
     &                          '     Re(u)','     Im(u)',
     &                          '     Re(w)','     Im(w)'
           ffx=2*pi/deltz
           do jr=1,ir
            write(16,'(f10.3,4G14.6)') rdc(jr),ffx*uh(jr),ffx*uv(jr)
           end do
          end if
         else
c 
c >>> set physical source
c
          srctyp=insrct
          ls=inls
          do is=1,inls
           sdc(is)=sdcsav(is)
           CALL SOURCE(V,NUML,SDC(is),LAYS(is),ZUS(is),ZLS(is))
          end do
         end if
c >>> re-initialize for new sector
         call pinit1()
        end if
       else
c >>> backward propagation
c >>> read in environmental parameters for sector
        if (isect.ne.isfst) then
         do ila=1,numl
          do j=1,6
           vold(ila,j)=v(ila,j)
          end do
          laytold(ila)=laytyp(ila)
         end do
         do ird=1,ir
          laysold(ird)=lay(ird)
         end do
         rangs=1e3*sectl
        end if
c
c >>> read environmental data
c
        rewind(81)
        sleft=0
        sectl=0
        do isct=1,isect
         sleft=sleft+sectl
         call rclenv()
         sright=sleft+sectl
        end do
c >>> Initialize receiver array
        call recarr()
c >>> If not first sector, set source array
        call souarr()
        deltz=(rdc(2)-rdc(1))
c >>> read source strengths
c >>> computation of source strengths
        if (isect.ne.isfst) then
         rewind(18)
         call csarr(nsym,deltz,rangs,vold,laytold,laysold,
     &              uh,uv,uhb,uvb)
         read(9,rec=isect) (cbuf(jr),jr=1,ir)
         call vadd(uh,1,cbuf,1,uh,1,2*ir)
         read(9,rec=isect+(nsect-1)*ir) (cbuf(jr),jr=1,ir)
         call vadd(uv,1,cbuf,1,uv,1,2*ir)
        else
         read(9,rec=isect) (uh(jr),jr=1,ir)
         read(9,rec=isect+(nsect-1)*ir) (uv(jr),jr=1,ir)
        end if
        if (v_field) then
         write(16,*)
         write(16,'(i10,a)') idir,  '   # Prop. dir.: 1 forw; 2 backw'
         write(16,'(i10,a)') isect,  '   # Sector #'
         write(16,'(f10.3,a)') sleft,'   # Range in km '
         write(16,'(a10,4a14)') ' Depth (m)',
     &                          '     Re(u)','     Im(u)',
     &                          '     Re(w)','     Im(w)'
         ffx=2*pi/deltz
         do jr=1,ir
          write(16,'(f10.3,4G14.6)') rdc(jr),ffx*uh(jr),ffx*uv(jr)
         end do
        end if
c >>> re-initialize for new sector
        call pinit1()
       end if
c        if (debug) pause

       if (sleft.lt.xright) then
c 
c >>> make tables for dispersive media
c
      if (idlmax.gt.0) then
       call dltable(2,2,freq)
      end if 

      nactf=2
      IF (NWVNOin.LT.0) THEN
C ***  AUTOMATIC SAMPLING
       AUSAMP=.TRUE.
C ***  FORCE COMPLEX CONTOUR
       ICNTIN=1
       RMAXA=ABS(SLEFT-SRIGHT)
       RMINA=RMAXA-ABS(XRIGHT-XLEFT)
       OFFDB=0E0
       CALL AUTSAM(CMININ,CMAXIN,RMINA,RMAXA,CMIN,CMAX,NWVNO,ICW1,ICW2)
      ELSE
       nwvno=nwvnoin
       CMIN=CMININ
       CMAX=CMAXIN
       NWVNO=MIN0(NWVNO,NP)
       ICW2=MIN0(NWVNO,ICW2)
       ICW1=MAX0(1,ICW1)
      END IF
c >>> Hanning tapering
      icut1=1
      icut2=nwvno
C *** FOR TL TAPERING BY CHERMIT 
c      ICUT1=ICW1
c      ICUT2=ICW2
C
C     CHECK THAT NWVNO IS INTEGER POWER OF 2
C
      IPOW=0
 991  IPOW=IPOW+1
      II=2**IPOW
      IF (II.LT.NWVNO) GO TO 991
      IF (II.NE.NWVNO) THEN
        NWVNO=II
        WRITE(6,*)
        WRITE(6,*) '>>>> NW MUST BE INTEGER POWER OF 2, CHANGED TO',
     &             NWVNO,' <<<<'      
      END IF
      IF (NWVNO.GT.NP) STOP '>>> TOO MANY WAVENUMBERS <<<'
      IF (.NOT.NFLAG) THEN
       WK0 = 2*PI*FREQ / CMAX
       WKMAX = 2*PI*FREQ / CMIN               
      ELSE
       WKMAX=2*PI*FREQ/CMIN
       dlwvno=wkmax/(nwvno-0.5)
       WK0=0.5*dlwvno
       CMAX=2*PI*FREQ/WK0
       ICUT1=1
       ICW1=1
      END IF
      WRITE(6,600)CMIN,CMAX       
      WRITE(6,550)NWVNO,ICW1,ICW2                      
      IF (NFLAG) WRITE(6,*) '*** NEGATIVE SPECTRUM BY SYMMETRY ***'

      DLWVNO = ( WKMAX - WK0 ) / ( FLOAT(NWVNO-1) )      
      DLRAN=2E0*PI/(NWVNO*DLWVNO)
      dlrfft=dlran 
      RSTEP=DLRAN*1.0E-3          
      if (idir.eq.1) then
c >>> find minimum speed
       sspmin=1e20
       do ii=1,numl
        if (v(ii,2).gt.100) then
         if (laytyp(ii).gt.2.and.v(ii,3).gt.0.1*v(ii,2)) then
          sspmin=min(sspmin,v(ii,3))
         else
          sspmin=min(sspmin,0.5*v(ii,2))
         end if
        end if
       end do
       r0=max(rstep,0.5*1E-3*sspmin/freq)
       r0=min(r0,(sright-sleft)*0.5)
       if (backscat.and.isect.lt.nsect) then              
        lf=nint((sright-sleft-2*r0)/rstep)+1        
        lf=max(lf,1)
       else
        lf=nint((sright-sleft-r0)/rstep)+2        
        lf=max(lf,1)
       end if
c >>> save sector variables
       slf(isect)=sleft+r0
       dlr(isect)=rstep
       lfs(isect)=lf
      else
       rstep=dlr(isect)
       dlran=rstep*1E3
       lf=lfs(isect)
       r0=sright-(slf(isect)+(lf-1)*dlr(isect))
      end if 
c >>> determine whether to force FFT integration
      fftsec=fftint.or.((lf*ir).gt.np)
c >>> Determine integration contour
      cont_45 = con_int .and. (.not.fftsec)
      R1=R0*1E3
      RANMAX=NWVNO*DLRFFT
       write(6,*)
       write(6,*) ' >>> Frequency    ',freq,' Hz <<<'
       write(6,*) ' >>> Sector no.   ',isect,' <<<'
       write(6,*) ' >>> Left border: ', sleft,' km <<<'
       write(6,*) ' >>> Right border:', sright,' km <<<'
       write(6,*) ' >>> Samples:     ',lf,' <<<'
       write(6,*) ' >>> Min range:   ',r0*1e3,' m  <<<'
       write(6,*) ' >>> Range step:  ',dlran,' m  <<<'

      WRITE(6,360) DLRFFT,RANMAX*1E-3
 360  FORMAT(1H ,' ',/1H ,'FFT RANGE STEP:   ',F12.3,' m',
     &               /1H ,'MAX FFT RANGE:',F12.3,' km')

C
C     IF OPTION 'J', SET DEFAULT OFFDB
C
      IF (ICNTIN.GT.0) THEN
       IF (OFFDB.LT.1E-10) THEN
        OFFDB=60.0*V(LAYS((LS-1)/2+1),2)/(FREQ*RANMAX)
        WRITE(6,*)
        WRITE(6,*) 'THE DEFAULT CONTOUR OFFSET IS APPLIED'
       ELSE
        WRITE(6,*)
        WRITE(6,*) 'THE USER DEFINED CONTOUR OFFSET IS APPLIED'
       END IF
       ATTMAX=OFFDB*FREQ*RANMAX/V(LAYS((LS-1)/2+1),2)
       WRITE(6,361) OFFDB,ATTMAX,ATTMAX*XRIGHT*1E3/RANMAX
 361   FORMAT(1H ,'CONTOUR OFFSET:         ',F12.6,' dB/wavelength',
     &       /1H ,'AT MAX FFT RANGE:       ',F12.2,' dB',
     &       /1H ,'AT MAX WINDOW RANGE:    ',F12.2,' dB',
     &       /1H ,'>> NOTE THAT COMPENSATION IS AUTOMATIC <<')
      END IF

C
C     OPEN SCRATCH FILE FOR INTEGRAND CONTOUR DATA
C
       IF (ICONTU) THEN
         NN=ICUT2-ICUT1+1
         NDEC=NN/NCONM
         IF (NDEC.GT.0) THEN
           NCON=(NN-1)/NDEC+1
         ELSE
           NDEC=1
           NCON=NN
         END IF
        CALL OPNBUF(27,NCON,IR*NOUT,100)
       END IF

c
c >>> Compute Green's functions
c

      CALL PINIT2 
c >>> nsym = 1 for only symmetric sources
c >>> nsym =2 for sym. and antisym. sources.
       nsym=1
      if (isect.ne.1.or.idir.gt.1) then
       srctyp=98
       do ila=1,nla
        if (laytyp(ila).gt.2) then
         nsym=2
        end if
       end do
      end if

      CALL CLTIME

      rewind(18)
c *** kernels
      CALL RD_CALIN(nsym,uh,uv)
c ***
      CALL RDTIME(T1)
      TOTTIM=TOTTIM+T1
      CALL CHKSOL
      WRITE(6,310) T1
 310  FORMAT(//1H ,'Integration kernels, CPU= ',F12.3,' SECS.')
C
C     WAVENUMBER INTEGRATION LOOP
C
      ftime=0e0
      if (fftsec) then
       write(6,*) '>>> Using FFT integration <<<'
      else if (cont_45) then
       write(6,*) '>>> 45 deg contour integration <<<'
      end if
c
      DO 20 NREC=1,IR
       CALL CLTIME
       if (fftsec.or.icontu) then
        CALL RD_PHINT(nsym)
       end if

       if ((plkern.or.anspec).and.
     &     (mod(isect-1,incsec).eq.0.or.isect.eq.nsect)) then
        do jjj=1,nrtl
         IF (NREC.EQ.irtl(jjj)) THEN
          if (.not.(icontu.or.fftsec)) then
           CALL RD_PHINT(nsym)
          end if
          do isym=1,nsym
Cs3
C     PLOT OF HANKEL OR FOURIER TRANSFORMS
C
           ltit=lenstr(atitle)
           write(stitle,'(a,a,i3,2(a,i2))') atitle(1:ltit),
     &                ' - Sect=',isect,' Dir=',idir,' Sym=',isym
           if (isym.eq.2) then
c >>> switch symmetric and anti-symmetric
            do i=1,npar
             if (iout(i).gt.0) then
              call vmov(cff(1,i),1,cffs(1),1,2*nwvno)
              call vmov(cff(nwvno+1,i),1,cff(1,i),1,2*nwvno)
              call vmov(cffs(1),1,cff(nwvno+1,i),1,2*nwvno)
             end if
            end do
           end if
           IF (PLKERN) THEN
            CALL PLINTGR(DLWVNO,WK0,SD,RDC(NREC),stitle,XAXIS,YIAXIS)
           END IF
C
C     PLOT OF ANGULAR SPECTRA
C
           IF (ANSPEC) THEN
            CALL PLSPECT(DLWVNO,WK0,SD,RDC(NREC),stitle,
     1                   XAXIS,YIAXIS,LAY(NREC))
           END IF
           if (isym.eq.2) then
c >>> switch symmetric and anti-symmetric
            do i=1,npar
             if (iout(i).gt.0) then
              call vmov(cff(1,i),1,cffs(1),1,2*nwvno)
              call vmov(cff(nwvno+1,i),1,cff(1,i),1,2*nwvno)
              call vmov(cffs(1),1,cff(nwvno+1,i),1,2*nwvno)
             end if
            end do
           end if
          end do
         END IF
        end do
       end if

C
C     DECIMATE INTEGRAND FOR DEPTH CONTOUR PLOT
C
       IF (ICONTU.and.
     &     (mod(isect-1,incsec).eq.0.or.isect.eq.nsect)) THEN
        NN=ICUT2-ICUT1+1
        do i=1,npar
         IF (IOUT(I).GT.0) THEN
          CALL CVMAGS(CFF(ICUT1,I),2,FFS,1,NN)
          CALL VCLIP(FFS,1,1E-20,1E20,FFS,1,NN)
          CALL VALG10(FFS,1,FFS,1,NN)
          CALL VSMUL(FFS,1,1E1,FFS,1,NN)
          CALL VDECIM(FFS,1,FFS,1,NN,NDEC,NC)
          CALL WRBUF(27,FFS,NC)
          CALL VMAX(FFS,1,AAA,NC)
          CONMAX(I)=MAX(CONMAX(I),AAA)
         END IF
        end do
       END IF

       if (fftsec) then
C
C     TRANSMISSION LOSS CALCULATION
C
        CALL RD_TRF(nsym)
C
C     WRITE TRANSFER FUNCTIONS TO SCRATCH FILES
C
        do i=1,npar
         if (iout(i).gt.0) then
          lognum=30+I
          if (mod(idir,2).eq.1) then
           write(lognum) (cff(j,I),j=1,LF)
          else
           if (i.eq.3.or.i.eq.7) then
            write(lognum) (-cff(j,I),j=LF,1,-1)
           else
            write(lognum) (cff(j,I),j=LF,1,-1)
c           write(6,*) 'cff=',cff(lf,i),cff(lf-9,i)
           end if
          end if
         end if
        end do
       end if
       CALL RDTIME(T2)
       ftime=ftime+T2
c      WRITE(6,320) NREC,T2
 320   FORMAT(1H ,'RECEIVER NO. ',I4,' TIME=',F8.3,' SECONDS')
 20   CONTINUE
  
      if (.not.fftsec) then
       ftime=0e0
       call cltime()
C
C     TRANSMISSION LOSS CALCULATION
C
       CALL rd_di_trf(nsym,lf,ir,r1,dlran,cbuf)
C
C     WRITE TRANSFER FUNCTIONS TO SCRATCH FILES
C
       do nrec=1,ir
        ioff=(nrec-1)*lf
        DO I=1,npar
         IF (IOUT(I).gt.0) then
          lognum=30+I
          if (mod(idir,2).eq.1) then
           write(lognum) (cff(j+ioff,I),j=1,LF)
          else
           if (i.eq.3.or.i.eq.7) then
            write(lognum) (-cff(j+ioff,I),j=LF,1,-1)
           else
            write(lognum) (cff(j+ioff,I),j=LF,1,-1)
c           write(6,*) 'cff=',cff(lf+ioff,i),cff(lf-9+ioff,i)
           end if
          end if
         end if
        end do
       end do
       CALL RDTIME(T2)
       ftime=ftime+T2
      end if

      write(6,'(1h ,a,f8.3,a)') 'Integration done, CPU=',ftime,' s'
      tottim=tottim+ftime

C*****  CONTOUR PLOTS OF INTEGRANDS VERSUS depth

      IF (ICONTU.and.
     &    (mod(isect-1,incsec).eq.0.or.isect.eq.nsect)) THEN
       CALL ENFBUF(27)
       ltit=lenstr(atitle)
       write(stitle,'(a,a,i3,a,i2)') atitle(1:ltit),
     &                ' - Sect=',isect,' Dir=',idir
        WN1=WK0+(ICUT1-1)*DLWVNO
        WN2=WK0+(NCON-1)*NDEC*DLWVNO
C
C     SLOWNESS AXIS
C
        WN1=WN1/(2E0*PI*FREQ)
        WN2=WN2/(2E0*PI*FREQ)      
        IF (WN1.LT.WN2*1E-1) THEN
          XLFTI=0
        ELSE
          XLFTI=INT(WN1*1E4)*1E-4
        END IF
        XRGHTI=WN2
        XINCI=INT(XRGHTI*1E4/5E0)*1E-4
        XSCALI=(XRGHTI-XLFTI)/20.0
        DO 710 I=1,npar
         IF (IOUT(I).LE.0) GO TO 710
          CALL RWDBUF(27)
          ZMINI=INT(CONMAX(I))
          ZINCI=10E0
          ZMAXI=ZMINI-ZINCI*10E0
          CALL INTCON(stitle,NCON,IR,NCON,IR,
     1                XLFTI,XRGHTI,XSCALI,XINCI,
     2                RDUP(1),RDDOWN(1),YSCALE,RDINC(1),
     3                ZMINI,ZMAXI,ZINCI,FREQ,SD,
     4                RD,RDLOW,WN1,WN2,I,1,0)
          DO 709 JR=1,IR
            DO 708 J=1,npar
              IF (IOUT(J).LE.0) GO TO 708
              CALL RDBUF(27,FFS,NCON)
              IF (J.EQ.I) THEN
                WRITE(29,'(1X,6G13.5)') (XS(JJJ),JJJ=1,NCON)
              END IF
 708        CONTINUE
 709      CONTINUE
 710    CONTINUE
       CALL CLSBUF(27)
      END IF

C
C     close and Rewind SCRATCH FILES
C
      CALL CLSBUF(30)


c
c >>> End of sector loop
c
       sleft=sright
       end if
       if (debug) then
         write(6,*) '>>> end of sector',isect
c         pause
       end if
      end do
c
c >>> make file with transfer functions
c
      if (idir.eq.1) then
       islst=ncsect
       lfsum=0
       indrng(1)=1
       do isect=1,ncsect
        lfsum=lfsum+lfs(isect)
        indrng(isect+1)=indrng(isect)+lfs(isect)
       end do
       XXL=slf(1)
       XXM=xright
       dxx=(xxm-xxl)/(lfsum-1)
       WDEPTH=0.0
       do i=1,npar
        if (iout(i).gt.0) then
         logtrf=30+npar+i
         open(logtrf,status='scratch',form='unformatted',
     &               access='direct',recl=2*ldaun*lfsum)
        end if
       end do
      end if       
c >>> write transfer functions to file
      do i=1,npar
       if (iout(i).gt.0) then
        do nrec=1,ir 
         call vclr(cff(1,i),1,2*lfsum)
         lognum=30+I
         rewind(lognum)
         do isect=isfst,islst,istp
          indran=indrng(isect)
c          write(6,*) 'nrec,isect,lf=',nrec,isect,lfs(isect)
          do irv=1,ir
           read(lognum) (cbuf(j),j=1,LFs(isect))
           if (irv.eq.nrec) then
c            write(6,*) 'cff=',cbuf(1),cbuf(10)
            call vmov(cbuf,1,cff(indran,i),1,2*lfs(isect))
           end if
          end do
         end do
         logtrf=30+npar+i
         irecord=nrec+(idir-1)*ir
         write(logtrf,rec=irecord)(cff(j,i),j=1,lfsum)
        end do
       end if
      end do

c >>> Close unassembled transfer function file

      DO I=npar,1,-1
       IF (IOUT(I).NE.0) THEN
        lognum=30+I
        close(lognum,status='delete')
       END IF
      end do
c 
c >>> end of forward/backward loop
c
      end do
c
c >>> Total field is idir=3 for backscattering case
      if (ndir.eq.1) then
       npdir=1
      else
       npdir=3
      end if

c >>> start of second forward/backward loop
      do idir=1,npdir

       ltit=lenstr(atitle)
       if (idir.eq.1) then
        write(6,*) '>>> Plots of forward propagated field <<<'
        write(ctitle,'(a,a)') atitle(1:ltit),' - Forward'
       else if (idir.eq.2) then
        write(6,*) '>>> Plots of backscattered field <<<'
        write(ctitle,'(a,a)') atitle(1:ltit),' - Backward'
       else
        write(6,*) '>>> Plots of total field <<<'
        write(ctitle,'(a,a)') atitle(1:ltit),' - Total'
       end if

       if (idir.eq.1) then
        isfst=1
        islst=ncsect
        istp=1
       else
        isfst=ncsect-1
        islst=1
        istp=-1
       end if
c
c >>> gridding of field magnitude
c 
       do i=1,npar
        if (iout(i).gt.0) then
         logtl=30+i
         open(unit=logtl,status='scratch',form='unformatted')
         call vclr(fac,1,lfsum)
         logtrf=30+npar+I
c Range array for computed transfer functions
         do nrec=1,ir 
          do isect=1,ncsect
           indran=indrng(isect)
           do irh=1,lfs(isect)
            xs(indran+irh-1)=slf(isect)+(irh-1)*dlr(isect)
            if (rdtest.and.nrec.eq.irtl(1)) then
              write(6,*) isect,irh,xs(indran+irh-1),x(indran+irh-1,i)
            end if
           end do
          end do
c >>> read transfer functions and compute magnitude squared
          if (idir.eq.3) then
           irecord=nrec
           read(logtrf,rec=irecord)(cbuf(j),j=1,lfsum)
           irecord=nrec+ir
           read(logtrf,rec=irecord)(cbuf(j+lfsum),j=1,lfsum)
           call vadd(cbuf(1),1,cbuf(1+lfsum),1,cbuf(1),1,2*lfsum)
          else
           irecord=nrec+(idir-1)*ir
           read(logtrf,rec=irecord)(cbuf(j),j=1,lfsum)
          end if
          call cvmags(cbuf,2,x(1,i),1,lfsum)
          if (cylgeo) then
           do irh=1,lfsum
            rangem=max(1e3*xs(irh),1e0)
            fcr=1e0/rangem
            x(irh,i)=fcr*x(irh,i)
           end do
          end if
c >>> now interpolate
          do 167 ic=1,lfsum
           xx=xxl+(ic-1)*dxx
           arg(ic)=0
           do ii=1,lfsum-1
            if (xs(ii).le.xx.and.xs(ii+1).gt.xx) then
c             arg(ic)=((xx-xs(ii))*x(ii+1,i)
c     &              + (xs(ii+1)-xx)*x(ii,i))/(xs(ii+1)-xs(ii))
             arg(ic)=x(ii,i)+ (x(ii+1,i)-x(ii,i))
     &                        *(xx-xs(ii))/(xs(ii+1)-xs(ii))
             go to 167
            end if
           end do
 167      continue
          write(logtl) (arg(ic),ic=1,lfsum)
          call vadd(arg(1),1,fac(1),1,fac(1),1,lfsum)
         end do
c >>> depth averaged loss
         dfact=1e0/ir
         call vsmul(fac(1),1,dfact,fac(1),1,lfsum)
         write(logtl)  (fac(ic),ic=1,lfsum)
        end if
       end do


       if (pltl.or.deptav) then
        do i=1,npar
         if (iout(i).gt.0) then
          logtl=30+i
          rewind(logtl)
          do nrec=1,ir 
           read(logtl) (arg(ic),ic=1,lfsum)
           do jjj=1,nrtl
            if (nrec.eq.irtl(jjj).and.pltl) then
c >>> convert to tl
             CALL VCLIP(arg,1,1E-20,1E20,arg,1,lfsum)
             call valg10(arg,1,arg,1,lfsum)
             call vsmul(arg,1,-10.0,x(1,i),1,lfsum)
             stitle=ctitle(1:lenstr(ctitle))//pident(i)             
             CALL PLTLOS(lfsum,xxl,dxx,stitle,I,XAXIS,YAXIS(I),
     1            XLEFT,XRIGHT,XINC,YUP(I),YDOWN(I),YINC(I),
     2            SD,RDC(NREC))
            end if
           end do
          end do
          if (deptav) then
           read(logtl) (arg(ic),ic=1,lfsum)
           if (.not.frcont) then
            CALL VCLIP(arg,1,1E-20,1E20,arg,1,lfsum)
            call valg10(arg,1,arg,1,lfsum)
            call vsmul(arg,1,-10.0,x(1,i),1,lfsum)             
            stitle=ctitle(1:lenstr(ctitle))//pident(i)             
            CALL PLDAV(lfsum,xxl,dxx,stitle,I,XAXIS,YAXIS(I),
     1           XLEFT,XRIGHT,XINC,YUP(I),YDOWN(I),YINC(I),SD)
           else
            dnew=(xright-xleft)/(nconm-1)
            lfn=nconm
            CALL Vfltsm(arg,Xs,lfsum,dxx,lfn,dnew)
            write(73) (xs(l),l=1,lfn)
           end if   
          end if
         end if
        end do
       end if
C
C     PLOT TRANSMISSION LOSS OVER DEPTH
C     INTERPOLATION PERFORMED BETWEEN DATA POINTS
C
      IF (TLDEP) THEN
       do i=1,npar
        if (iout(i).gt.0) then
         DO L=1,NTLDEP
          RTLDEP=XLEFT+(L-1)*XINC
          ITL=INT((RTLDEP-xxl)/dxx)+1
          IF (ITL.GT.0.AND.ITL.LT.lfsum) THEN
           RTL=xxl+(ITL-1)*dxx
           RAT=(RTLDEP-RTL)/dxx
           logtl=30+i
           rewind(logtl)
           DO JR=1,IR
            read(logtl) (XS(jrh),jrh=1,lfsum)
            X(JR,I)= XS(ITL)+RAT*(XS(ITL+1)-XS(ITL))
            x(jr,I)=-10.0*log10(x(jr,i))
           end do
           stitle=ctitle(1:lenstr(ctitle))//pident(i)             
           CALL PTLDEP(IR,RDC(1),RDSTEP,stitle,I,XAXIS,CYAXIS(1),
     1            YDOWN(I),YUP(I),YINC(I),RDUP(1),RDDOWN(1),RDINC(1),
     2            SD,RTLDEP)
          end if
         end do
        end if
       end do
      end if
C
C     GENERATE TL-CONTOUR PLOT FILES
C
      IF (DRCONT) THEN
        NDECc=LFsum/NCONM
        IF (NDECc.GT.0) THEN
          NCON=(LFsum-1)/NDECc+1
        ELSE
          NDECc=1
          NCON=LFsum
        END IF
        cXXL=slf(1)
        cXXM=xright
        cdxx=(cxxm-cxxl)/(ncon-1)
        WDEPTH=0.0

        DO I=1,npar
         IF (IOUT(I).NE.0) THEN
          stitle=ctitle(1:lenstr(ctitle))//pident(i)             
          CALL CONDRW(stitle,NCON,IR,NCON,IR,XLEFT,XRIGHT,XSCALE,XINC,
     1            RDUP(1),RDDOWN(1),YSCALE,RDINC(1),ZMIN(I),ZMAX(I),
     2            ZSTEP(I),FREQ,SD,RD,RDLOW,cXXL,cXXM,PX,icdr)
c
c >>> add bottom shading
c
          rewind(22)
 566      read(22,'(a)',end=567,err=567) botbuf
          write(28,'(a)') botbuf
          go to 566
 567      continue
          logtl=30+I
          rewind(logtl)
          DO JR=1,IR
           read(logtl) (arg(j),j=1,lfsum)
           CALL VDECIM(arg,1,XS,1,lfsum,NDECC,NCon)
           CALL VCLIP(xs,1,1E-20,1E20,xs,1,ncon)
           call valg10(xs,1,xs,1,ncon)
           call vsmul(xs,1,-10.0,xs,1,ncon)
           CALL CONDRB(1,NCON,NCON,XS)
          end do
         end if
        end do
       end if
c
c >>> Write TL to file 72 for option frcont
c
       if (frcont) then
        dnew=(xright-xleft)/(nconm-1)
        lfn=nconm
        DO 79 I=1,npar
          IF (IOUT(I).NE.0) THEN
            logtl=30+i
            rewind(logtl)
            DO 78 JR=1,IR
c >>> read linear field and sample
              read(logtl) (xs(ic),ic=1,lfsum)
              CALL Vfltsm(XS,X,lfsum,dxx,lfn,dnew)
c >>> convert to loss
              CALL VCLIP(x,1,1E-20,1E20,x,1,lfn)
              CALL VALG10(x,1,x,1,lfn)
              CALL VSMUL(x,1,10.0,x,1,lfn)
              CALL VNEG(x,1,x,1,lfn)
              write(72) (x(j,1),j=1,lfn)
 78         CONTINUE
          END IF
 79     CONTINUE
       end if

       DO I=npar,1,-1
        IF (IOUT(I).NE.0) THEN
          logtl=30+I
          close(logtl,status='delete')
        END IF
       end do

c >>> end of second idir loop
      end do

c
c >>> Close assembled transfer funtion file
c
       DO I=npar,1,-1
        IF (IOUT(I).NE.0) THEN
          logtrf=30+npar+I
          close(logtrf,status='delete')
        END IF
       end do


c >>>
c >>> End of frequency loop
c
 15   continue
c
c >>> close kernel file
c
      close(18,status='delete')
c
c >>> Contours vs frequency and range
c
      if (frcont) then
       dnew=(xright-xleft)/(nconm-1)
       lfn=nconm
       if (.not.deptav) then
       do 18 jr=1,ir
       do 18 j=1,npar
        if (iout(j).ne.0) then
         CALL CONFR(j,nconm,NFREQ,xleft,xright,
     1              XLEFT,XRIGHT,XINC,XAXIS,FREQ1,FREQ2,
     2              freq1,freq2,9.0,2.0,ZMIN(j),ZMAX(j),Zstep(j),
     3              TITLE,dnew,sd,rdc(jr))
         rewind(72)
         do 17 jf=1,nfreq
          FREQ=EXP(F1LOG+(Jf-1)*DFLOG)
          DO 17 I=1,npar
           IF (iout(i).ne.0) THEN
            DO 16 LR=1,IR
             read(72) (x(jj,2),jj=1,lfn)
             if (i.eq.j.and.lr.eq.jr) then
              call confab(lfn,freq)
             end if
 16         CONTINUE
           END IF
 17       CONTINUE
        end if
 18    continue
       end if
c
c >>> Frequency-range contours for depth average
c
       if (deptav) then
        do 28 j=1,npar
        if (iout(j).ne.0) then
         CALL CONFR(j,nconm,NFREQ,xleft,xright,
     1              XLEFT,XRIGHT,XINC,XAXIS,FREQ1,FREQ2,
     2              freq1,freq2,1.0,2.0,ZMIN(j),ZMAX(j),Zstep(j),
     3              TITLE,dnew,sd,0.0)
         rewind(73)
         do 27 jf=1,nfreq
          FREQ=EXP(F1LOG+(Jf-1)*DFLOG)
          DO 27 I=1,npar
           IF (iout(i).ne.0) THEN
             read(73) (x(jj,2),jj=1,lfn)
             if (i.eq.j) then
              call confab(lfn,freq)
             end if
           END IF
 27       CONTINUE
         end if
 28     continue
       end if
      end if
      OPTION(2)='PLTEND'
      WRITE(19,777) OPTION
 777  FORMAT(1H ,2A6)
C          
C          
      WRITE(6,9960)
 9960 FORMAT(//1H ,'*** OASTL FINISHED ***')
C
C     CLOSE PLOT FILES
C
      CLOSE(UNIT=19,STATUS='KEEP')
      CLOSE(UNIT=20,STATUS='KEEP')
      CLOSE(UNIT=21,STATUS='KEEP')
      IF (DRCONT) THEN
        CLOSE(UNIT=22,STATUS='KEEP')
        write(28,'(a)') '@EOF'
        CLOSE(UNIT=28,STATUS='KEEP')
        CLOSE(UNIT=29,STATUS='KEEP')
      END IF
      WRITE(6,9962) TOTTIM
 9962 FORMAT(//1H ,'*** TOTAL TIME: ',F10.3,' SECONDS ***')
      END  
C          
      SUBROUTINE GETOPT(IPROF,ICNTIN,ICONTU,cylgeo,backscat,fftint,
     &                  v_field,con_int)

C          
C     INPUT OF OPTIONS            
C          
      INCLUDE 'compar.f'
      INCLUDE 'comfip.f'
      CHARACTER*1 OPT(40)
      LOGICAL ICONTU
      logical cylgeo,backscat,fftint,oneway,v_field,con_int
      common /mconew/ oneway
      WRITE(6,300)                
 300  FORMAT(//1H ,'OPTIONS:',/)    
      DEBUG=.FALSE.
      DEPTAV=.FALSE.
      DRCONT=.FALSE.
      FRCONT=.false.
      PLTL=.FALSE.
      PLKERN=.FALSE.
      ANSPEC=.FALSE.
      TLDEP=.FALSE.
      SHEAR=.FALSE.
      mom_sou=.FALSE.
      ICONTU=.FALSE.
      SCTOUT=.FALSE.
      cylgeo=.true.
      oneway=.false.
      IREF=0                      
      ISTYP=0                     
c >>> force plane geometry
      ICDR=1
      IPROF=0
      ICNTIN=0
      srctyp=1
      nout=0
      DO 10 I=1,npar                 
 10   IOUT(I)=0                   
c >>>
      READ(1,200) OPT             
 200  FORMAT(40A1)                
      DO 50 I=1,40                
      IF (OPT(I).EQ.'N') THEN 
       IF (IOUT(1).GT.0) GO TO 50  
       NOUT=NOUT+1                 
       IOUT(1)=1                   
       WRITE(6,301)                
 301   FORMAT(1H ,'NORMAL STRESS') 
      ELSE IF (OPT(I).EQ.'V') THEN 
       IF (IOUT(2).GT.0) GO TO 50  
       NOUT=NOUT+1                 
       IOUT(2)=1                   
       WRITE(6,302)                
 302   FORMAT(1H ,'VERTICAL VELOCITY')                   
      ELSE IF (OPT(I).EQ.'H') THEN 
       IF (IOUT(3).GT.0) GO TO 50  
       NOUT=NOUT+1                 
       IOUT(3)=1                   
       WRITE(6,303)                
 303   FORMAT(1H ,'HORIZONTAL VELOCITY')                 
      ELSE IF (OPT(I).EQ.'R') THEN
       IF (IOUT(5).GT.0) GO TO 50              
       NOUT=NOUT+1            
       IOUT(5)=1              
       WRITE(6,3031)           
 3031  FORMAT(1H ,'RADIAL STRESS')             
      ELSE IF (OPT(I).EQ.'K') THEN
       IF (IOUT(6).GT.0) GO TO 50              
       NOUT=NOUT+1            
       IOUT(6)=1              
       WRITE(6,3032)           
 3032  FORMAT(1H ,'BULK STRESS')             
      ELSE IF (OPT(I).EQ.'S') THEN
       IF (IOUT(7).GT.0) GO TO 50              
       NOUT=NOUT+1            
       IOUT(7)=1              
       WRITE(6,'(1h ,a)') 'SHEAR STRESS'             
      ELSE IF (OPT(I).EQ.'A') THEN 
       IF (DEPTAV) GO TO 50     
       DEPTAV=.TRUE.
       WRITE(6,304)                
 304   FORMAT(1H ,'DEPTH AVERAGE')
      ELSE IF (OPT(I).EQ.'C') THEN 
       IF (DRCONT) GO TO 50    
       DRCONT=.TRUE.
       WRITE(6,305)           
 305   FORMAT(1H ,'DEPTH-RANGE CONTOURS')
      ELSE IF (OPT(I).EQ.'f') THEN 
       IF (FRCONT) GO TO 50    
       FRCONT=.TRUE.
       WRITE(6,3051)           
 3051  FORMAT(1H ,'FREQUENCY-RANGE CONTOURS')
      ELSE IF (OPT(I).EQ.'c') THEN 
       IF (ICONTU) GO TO 50    
       ICONTU=.TRUE.
       WRITE(6,306)           
 306   FORMAT(1H ,'DEPTH INTEGRAND CONTOURS')
      ELSE IF (OPT(I).EQ.'T') THEN 
       IF (PLTL) GO TO 50    
       PLTL=.TRUE.
       WRITE(6,307)           
 307   FORMAT(1H ,'TRANSMISSION LOSS')
      ELSE IF (OPT(I).EQ.'I') THEN
       IF (PLKERN) GO TO 50
       PLKERN=.TRUE.
       WRITE(6,309)
 309   FORMAT(1H ,'HANKEL TRANSFORM INTEGRANDS')
      ELSE IF (OPT(I).EQ.'P') THEN
       IF (.not.cylgeo) GO TO 50
       cylgeo=.false.
       WRITE(6,313)
 313   FORMAT(1H ,'PLANE GEOMETRY')
      ELSE IF (OPT(I).EQ.'L') THEN
       IF (LINA.GT.0) GO TO 50
       LINA=1
       extlar=.false.
       WRITE(6,'(1H ,a)') 'VERTICAL SOURCE ARRAY - Internal'
      ELSE IF (OPT(I).EQ.'l') THEN
       IF (LINA.GT.0) GO TO 50
       LINA=1
       extlar=.true.
       WRITE(6,'(1H ,a)') 'VERTICAL SOURCE ARRAY - External'
      ELSE IF (OPT(I).EQ.'a') THEN
       IF (ANSPEC) GO TO 50
       ANSPEC=.TRUE.
       WRITE(6,3091)
 3091  FORMAT(1H ,'ANGULAR SPECTRA')
      ELSE IF (OPT(I).EQ.'s') THEN
       IF (SCTOUT) GO TO 50
       SCTOUT=.TRUE.
       WRITE(6,3092)
 3092  FORMAT(1H ,'OUTPUT OF SCATTERING DISCONTINUITIES')
      ELSE IF (OPT(I).EQ.'Z') THEN
       IF (IPROF.GT.0) GO TO 50
       IPROF=1
       WRITE(6,314)
 314   FORMAT(1H ,'PLOT OF VELOCITY PROFILES')
      ELSE IF (OPT(I).EQ.'J') THEN
        IF (ICNTIN.GT.0) GO TO 50
        ICNTIN=1
        WRITE(6,315)
 315    FORMAT(1H ,'COMPLEX INTEGRATION CONTOUR')
      ELSE IF (OPT(I).EQ.'D') THEN
        IF (TLDEP) GO TO 50
        TLDEP=.TRUE.
        WRITE(6,316)
 316    FORMAT(1H ,'TRANSMISSION LOSS AS FUNCTION OF DEPTH')
      else if (opt(i).eq.'d') then
        if (doppler) go to 50
        doppler=.true.
        write(6,'(1h ,a,f6.1,a)') 'Doppler compensation, speed',
     &                            vrec,' m/s'
      ELSE IF (OPT(I).EQ.'X'.or.opt(i).eq.'2') THEN
        IF (SHEAR) GO TO 50
        srctyp=2
        SHEAR=.TRUE.
        ver_for=.true.
        WRITE(6,317)
 317    FORMAT(1H ,'VERTICAL POINT FORCE IN SOLID MEDIA')
      ELSE IF (opt(i).eq.'h'.or.opt(i).eq.'3') THEN
        IF (hor_for) GO TO 50
        srctyp=3
        hor_for=.true.
        WRITE(6,3170)
 3170   FORMAT(1H ,'HORIZONTAL POINT FORCE IN SOLID MEDIA')
      ELSE IF (opt(i).eq.'4') THEN
        IF (dip_sou) GO TO 50
        srctyp=4
        dip_sou=.TRUE.
        WRITE(6,3171)
 3171   FORMAT(1H ,'DIP-SLIP SOURCE IN SOLID MEDIA')
      ELSE IF (OPT(I).EQ.'m'.or.opt(i).eq.'5') THEN
        IF (mom_sou) GO TO 50
        srctyp=5
        mom_sou=.TRUE.
        WRITE(6,3172)
 3172   FORMAT(1H ,'MOMENT SOURCE IN SOLID MEDIA')
      ELSE IF (OPT(I).EQ.'Q') THEN
        IF (DEBUG) GO TO 50
        DEBUG=.TRUE.
        WRITE(6,318)
 318    FORMAT(1H ,'***** DEBUGGING *****')
      ELSE IF (OPT(I).EQ.'F') THEN
        fftint=.true.
        WRITE(6,319)
 319    FORMAT(1H ,'FFP INTEGRATION SCHEME')
      ELSE IF (OPT(I).EQ.'i') THEN
        icerc=.true.
        WRITE(6,320)
 320    FORMAT(1H ,'LePage ice reflection coefficient')
      ELSE IF (OPT(I).EQ.'r') THEN
        freerc=.true.
        WRITE(6,321)
 321    FORMAT(1H ,'LePage rough surface reflection coefficient')
      ELSE IF (OPT(I).EQ.'t') THEN
        tablerc=.true.
        WRITE(6,322)
 322    FORMAT(1H ,'TABULATED SURFACE REFLECTION COEFFICIENT')
      ELSE IF (OPT(I).EQ.'b') THEN
        bottomrc=.true.
        WRITE(6,323)
 323    FORMAT(1H ,'TABULATED BOTTOM REFLECTION COEFFICIENT')
      ELSE IF (OPT(I).EQ.'g') THEN
       IF (goff) GO TO 50
       goff=.true.
       WRITE(6,'(a)') 'Goff-Jordan power spectrum'
      else if (opt(i).eq.'x') then
        extrap=.true.
        write(6,'(A)') 'Using kernel extrapolation'
      ELSE IF (OPT(I).EQ.'v') THEN
       IF (v_field) GO TO 50
       v_field=.true.
       WRITE(6,'(a)') 'Virtual source strengths'
      ELSE IF (OPT(I).EQ.'B'.or.opt(i).eq.'b') THEN
       IF (backscat) GO TO 50
       if (oneway) then
        stop '>>> Option b and o mutually exclusive <<<'
       end if
       backscat=.true.
       WRITE(6,'(a)') 'Backscattering computed'
      else if (opt(i).eq.'o') then
        if (oneway) go to 50
        if (backscat) then
         stop '>>> Option b and o mutually exclusive <<<'
        end if
        oneway=.true.
        write(6,'(1h ,a)') 'Oneway Solution'
      else if (opt(i).eq.'q') then
        if (con_int) go to 50
        con_int=.true.
        write(6,'(1h ,a)') '45 deg contour integration'
      ELSE IF (OPT(I).NE.' ') THEN
        WRITE(6,399) OPT(I)
 399    FORMAT(1H ,'>>>> UNKNOWN OPTION: ',A1,' <<<<')
      ELSE
      END IF
 50   CONTINUE                    
      IF (NOUT.NE.0) RETURN       
      IOUT(1)=1                   
      NOUT=1                      
      WRITE(6,301)                
      RETURN                      
      END  
      SUBROUTINE AUTSAM(C1,C2,RMIN,RMAX,CMIN,CMAX,NW,IC1,IC2)
      PARAMETER (RFAC=3.0)
      PARAMETER (NR=100,NKMIN=2**8)
      common /refnum/ csref,isref
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnrd.f'
C *** DETERMINE WAVENUMBER SAMPLING INTERVAL
      DK=2*PI/(RFAC*RMAX*1E3)
      CREFA=min(csref,1500.0)
      DK=MIN(DK,2*PI*FREQ/(CREFA*NKMIN))
C *** MINIMUM WAVENUMBER INTERVAL
      DR=1E3*(RMAX-RMIN)/(NR-1)
c >>> The following was commented out for while. Re-enacted 960717
      dr=min(dr,rmax*1E3)
c >>> Make at least 1/5 wavelength
      dr=max(0.2*csref/freq,dr)
      WM=2*PI/DR
C *** INTEGRATION LIMITS
      WN1=2*PI*FREQ/C2
      WN2=2*PI*FREQ/C1
      WNMAX=MAX(WM,1.1*(WN2-WN1)+WN1)
      WNMIN=MIN(WNMAX-WM,WN1-0.1*(WN2-WN1))
      WNMIN=MAX(WNMIN,DK*0.5)
C *** NUMBER OF WAVENUMBER SAMPLING POINTS
       NW1=(WNMAX-WNMIN)/DK+1
       NW=NKMIN/2
 1     NW=NW*2
       IF (NW.LT.NW1) GO TO 1
       WNMAX=WNMIN+(NW-1)*DK
       IC1=(WN1-WNMIN)/(WNMAX-WNMIN)*(NW-1)+1
       IC2=(WN2-WNMIN)/(WNMAX-WNMIN)*(NW-1)+1
       IC1=MAX(1,IC1)
       IC2=MIN(NW,IC2)
       CMIN=2*PI*FREQ/WNMAX
       CMAX=2*PI*FREQ/WNMIN
       RETURN
       END
      BLOCK DATA SAFBK1
      INCLUDE 'compar.f'      
      DATA PROGNM /'OASTRD'/
      DATA OMEGIM /0.0/
      DATA LUGRN,LUTRF,LUTGRN,LUTTRF /30,35,30,35/
      DATA SHEAR,DECOMP,SCTOUT,NFLAG,PADE 
     &     /.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE./
      DATA MSUFT,MBMAX,MBMAXI,SRCTYP,ISROW,ISINC /1,1,2,1,1,0/
c      data rdtest /.true./
      END

