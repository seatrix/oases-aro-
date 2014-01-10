      PROGRAM OASPRD
c ********************************************************
c *                       OASES                          *
c *  Ocean Acoustic and Seismic Exploration Synthetics   *
c *                   Copyright (C)                      *
c *                  Henrik Schmidt                      *
c *       Massachusetts Institute of Technology          *
c *               Cambridge, MA 02139                    *
c ********************************************************
c
C     2-d - Range Dependent Environments - PULSE VERSION          
C     inttyp=-1 : tau-p seismograms
      parameter (maxsect = 1000)
      
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comnrd.f'

c >>> arrays for saving source parameters
      real vold(nla,6), sdcsav(nrd)
      integer laytold(nla),laysold(nrd)
      real vnxt(nla,6)
      integer laytnxt(nla),laysnxt(nrd)
c >>> Source strength array for backscattering
      complex uh(nrd),uv(nrd),uhb(nrd),uvb(nrd)
c >>> Sector variables
c >>> forward direction
      real slf(maxsect),dlr(maxsect)
      integer lfs(maxsect),indrng(maxsect),nvs(maxsect)
c >>> backward direction
      integer lfs_b(maxsect),indrng_b(maxsect)
c >>> vertical receiver pointers
      integer irtl(nrd)

      COMPLEX SLOW
      LOGICAL CFRFL,ICONTU,PLPOUT,CDROUT,GENTS,AUSAMP
      character*80 atitle,ctitle,stitle
      character*4 pident(npar)
      LOGICAL cylgeo,backscat,fftint,fftsec,v_field,trfs_save
      common /refnum/ csref,isref


      CHARACTER*50 FILENM
      DIMENSION X(NP2,NPAR),FF(2,NP3),PX(MODULO)   
      DIMENSION CONMAX(NPAR)
      COMPLEX CTRF(NPAR)
      DIMENSION FFS(2,NP),XS(NP2),AKM(100),AKM1(100),AKM2(100)           
      DIMENSION GRPVEL(NMOD),PHVEL(NMOD)
      CHARACTER*16 TITLEY(12)
      CHARACTER*6  OPTION(2)
      CHARACTER*80 TITLE

      EQUIVALENCE (X(1,1),CFF(1,1)),(XS(1),CFFS(1))              
      EQUIVALENCE (FF(1,1),CFF(1,1)),(FFS(1,1),CFFS(1))              
      equivalence (atitle,title)

C
C     COMPLEX FREQUENCY INTEGRATION ADDED AS OPTION 'O' 880414
C
      COMPLEX CCSQ,CDSQ,CC

      CHARACTER*6 TRFFILE(2)
      DATA TRFFILE /'trf   ',
     &              'trfb  '/
      DATA OPTION /'RDOASP','      '/         
      DATA TITLEY(1) /'NORMAL STRESS  ('/,        
     1     TITLEY(2) /'VERTICAL VEL.  ('/,        
     2     TITLEY(3) /'RADIAL VEL.    ('/,        
     3     TITLEY(4) /'Um - Vm        ('/,        
     4     TITLEY(5) /'RADIAL STRESS  ('/,        
     &     TITLEY(6) /'BULK STRESS    ('/,        
     &     TITLEY(7) /'SHEAR STRESS   ('/        
      DATA XAXI,YAXI /20.0,12.0/
      DATA CONMAX /NPAR*-1E20/
      data pident /' Szz',' w  ',' u  ',' v  ',
     &             ' Sxx',' P  ',' S  '/                     
C           
C ********************FORMATS*******************               
C           
C           
 209  FORMAT(1H ,'**********************************',
     &      /1H ,'*       OASES Version ',f2.2,'       *',
     &      /1H ,'**********************************',
     &     //1H ,'Range-Dependent Pulse Module RDOASP' )
 210  FORMAT(//1H ,A)         
500   FORMAT(/1H ,'SOURCE CENTRE FREQUENCY:',F10.2,' HZ ')             
550   FORMAT(/1H ,3X,'NW',5X,' ICW1',5X,' ICW2',
     &       5X,'ICUT1',5X,'ICUT2',/(1X,I5,4I10))  
  600 FORMAT(//,'  CMIN = ',G16.6,' M/S ',/,'  CMAX = ',G16.6,' M/S ')          
C           
C           
C **********************************************               
C           
C
c >>> default real axis contour
      cont_45=.false.
 
c
c Direct access file record length
c
      call dasrec(ldaun)

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
     &   '>>>>        CHANGE IN ALL SOURCE FILES AND RELINK   <<<<'      
      END IF
C           
c >>> flag range dependent version
      rdoas=.true.
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
      srctyp=1
      MSUFT=1
      MBMAX=1
      MBMAXI=2
      ISROW=1
c >>> force plane geometry
      icdr=1
C
C     DEFAULT: ISOTROPIC LAYERS
C
      NTISOL=0
C           
C           
C           
      CALL OPFILR(1,IOER)
      IF (IOER.GT.0) STOP '>>>> ERROR: INPUT FILE NOT FOUND <<<<'
      READ(1,'(A)')TITLE       
      write(6,209) version
      WRITE(6,210)TITLE      
      CALL GETOPT(IGRP,ISTACK,IBODY,ISDEP,IPROF,ICNTIN,CFRFL,ICONTU,
     &            PLPOUT,CDROUT,GENTS,cylgeo,backscat)
      insrct=srctyp
      trfs_save=trfsou

      if (doppler) then
        read(1,*) freqs,offdbin,istyp,vsou,vrec
      else IF (ICNTIN.GT.0) THEN
        READ(1,*) FREQS,OFFDBIN
      ELSE
        OFFDBIN=0E0
        READ(1,*)FREQS       
      END IF
      TOTTIM=0.
      WRITE(6,500)FREQS      
      if (doppler) then
       write(6,*) 'Source type:            ',istyp
       write(6,*) 'Source velocity:   ',vsou
       write(6,*) 'Receiver velocity: ',vrec
      end if
C           
C     READ IN ENVIRONMENTAL DATA

c >>> open file for sector environmental parameters
      call opnenv()

      if (iprof.ne.0.or.INTF.gt.0) then
       read(1,*,err=771) nsect,incsec
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

      call opfilw(22,ioer)
      write(22,'(a)') 'BOTTOM         1'
      rlft=0.0
      do isect=1,nsect
       read(1,*) numl,sectl
       backspace(1)
       CALL INENVI
       call savenv()
        rght=rlft+sectl*1000.0
        dbottom=fndbot()
        write(22,'(2F12.3)') rlft,dbottom
        write(22,'(2F12.3)') rght,dbottom
        rlft=rght
      end do      
c >>> initialize first sector
       call rwdenv()
       call rclenv()
       call rwdenv()

C          
C     SOURCE AND RECEIVER DATA
C
      CALL INSRC(SD)
      insrct=srctyp
      write(6,'(a,i3)') ' Source type:',insrct


      isref=lays((ls-1)/2+1)
      csref=v(isref,2)
      inls=ls
      do is=1,inls
       sdcsav(is)=sdc(is)
      end do
      CALL INREC(RD,RDLOW,idummy)
c
c >>> input receiver depths for TRF file
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
c     Wavenumber sampling parameters
c
      READ(1,*)CMININ,CMAXIN   
      READ(1,*)NWVNOIN,ICW1,ICW2,INTF

      AUSAMP=(NWVNOIN.LT.1)  
C ***
      IF (CMININ.EQ.0) STOP '*** CMIN MUST BE NON-ZERO ***'
      IF (((CMAXIN.GT.0).AND.(CMININ.GT.0)).OR.
     1    ((CMAXIN.LT.0).AND.(CMININ.LT.0))) THEN
       NFLAG=.FALSE.
      ELSE
       IF (CMININ.LE.0) STOP '*** CMIN/CMAX CONFLICT ***'
       CMAXIN=1E8
       NFLAG=.TRUE.
      END IF

c >>> Force full spectrum
c       CMAXIN=1E8
c       NFLAG=.TRUE.

      PLPOUT=PLPOUT.OR.((INTF.GT.0).AND.(.NOT.ICONTU))

c *** TIME-FREQUENCY-RANGE PARAMETERS
      READ(1,*) NX,FR1,FR2,DT,R0,RSPACE,NPLOTS

c >>> save range sampling
      r0trf=r0
      drtrf=rspace
      nrtrf=nplots
  
c >>> finish bottom shading file
       write(22,'(2F12.3)') (r0+(nplots-1)*rspace)*1000.,dbottom
       write(22,'(2F12.3)') (r0+(nplots-1)*rspace)*1000.,rdtllo
       write(22,'(2F12.3)')  0E0,rdtllo
       write(22,'(a)') '%EOF'
       endfile(22)
       close(22,status='keep')

C *** MODIFY FREQUENCY SAMPLING FOR AUTOMATIC SAMPLING
      IF (AUSAMP) THEN
C     **** commented out by bjs 2/26/97 ****
c      CALL VMAX(V(1,2),1,CREF,NUML)
c      RM=MAX(ABS(R0),ABS(R0+(NPLOTS-1)*RSPACE))*1E3
c      TREQ=RM/CREF
c989   NX=MAX(NX,nint(TREQ/DT))
c      IPOW=0
c990   IPOW=IPOW+1
c      II=2**IPOW
c      IF (II.LT.NX) GO TO 990
c      NX=II
c      IF (NX.GT.NP*2) THEN
c       DT=DT*2
c       GO TO 989
c      END IF
       WRITE(6,*)
c      WRITE(6,*) '>>> AUTOMATIC FREQUENCY SAMPLING PARAMETERS'
       WRITE(6,*) '>>> NON-AUTOMATIC FREQUENCY SAMPLING PARAMETERS'
       WRITE(6,*) '     NX=',NX
       write(6,*) '     DT=',DT
      ELSE
      END IF
C
C     CHECK THAT NX IS INTEGER POWER OF 2
C
       IPOW=0
  991  IPOW=IPOW+1
       II=2**IPOW
       IF (II.LT.NX) GO TO 991
       IF (II.NE.NX) THEN
        NX=II
        WRITE(6,*)
        WRITE(6,*) '>>>> NT MUST BE INTEGER POWER OF 2, CHANGED TO',
     &             NX,' <<<<'      
       END IF
       NX=MIN0(NX,2*NP)

      R1=R0*1000.            

      DLFREQ=1E0/(DT*NX)     
      MX=FR2/DLFREQ+2
      LX=FR1/DLFREQ+1
      LX=MAX0(LX,1)
      MX=MIN0(MX,NX/2)
      LXP1=MAX(2,LX)        
      if ((fr2-fr1).lt.dlfreq) then
       if (dlfreq*(mx-1)-fr1 .lt. fr1-dlfreq*(lxp1-1)) then
        lxp1=mx
       else
        mx=lxp1
       end if
      end if        
      FREQM=DLFREQ*(MX-1)    
      FREQ0=DLFREQ*(LXP1-1)    
      freq1=freq0
      freq2=freqm

      xright = r0trf + (nrtrf-1)*drtrf
c
c >>> Read sources from strf file and make frequency sampling 
c     compatible
c
      if (trfsou) then
       call rdstrf(freq0,freqm,nx,lxp1,mx,dt,dlfreq)
       write(6,*) 'Sources read from strf file:', ls
       write(6,*) 'LXP1 =', lxp1
       write(6,*) 'MX =', mx
       sd=sdc(1)
       inls=ls
C
C     DETERMINATION OF SOURCE LAYERS
C
       WRITE(6,908)
 908   FORMAT(/1H ,'SOURCE DATA:',//1H ,'  N0. ','   DEPTH  ',
     1       'LAYER','      ZU        ZL')
       DO 906 I=1,LS
 906    CALL SOURCE(V,NUML,SDC(I),LAYS(I),ZUS(I),ZLS(I))
       WRITE(6,907) 1,SDC(1),LAYS(1),ZUS(1),ZLS(1)
       IF (LS.GT.1) THEN
        WRITE(6,907) LS,SDC(LS),LAYS(LS),ZUS(LS),ZLS(LS)
       END IF
 907   FORMAT(1H ,I6,F10.1,I6,2F10.1)
       isref=lays((ls-1)/2+1)
       csref=v(isref,2)
       inls=ls
       do is=1,inls
        sdcsav(is)=sdc(is)
       end do
      end if
c 
c >>> make tables for dispersive media
c
      if (idlmax.gt.0) then
       call dltable(lxp1,mx,dlfreq)
      end if 

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

      IF (AUSAMP) THEN
C *** AUTOMATIC SAMPLING
       AUSAMP=.TRUE.
c >>> If option J, use complex wn, otherwize complex frequency.
       if (ICNTIN.eq.1) then
        CFRFL=.FALSE.
       else 
        CFRFL=.TRUE.
        ICNTIN=0
       end if
       cref=0e0
       do 981 l=1,numl
        if (v(L,2).lt.2E4) cref=max(cref,v(L,2))
 981    continue
       cref=csref
       RANREF=CREF*(NX*DT)
c >>> max and min ranges
       rm =0e0
       rmi=1e20
       do ii=1,nplots
        rr=abs(1E3*(r0+(ii-1)*rspace))
        rm=max(rm,rr)
        rmi=min(rmi,rr)
       end do
       rmi=max(rmi,rm*0.01)
       if (cfrfl) then
        RANM=1.5*max(RANREF,RM)
       else
        ranm=1.5*RM
       end if
       OFFDBIN=0E0
       WRITE(6,*)
       WRITE(6,*) '>>> AUTOMATIC SAMPLING '
       write(6,*) '    REFERENCE SPEED:',CREF
       if (ICNTIN.eq.1) then
        write(6,*) '    Real Frequency Integration <<<' 
        write(6,*) '    Complex Wavenumber Integration <<<'
       else 
        WRITE(6,*) '    Complex Frequency Integration'
        WRITE(6,*) '    Real Wavenumber Integration'
        WRITE(6,*) '    REFERENCE RANGE:',RANREF
       end if
       FREQ=FREQM
       CALL AUTSAM(CMININ,CMAXIN,rmi,RANM,CMIN,CMAX,
     &             NWVNO,ICW1,ICW2,ibody)
       WRITE(6,*) '    MAX NO. OF WAVENUMBERS:',NWVNO
       ICUT1=1
       ICUT2=NWVNO
      ELSE
       CMIN=CMININ
       CMAX=CMAXIN
       NWVNO=NWVNOIN
      END IF 

C
C     CHECK WHETHER TO INVOKE DEFAULT CONTOUR OFFSET
C
      IF (ICNTIN.GT.0.AND.OFFDBIN.LT.1E-10) THEN
        OFFDB=60.0*cref*(1E0/CMIN-1E0/CMAX)/NWVNO
        WRITE(6,*) 
        WRITE(6,*) 'DEFAULT CONTOUR OFFSET APPLIED,',OFFDB,
     &             ' dB/wavelength'
      ELSE
        OFFDB=OFFDBIN
      END IF
C
C     MOVE FREQUENCY INTEGRATION CONTOUR TO ATTENUATE
C     WRAP-AROUND BY FACTOR 50
C
      if (.not.trfsou) then
       IF (CFRFL) THEN
        OMEGIM=-LOG(50.0)*DLFREQ
       ELSE
        OMEGIM=0.0
       END IF
      else
       write(6,*)
       write(6,*) '>>> Imag(omega) = ', omegim,' from source TRF file'
      end if
      NIPLOT=0       
      IF (INTF.GT.0) NIPLOT=(MX-LXP1)/INTF+1  

C >>> Wavenumber sampling
      IF (.NOT.NFLAG) THEN
       WK0 = 2*PI*FREQM / CMAX
       WKMAX = 2*PI*FREQM / CMIN               
      ELSE
       WKMAX=2*PI*FREQM/CMIN
       dlwvno=wkmax/(nwvno-0.5)
       WK0=0.5*dlwvno
       CMAX=2*PI*FREQM/WK0
       ICUT1=1
       ICW1=1
      END IF
      IF (NWVNO.GT.1) THEN
        DLWVNO = ( WKMAX - WK0 ) / ( FLOAT(NWVNO-1) )            
      ELSE
        DLWVNO=1.0
      END IF
      DLRAN=1E3*RSPACE
      RRMAX=R1+(NPLOTS-1)*DLRAN
      if (tilt) rrmax=max(rrmax+abs(ofstar(1)),rrmax+abs(ofstar(ir)))

      IF (.NOT.AUSAMP) THEN
C           
C *** EXTEND WAVENUMBER INTERVAL FOR TAPERING
C
       IF (ICW1.LT.2) THEN
         ICUT1=ICW1
       ELSE
c         ICUT1=MAX(1,ICW1-NINT(0.05*NWVNO))
         icut1=1 
       END IF
       IF (ICW2.EQ.NWVNO) THEN
         ICUT2=ICW2
       ELSE
c         ICUT2=MIN(NWVNO,ICW2+NINT(0.05*NWVNO))
         icut2=nwvno
       END IF
       WRITE(6,600)CMIN,CMAX  
       WRITE(6,550)NWVNO,ICW1,ICW2,ICUT1,ICUT2           
      END IF
      IF (NFLAG) THEN
        WRITE(6,*) 'NEGATIVE SPECTRUM BY SYMMETRY'
      END IF
 
      WRITE(6,9980) DT,TMIN,DLFREQ,FREQ0,FREQM,NIPLOT,NOUT,NPLOTS,
     1            NOUT       
 9980 FORMAT(/1H ,'TIME STEP:        ',F12.6,' SECS',           
     1      /1H ,'MIN. TIME:        ',F12.6,' SECS',           
     1      /1H ,'FREQUENCY STEP:   ',F12.6,' HZ',             
     2      /1H ,'MIN. FREQUENCY:   ',F12.6,' HZ',             
     3      /1H ,'MAX. FREQUENCY:   ',F12.6,' HZ',             
     4      //1H ,'INTEGRAND PLOTS:  ',I5,' (*',I1,')',         
     5      /1H ,'PULSE PLOTS:      ',I5,' (*',I1,')')         
      WRITE(6,9979) R0,RSPACE         
 9979 FORMAT(/1H ,'MINIMUM RANGE:    ',F12.6,' KM',             
     1      /1H ,'RANGE STEP:       ',F12.6,' KM')             
c      LXP1=LX                
c      IF (LX.LE.1) LXP1=2
      NUMFR=MX-LXP1+1
C
C     OPEN FILES
C
      IF (PLPOUT) THEN
       CALL OPFILW(19,IOER)
       CALL OPFILW(20,IOER)
       WRITE(19,'(I8,10X,A40)') MODULO,'MODULO'
      END IF
      IF (DEBUG) CALL OPFILW(21,IOER)
C
C     OPEN TRANSFER FUNCTION SCRATCH FILES
C
c      IF (DECOMP) THEN
c        NCOMPO=5
c      ELSE
        NCOMPO=1
c      END IF

      LUGRN=30
      LUTRF=35

      CALL MORDER()

      write(6,*)
      write(6,*) 'Source type:  ',SRCTYP
      write(6,*) 'Fourier terms:',msuft
      write(6,*)
C
C     GENERATE TRANSFER FUNCTION FILE FOR POST-PROCESSOR
C
      if (nrtrf.gt.0) then
       call vmov(rdc,1,arg,1,abs(ir))
       irs=ir
       do ii=1,nrtl
        rdc(ii)=arg(irtl(ii))
       end do
       ir=-nrtl
       if (backscat) then
        ndir=2
       else
        ndir=1
       end if
       do I=1,ndir
        LUTTRF=LUTRF+I-1
        CALL TRFHEAD(TRFFILE(I),TITLE,rdtlup,rdtllo,R0,RSPACE,
     &              NX,LXP1,MX,DT,FREQS,SD)
       end do
       ir=irs
       call vmov(arg,1,rdc,1,abs(ir))
      end if

C *** INITIALIZE POINTERS AND PARAMETERS          
      CALL PINIT1
      IF (DEBUG) CALL PREQV(NUML,NUMI)

C *** PREPARE BESSEL FUNCTIONS FOR FULL INTEGRATION
c      IF (INTTYP.EQ.2) THEN
c       CALL PREPBF(RRMAX,WKMAX*1.3)
c       ICDR=1
c      END IF
c
c >>> source spectrum for doppler compensation
c
      if (doppler) then
       call spulse(freqs,dt,nx,lxp1,mx)
      end if

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

       call dasrec(ldaun)
       luntrf=8
       lunpnt=7
c >>> open direct access file or virtual source transfer functions
       open(luntrf,status='scratch',access='direct',form='unformatted',
     &       recl=2*ldaun*ir)
c >>> open direct access file or virtual source pointers
       open(lunpnt,status='scratch',access='direct',form='unformatted',
     &       recl=2*ldaun*ir)

c >>>
c >>> Frequency loop
c >>>
      DO 15 JJ=LXP1,MX       
       srctyp=insrct
       trfsou=trfs_save
       ftime=0e0
       CALL CLTIME
       nactf=jj
       KPLOT=0                
       IF (NIPLOT.LE.0) GO TO 14               
       IF (MOD(JJ-LXP1,INTF).EQ.0) KPLOT=1        
 14    CONTINUE               
       FREQ=(JJ-1)*DLFREQ     
       write(6,309) freq
 309   FORMAT(//1H ,'Frequency: ',F10.2,' Hz')
       write(6,'(a,i3)') ' Source type:',insrct

       DSQ=2E0*PI*FREQ+CMPLX(0E0,OMEGIM)
       RDSQ=DSQ
       CSQ=DSQ*DSQ            

c
c Here the computation of sector interface transfer functions should be done
c
       cosmic=1.0
       eps_r = cosmic*(rdlow-rd)/(ir-1)
       luntrf=8
       lunpnt=7
       call rwdenv()

       do isect=1,nsect
        call rclenv()   
        call vs_trf(luntrf,isect,eps_r,cminin,cmaxin)    
        nvs(isect)=isrow
       end do 
       endfile(8)
       endfile(7)
c
c >>> forward/backward loop
c
       srctyp=insrct
       do idir=1,ndir
C     
C     OPEN SCRATCH FILES FOR unassembled transfer functions
C
        LOGNUM=30+idir
        open(unit=lognum,status='scratch',form='unformatted')

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
           numlold=numl
           do ila=1,numl
            do j=1,6
             vold(ila,j)=v(ila,j)
            end do
            laytold(ila)=laytyp(ila)
           end do
           do ird=1,ir
            laysold(ird)=lay(ird)
           end do
c add virtual source offset to sector length
           rangs=1e3*sectl + eps_r
          end if

          if (sleft.lt.xright) then
           ncsect=ncsect+1
c >>> read in environmental parameters for next sector
           call rclenv()
           numlnxt=numl
           do ila=1,numl
            do j=1,6
             vnxt(ila,j)=v(ila,j)
            end do
            laytnxt(ila)=laytyp(ila)
           end do
           do ird=1,ir
            laysnxt(ird)=lay(ird)
           end do
c restore old environment
           numl=numlold
           do ila=1,numl
            do j=1,6
             v(ila,j)=vold(ila,j)
            end do
            laytyp(ila)=laytold(ila)
           end do
           do ird=1,ir
            lay(ird)=laysold(ird)
           end do
           if (isect.eq.nsect) then
            sectl=max(sectl,xright-sleft)
           end if
           sright=min(sleft+sectl,xright)

c >>> If not first sector, set source array
           if (isect.ne.isfst) then
c >>> insert here computation of source strengths
            rewind(18)
            call vs_arr(nsym,rangs,vnxt,laytnxt,laysnxt,
     &               isect,nvs(isect),isect-1,nvs(isect-1),
     &               uh,uv,uhb,uvb)
            if (backscat) then
             write(6,*) 'Writing sector ',isect-1
             write(9,rec=isect-1) (uhb(jr),jr=1,ir)
             write(9,rec=isect-1+(nsect-1)*ir) (uvb(jr),jr=1,ir)
            end if
           end if
c >>> re-initialize for next sector
           numl=numlnxt
           do ila=1,numl
            do j=1,6
             v(ila,j)=vnxt(ila,j)
            end do
            laytyp(ila)=laytnxt(ila)
           end do
           do ird=1,ir
            lay(ird)=laysnxt(ird)
           end do
c >>> Initialize receiver array
           call recarr()
c >>> If not first sector, set source array
           if (isect.ne.isfst) then
            srctyp=98
            trfsou=.false.
            call souarr()
           else
c 
c >>> set physical source
c
            srctyp=insrct
            ls=inls
            trfsou=trfs_save
            do is=1,inls
             sdc(is)=sdcsav(is)
             CALL SOURCE(V,NUML,SDC(is),LAYS(is),ZUS(is),ZLS(is))
            end do
            WRITE(6,908)
            write(6,'(a,i3)') ' Source type:',srctyp 
            WRITE(6,907) 1,SDC(1),LAYS(1),ZUS(1),ZLS(1)
            IF (LS.GT.1) THEN
             WRITE(6,907) LS,SDC(LS),LAYS(LS),ZUS(LS),ZLS(LS)
            END IF
           end if

           call pinit1()
          end if
         else
c >>> backward propagation
c >>> read in environmental parameters for sector
          if (isect.ne.isfst) then
           numlold=numl
           do ila=1,numl
            do j=1,6
             vold(ila,j)=v(ila,j)
            end do
            laytold(ila)=laytyp(ila)
           end do
           do ird=1,ir
            laysold(ird)=lay(ird)
           end do
           rangs=1e3*sectl + eps_r
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

          numlnxt=numl
          do ila=1,numl
           do j=1,6
            vnxt(ila,j)=v(ila,j)
           end do
           laytnxt(ila)=laytyp(ila)
          end do
          do ird=1,ir
           laysnxt(ird)=lay(ird)
          end do

          numl=numlold
          do ila=1,numl
           do j=1,6
            v(ila,j)=vold(ila,j)
           end do
           laytyp(ila)=laytold(ila)
          end do
          do ird=1,ir
           lay(ird)=laysold(ird)
          end do

c >>> read source strengths
c >>> computation of source strengths
          if (isect.ne.isfst) then
           rewind(18)
           call vs_arr(nsym,rangs,vnxt,laytnxt,laysnxt,
     &               isect,nvs(isect),isect+1,nvs(isect+1),
     &               uh,uv,uhb,uvb)
           read(9,rec=isect) (cbuf(jr),jr=1,ir)
           call vadd(uh,1,cbuf,1,uh,1,2*ir)
           read(9,rec=isect+(nsect-1)*ir) (cbuf(jr),jr=1,ir)
           call vadd(uv,1,cbuf,1,uv,1,2*ir)
          else
           read(9,rec=isect) (uh(jr),jr=1,ir)
           read(9,rec=isect+(nsect-1)*ir) (uv(jr),jr=1,ir)
          end if
c >>> re-initialize for new sector
          numl=numlnxt
          do ila=1,numl
           do j=1,6
            v(ila,j)=vnxt(ila,j)
           end do
           laytyp(ila)=laytnxt(ila)
          end do
          do ird=1,ir
           lay(ird)=laysnxt(ird)
          end do
c >>> Initialize receiver array
          call recarr()
c >>> If not first sector, set source array
          call souarr()
          call pinit1()
         end if

         if (sleft.lt.xright) then
c 
c >>> make tables for dispersive media
c
          if (idlmax.gt.0) then
           call dltable(2,2,freq)
          end if 

          if (idir.eq.1) then
           lf=0
           indrng(isect)=0
           r0=sright-sleft
           rstep=drtrf
           do jran=1,nrtrf
            rng=r0trf+(jran-1)*drtrf
            if (rng.gt.(sleft*1.0001).or.isect.eq.1) then
c >>> range in sector
             if (isect.eq.1) then
              r0=rng-sleft
             else
              r0=rng-sleft + 1e-3*eps_r
             end if
             indrng(isect)=jran
             do kran=jran,nrtrf
              rng=r0trf+(kran-1)*drtrf
              if (rng.le.(sright*1.0001).or.isect.eq.nsect) then
               lf=lf+1
              else
               go to 666
              end if 
             end do
             go to 666
            end if
           end do
c >>> save sector variables
 666       slf(isect)=sleft+r0
           dlr(isect)=rstep
           lfs(isect)=lf
          else
c >>> backward propagation
           lf=0
           indrng_b(isect)=0
           r0=sright-sleft
           rstep=drtrf
           do jran=nrtrf,1,-1
            rng=r0trf+(jran-1)*drtrf
            if (rng.lt.(sright*0.9999).or.isect.eq.nsect) then
c >>> range in sector
             r0=sright - rng 
             indrng_b(isect)=jran
             do kran=jran,1,-1
              rng=r0trf+(kran-1)*drtrf
              if (rng.ge.(sleft*0.9999).or.isect.eq.1) then
               lf=lf+1
              else
               go to 777
              end if 
             end do
             go to 777
            end if
           end do
c >>> save sector variables
 777       lfs_b(isect)=lf
          end if 
          R1=R0*1E3 + eps_r
          dlran=rstep*1e3

          write(6,*)
          write(6,*) ' >>> Frequency    ',freq,' Hz <<<'
          if (mod(idir,2).eq.1) then
           write(6,*)' >>> Forward Propagation  <<<'
          else
           write(6,*)' >>> Backward Propagation <<<'
          end if
          write(6,*) ' >>> Sector no.   ',isect,' <<<'
          write(6,*) ' >>> Left border: ', sleft,' km <<<'
          write(6,*) ' >>> Right border:', sright,' km <<<'

          write(6,*)
          write(6,*) ' >>> Range samples:',lf,' <<<'
          if (mod(idir,2).eq.1) then
           do kran=indrng(isect),indrng(isect)+lfs(isect)-1
            write(6,*) '    Range:       ', r0trf+(kran-1)*drtrf
           end do
          else
           do kran=indrng_b(isect),indrng_b(isect)-lfs_b(isect)+1,-1
            write(6,*) '    Range:       ', r0trf+(kran-1)*drtrf
           end do
          end if
          nactf=2
          IF (ausamp) THEN
C ***  AUTOMATIC SAMPLING
           AUSAMP=.TRUE.
c >>> max and min ranges
           rm =0e0
           rmi=1e3*abs(sleft-sright)
           do ii=1,lf
            rr=abs(1E3*(r0+(ii-1)*rstep))
            rm=max(rm,rr)
            rmi=min(rmi,rr)
           end do
           rmi=max(rmi,0.05*1E3*abs(sleft-sright))
           if (cfrfl) then
            RANM=1.5*max(RANREF,RM)
           else
            ranm=1.5*RM
           end if
           RMAXA=1E3*ABS(SLEFT-SRIGHT)
           ranm=max(ranm,2E0*rmaxa)
c           OFFDB=0e0
           CALL AUTSAM(CMININ,CMAXIN,RMI,RANM,CMIN,CMAX,
     &               NWVNO,ICW1,ICW2,ibody)
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
c         ICUT1=ICW1
c         ICUT2=ICW2

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
          WRITE(6,550)NWVNO,ICW1,ICW2,icut1,icut2
                      
          IF (NFLAG) WRITE(6,*) '*** NEGATIVE SPECTRUM BY SYMMETRY ***'

          DLWVNO = ( WKMAX - WK0 ) / ( FLOAT(NWVNO-1) )      
          rstfft=2e0*pi/(wkmax-wk0)
          RANMAX=NWVNO*rstfft

C
C     IF OPTION 'J', SET DEFAULT OFFDB
C
          IF (ICNTIN.GT.0) THEN
           IF (OFFDBIN.LT.1E-10) THEN
            OFFDB=60.0*cref*(1E0/CMIN-1E0/CMAX)/NWVNO
c            OFFDB=60.0*V(LAYS((LS-1)/2+1),2)/(FREQ*RANMAX)
            WRITE(6,*)
            WRITE(6,*) 'THE DEFAULT CONTOUR OFFSET IS APPLIED'
           ELSE
            WRITE(6,*)
            WRITE(6,*) 'THE USER DEFINED CONTOUR OFFSET IS APPLIED'
            offdb=offdbin
           END IF
           ATTMAX=OFFDB*FREQ*RANMAX/V(LAYS((LS-1)/2+1),2)
           WRITE(6,361) OFFDB,ATTMAX,ATTMAX*XRIGHT*1E3/RANMAX
 361       FORMAT(1H ,'CONTOUR OFFSET:         ',F12.6,' dB/wavelength',
     &           /1H ,'AT MAX FFT RANGE:       ',F12.2,' dB',
     &           /1H ,'AT MAX WINDOW RANGE:    ',F12.2,' dB',
     &           /1H ,'>> NOTE THAT COMPENSATION IS AUTOMATIC <<')
          END IF

c
c >>> Compute Green's functions
c

          CALL PINIT2 
c          write(6,*) 'omegim,offima=',omegim,offima
c >>> nsym = 1 for only symmetric sources
c >>> nsym =2 for sym. and antisym. sources.
          nsym=1
          if (isect.ne.1.or.idir.gt.1) then
           srctyp=98
           trfsou=.false.
           do ila=1,numl
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
          ftime=ftime+t1
          CALL CHKSOL
          WRITE(6,310) T1
 310      FORMAT(//1H ,'Integration kernels, CPU= ',F12.3,' SECS.')

C
C     WAVENUMBER INTEGRATION LOOP
C
c
          CALL CLTIME
          if ((kplot.eq.1).and.
     &       (mod(isect-1,incsec).eq.0.or.isect.eq.nsect)) then
           nrec=irtl(nrtl/2+1)
           write(6,*) 'PLINTGR, nrec=',nrec
           CALL RD_PHINT(nsym,nrec)
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
            CALL PLINTGR(DLWVNO,WK0,SD,RDC(NREC),stitle,20.0,12.0)
           end do
          end if

          call cltime()
C
C     TRANSFER FUNCTION CALCULATION
C
c          write(6,*) '>>> Calling rd_di_trf'
c          write(6,*) '>>> omegim,offima=',omegim,offima
          CALL rd_di_trf(nsym,lf,ir,r1,dlran,cbuf)
C
C     WRITE TRANSFER FUNCTIONS TO SCRATCH FILES
C
          lognum=30+idir
          do jran=1,lf
           do jdep=1,nrtl
            jout=0
            indx=jran+(irtl(jdep)-1)*lf
            do I=1,npar
             if (IOUT(I).gt.0) then
              jout=jout+1
              if (mod(idir,2).eq.1) then
               cfile(jout)=cff(indx,i)
              else
               if (i.eq.3.or.i.eq.7) then
                cfile(jout)=-cff(indx,i)
               else
                cfile(jout)=cff(indx,i)
               end if
              end if
             end if
            end do
            write(lognum) (cfile(jout),jout=1,nout)
           end do
          end do

          CALL RDTIME(T2)
          ftime=ftime+T2
          write(6,'(1h ,a,f8.3,a)') 'Integration done, CPU=',t2,' s'
C
C     close and Rewind SCRATCH FILES
C
          CALL CLSBUF(30)
c
c >>> End of sector loop
c
          sleft=sright
          if (debug) then
           write(6,*) '>>> end of sector',isect
           pause
          end if
         end if
c >>> end of sector loop
        end do
c >>> end of direction loop
       end do
       call cltime()
   
c >>> assemple trf file
       idir=1

       do i=1,nout
        call vclr(cff(1,i),1,2*nrtrf*nrtl)
       end do

       lognum=30+idir
       rewind(lognum)

       do isect=1,ncsect
        do jran=indrng(isect),indrng(isect)+lfs(isect)-1
c         write(6,*) 'jran=',jran, 'nrtl=',nrtl
         do jdep=1,nrtl
          indx=jdep+(jran-1)*nrtl
          read(lognum) (cff(indx,i),i=1,nout)
         end do
        end do
       end do
       close(lognum)

c >>> Write forward trf file

       luttrf=lutrf+idir-1
       do jran=1,nrtrf
        if (cylgeo) then
         rangem=max(1e3*(r0trf+(jran-1)*drtrf),1e0)
         fcr=1e0/sqrt(rangem)
        else
         fcr=1e0
        end if
        do jdep=1,nrtl
         indx=jdep+(jran-1)*nrtl
         write(LUTTRF) (cff(indx,i)*fcr,i=1,nout)
        end do
       end do
       
c >>> Compute backscatter

       if (ndir.gt.1) then
        idir=2
        do i=1,nout
         call vclr(cff(1,i),1,2*nrtrf*nrtl)
        end do
        lognum=30+idir
        rewind(lognum)
        do isect=ncsect-1,1,-1
         do jran=indrng_b(isect),indrng_b(isect)-lfs_b(isect)+1,-1
          do jdep=1,nrtl
           indx=jdep+(jran-1)*nrtl
           read(lognum) (cfile(i),i=1,nout)
           do i=1,nout
            cff(indx,i)=cfile(i)
           end do
          end do
         end do
        end do
        close(lognum)

c >>> Write backward trf file

        luttrf=lutrf+idir-1
        do jran=1,nrtrf
         if (cylgeo) then
          rangem=max(1e3*(r0trf+(jran-1)*drtrf),1e0)
          fcr=1e0/sqrt(rangem)
         else
          fcr=1e0
         end if
         do jdep=1,nrtl
          indx=jdep+(jran-1)*nrtl
          write(LUTTRF) (cff(indx,i)*fcr,i=1,nout)
         end do
        end do
       end if

       CALL RDTIME(T1)
       ftime=ftime+t1
       TOTTIM=TOTTIM+ftime
       WRITE(6,9988) JJ,FREQ,ftime
 9988 FORMAT(1H ,'FREQ. NO.',I4,' : ',F10.3,' HZ',             
     1       F10.3,' SECS')

 15   CONTINUE               

      IF (NPLOTS.GT.0) THEN
       DO 352 I=0,Ndir-1
        LUTTRF=LUTRF+I
        CLOSE(LUTTRF,STATUS='KEEP')
 352   CONTINUE
      END IF


C
C     THAT'S IT FOLKS
C
      IF (PLPOUT) THEN
       OPTION(2)='PLTEND'
       WRITE(19,'(1H ,2A6)') OPTION
      END IF
C           
C           
      CALL RDTIME(T1)
      TOTTIM=TOTTIM+T1
      WRITE(6,9960) TOTTIM
 9960 FORMAT(/1H ,'*** OASES PULSE FINISHED ***',
     1      //1H ,'    TOTAL TIME: ',F10.3,' SECS')
      END   
C           
      SUBROUTINE GETOPT(IGRP,ISTACK,IBODY,ISDEP,IPROF,ICNTIN,CFRFL,
     1                  ICONTU,PLPOUT,CDROUT,GENTS,cylgeo,backscat)      
c ********************************************************
c *                       OASES                          *
c *  Ocean Acoustic and Seismic Exploration Synthetics   *
c *                   Copyright (C)                      *
c *                  Henrik Schmidt                      *
c *       Massachusetts Institute of Technology          *
c *               Cambridge, MA 02139                    *
c ********************************************************
c
C     
C           
C     INPUT OF OPTIONS       
C           
      INCLUDE 'compar.f'
      LOGICAL cylgeo,backscat,oneway
      LOGICAL CFRFL,ICONTU,PLPOUT,CDROUT,GENTS
      CHARACTER*1 OPT(40)
      common /mconew/ oneway
      WRITE(6,300)           
 300  FORMAT(/1H ,'OPTIONS:',/)                
      INTTYP=0
      NOUT=0
      IREF=0
      ISTYP=-1           
      ICDR=0
      ISTACK=-1
      IBODY=0
      ISDEP=0
      IGRP=0
      IPROF=0
      ICNTIN=0
      CFRFL=.FALSE.
      SHEAR=.FALSE.
      mom_sou=.false.
      ICONTU=.FALSE.
      DECOMP=.FALSE.
      TRFOUT=.TRUE.
      PLPOUT=.FALSE.
      CDROUT=.FALSE.
      GENTS=.FALSE.
      DETERM=.FALSE.
      sctout=.false.
c >>> force plane geometry for internal computations
      ICDR=1
c >>> but default cylindrical geometry for spreading law
      cylgeo=.true.

      DO 10 I=1,NPAR            
 10   IOUT(I)=0              
      READ(1,200) OPT        
 200  FORMAT(40A1)           
      DO 50 I=1,40           
      IF (OPT(I).EQ.'B') THEN
       IF (IBODY.GT.0) GO TO 50              
       IBODY=1
       WRITE(6,202)           
 202   FORMAT(1H ,'SLOWNESS INTEGRATION')
      ELSE IF (OPT(I).EQ.'A') THEN
       IF (IBODY.GT.0) GO TO 50              
       IBODY=2
       WRITE(6,201)           
 201   FORMAT(1H ,'ANDYs INTEGRATION')
      ELSE IF (OPT(I).EQ.'N') THEN
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
c       IF (IOUT(3).GT.0.or.IOUT(NPAR).GT.0) GO TO 50              
       NOUT=NOUT+1            
       IOUT(3)=1              
       WRITE(6,303)           
 303   FORMAT(1H ,'HORIZONTAL VELOCITY') 
      ELSE IF (OPT(I).EQ.'R') THEN
       IF (IOUT(5).GT.0) GO TO 50              
       NOUT=NOUT+1            
       IOUT(5)=1              
       WRITE(6,306)           
 306   FORMAT(1H ,'RADIAL STRESS')             
      ELSE IF (OPT(I).EQ.'K') THEN
       IF (IOUT(6).GT.0) GO TO 50              
       NOUT=NOUT+1            
       IOUT(6)=1              
       WRITE(6,3061)           
 3061  FORMAT(1H ,'BULK STRESS')             
      ELSE IF (OPT(I).EQ.'S') THEN
       IF (IOUT(7).GT.0) GO TO 50              
       NOUT=NOUT+1            
       IOUT(7)=1              
       WRITE(6,'(1h ,a)') 'SHEAR STRESS'             
      ELSE IF (OPT(I).EQ.'L') THEN
       IF (LINA.GT.0) GO TO 50
       LINA=1
       extlar=.false.
       WRITE(6,308)
 308   FORMAT(1H ,'VERTICAL SOURCE ARRAY - Internal')
      ELSE IF (OPT(I).EQ.'l') THEN
       IF (LINA.GT.0) GO TO 50
       LINA=1
       extlar=.true.
       WRITE(6,309)
 309   FORMAT(1H ,'VERTICAL SOURCE ARRAY - External')
      ELSE IF (OPT(I).EQ.'v') THEN
       IF (LINA.GT.0) GO TO 50
       LINA=1
       trfsou=.true.
       WRITE(6,'(a)') 'VERTICAL SOURCE ARRAY - trf-file'
      ELSE IF (OPT(I).EQ.'P') THEN
       IF (.not.cylgeo) GO TO 50
       cylgeo=.false.
       WRITE(6,313)
 313   FORMAT(1H ,'PLANE GEOMETRY')
      ELSE IF (OPT(I).EQ.'Z') THEN
       IF (IPROF.GT.0) GO TO 50
       IPROF=1
       PLPOUT=.TRUE.
       WRITE(6,314)
 314   FORMAT(1H ,'PLOT OF VELOCITY PROFILES')
      ELSE IF (OPT(I).EQ.'J') THEN
        ICNTIN=1
        WRITE(6,315)
 315    FORMAT(1H ,'COMPLEX INTEGRATION CONTOUR')
      ELSE IF (OPT(I).EQ.'F') THEN
        INTTYP=1
        WRITE(6,316)
 316    FORMAT(1H ,'FILON INTEGRATION SCHEME')
c      ELSE IF (OPT(I).EQ.'f') THEN
c        INTTYP=2
c        WRITE(6,3161)
c 3161    FORMAT(1H ,'FULL BESSEL INTEGRATION SCHEME')
      ELSE IF (OPT(I).EQ.'O') THEN
        CFRFL=.TRUE.
        WRITE(6,317)
 317    FORMAT(1H ,'COMPLEX FREQUENCY INTEGRATION')
      else if (opt(i).eq.'8') then
        if (.not.double_trf) then
         double_trf=.true.
         write(6,'(1h ,a)') 'DOUBLE PRECISION TRF FILE'
        end if
      ELSE IF (OPT(I).EQ.'X'.or.opt(i).eq.'2') THEN
        IF (SHEAR) GO TO 50
        srctyp = 2
        SHEAR=.TRUE.
        ver_for=.true.
        WRITE(6,318)
 318    FORMAT(1H ,'VERTICAL POINT FORCE IN  SOLID MEDIA')
      ELSE IF (opt(i).eq.'h'.or.opt(i).eq.'3') THEN
        IF (hor_for) GO TO 50
        srctyp=3
        hor_for=.true.
        WRITE(6,3180)
 3180   FORMAT(1H ,'HORIZONTAL POINT FORCE IN SOLID MEDIA')
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
      ELSE IF (OPT(I).EQ.'C') THEN
        IF (ICONTU) GO TO 50
        ICONTU=.TRUE.
        CDROUT=.TRUE.
        WRITE(6,319)
 319    FORMAT(1H ,'FREQUENCY/SLOWNESS CONTOURS OF INTEGRANDS')
      ELSE IF (OPT(I).EQ.'U') THEN
        IF (DECOMP) GO TO 50
        DECOMP=.TRUE.
        WRITE(6,320)
 320    FORMAT(1H ,'WAVE FIELD DECOMPOSITION')
        IF (.NOT.TRFOUT) THEN
         TRFOUT=.TRUE.
         WRITE(6,321)
        END IF
      ELSE IF (OPT(I).EQ.'T') THEN
        IF (TILT) GO TO 50
        TILT=.TRUE.
        WRITE(6,321)
 321    FORMAT(1H ,'Tilted receiver arrays')
      ELSE IF (OPT(I).EQ.'Q') THEN
        IF (debug) GO TO 50
        debug=.TRUE.
        WRITE(6,322)
 322    FORMAT(1H ,'>>> debugging <<<')
      else if (opt(i).eq.'t') then
        inttyp=-1
        write(6,'(A)') 'tau-p SEISMOGRAMS'
      else if (opt(i).eq.'d') then
        doppler=.true.
        write(6,'(A)') 'Moving source and receivers'
      else if (opt(i).eq.'x') then
        extrap=.true.
        write(6,'(A)') 'Using kernel extrapolation'
      ELSE IF (OPT(I).EQ.'s') THEN
       IF (SCTOUT) GO TO 50
       SCTOUT=.TRUE.
       WRITE(6,3092)
 3092  FORMAT(1H ,'OUTPUT OF SCATTERING DISCONTINUITIES')
      ELSE IF (OPT(I).EQ.'g') THEN
       IF (goff) GO TO 50
       goff=.true.
       WRITE(6,'(a)') 'Goff-Jordan power spectrum'
      ELSE IF (OPT(I).EQ.'b') THEN
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
        if (cont_45) go to 50
        cont_45=.true.
        write(6,'(1h ,a)') '45 deg contour integration'
      ELSE
      END IF
 50   CONTINUE               
      GENTS=(ISTACK.GT.-1)
      IF (ISTYP.LT.0) THEN   
      ISTYP=2                
      WRITE(6,305) ISTYP     
      END IF
 305  FORMAT(1H ,'SOURCE TYPE:',I2)           
      IF (INTTYP.EQ.2) THEN
       ICNTIN=0
      END IF
      IF (NOUT.NE.0) RETURN  
      IOUT(1)=1              
      NOUT=1
      WRITE(6,301)           
      RETURN
      END


      BLOCK DATA BLKPUL
c ********************************************************
c *                       OASES                          *
c *  Ocean Acoustic and Seismic Exploration Synthetics   *
c *                   Copyright (C)                      *
c *                  Henrik Schmidt                      *
c *       Massachusetts Institute of Technology          *
c *               Cambridge, MA 02139                    *
c ********************************************************
c
      INCLUDE 'compar.f'
C
C**** DEFINITION OF MAX REAL ARGUMENT TO THE EXPONENTIAL FUNCTION
      COMMON /ARGMAX/ AM
C**** THE FOLLOWING DEFINITION SHOULD BE USED FOR THE FPS164
CFPS  DATA AM /300./
C**** THE FOLLOWING DEFINITION SHOULD BE USED FOR THE VAX
      DATA AM /65./     
      DATA OMEGIM /0.0/
      DATA PROGNM /'OASPRD'/
C3D   DATA PROGNM /'OASP3 '/
      DATA LUGRN,LUTRF,LUTGRN,LUTTRF /30,35,30,35/
      DATA SHEAR,DECOMP,SCTOUT,NFLAG /.FALSE.,.FALSE.,.FALSE.,.FALSE./
      DATA MSUFT,MBMAX,MBMAXI,SRCTYP,ISROW,ISINC /1,1,2,1,1,0/
      data bintrf /.true./
      END

      SUBROUTINE AUTSAM(C1,C2,rmin,RMAX,CMIN,CMAX,NW,IC1,IC2,ibody)
c ********************************************************
c *                       OASES                          *
c *  Ocean Acoustic and Seismic Exploration Synthetics   *
c *                   Copyright (C)                      *
c *                  Henrik Schmidt                      *
c *       Massachusetts Institute of Technology          *
c *               Cambridge, MA 02139                    *
c ********************************************************
c
c      PARAMETER (RFAC=1.5)
      PARAMETER (NR=150)
c >>> NOTE: ranges in meter
c >>> Minimum # of k, # of periods of exponential over taper interval
      PARAMETER (NKTMIN=2**9,WTAPER=10.0,NKDMIN=2**7)
      common /refnum/ csref,isref
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnrd.f'
      write(6,*) '>>> Reference range for sampling:',rmax,' m'
      if (inttyp.eq.1) then
        RFAC=2E0
      else
        RFAC=3E0
      end if
C *** DETERMINE WAVENUMBER SAMPLING INTERVAL
      DK=2*PI/(RFAC*RMAX)
      CREFA=min(csref,1500.0)
      DK=MIN(DK,2*PI*FREQ/(CREFA*NKDMIN))
c >>> Determine tapering interval from 
c     minimum source receiver separations
c >>> Depth-separation
      zsep=1E20
      zav=0e0
      do 200 isou=1,ls
       do 200 ircv=1,ir
        zsep=min(zsep,abs(rdc(ircv)-sdc(isou)))
 200    continue
      zsep=max(zsep,abs(sdc(ls)-sdc(1))/(10*ls))
      write(6,*) '>>> Reference depth for tapering:', zsep,' m'
      rref=rmin+wtaper*zsep
c >>> set to minimum one wavelength
      rref= max(rref,csref/freq)
      
      write(6,*) '>>> Reference range for tapering:', rref,' m'
      taperk=2*pi*wtaper/rref
      taperh=2*pi*wtaper/(max(rmin,v(lays(1),2)/freq))
      write(6,*) 'taperh,taperk=',taperh,taperk
C *** INTEGRATION LIMITS
      WN1=2*PI*FREQ/C2
      if (ibody.eq.2) then
       WN2=2*PI*(FREQ+0.2*freq2)/C1
      else
       WN2=2*PI*FREQ/C1
      end if
      wnmax=wn2+taperk
      WNMIN=max(dk,wn1-0.1*(wn2-wn1))
      DK=MIN(DK,(WNmax-WNmin)/(NKTMIN-1))
C *** NUMBER OF WAVENUMBER SAMPLING POINTS
       NW=(WNMAX-WNMIN)/DK+1
       WNMAX=WNMIN+(NW-1)*DK
       IC1=(WN1-WNMIN)/(WNMAX-WNMIN)*(NW-1)+1
       IC1=MAX(IC1,1)
       IC2=(WN2-WNMIN)/(WNMAX-WNMIN)*(NW-1)+1
       IC2=MIN(IC2,NWVNO)
       CMIN=2*PI*FREQ/WNMAX
       CMAX=2*PI*FREQ/WNMIN
       RETURN
       END

