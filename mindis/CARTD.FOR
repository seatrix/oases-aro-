C
C     ***DATA AREA HOLDING CARTESIAN COORDINATE INFORMATION
C     *** JOHN STEIERT   3-DEC-81
C
C          CTXORG        VALUE OF X AT AXIS ORIGIN (CARTESIAN)
C          CTYORG        VALUE OF Y AT AXIS ORIGIN (CARTESIAN)
C          CTXSTP        X STEP INTERVAL IN USER'S UNITS
C          CTYSTP        Y STEP INTERVAL IN USER'S UNITS
C          CTXMAX        VALUE OF X AT END OF X-AXIS
C          CTYMAX        VALUE OF Y AT END OF Y-AXIS
C          CTXPST        X STEP INTERVAL IN PLOTTER UNITS
C          CTYPST        Y STEP INTERVAL IN PLOTTER UNITS
C
      REAL CTXORG, CTXSTP, CTXMAX, CTYORG, CTYSTP, CTYMAX
      REAL CTXPST, CTYPST
      COMMON /CARTD/ CTXORG, CTXSTP, CTXMAX, CTYORG, CTYSTP,
     $                CTYMAX, CTXPST, CTYPST
