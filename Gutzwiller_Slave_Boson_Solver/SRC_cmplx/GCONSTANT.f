!****************************************************
      MODULE GCONSTANT
        USE gprec
        REAL(gq)   ,PARAMETER :: PI =3.141592653589793238_gq,TPI=2*PI
        COMPLEX(gq),PARAMETER :: ZI=(0._gq,1._gq)
        COMPLEX(gq),PARAMETER :: ZITPI = (0._gq,1._gq)*TPI
        REAL(gq)   ,PARAMETER :: RYTOEV=13.60569193_gq
        COMPLEX(gq),PARAMETER :: Z1=(1._gq,0._gq),Z0=(0._gq,0._gq)
        REAL(gq)   ,PARAMETER :: D1=1._gq,D0=0._gq
        REAL(gq)   ,PARAMETER :: SMALL=1.E-10_gq
        REAL(gq)   ,PARAMETER :: KTOEV=8.617E-5_gq ! eV/K
        REAL(gq)   ,PARAMETER :: KTORY=0.63333787391E-5_gq ! Rydberg/K
        INTEGER    ,PARAMETER :: MAXNNZ=100000000
        INTEGER    ,PARAMETER :: GIU=77 ! Temporary file ID
        REAL(gq)   ,PARAMETER :: RUBOUND=1.E30_gq, RLBOUND=1.E-30_gq
        INTEGER    ,PARAMETER :: IU_KGEN=14
!
      END MODULE GCONSTANT
!
