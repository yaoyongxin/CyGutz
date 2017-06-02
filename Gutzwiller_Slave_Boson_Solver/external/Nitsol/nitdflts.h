
c  This include file contains defaults for various parameters that
c  control the nonlinear iterations.

c  Note that use of parameter statements here is not necessary, it
c  merely facilitates conversion to single precision:  merely change
c  each double precision declaration to real, and edit the definitions
c  of the constants in the first two parameter statements.

      double precision one,       two,       rtfiv
      parameter      ( one=1.0d0, two=2.0d0, rtfiv=2.23606797749978981 )
  
      double precision tenth,        half,        fournines
      parameter      ( tenth=0.10d0, half=0.50d0, fournines=one-1.0d-4 )
  
      double precision DFLT_CHOICE1_EXP
      parameter      ( DFLT_CHOICE1_EXP=(one+rtfiv)*half )
  
      double precision DFLT_CHOICE2_EXP
      parameter      ( DFLT_CHOICE2_EXP=two )
  
      double precision DFLT_CHOICE2_COEF
      parameter      ( DFLT_CHOICE2_COEF=one )
  
      double precision DFLT_ETA_CUTOFF
      parameter      ( DFLT_ETA_CUTOFF=tenth )
  
      double precision DFLT_ETA_MAX
      parameter      ( DFLT_ETA_MAX=fournines )
  
      double precision DFLT_THMIN
      parameter      ( DFLT_THMIN=tenth )
  
      double precision DFLT_THMAX
      parameter      ( DFLT_THMAX=half )
  
      double precision DFLT_ETA_FIXED
      parameter      ( DFLT_ETA_FIXED=tenth )

      integer     DFLT_PRLVL
      parameter ( DFLT_PRLVL=0 )

      integer     STDOUT
      parameter ( STDOUT=6 )
