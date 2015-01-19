#WFOPT = -O3  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback #-check # -O -p -g
#WLDFLAGS = $(WFOPT) -L/opt/intel/mkl/10.0.011/lib/32 -lpthread
#WLIBS = -lmkl_lapack -lmkl -lguide -pthread
#F90 = ifort

include ../Makefile.in

EXEC = dmft

f2pyopts = --fcompiler=intel 

###############################################################################
###############################################################################

FFLAGS = $(WFOPT)

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS    = reallocate.o modules.o mod_mpi.o gtfnam2.o w_atpar.o angle.o abcd.o atpar.o dvbes1.o \
          $(SRC_GUTZ) \
          errclr.o errflg.o harmon.o inverss.o l2main.o eigsys.o \
          dmftmain.o  latgen.o outerr.o outwin.o rint13.o \
	  rotate.o sphbes.o symop.o ylm.o cputim.o mknam.o diracout.o \
          dergl.o sym.o radint.o rout.o lomain.o inth.o inouh.o rotdef.o \
          cmpxqtl.o testxqtldiff.o testtransfo.o printgloc.o  \
	  cmpdmftrans.o findex.o zinv.o overlap1.o fermi.o  GUTZ_LINK.o trafo.o spin_rotate.o

pylibs = trafo.so trafoso.so


#..............................................................................

all : $(EXEC) 
# $(pylibs)

$(EXEC): $(OBJS) $(pylibs)
	$(PF90) -o $(EXEC) $(OBJS) $(WLDFLAGS) $(WLIBS) $(LIBGUTZ)

trafo.so : trafo.f
	$(CMP) -c trafo.f -m trafo
trafoso.so : trafoso.f
	$(CMP) -c trafoso.f -m trafoso
occupr.so : occupr.f
	$(CMP) -c occupr.f -m occupr
occupi.so : occupi.f
	$(CMP) -c occupi.f -m occupi
rdU.so : rdU.f90
	$(CMP) -c rdU.f90 -m rdU

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
#$(OBJS):	Makefile

#..............................................................................
#  remove unnecessary files (executable(s) are not removed)
#
clean:
	rm -f *.o *_tmp_.* *.P .sequential .parallel *.mod *.so $(EXEC) mod_mpi.f

#..............................................................................
#
#  define inference rules to generate object files from source files
#
.f.o:
	$(PF90) $(FFLAGS) -c $<
.c.o:
	$(CC)  -c $<


.f90.o:
	$(PF90) $(FFLAGS) -c $<

#.SUFFIXES : .F
.F.o:
	$(preproc) -P  -C -traditional $(Fmpi_define) $< $*.f
	$(PF90) $(FFLAGS) -c $*.f
	rm $*.f

