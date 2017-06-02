include ../../../Makefile.in

EXEC = dmft

FORT = $(PF90) ${FREE_FORM}

###############################################################################

#WFOPT = -O0  -check bound -check pointer -traceback

FFLAGS = $(WFOPT)

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS    = reallocate.o modules.o mod_mpi.o gtfnam2.o w_atpar.o angle.o abcd.o atpar.o dvbes1.o \
          errclr.o errflg.o harmon.o inverss.o l2main.o eigsys.o \
          dmftmain.o  latgen.o outerr.o outwin.o rint13.o \
	  rotate.o sphbes.o symop.o ylm.o cputim.o mknam.o diracout.o \
          dergl.o sym.o radint.o rout.o lomain.o inth.o inouh.o rotdef.o \
          cmpxqtl.o testxqtldiff.o testtransfo.o printgloc.o  \
	  cmpdmftrans.o findex.o zinv.o overlap1.o fermi.o  GUTZ_LINK.o spin_rotate.o

#..............................................................................

all : $(EXEC) 

$(EXEC): $(OBJS)
	$(FORT) -o $(EXEC) $(OBJS) $(MKL_LIB) 

#  remove unnecessary files (executable(s) are not removed)
clean:
	rm -f *.o *_tmp_.* *.P .sequential .parallel *.mod $(EXEC) mod_mpi.f GUTZ_LINK.f

#..............................................................................
#
#  define inference rules to generate object files from source files
#
.f.o:
	$(FORT) $(FFLAGS) -c $<
.c.o:
	$(CC)  -c $<


.f90.o:
	$(FORT) $(FFLAGS) -c $<

#.SUFFIXES : .F
.F.o:
	$(preproc) -P  -traditional $(Fmpi_define) $< $*.f
	$(FORT) $(FFLAGS) -c $*.f
	rm $*.f

