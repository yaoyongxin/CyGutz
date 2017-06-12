include ../../../Makefile.in

EXEC = dmft

###############################################################################

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS    = reallocate.o modules.o mod_mpi.o gtfnam2.o \
		  w_atpar.o angle.o abcd.o atpar.o dvbes1.o \
          errclr.o errflg.o harmon.o inverss.o l2main.o eigsys.o \
          dmftmain.o  latgen.o outerr.o outwin.o rint13.o \
	  	  rotate.o sphbes.o symop.o ylm.o cputim.o mknam.o diracout.o \
          dergl.o sym.o radint.o rout.o lomain.o inth.o inouh.o rotdef.o \
          cmpxqtl.o testxqtldiff.o testtransfo.o printgloc.o  \
	      cmpdmftrans.o findex.o zinv.o overlap1.o fermi.o  \
		  GUTZ_LINK.o spin_rotate.o


all : $(EXEC) 

$(EXEC): $(OBJS)
	$(PF90) ${FREE_FORM} -o $(EXEC) $(OBJS) $(MKL_LIB) 

#  remove unnecessary files (executable(s) are not removed)
clean:
	rm -f *.o *_tmp_.* *.P .sequential .parallel *.mod $(EXEC) 

#
#  define inference rules to generate object files from source files
#
.f.o:
	$(PF90) ${FREE_FORM} $(FFLAGS) -c $<

.c.o:
	$(CC)  -c $<

.F.o:
	$(PF90) ${FREE_FORM} $(FFLAGS) -cpp $(Fmpi_define) -c $<

