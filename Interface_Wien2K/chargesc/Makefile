DIR_GUTZ = /home/ykent/GUTZ/Gutz.MPI.5
SRC_GUTZ = gbase.o GMPI.o 

include ../Makefile.in

#.SUFFIXES:      .F
.SUFFIXES:      .F90
SHELL = /bin/bash 
FORT = $(PF90)
CC = gcc


###############################################################################

#..............................................................................
#
#  Object files 
#
OBJS    =  modules.o mod_mpi.o cputim.o zeroin.o reallocate.o w_atpar.o abc.o  angle.o atpar.o c3fft.o charge.o \
           $(SRC_GUTZ) \
           csplit.o d5splt.o dergl.o deter.o diracout.o dvbes1.o \
           errclr.o errflg.o f7splt.o fermi.o \
           gaunt.o getfft.o gtfnam.o \
           harmon.o ifflim.o inouh.o inth.o kdelta.o l2mainSubs.o l2main.o dmft2.o \
           latgen.o lomain.o mknam.o notri.o open_helpfile.o outerr.o outwin.o \
           p3splt.o ph.o psplit.o recfil.o recpr.o rint13.o \
           rotate.o rotdef.o setfft1.o sphbes.o stern.o \
           sumupfft.o t3j.o t3j0.o ylm.o  sortag.o eigsys.o


OBJSMOD = abc.o atpar.o charge.o csplit.o d5splt.o dergl.o diracout.o \
          f7splt.o fermi.o \
          getfft.o gtfnam.o harmon.o inouh.o inth.o kdelta.o \
          l2mainSubs.o l2main.o dmft2.o latgen.o lomain.o outwin.o p3splt.o \
          psplit.o recfil.o recpr.o rint13.o rotdef.o setfft1.o \
          stern.o ylm.o

RC_OBJS = 
SP_OBJS =

#..............................................................................
#
#  Build executable 
#

default: dmft2 combineud
all : dmft2 combineud

dmft2: $(OBJS)
	$(FORT) $(WLDFLAGS) -o dmft2 $(OBJS) $(WLIBS)

combineud : combineud.f90
	$(F90) -o combineud combineud.f90

gmpi:
	for file in $(DIR_GUTZ)/gbase.f $(DIR_GUTZ)/GMPI.f ; \
          do \
            fname=$$(basename $$file); \
            if [ -e $$fname ]; then\
              if ! diff $$fname  $$file  >/dev/null ; then cp $$file . ; fi; \
            else\
              cp $$file .; \
            fi; \
          done

#  All routines depend upon an include file (contains common PARAMETERS)
#
$(OBJS):		Makefile
$(OBJSMOD):		modules.o w_atpar.o

mod_mpi.o : mod_mpi.F
	$(preproc) -P  -C -traditional $(Fmpi_define) -o mod_mpi.f mod_mpi.F
	$(FORT) $(WFOPT) -c mod_mpi.f
	rm mod_mpi.f


#..............................................................................
#
#  remove unnecessary files (executable(s) are not removed)
#
clean:
	rm -f *.o *_tmp_.* *.P .real .complex .sequential .parallel *.mod dmft2 combineud

#..............................................................................
#
#  define inference rules to generate object files from source files
#  (From .F files the appropriatly typed subroutine (REAL/COMPLEX) has to
#  be chosen.)
#

.f.o:
	$(FORT) $(WFOPT) -c $<
.c.o:
	$(CC) -c $<
.f90.o:
	$(FORT) $(WFOPT) -c $<
.F90.o:
	@sed -e 's/^$(TYPE_COMMENT)//'   $< >$*_tmp_.F90
	@echo "$<: $(TYPE) version extracted"
	$(FORT) $(WFOPT) -c $*_tmp_.F90
	mv $*_tmp_.o $*.o
	rm $*_tmp_.F90

.SUFFIXES : .F
.F.o:
	$(preproc) -P  -C -traditional $(Fmpi_define) $< $*.f
	$(PF90) $(WFOPT) -c $*.f
	rm $*.f


#.F.o:
#	@sed -e 's/^$(TYPE_COMMENT)//'   $< >$*_tmp_.F
#	@echo "$<: $(TYPE) version extracted"
#	$(FORT) $(WFOPT) -c $*_tmp_.F
#	mv $*_tmp_.o $*.o
#	rm $*_tmp_.F
