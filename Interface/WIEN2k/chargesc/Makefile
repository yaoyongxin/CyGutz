SRC_GUTZ = gbase.o GMPI.o 

include ../../../Makefile.in

#.SUFFIXES:      .F
.SUFFIXES:      .F90
SHELL = /bin/bash 
CC = gcc


###############################################################################

#..............................................................................
#
#  Object files 
#
OBJS    =  modules.o mod_mpi.o cputim.o zeroin.o reallocate.o \
		   w_atpar.o abc.o  angle.o atpar.o c3fft.o charge.o \
           $(SRC_GUTZ) \
           csplit.o d5splt.o dergl.o deter.o diracout.o dvbes1.o \
           errclr.o errflg.o f7splt.o fermi.o \
           gaunt.o getfft.o gtfnam.o \
           harmon.o ifflim.o inouh.o inth.o kdelta.o \
		   l2mainSubs.o l2main.o dmft2.o \
           latgen.o lomain.o mknam.o notri.o open_helpfile.o \
		   outerr.o outwin.o \
           p3splt.o ph.o psplit.o recfil.o recpr.o rint13.o \
           rotate.o rotdef.o setfft1.o sphbes.o stern.o \
           sumupfft.o t3j.o t3j0.o ylm.o  sortag.o eigsys.o


OBJSMOD = abc.o atpar.o charge.o csplit.o d5splt.o dergl.o diracout.o \
          f7splt.o fermi.o \
          getfft.o gtfnam.o harmon.o inouh.o inth.o kdelta.o \
          l2mainSubs.o l2main.o dmft2.o latgen.o lomain.o outwin.o p3splt.o \
          psplit.o recfil.o recpr.o rint13.o rotdef.o setfft1.o \
          stern.o ylm.o


#  Build executable 

default: dmft2 
all : dmft2 

dmft2: $(OBJS)
	$(PF90) ${FREE_FORM} $(WLDFLAGS) -o dmft2 $(OBJS) $(MKL_LIB)

#  All routines depend upon an include file (contains common PARAMETERS)
#
$(OBJS):		Makefile
$(OBJSMOD):		modules.o w_atpar.o

#
#  remove unnecessary files (executable(s) are not removed)
#
clean:
	rm -f *.o *_tmp_.* *.P .real .complex .sequential .parallel \
		*.mod dmft2 combineud

.f.o:
	$(PF90) ${FREE_FORM} $(FFLAGS) -c $<

.c.o:
	$(CC) -c $<

.F.o:
	$(PF90) ${FREE_FORM} $(FFLAGS) -cpp $(Fmpi_define) -c $<
