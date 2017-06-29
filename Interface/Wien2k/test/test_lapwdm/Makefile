###############################################################################
#                                                                             #
#                          IBM RS6000  AIX 3.2.x                              #
#                                                                             #
#  make           ... generate executable for the REAL-version                #
#  make complex   ... generate executable for the COMPLEX-version             #
#  make clean     ... delete unnecessary files                                #
#                                                                             #
###############################################################################
#
# FC ........... compiler name
# FOPT ......... optimization flags
# FGEN ......... code generation flags (flags not related to optimization)
# LDFLAGS ...... linker flags
# R_LIBS ....... libraries needed to build the REAL executable
# C_LIBS ....... libraries needed to build the COMPLEX executable
# DESTDIR ...... dir. where the executable should go 
# R_EXECNAME ... name of the resulting REAL executable (without prefixed path!)
# C_EXECNAME ... name of the resulting COMPLEX executable ( -"- )
#
SHELL = /bin/sh
FC = mpiifort
CC = gcc
FOPT = -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include
FGEN = 
LDFLAGS = $(FOPT) -L$(MKLROOT)/lib/$(MKL_TARGET_ARCH) -pthread
R_LIBS = -lfftw3 -L/home/ykent/OPT/lib/fftw-3.3.3-intel/lib -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread
C_LIBS = $(R_LIBS)
DESTDIR	= ./
R_EXECNAME = lapwdm
C_EXECNAME = lapwdmc

###############################################################################

TYPE='REAL'
TYPE_COMMENT='!_REAL'
FFLAGS = $(FGEN) $(FOPT)
R_EXEC = $(DESTDIR)$(R_EXECNAME)
C_EXEC = $(DESTDIR)$(C_EXECNAME)

#..............................................................................
#
#  Object files common to both REAL and COMPLEX type subroutines
#
#
OBJS    = reallocate.o modules.o abc.o angle.o atpar.o dvbes1.o dmat.o\
          errclr.o errflg.o euler.o diracout.o dergl.o  \
          gtfnam2.o harmon.o inispl.o \
          inverss.o l2main.o lmtmat.o lapwdm.o  latgen.o outerr.o \
          outwin.o orb.o ph.o rint13.o inouh.o inth.o \
          rotate.o rotdef.o sphbes.o  symop.o t3j.o\
          SymmRot.o symmetrization.o ylm.o cputim.o mknam.o \
          xsplt.o locmat.o timeinv.o sym.o readw.o\
          radint.o rout.o lomain.o spin.o couplx.o printx.o xoper.o\
          xqproduct.o symmso.o addtinvso.o addtinv.o output.o 


#..............................................................................
#
#  Object files for REAL and COMPLEX versions
#
R_OBJS = $(OBJS) 
C_OBJS = $(OBJS)

#..............................................................................
#
#  Some source files contain both REAL and COMPLEX versions. These files
#  have the file extension '.frc' and must be treated differently.
#
.SUFFIXES:	.frc

#..............................................................................
#
#  Build executable (either REAL or COMPLEX versions)
#
real: rm_c_files
	$(MAKE) TYPE='REAL' TYPE_COMMENT='!_REAL' $(R_EXEC)

complex: rm_r_files
	$(MAKE) TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' $(C_EXEC)

#..............................................................................

$(R_EXEC): $(R_OBJS)
	$(FC) -o $(R_EXEC) $(R_OBJS) $(LDFLAGS) $(R_LIBS)

$(C_EXEC): $(C_OBJS)
	$(FC) -o $(C_EXEC) $(C_OBJS) $(LDFLAGS) $(C_LIBS)

rm_c_files:
	if [ -f .complex ]; then \
	   $(MAKE) clean; \
	fi
	touch .real

rm_r_files:
	if [ -f .real ]; then \
	   $(MAKE) clean; \
	fi
	touch .complex

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(R_OBJS):		
$(C_OBJS):		

#..............................................................................
#
#  remove unnecessary files (executable(s) are not removed)
#
# Dependancies
$(OBJS):		Makefile

clean:
	-rm -f *.o _tmp_.* *.P .real .complex *.mod

#..............................................................................
#
#  define inference rules to generate object files from source files
#  (From .frc files the appropriatly typed subroutine (REAL/COMPLEX) has to
#  be extracted first.)
#
.f.o:
	$(FC) $(FFLAGS) -c $<
.frc.o:
	@sed -e 's/^$(TYPE_COMMENT)//'   $< >$*_tmp.f
	@echo "$<: $(TYPE) version extracted"
	$(FC) $(FFLAGS) -c $*_tmp.f
	mv $*_tmp.o $*.o
	rm $*_tmp.f
