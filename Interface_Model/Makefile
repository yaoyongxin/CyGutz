# Gutz set up
#-----------------------------------------------------------------------

# fortran compiler and linker
#-----------------------------------------------------------------------
FC=ifort
CC=gcc

#-----------------------------------------------------------------------
# general fortran flags  (there must a trailing blank on this line)
#-----------------------------------------------------------------------

FFLAGS =  -FR -assume byterecl 

#-----------------------------------------------------------------------
# optimization
#-----------------------------------------------------------------------

OFLAG= -O3 
OFLAG= -FR -O0 -check bounds -check pointer -traceback

LAPACK= -mkl ./libdmy.a

LIB  = $(LAPACK) $(BLAS)


#-----------------------------------------------------------------------
# general rules and compile lines
#-----------------------------------------------------------------------

SOURCE= prec.o  TB.o GUTZ_LINK.o

TBY: $(SOURCE) $(FFT3D) MAIN.o 
	rm -f TBY
	$(FC) -o TBY $(LINK) $(SOURCE) MAIN.o  $(FFT3D) $(LIB) 

clean:	
	-rm -f  *.o *.mod ; touch *.f

.f.o:
	$(FC) $(FFLAGS) $(OFLAG)  -c $<
.c.o:
	$(CC)  -c $<
