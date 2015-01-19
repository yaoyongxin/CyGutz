FC= mpif90

OFLAG= -O2 -traceback -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML
#OFLAG= -O0 -check bounds -check pointer -traceback

LIB= -mkl=sequential

OBJ  = GPREC.o GCONSTANT.o gasa.o GMPI.o UTIL.o SPARSE.o WAREHOUSE.o GGREENFUN.o BANDSTRU.o FOCKSTATE.o CORRORB.o \
       LOCALHAMIL.o PROJECTOR.o GUTZ.o garpack.o FERMI.o arpack_lib.o minpack_lib.o sparsekit_lib.o GAGENT.o TEST.o

all: $(OBJ) MAIN.o
	rm -f CyGutz
	$(FC) -o CyGutz $(OBJ) MAIN.o $(LIB)

touch:
	touch *f

.f.o:
	$(FC) $(OFLAG) -FR  -c $*.f

clean:
	rm -f *o *mod

install:
	cp CyGutz ~/WIEN_GUTZ/bin

arpack_lib.o:
	$(FC) $(OFLAG) -FI -c arpack_lib.f

sparsekit_lib.o:
	$(FC) $(OFLAG) -FI -c sparsekit_lib.f
