include ../../Makefile.in

INCLUDE=-I../include

#WFOPT =  -O0 -CB -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback

LIB= ../lib/libcubpack.a \
		 ../lib/libprimme.a  \
		 ../lib/libnitsol.a \
		 ${MKL_LIB}

OBJ  = GPREC.o GINFO.o GCUBPACK.o GCONSTANT.o gasa.o GMPI.o UTIL.o GBROYDEN.o \
	   GHDF5_BASE.o SPARSE.o GHDF5.o WAREHOUSE.o \
	   GGREENFUN.o GGREENFUN2.o BANDSTRU.o FOCKSTATE.o CORRORB.o \
       LOCALHAMIL.o garpack.o gprimme.o PROJECTOR.o GUTZ.o FERMI.o \
	   arpack_lib.o minpack_lib.o sparsekit_lib.o GAGENT.o ghybrd.o gnitsol.o

all: $(OBJ) MAIN.o
	$(PF90) -openmp -o CyGutz $(OBJ) MAIN.o $(LIB) ${HDF5_LIB}

touch:
	touch *f

%.o: %.f
	$(PF90) -openmp $(FREE_FORM) $(WFOPT) -c $<

clean:
	rm -f *o *mod

install:
	cp CyGutz ${DESTDIR}/.

arpack_lib.o: arpack_lib.f
	$(PF90) $(FIX_FORM) $(WFOPT) -c $*.f

sparsekit_lib.o: sparsekit_lib.f
	$(PF90) $(FIX_FORM) $(WFOPT) -c $*.f

GHDF5_BASE.o: GHDF5_BASE.f
	$(PF90) $(FREE_FORM) $(WFOPT) -c $*.f ${HDF5_INCLUDE}

GHDF5.o: GHDF5.f
	$(PF90) $(FREE_FORM) $(WFOPT) -c $*.f ${HDF5_INCLUDE} 

GCUBPACK.o: GCUBPACK.f
	$(PF90) $(FREE_FORM) $(WFOPT) -c $*.f $(INCLUDE)

gprimme.o: gprimme.f
	$(PF90) $(FIX_FORM) $(WFOPT) -c $*.f $(INCLUDE)
