F90 = ifort
FFLAGS = -O0 -check bound -check pointer -check arg_temp_created \
	    -traceback -check all

obj = GPREC.o GHDF5_BASE.o 

HDF5_BASE = /home/ykent/OPT/lib/hdf5-1.8.15-intel/
HDF5_LIB = -L$(HDF5_BASE)/lib/ -lhdf5_fortran -lhdf5
HDF5_INCLUDE = -I$(HDF5_BASE)/include


all: ${obj} vectorso2hdf5.o hdf52vectorso.o
	${F90} -o vectorso2hdf5 ${FFLAGS} vectorso2hdf5.o ${obj} \
		${HDF5_LIB} 
	${F90} -o hdf52vectorso ${FFLAGS} hdf52vectorso.o ${obj} \
		${HDF5_LIB} 


%.o: %.f90
	${F90} ${FFLAGS} -cpp -Dcmplx_mode -c $< ${HDF5_INCLUDE}

clean:
	rm -f *o *mod *h5 vectorso2hdf5
