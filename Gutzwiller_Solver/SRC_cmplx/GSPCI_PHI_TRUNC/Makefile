include ./Makefile.in


all: ${obj} MAIN.o
	${PF90} ${FFLAGS} -o exe_spci_phi_trunc MAIN.o \
		${obj} ${MKL_LIB} ${LIB}

GPRIMME.o: ../GPRIMME.f90
	${F90} ${FIX_FORM} ${FFLAGS} -c $< $(INCLUDE)

%.o: ./%.f90
	${PF90} ${FREE_FORM} ${FFLAGS} -cpp ${DCPP} \
			-c $< $(INCLUDE)

%.o: ../%.f90
	${PF90} ${FREE_FORM} ${FFLAGS} -cpp  ${DCPP} \
			-c $< $(INCLUDE)

clean:
	rm -f *o *mod *out exe_spci_phi_trunc

install:
	cp exe_spci_phi_trunc ${DESTDIR}
