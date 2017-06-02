include ../../../Makefile.in

OBJ  = buckley.o internal_types.o ds_routines.o error_handling.o rule_1.o rule_c2.o \
			 rule_c3.o rule_cn.o rule_t2.o rule_t3.o rule_tn.o rule_general.o volume.o \
			 divide.o region_processor.o global_all.o check.o cui.o 

libcubpack.a: $(OBJ)
	ar qc $@ *.o *.mod
	cp $@ ../../lib/.
	cp cui.mod ../../include/.

clean:
	rm -f *.o *.mod *.a ../../lib/libcubpack.a ../../include/cui.mod

%.o: %.f90
	$(F90) $(WFOPT) -c $<
