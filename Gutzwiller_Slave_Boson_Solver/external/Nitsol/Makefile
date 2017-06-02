include ../../../Makefile.in

OBJ = nitsol.o nitdrv.o nitjv.o nitbt.o nitgm.o nitstb.o nittfq.o nitbd.o

libnitsol.a:  $(OBJ)
	ar rs $@ $?
	cp $@ ../../lib/.

clean:
	rm -f *.o *.a ../../lib/libnitsol.a

# Dependencies.
nitsol.o:	nitsol.f nitprint.h nitinfo.h nitparam.h nitdflts.h
nitdrv.o:	nitdrv.f nitprint.h nitinfo.h nitparam.h
nitjv.o:	nitjv.f nitprint.h
nitbt.o:	nitbt.f nitprint.h nitparam.h
nitgm.o:	nitgm.f nitprint.h
nitstb.o:	nitstb.f nitprint.h
nittfq.o: 	nittfq.f nitprint.h
nitbd.o:	nitbd.f nitprint.h nitparam.h nitdflts.h


%.o: %.f
	$(F90) -c $(WFOPT) $<
