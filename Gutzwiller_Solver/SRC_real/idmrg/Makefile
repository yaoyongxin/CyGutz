APP=idmrg

HEADERS=

CCFILES=$(APP).cc 

include ../../../Makefile.in

# Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))

# Rules ------------------
%.o: %.cc $(HEADERS)
	$(CCCOM) -c ${itensor_includes} -o $@ $<

$(APP): $(OBJECTS)
	$(CCCOM) ${itensor_includes} $(OBJECTS) -o exe_idmrg $(itensor_libs)

clean:
	rm -f  *.o $(APP) 

install:
	cp exe_idmrg ${DESTDIR}
