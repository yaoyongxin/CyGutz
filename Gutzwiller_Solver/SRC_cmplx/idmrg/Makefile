APP=idmrg

HEADERS=

CCFILES=$(APP).cc 

include ../../../Makefile.in

CCCOM=g++ -m64 -std=c++11 -fPIC -O3

# Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))

# Rules ------------------
%.o: %.cc $(HEADERS)
	$(CCCOM) -c ${itensor_includes} -o $@ $<

$(APP): $(OBJECTS)
	$(CCCOM) ${itensor_includes} $(OBJECTS) -o exe_idmrg $(itensor_libs)

clean:
	rm -f  *.o exe_idmrg

install:
	cp exe_idmrg ${DESTDIR}
