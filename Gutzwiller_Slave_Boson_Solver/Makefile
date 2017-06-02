SRC = external/cubpack external/Nitsol external/primme SRC_cmplx

all:
	mkdir -p lib include
	for dir in ${SRC} - ; do \
    if [ "$$dir" = "-" ] ; then break ; fi ; \
    (cd $$dir && ${MAKE} ; ) ; done

clean:
	for dir in ${SRC} - ; do \
    if [ "$$dir" = "-" ] ; then break ; fi ; \
    (cd $$dir && ${MAKE} clean ; ) ; done

install :
	cd SRC_cmplx && ${MAKE} install

