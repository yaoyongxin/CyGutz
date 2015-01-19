include Makefile.in

SRC = ksum chargesc CyGutz

all:
	for dir in ${SRC} - ; do \
		if [ "$$dir" = "-" ] ; then break ; fi ; \
		(cd $$dir && ${MAKE} all ; ) ; done

clean:
	for dir in ${SRC} - ; do \
		if [ "$$dir" = "-" ] ; then break ; fi ; \
		(cd $$dir && ${MAKE} clean ; ) ; done
	rm -f *.pyc


install :
	if [ ! -d ${DESTDIR} ] ; then mkdir -p ${DESTDIR} ; fi
	if [ ! -d ${DESTDIR}/database ] ; then mkdir -p ${DESTDIR}/database ; fi
	cp ksum/dmft ${DESTDIR}/
	cp ksum/*.so ${DESTDIR}/
	cp chargesc/dmft2 ${DESTDIR}/
	cp chargesc/combineud ${DESTDIR}/
	cp python/* ${DESTDIR}/
	cp CyGutz/CyGutz ${DESTDIR}/
