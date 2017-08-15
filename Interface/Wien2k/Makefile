include ../../Makefile.in

SRC = ksum chargesc

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
	cp ksum/gwien1 ${DESTDIR}/
	cp chargesc/gwien2 ${DESTDIR}/

