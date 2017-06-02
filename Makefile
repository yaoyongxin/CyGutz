GUTZ_ROOT = ${HOME}/WIEN_GUTZ/bin

install:
	#--------------------------------------------------------
	# Build CyGutz and install to ~/WIEN_GUTZ/bin
	# Build Gutzwiller_Slave_Boson_Solver
	mkdir -p ${HOME}/WIEN_GUTZ
	mkdir -p ${HOME}/WIEN_GUTZ/bin
	rm -rf  ${HOME}/WIEN_GUTZ/bin/*
	if [ "${GUTZ_ROOT}" != "${WIEN_GUTZ_ROOT}" ]; \
		then \
	  echo 'export WIEN_GUTZ_ROOT=${GUTZ_ROOT}' >> ~/.bashrc; \
		export WIEN_GUTZ_ROOT=${GUTZ_ROOT}; \
	fi
	cd Gutzwiller_Slave_Boson_Solver && make && make install && cd ..
	# Copy tools, examples to ~/WIEN_GUTZ/bin/
	cp -r tools examples ${GUTZ_ROOT}
	# Installation finished gracefully! Enjoy!


pyupdate:
	rm -r ${GUTZ_ROOT}/tools ${GUTZ_ROOT}/examples 
	cp -r tools examples ${GUTZ_ROOT}

clean:
	cd Gutzwiller_Slave_Boson_Solver && make clean && cd ..

