GUTZ_ROOT2 = ${HOME}/WIEN_GUTZ/bin2


install:
	#--------------------------------------------------------
	# Build CyGutz and install to ~/WIEN_GUTZ/bin2
	# Build Gutzwiller_Solver
	mkdir -p ${GUTZ_ROOT2}
	if [ "${GUTZ_ROOT2}" != "${WIEN_GUTZ_ROOT2}" ]; \
		then \
	  echo 'export WIEN_GUTZ_ROOT2=${GUTZ_ROOT2}' >> ~/.bashrc; \
		export WIEN_GUTZ_ROOT2=${GUTZ_ROOT2}; \
	fi
	# Build Gutzwiller solver 
	cd Gutzwiller_Solver && make && make install && cd ..
	# Build WIEN2k interface
	cd Interface/Wien2k &&  make && make install && cd ../..
	# Copy pygtool to ~/WIEN_GUTZ/bin2/
	cp -r pygtool/* ${GUTZ_ROOT2}
	# Install pyglib
	pip install -e ./pyglib --user
	# Installation finished gracefully! Enjoy!

gtool:
	cp -r pygtool/* ${GUTZ_ROOT2}

tests:
	cd test && ./tscan.py && cd ..

clean:
	cd Gutzwiller_Solver && make clean && cd ..
	cd Interface/Wien2k && make clean && cd ../..

clean_bin:
	rm -rf ${GUTZ_ROOT2}/*
