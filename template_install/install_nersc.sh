#!/bin/bash

echo
# Default WIEN_GUTZ_ROOT2.
GUTZ_ROOT2=${HOME}/WIEN_GUTZ/bin2
echo "Build CyGutz and install to ${GUTZ_ROOT2}"

source ${HOME}/.bashrc.ext
if [ "${GUTZ_ROOT2}" != "${WIEN_GUTZ_ROOT2}" ]; then 
  echo "# CyGutz" >> ~/.bashrc.ext;
  echo "export WIEN_GUTZ_ROOT2=${GUTZ_ROOT2}" >> ~/.bashrc.ext; 
  echo "export PATH=\${WIEN_GUTZ_ROOT2}:\${PATH}" >> ~/.bashrc.ext;
  source ${HOME}/.bashrc.ext
fi

make install
