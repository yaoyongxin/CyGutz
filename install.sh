#!/bin/bash

echo

# Default WIEN_GUTZ_ROOT2.
GUTZ_ROOT2=${HOME}/WIEN_GUTZ/bin2
echo "Build CyGutz and install to ${GUTZ_ROOT2}"

source ${HOME}/.bashrc
if [ "${GUTZ_ROOT2}" != "${WIEN_GUTZ_ROOT2}" ]; then 
  echo "export WIEN_GUTZ_ROOT2=${GUTZ_ROOT2}" >> ~/.bashrc; 
  echo "export PATH=\${WIEN_GUTZ_ROOT2}:\${PATH}" >> ~/.bashrc;
  source ${HOME}/.bashrc
fi

make install
