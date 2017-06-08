#!/bin/bash

echo

# Default WIEN_GUTZ_ROOT2.
GUTZ_ROOT2=${HOME}/WIEN_GUTZ/bin2
echo "Build CyGutz and install to ${GUTZ_ROOT2}"

source ${HOME}/.bashrc
if [ "${GUTZ_ROOT2}" != "${WIEN_GUTZ_ROOT2}" ]; then 
  echo "export WIEN_GUTZ_ROOT2=${GUTZ_ROOT2}" >> ~/.bashrc; 
export WIEN_GUTZ_ROOT2=${GUTZ_ROOT2}
fi

make install
