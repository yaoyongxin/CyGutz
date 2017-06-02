#!/bin/bash

echo

# Default WIEN_GUTZ_ROOT.
GUTZ_ROOT=${HOME}/WIEN_GUTZ/bin
echo "Build CyGutz and install to ${GUTZ_ROOT}"

source ${HOME}/.bashrc
if [ "${GUTZ_ROOT}" != "${WIEN_GUTZ_ROOT}" ]; then 
  echo "export WIEN_GUTZ_ROOT=${GUTZ_ROOT}" >> ~/.bashrc; 
  echo "export PYTHONPATH=${GUTZ_ROOT}/tools" >> ~/.bashrc;
export WIEN_GUTZ_ROOT=${GUTZ_ROOT}
fi

make install
