#!/bin/bash

case=Am_fcc

for i in {0..23}
do
  echo ${case}_${i}
  cd ${case}_${i}
  if [ -d ykent* ]; then
    mv ykent*/* .
    echo 'ykent* exist!'
  fi
  cd ..
done

