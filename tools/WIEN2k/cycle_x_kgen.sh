#!/bin/bash

case=UO2

for i in {0..15}
do 
    cd ${case}_${i}
    printf "4096 \n 0 \n" > k.inp
    x kgen < k.inp
    rm k.inp
    cd ..
done
