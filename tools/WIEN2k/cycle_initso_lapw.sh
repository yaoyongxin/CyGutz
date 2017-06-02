#!/bin/bash

case=UO2

for i in {0..15}
do 
    cd ${case}_${i}
    initso_lapw
    cd ..
done
