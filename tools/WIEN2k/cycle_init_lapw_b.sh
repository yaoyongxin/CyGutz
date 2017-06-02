#!/bin/bash

case=UO2

for i in {0..15}
do 
    cd ${case}_${i}
    init_lapw -b -vxc 5 -rkmax 8.5 -numk 4096
    cd ..
done
