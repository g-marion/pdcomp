#!/bin/bash

st=1  # Starting CM1 output number for looping
en=10  # Ending CM1 output number for looping

for(( i=$st; i<=$en; i++ )); do
  time

  cp /path/to/pdcomp/directory/run/def.pdcomp.input /path/to/pdcomp/directory/pdcomp/run/pdcomp.input

  sed -i "s/cm1out_000001/cm1out_$i/g" /path/to/pdcomp/directory/pdcomp/run/pdcomp.input 
  sed -i "s/pdcomp_000001/pdcomp_$i/g" /path/to/pdcomp/directory/pdcomp/run/pdcomp.input

  ./pdcomp.exe
  wait
done

time
