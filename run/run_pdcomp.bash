#!/bin/bash

st=101
en=170

for(( i=$st; i<=$en; i++ )); do
  time

  cp /home/geoffrey.marion/pdcomp/run/def.pdcomp.input /home/geoffrey.marion/pdcomp/run/pdcomp.input

  sed -i "s/cm1out_101/cm1out_$i/g" /home/geoffrey.marion/pdcomp/run/pdcomp.input
  sed -i "s/pdcomp_101/pdcomp_$i/g" /home/geoffrey.marion/pdcomp/run/pdcomp.input

  ./pdcomp.exe
  wait
done

time
