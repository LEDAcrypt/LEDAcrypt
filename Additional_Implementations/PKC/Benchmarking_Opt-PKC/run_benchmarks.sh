#!/bin/bash
for l in `seq 128 64 256`
do 
  for blocks in SL 64
    do 
     bin/$l\SL_\DFR$blocks/test_$l\SL_\DFR$blocks 
  done 
done
