#!/bin/bash
EXPECTED_ARGS=2
ARGC=$#
if [ $ARGC -ne $EXPECTED_ARGS ]; then
  echo "Benchmarking script for LEDACrypt"
  echo "first parameter: primitive to benchmark among KEM-CPA,KEM,PKC"
  echo "second parameter: 0 for human readable output 1 for machine readable output"
else
  CIPHER=$1
  PREFIX="./build/bin/$CIPHER/LEDAcrypt_$CIPHER"
  if [ $CIPHER != "KEM-CPA" ]; then
    for DFR in {64,SL}
    do
       for SL in {1,3,5}
       do
          for NZERO in {2,3,4}
          do
                 COMMAND=$PREFIX"_SL_$SL""_N0_$NZERO""_DFR_$DFR $2"
                 $COMMAND
             done
          done
    done
  else
  for SL in {1,3,5}
  do
     for NZERO in {2,3,4}
     do
     COMMAND=$PREFIX"_SL_$SL""_N0_$NZERO $2"
     $COMMAND
     done
  done
  fi;
fi;
