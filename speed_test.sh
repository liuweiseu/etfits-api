#!/bin/bash
# Get the speed (time it took for each program to finish 
# on every file in a given path, using s6fits.c, s6fits.so (python module)
# and Kyle's program that does basically the same thing with fitsio.


echo "file sizes in nearest MB (file sizes below 1MB = 1MB)" 
while read p; do
  file=$p
  ls -s --block-size=M $file 

  echo "number of hits: $(./test $file 0)" 
  echo "speed for s6fits.c:"
  { time ./test $file ; } 2> temp.txt
  cat temp.txt
  rm temp.txt
  echo "speed for s6python:" 2> temp.txt
  { time python s6python.py $file ; } 2> temp.txt
  cat temp.txt
  rm temp.txt
  echo "speed for fitsio_fitsdata:" 
  { time python fitsio_fitsdata.py $file; } 2> temp.txt
  cat temp.txt
  rm temp.txt
done <files_to_be_timed.txt
