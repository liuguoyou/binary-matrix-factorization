#!/bin/bash
for i in *.pbm
do 
  time pnmtopng -compression 9 ${i} > ${i/pbm/png} 
  time 7zr a ${i/pbm/7z} $i
  time minidjvu $i ${i/pbm/djvu} 
done
