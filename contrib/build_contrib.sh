#!/bin/bash
export CONTRIB=`pwd`
#
# GSL
#
if [[ ! -f lib/libgsl.a ]]
then
  cd gsl-2.4
  CFLAGS="-mtune=native -march=native -std=c99" ./configure --prefix=${CONTRIB}
  make -j 4
  make install
  cd ..
fi
#
# IRP
#
#if [[ ! -f lib/libirp.a ]]
#then
#  cd irp
#  git pull
#  cp Makefile.manual Makefile
#  CONTRIB=${CONTRIB} PREFIX=${CONTRIB} make 
#  CONTRIB=${CONTRIB} PREFIX=${CONTRIB} make  install
#  cd ..
#fi
#
# Image Quality Assesment library 
#
#if [[ ! -f lib/libiqa.a ]]
#then
# cd iqa
#  make clean
#  RELEASE=1make 
#  PREFIX=${CONTRIB} make install
#cd ..
#fi
#
# Audio File Signal Processing (AFsp)
#
#if [[ ! -f lib/libAO.a ]]
#then
#cd AFsp*
#CFLAGS="-O3 -mtune=native -march=native" make 
# merge libraries into one; I don't like the name inconsistencies ...
#cd lib
#ar x libAO.a
#ar x libtsplite.a
#ar rcs libafsp.a *.o
#rm *.o
#cd ..
#cp lib/libafsp.a ../lib/
#mkdir -p ../include/afsp
#cp -r include/* ../include/afsp/
#cd ..
#fi
#
# zlib and png
#

#
# fftw
#
#if [[ ! -f lib/libfftw3.a ]]
#then
#  cd fftw-3.3.7
#  ./configure --prefix=${CONTRIB} --enable-sse2 --enable-avx2 --enable-generic-simd256 --enable-fma --enable-openmp
#  make install
#fi
