#!/bin/bash
#
# IRP
#
#echo "-------- GET IRP (git) -------- "
#if [[ -d irp ]]
#then
#  cd irp
#  git pull
#  cd ..
#else
#  git clone https://@gitlab.fing.edu.uy/nacho/irp.git
#fi
#
# GSL
#
echo "-------- GET GSL (http) -------- "
[[ -d gsl ]] || wget -O- https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz | tar xzf -
#
# AFsp
#
#echo "-------- GET  AFSP (http) -------- "
#[[ -d AFsp-v10r1a ]] || wget -O- http://www-mmsp.ece.mcgill.ca/Documents/Downloads/AFsp/AFsp-v10r1a.tar.gz | tar xzf - 
#
# zlib
# 
#echo "-------- GET  ZLIB (http) -------- "
#[[ -d zlib ]] || git clone https://github.com/madler/zlib.git zlib
#
# libpng
#
#echo "-------- GET  PNG (http) -------- "
#[[ -d png ]] || git clone git://git.code.sf.net/p/libpng/code png

#echo "-------- FFTW 3.3.7 ------------- "
#[[ -d fftw-3.3.7 ]] || wget -O- http://www.fftw.org/fftw-3.3.7.tar.gz | tar xzf - 
