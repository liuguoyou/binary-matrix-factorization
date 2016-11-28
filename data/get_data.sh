#!/bin/bash
for file in mnist usps
do
  wget -c http://iie.fing.edu.uy/~nacho/data/bmf/${file}.7z
  7zr x ${file}.7z
done
