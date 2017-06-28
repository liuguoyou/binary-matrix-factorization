#!/bin/bash
for file in billie billie2 camera cuadro cuadro2 einstein
do
  wget -c http://iie.fing.edu.uy/~nacho/data/bmf/${file}.7z
  7zr x ${file}.7z
done
