#!/bin/bash
#PBS -N scale
#PBS -l nodes=1:ppn=1:gpus=1

cd $PBS_O_WORKDIR

./a.out lena.png salida.png
