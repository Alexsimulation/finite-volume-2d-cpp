#! /bin/bash

mpirun -n 6 ./main > log.txt &
sleep 1
gnuplot plot.gpl
