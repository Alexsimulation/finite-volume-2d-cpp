#! /bin/bash

while getopts n: flag
do
    case "${flag}" in
        n) nprocesses=${OPTARG};;
    esac
done

mpirun -n ${nprocesses} ./main > log.txt &
sleep 1
gnuplot plot.gpl
