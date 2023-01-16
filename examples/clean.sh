#!/bin/bash

names='advection eulerTube heatTransfer eulerNaca navier'

for name in $names
do
    rm $name/*.vtu
    rm $name/*.pvtu
    rm $name/*.msh
    rm $name/log.txt
    rm -r $name/logs
    rm -r $name/results
    ls */** | grep -v -e "\." -e "\makefile" | xargs -r rm
done

