#!/bin/bash

rm *.vtu
rm *.pvtu
rm *.msh
ls . | grep -v -e "\." -e "\makefile" | xargs -r rm
