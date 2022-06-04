#!/usr/bin/env bash
gfortran -c load_config.F90 -o load_config.o -I$HOME/.local/jsonfortran-gnu-8.3.0/lib/ 

gfortran -o load_config load_config.o  -L$HOME/.local/jsonfortran-gnu-8.3.0/lib/ -ljsonfortran 


./load_config
