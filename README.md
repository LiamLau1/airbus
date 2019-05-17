# airbus
Simulations for Airbus FYI

gcc -lm -c -Wall -fpic -o model.o radbelt_trmfun.c
gcc -shared -o libmodel.so model.o
