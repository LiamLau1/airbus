# airbus
Simulations for Airbus FYI

## Compiling to object file
gcc -lm -c -Wall -fpic -o model.o radbelt_trmfun.c
## Compiling to shared lib
gcc -shared -o libmodel.so model.o
