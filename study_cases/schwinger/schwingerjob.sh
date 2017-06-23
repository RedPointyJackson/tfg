#!/bin/sh

cd ../../
make annealer > /dev/null
make measure > /dev/null

cd study_cases/schwinger/

T1=1.5
T2=1
T3=0.5
T4=2.0

../../exe/annealer -p -l 8 -i 1e6 -m 50 -n 2 ${T1} > T1.net
../../exe/annealer -p -l 8 -i 1e6 -m 50 -n 2 ${T3} > T2.net
../../exe/annealer -p -l 8 -i 1e6 -m 50 -n 2 ${T2} > T3.net
../../exe/annealer -p -l 8 -i 1e6 -m 50 -n 2 ${T4} > T4.net

../../exe/measure -pfs *.net > meas.dat

./plotter
