#!/bin/sh



cd ../../
make annealer > /dev/null
make measure > /dev/null

cd study_cases/ising/

for T in $(seq 10 -0.1 0); do
    echo "Doing T = $T"
    ../../exe/annealer -l 16 -f -i 1e3 -m 5 -t 10 -T 100 -n 5 ${T} > thermalized.net
    cat thermalized.net | ../../exe/annealer -c 1e3 -l 16 -i 100 -m 5 -t 10 -T 1000 -n 5 ${T} > net_${T}.net
done

../../exe/measure -ec net_*.net > ising.data

./plotter

rm *.net
