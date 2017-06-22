#!/bin/sh

./exe/annealer -i 10000 -m 1000 -n 2 -t 1 -T 2 1 A.net
./exe/measure -ce A.net  > /tmp/A.csv

cat A.net | ./exe/annealer -c 10 -i 100 -m 100 -n 2 -t 1 -T 2 1 B.net
./exe/measure -ce B.net  > /tmp/B.csv

cat B.net | ./exe/annealer -c 100 -i 10000 -m 1000 -n 2 -t 1 -T 2 1 C.net
./exe/measure -ce C.net > /tmp/C.csv

python test/testcontinuation.py

echo "Line should be continuous."

rm A.net B.net C.net
