#!/bin/bash
for line in {1..4}
do
    for stride in {16..25}
    do
        ./sim -i 01-mcf-gem5-xcg.trace -r 50 -p$line -s$stride >> 1.txt
    done
done