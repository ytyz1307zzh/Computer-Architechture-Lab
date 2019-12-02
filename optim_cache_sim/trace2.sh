#!/bin/bash
for line in {1..4}
do
    for stride in {1..15}
    do
        ./sim -i 02-stream-gem5-xaa.trace -r 50 -p$line -s$stride >> 2.txt
    done
done