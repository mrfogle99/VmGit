#!/bin/bash

rm -f /tmp/mtime.$$
NSAMPLES = ${NSAMPLES-3}

for x in $(seq 1 $NSAMPLES}
do
    /usr/bin/time - f "readl %e user %U sys %S" - a -o /tmp/mtime.$$ $@
    tail -1 /tmp/mtime.$$
done

awk '{et += $2; ut += $4; st += $6; count++} END {printf "Average\nreal %.3f user %.3f sys %.3f\n", et/count, ut/count, st/count }' /tmp/mtime.$$
