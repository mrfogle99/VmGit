#!/bin/bash
# Author: David Holmqvist <daae19@student.bth.se>

echo "NOTE: this script relies on the binaries pearson and pearson_par to exist"

status=0
temp_dir=$(mktemp -d 2> /dev/null || mktemp -d -t 'mytmpdir')
red=$(tput setaf 1)
reset=$(tput sgr0)

for thread in 1 2 4 8 16 32
do
    for size in 128 256 512 1024
    do
        ./pearson "data/$size.data" "$temp_dir/${size}_seq.data"
        ./pearson_par "data/$size.data" "$temp_dir/${size}_par.data" $thread

        if ! cmp -s "$temp_dir/${size}_seq.data" "$temp_dir/${size}_par.data"
        then
            echo "${red}Incongruent output data detected with size $size and $thread thread(s)${reset}"
            status=1
        fi

        rm "$temp_dir/${size}_seq.data" "$temp_dir/${size}_par.data"
    done
done

exit $status