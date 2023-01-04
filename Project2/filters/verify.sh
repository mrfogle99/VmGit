#!/bin/bash

success=1
red=$(tput setaf 1)
green=$(tput setaf 2)
reset=$(tput sgr0)

if ! [ $# -eq 2 ]
then
    echo "${red}Usage: $0 [optimized_blur_executable] [optimized_threshold_executable]$reset" && exit 1
fi

echo "NOTE: The provided optimized executables must be named differently compared to the unoptimised 'blur' and 'threshold' executables."

function check_equal() {
    if ! cmp -s $3 $4
    then
        echo "${red}Incorrect $2 transformation ($1)${reset}"
        success=0
    else
        echo "${green}Correct $2 transformation ($1)$reset"
    fi
}

files=(im1.ppm im2.ppm im3.ppm im4.ppm)
filters=(blur threshold)
user_filters=($1 $2)
radius=15

for filter in ${filters[@]} ${user_filters[@]}
do
    ! test -x $filter && echo "$red'$filter' is not an executable file$reset" && exit 1
done

for file in ${files[@]}
do
    ./blur $radius data/$file data/blurred_$file
    $1 $radius data/$file data/user_blurred_$file
    check_equal $file "blur" data/blurred_$file data/user_blurred_$file
    rm data/blurred_$file data/user_blurred_$file
    ./threshold data/$file data/thresholded_$file
    $2 data/$file data/user_thresholded_$file
    check_equal $file "threshold" data/thresholded_$file data/user_thresholded_$file
    rm data/thresholded_$file data/user_thresholded_$file
done

echo
echo "======= Verdict ======="

if [ $success -eq 1 ]
then
    echo "${green}All good$reset"
    exit 0
else
    echo "${red}Incorrect transformations encountered${reset}"
    exit 1
fi
