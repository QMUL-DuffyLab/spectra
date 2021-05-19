#!/bin/zsh

# first argument: directory to look at
# second argument: frame
# third: output prefix i.e. what protein is this

mkdir -p out/long_short_frame_comp/"${3}"/"${2}"

for file in "${1}"*/frame"${2}".csv; do
  pigment=$(basename $(dirname $file))
  print $file, out/long_short_frame_comp/"${3}"/"${2}"/$pigment.csv
  cp $file out/long_short_frame_comp/"${3}"/"${2}"/$pigment.csv
done
