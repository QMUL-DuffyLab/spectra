#!/bin/bash

# sed -e 's/[\[\]\"]//g' -e 's/,/\\n/g' -e 's/\s\+/,/g' $1
for file in LHCII/**/*; do
  # cp $file ${file}_temp
  # awk '{ print $7, $8, $9 }' ${file}_temp > $file

  dir=${file%/*} # strip trailing filename frame.x
  dir=${dir##*/} # strip leading LHCII directory - get e.g. CLA602
  case $dir in
    *CLA*)
      tresp_file=LHCII/cla_tresp.dat
      ;;
    *CHL*)
      tresp_file=LHCII/chl_tresp.dat
      ;;
    *LUT*)
      tresp_file="NO LUT TRESP FILE"
      ;;
  esac
  echo $dir, $tresp_file
  # rm ${file}_temp
done
