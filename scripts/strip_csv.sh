#!/bin/bash

# sed -e 's/[\[\]\"]//g' -e 's/,/\\n/g' -e 's/\s\+/,/g' $1
for file in LHCII/**/*; do
  name1=${file/_test/}
  name2=${name1/./}.csv
  mv $file $name2
  # tr -d '[]"' < ${file}_temp | tr ',' '\n' > $file;
  # rm ${file}_temp
done
