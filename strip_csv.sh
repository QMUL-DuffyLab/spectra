#!/bin/bash

# sed -e 's/[\[\]\"]//g' -e 's/,/\\n/g' -e 's/\s\+/,/g' $1
for file in LHCII/**/*; do
  cp $file ${file}_temp
  tr -d '[]"' < ${file}_temp | tr ',' '\n' > $file;
  rm ${file}_temp
done
