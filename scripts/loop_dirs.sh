#!/bin/zsh

for i in structures/monomer/*/* ; do
  if [ -d "$i" ]; then
    ./scripts/run.py -T 300.0 -c 2 -ps 1 -pc 0 -f 0 -i $i
    out=${i/structures/out}
    echo $out
    python scripts/average_spectra.py -c 0 -d 1 -i $out
  fi
done
