#!/bin/zsh

for dir in *_CSV; do
  code=$(echo $dir | cut -f1 -d_)
  pdbfile=$(echo $dir/**.pdb)
  if [[ $code =~ "CL" || $code =~ "KC" ]]
  then
    s1file=$(echo $dir/**Qy*)
    s2file=$(echo $dir/**Qx*)
    # echo "Code = $code, PDB file = $pdbfile, S1 file = $s1file, S2 file = $s2file"
  else
    s1file=$(echo $dir/**S1*)
    s2file=$(echo $dir/**S2*)
  fi
  # awk adds a blank line at the start, not 100% sure why, so pipe through tail
  awk 'FNR==NR{a[NR]=$7 FS $8 FS $9; next} (FNR>1 && FNR-1 in a){ print a[FNR-1], $3 }' $pdbfile $s1file | tail -n +2 > $dir/frame1.csv
  awk 'FNR==NR{a[NR]=$7 FS $8 FS $9; next} (FNR>1 && FNR-1 in a){ print a[FNR-1], $3 }' $pdbfile $s2file | tail -n +2 > $dir/frame2.csv
done
