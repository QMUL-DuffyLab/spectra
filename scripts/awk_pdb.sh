#!/bin/bash
file=/home/callum/code/spectra/structures/LHCII_PDB/LHCII/LHCII_trimer/LHCII_trimer_pigments_no_conect.pdb
# awk -v FIELDWIDTHS="7 5 5 4 2 8 8 8 8 5 5" '
awk -v FIELDWIDTHS="7 5 5 5 10 8 7 9 5 6 5" '
{
if ($1 == "HETATM ") {
# if ($1 == "ATOM   ") {
  # the nitrogens are reported as N A, N B, etc.: the following two lines
  # delete external whitespace but keep the internal one
  # gsub(/^[ \t]+/, "", $3)
  # gsub(/[ \t]+$/, "", $3)
  # but we need to delete the internal space as well because of how the
  # TrESP data is spelled: just uses NA, NB, etc.
  # gsub(/ /, "", $3)

  # $4 and $6 used to construct filename
  # $7,8,9 are coordinates: remove all whitespace from these
  gsub(/ /, "", $4)
  gsub(/ /, "", $5)
  gsub(/ /, "", $6)
  gsub(/ /, "", $8)
  gsub(/ /, "", $9)
  # print $3, $7, $8, $9 > $4$6".csv"
  # file_files.py expects all the PDB columns
  # need to write fields manually because printing $0 would
  # ignore the gsub commands
  print $1, $2, $3, $4, $5, $6, $7, $8, $9 > $4$5".csv"
}
}' $file
