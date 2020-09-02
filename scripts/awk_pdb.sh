#!/bin/bash
file=/home/callum/code/spectra/PDB_structures/LHCII/LHCII_monomer/LHCII_monomer_chlorophyll.pdb
awk -v FIELDWIDTHS="7 5 5 4 2 9 7 9 8 5 5" '
{
if ($1 == "HETATM ") {
  # delete leading and trailing whitespace from the atom labels
  gsub(/^[ \t]+/, "", $3)
  gsub(/[ \t]+$/, "", $3)
  # $4 and $6 used to construct filename
  # $7,8,9 are coordinates: remove all whitespace from these
  gsub(/ /, "", $4)
  gsub(/ /, "", $6)
  gsub(/ /, "", $8)
  gsub(/ /, "", $9)
  print $3, $7, $8, $9 > $4$6".csv"
}
}' $file
