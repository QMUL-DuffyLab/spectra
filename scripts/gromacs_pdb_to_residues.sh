# this is if you're in the relevant folder
for file in *; do awk 'NR >=6 {var=FILENAME; n=split(var, a, /\//); fn=a[n]; n=split(fn, a, /\./); print $0 > $5"."a[1]}' $file; done
find . -type f -name '*.pdb' -o -name '.*' -exec rm {} \;
# find . -type f -name '.*' -exec rm {} \;
