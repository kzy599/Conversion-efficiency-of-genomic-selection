#!/bin/bash
maindir="/home/kangziyi/DRcand/checkwac"
gsmode=("pop")
nr=(0.1 0.3 0.5 0.7 0.9 1)

for g in "${gsmode[@]}"; do
  for l in "${nr[@]}"; do
          cp_path="$maindir/ProDT"
          
          find "$maindir/$g$l" \( -name '*[1-9].csv' -o -name '*10.csv' \) |xargs -i cp -vf {} $cp_path
        done
      done