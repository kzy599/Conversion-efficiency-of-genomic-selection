#!/bin/bash

maindir="/home/kangziyi/DRcand/checkwac"
gsmode=("pop")
nr=(1 0.1 0.3 0.5 0.7 0.9)

for g in "${gsmode[@]}"; do
  for l in "${nr[@]}"; do
          dir_path="$maindir/$g$l"
          cd "$dir_path"
          nohup nice -n 10 Rscript runme.R &
  done
done