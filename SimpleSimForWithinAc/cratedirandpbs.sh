#!/bin/bash
echo "Make sure you have already set the parameters in this sh and parameters.R"

echo "crate dir"

# Set absolute paths
maindir="/home/kangziyi/DRcand/checkwac"

# Define arrays
gsmode=("pop")
nr=(1 0.1 0.3 0.5 0.7 0.9)

basename="BasePoP"
dtname="ProDT"

# Create main directory
mkdir -p "$maindir"
mkdir -p "$maindir/$basename"
mkdir -p "$maindir/$dtname"
# Create subdirectories
for g in "${gsmode[@]}"; do
  for l in "${nr[@]}"; do
          dir_path="$maindir/$g$l"
          mkdir -p "$dir_path"
  done
done

echo "prepare the Rscript for each dir"

Rscript copyfile.R

echo "dir and script have been crated"

echo "Please run the population.R in BasePoP with nohup, because it will taking many time"
