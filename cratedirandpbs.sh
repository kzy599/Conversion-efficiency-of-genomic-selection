#!/bin/bash
echo "Make sure you have already set the parameters in this sh and parameters.R"

echo "crate dir"

ptwd="/home/mengxh/kangzy/WssvCode/WssvCodeFinal/" 

cd $ptwd

gsmode="pop"

lowdesenty=(0 3 12 23 46 114 228 455 682 1250)

reference=(10 30 50 70 100 150 200)

IM=("IF" "IT")

FileName="ge.pbs"

basename="BasePoP"

mkdir $basename

cd $basename

content="#!/bin/bash
#PBS -N basePoP
#PBS -l walltime=9999:00:00
#PBS -l nodes=1:ppn=20
#PBS -l mem=70gb
#PBS -q q_fat
#PBS -S /bin/bash
        
cd $ptwd$basename
        
Rscript population.R"

echo "$content" > "$FileName"

cd ..

for i in ${gsmode[@]}
do
  for g in ${lowdesenty[@]}
   do
    for w in ${IM[@]}
      do 
        for r in ${reference[@]}
          do
            if [ "$g" -eq 1250 ] && [ "$w" = "IT" ]; then
             continue
            fi
            if [ "$g" -eq 0 ] && [ "$w" = "IT" ]; then
             continue
            fi
            mkdir $i$g$w$r
            cd $i$g$w$r
            content="#!/bin/bash
#PBS -N $i$g$w$r
#PBS -l walltime=9999:00:00
#PBS -l nodes=1:ppn=20
#PBS -l mem=70gb
#PBS -q q_fat
#PBS -S /bin/bash
        
cd $ptwd$i$g$w$r
        
Rscript runme.R"
        echo "$content" > "$FileName"
        cd ..
        done
      done
    done
done

echo "prepare the Rscript for each dir"

Rscript copyfile.R

echo "dir and script have been crated"

echo "Please run the population.R in BasePoP with nohup, because it will taking many time"
