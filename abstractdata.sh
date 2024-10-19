#!/bin/bash
ptwd="/home/mengxh/kangzy/WssvCode/WssvCodeFinal/" 
cd $ptwd

result="output"

mkdir $result

gsmode="Pop"

lowdesenty=(0 3 12 23 46 114 228 455 682 1250)

reference=(10 30 50 70 100 150 200)

IM=("IF" "IT")

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
            
         newden=$((g * 10))

         find -name $i$newden$w$r*.csv |xargs -i cp -vf {} $ptwd$result
        
        done
      done
    done
done