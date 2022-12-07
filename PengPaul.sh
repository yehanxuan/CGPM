#!/bin/bash
id=2
mkdir -p Penglog

for n in 50
do 
    for m in LogDet EM REML
    do
        for d in easySin easy 
        do 
            for s in Gaussian
            do
                for I in EM LS LOC
                do
                    for e in Gaussian t uniform
                    do
                        nohup Rscript ./mainFPCA.R $m $d $id $n $s $I $e > ./Penglog/simu-$m-$d-$id-$n-$s-$I-$e.log 2>&1 < /dev/null &
			        wait
			        done
                done
            done
        done
    done
done
			
                    
