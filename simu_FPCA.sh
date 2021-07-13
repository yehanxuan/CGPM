#!/bin/bash


id=2

mkdir -p figs
mkdir -p log

for m in LogDet VNDiv frobDiverg
do
    for d in Fourier
    do 
        for n in 200 500
        do
            for s in Gaussian 
            do 
                for I in LS
                do
                    nohup Rscript ./mainFPCA.R $m $d $id $n $s $I > ./log/simu-$m-$d-$id-$n-$s-$I.log 2>&1 < /dev/null &
                wait
                done
            done
        done
    done
done
            

    
    