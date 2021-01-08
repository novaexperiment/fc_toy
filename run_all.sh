#!/bin/bash

for METH in wilks fc hc prof
do
    for TRUTH in nh ih
    do
        for FIT in nh ih either
        do
            for MOCK in nh ih
            do
                DESC=${METH}_${TRUTH}_${FIT}_${MOCK}
                echo $DESC
                ./fc $METH $TRUTH $FIT $MOCK > ${DESC}.txt
                echo Plotting...
                root -b -q -l plot.C'("'$DESC'.txt", "plots/cov_'$DESC'.pdf", "plots/crit_'$DESC'.pdf")'
            done
        done
    done
done
