#!/bin/bash

for METH in wilks fc
do
    for TRUTH in nh ih
    do
        for ASSUME in nh ih
        do
            DESC=${METH}_${TRUTH}_${ASSUME}
            echo $DESC
            ./fc $METH $TRUTH $ASSUME > ${DESC}.txt
            echo Plotting...
            root -b -q plot.C'("'$DESC'.txt", "plots/cov_'$DESC'.pdf", "plots/crit_'$DESC'.pdf")'
        done
    done
done
