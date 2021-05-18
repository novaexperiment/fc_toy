#!/bin/bash

mkdir -p plots

for METH in wilks fc hc prof post
do
    for TRUTH in nh ih
    do
        for FIT in nh ih either
        do
            # FIT is unused for the profile method
            if [[ $METH == prof && $FIT != either ]]; then continue; fi
            # Ditto post
            if [[ $METH == post && $FIT != either ]]; then continue; fi

            for MOCK in nh ih
            do
                # MOCK is only used for the FC method
                if [[ $METH != fc && $MOCK != nh ]]; then continue; fi

                DESC=${METH}_${TRUTH}
                if [[ $METH != prof && $METH != post ]]; then DESC=${DESC}_${FIT}; fi
                if [[ $METH == fc ]]; then DESC=${DESC}_${MOCK}; fi

                echo $DESC
                ./fc $METH $TRUTH $FIT $MOCK > ${DESC}.txt
                echo Plotting...
                root -b -q -l plot.C'("'$DESC'.txt", "plots/cov_'$DESC'.pdf", "plots/crit_'$DESC'.pdf")'
            done
        done
    done
done
