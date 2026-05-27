#!/usr/bin/env bash

cd output
grep -v chr *$1*assoc* | awk '$11 < 1e-06 {print $0}' > summary_$1.txt
cut -f 1,2,3 summary_more_$1.txt | cut -d ':' -f 2 | sort| uniq > intermediate_more_$1.txt
