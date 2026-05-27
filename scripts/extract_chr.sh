#!/usr/bin/env bash

cd output
n=$1
cut -f 1,2,3 conserved_sig.txt | awk '$1 == $n {print $0}' | cut -f 3 | sed -z 's/\n/,/g' > chr${n}_pos.txt
