#!/usr/bin/env bash

cd output
cut -f 1,2 conserved_sig.txt | awk '$1 == $1 {print $0}' | cut -f 2 | sed -z 's/\n/,/g' > chr$1_markers.txt
