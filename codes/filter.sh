#!/usr/bin/env bash

cd output
cut -f 1,2 summary_more_$1.txt | cut -d ':' -f 2 | sort| uniq > intermediate_more_$1.txt
