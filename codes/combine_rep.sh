#!/usr/bin/env bash

cd output
comm -12 intermediate_more_beta.txt intermediate_more_log2beta.txt | comm -12 - intermediate_more_quantile_beta.txt | comm -12 - intermediate_more_M.txt > intermediate_more.txt
