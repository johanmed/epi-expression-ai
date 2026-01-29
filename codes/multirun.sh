NUM_TRAITS=$(less $1_phenofile.txt | awk -F, '{print NF}' | tail -n 1)

for num in $(seq 1 "$NUM_TRAITS"); do
    output="$1_${num}"
    gemma -g new_BXD.8_geno.txt -p $1_phenofile.txt -n ${num} -a BXD.8_snps.txt -k output/rbxd.cXX.txt -lmm 9 -o ${output}
done
