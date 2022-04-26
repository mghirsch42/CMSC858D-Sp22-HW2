#/bin/bash

data_path="data/"
output_path="suffix_arrays/"
time_file="results/build_times.csv"

for r_len in 1000 10000 100000
do
    ref_file="${data_path}ref_${r_len}.fa"
    for k in 0 5 10
    do
        for mode in "naive" "simpaccel"
        do
            output_file="${output_path}rlen_${r_len}_k_${k}_mode_${mode}.dat"
            ./buildsa "--preftab" "${k}" "${ref_file}" "${output_file}" "--query_mode" "${mode}" "--time_fname" "${time_file}"
        done
    done
done