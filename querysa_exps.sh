#/bin/bash

data_path="data/"
sa_path="suffix_arrays/"
output_path="query_results/"
time_file_prefix="results/query_times"

for r_len in 1000 10000 100000
do
    for k in 0 5 10
    do
        for mode in "naive" "simpaccel"
        do
            for q_len in 10 100
            do
                for q_n in 100 1000 10000
                do
                    for q_type in "ref" "rand"
                    do
                        index_file="${sa_path}rlen_${r_len}_k_${k}_mode_${mode}.dat"

                        echo "r_len: ${r_len}, k: ${k}, mode: ${mode}, q_len: ${q_len}, q_n: ${q_n}, q_type: ${q_type}"                    
                        query_file="${data_path}query_rlen_${r_len}_qlen_${q_len}_qn_${q_n}_${q_type}.fa"
                        output_file="${output_path}rlen_${r_len}_k_${k}_mode_${mode}_qlen_${q_len}_qn_${q_n}_qtype_${q_type}.txt"
                        ./querysa "${index_file}" "${query_file}" "${mode}" "${output_file}" "--time_fname" "${time_file_prefix}_${q_type}.csv"
                    done
                done
            done
        done
    done
done