import random

output_folder = "data/"
alph = ["A", "C", "G", "T"]

reference_lengths = [1000, 10000, 100000]
query_lengths = [10, 100]
n_queries = [100, 1000, 10000]

for r_len in reference_lengths:
    reference = "".join(random.choice(alph) for i in range(r_len))
    ref_fname = output_folder+"ref_"+str(r_len)+".fa"
    with open(ref_fname, "w") as f:
        f.writelines([">ref_"+str(r_len)+"\n", reference+"\n"])
    for rand_query in [True, False]:
        for query_len in query_lengths:
            for n_query in n_queries:
                queries = []
                for n in range(n_query):
                    if rand_query:
                        query = "".join(random.choice(alph) for i in range(query_len))
                    else:
                        index = random.randint(0, len(reference)-1-query_len)
                        query = reference[index : index+query_len]
                    queries.append(query)
                if rand_query:
                    q_fname = output_folder+"query_rlen_"+str(r_len)+"_qlen_"+str(query_len)+"_qn_"+str(n_query)+"_rand.fa"
                else:
                    q_fname = output_folder+"query_rlen_"+str(r_len)+"_qlen_"+str(query_len)+"_qn_"+str(n_query)+"_ref.fa"
                with open(q_fname, "w") as f:
                    for n in range(len(queries)):
                        f.write(">" + str(n) + "\n")
                        f.write(queries[n] + "\n")