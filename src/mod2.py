from mod1 import PatternCount

#2.1
def symbol_array(genome, symbol):
    result = {}
    genome_len = len(genome)
    extended_genome = genome + genome[0:genome_len // 2]

    for i in range(genome_len):
        result[i] = PatternCount(extended_genome[i:i + (genome_len // 2)], symbol)
    
    return result

#2.1.1
def improved_symbol_array(genome, symbol):
    result = {}
    genome_len = len(genome)
    window_size = genome_len // 2
    extended_genome = genome + genome[0:window_size]

    result[0] = PatternCount(extended_genome[0:window_size], symbol)

    for i in range(1, genome_len):
        result[i] = result[i-1]

        if (extended_genome[i-1] == symbol):
            result[i] -= 1
        if (extended_genome[i + window_size - 1] == symbol):
            result[i] += 1
    
    return result


def test_functions():
    ecoli_genome = open("../decoded-genomes/escherichia_coli.txt", "r").read()
    #this function takes forever to execute due
    #to intentional bad optimization
    #print(symbol_array(ecoli_genome, "C")) 
    print(improved_symbol_array("AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", "CC"))
    #print (improved_symbol_array(ecoli_genome, "C"))

test_functions()