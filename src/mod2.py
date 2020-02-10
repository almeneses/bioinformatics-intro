from mod1 import PatternCount

def symbol_array(genome, symbol):
    result = {}
    genome_len = len(genome)
    extended_genome = genome + genome[0:genome_len // 2]

    for i in range(genome_len):
        result[i] = PatternCount(extended_genome[i:i + (genome_len // 2)], symbol)
    
    return result
