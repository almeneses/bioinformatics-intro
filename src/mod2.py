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

#2.4.1
def skew_array(genome):
    """
        Returns the skew array of the genome (string) 
    """

    result = [0]
    genome_len = len(genome)

    for i in range (genome_len):

        if genome[i] == "A" or genome[i] == "T":
            result.append(result[i])
        
        if genome[i] == "G":
            result.append(result[i] + 1)
        
        if genome[i] == "C":
            result.append(result[i] - 1)
    
    return result
    
#2.4.2
# Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.
# Input: A DNA string Genome.
# Output: All integer(s) i minimizing Skew[i] among all values of i (from 0 to len(Genome)).

def minimum_skew(genome):
    skew_arr = skew_array(genome)
    skew_min = min(skew_arr)
    return [ index for index, value in enumerate(skew_arr) if value == skew_min ]

#2.5.1
# Hamming Distance Problem: Compute the Hamming distance between two strings.
# Input: Two strings of equal length.
# Output: The Hamming distance between these strings.

def hamming_distance(p, q):
    """
        Returns the number of diferent characters (distance) between strings
        p and q
    """
    if len(p) != len(q):
        return

    distance = 0

    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance

# 2.5.2
# Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
# Input: Strings pattern and genome along with an integer d.
# Output: All starting positions where pattern appears as a substring of genome with at most d mismatches.

def aprox_pattern_matching(genome, pattern, d):
    """
        Returns a list of starting positions where pattern is a substring of genome
        with at most d mismatches
    """
    pattern_len = len(pattern)
    mismatch_count = 0
    positions = []

    for i in range(len(genome) - pattern_len + 1):
        
        mismatch_count = hamming_distance(genome[i : i + pattern_len], pattern)
        
        if mismatch_count <= d:
            positions.append(i)

    return positions

#2.5.3
# Approximate Pattern Count Problem: Count all the approximate ocurrences of a pattern in a string.
# Input: Strings genome and pattern along with an integer d.
# Output: The count of times the pattern is present in the genome with a difference of max d nucleotides.

def aprox_pattern_count(genome, pattern, d):
    """
        Returns the count of times pattern is present in genome
        with at most d mismatches
    """
    pattern_len = len(pattern)
    mismatch_count = 0
    count = 0

    for i in range(len(genome) - pattern_len + 1):
        
        mismatch_count = hamming_distance(genome[i : i + pattern_len], pattern)
        
        if mismatch_count <= d:
            count += 1

    return count    

def test_functions():
    with open("../decoded-genomes/escherichia_coli.txt", "r") as ecoli_file:
        ecoli_genome = ecoli_file.read()
    
    #this function takes forever to execute due
    #to intentional bad optimization
    
    #print(symbol_array(ecoli_genome, "C")) 
    
    #print (improved_symbol_array(ecoli_genome, "C"))
    print (skew_array("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))
    print (minimum_skew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))
    print(aprox_pattern_matching("CCAAATCCCCTCATGGCATGCATTCCCGCAGTATTTAATCCTTTCATTCTGCATATAAGTAGTGAAGGTATAGAAACCCGTTCAAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGCGCCATAATCCAAACA", "AATCCTTTCA", 3))
    print (aprox_pattern_count("TTTAGAGCCTTCAGAGG", "GAGG", 2))
    
test_functions()