import os
from typing import Dict, List
from mod1 import pattern_count


#2.1
def symbol_array(genome : str, symbol : str) -> List[str]:
    """
        Finds all the positions where the given symbol is present in the given genome.

        Parameters
        ---
            genome : str
                The genome string to search in.
            symbol : str
                The symbol string to search for.
        Returns
        ---
            List[str] : List of all positions where the symbol appears in the genome.
    """

    result = {}
    window_size = len(genome) // 2
    extended_genome = genome + genome[0:window_size]

    for i in range(len(genome)):
        result[i] = pattern_count(extended_genome[i:i + window_size], symbol)
    
    return result

#2.1.1
def improved_symbol_array(genome : str, symbol : str) -> List[str]:

    """
        Finds all the positions where the given symbol is present in the given genome.
        (Improved, faster version).

        Parameters
        ---
            genome : str
                The genome string to search in.
            symbol : str
                The symbol string to search for.
        Returns
        ---
            List[str] : List of all positions where the symbol appears in the genome.
    """
    result = {}
    window_size = len(genome) // 2
    extended_genome = genome + genome[0:window_size]
    result[0] = pattern_count(extended_genome[0:window_size], symbol)

    for i in range(1, len(genome)):

        result[i] = result[i-1]

        if (extended_genome[i-1] == symbol):
            result[i] -= 1
        if (extended_genome[i + window_size - 1] == symbol):
            result[i] += 1
    
    return result

#2.4.1
def skew_array(genome : str) -> List[int]:
    """
        Finds the skew array of the genome (string). The Skew array
        is an array where skew[i] = (ocurrences of G) - (ocurrences of C), 
        skew[0] is always 0. 

        Parameters
        ---
            genome : str
                the genome string to generate the skew array from.
        
        Returns
        ---
            List[int] : Skew array with the counts of # of G's minus # of C's.
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
def minimum_skew(genome : str) -> List[int]:
    """
        Finds a position in a genome where the skew diagram attains a minimum.

        Parameters
        ---
            genome : str
                A DNA string.
            
        Returns
        ---
            List[int] : All the positions in genome where skew is minimum.
    """

    skew_arr = skew_array(genome)
    skew_min = min(skew_arr)
    return [ index for index, value in enumerate(skew_arr) if value == skew_min ]

#2.5.1
def hamming_distance(p : str, q : str) -> int:
    """
        Calculates the number of different characters between two strings of equal length (Hamming distance).
        
        Parameters
        ---
            p : str
                String to be compared.
            q : str
                String to be compared.
                
        Returns
        ---
         int : The Hamming distance between the given strings.
    """
    if len(p) != len(q):
        return

    distance = 0

    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance

# 2.5.2
def aprox_pattern_matching(genome : str, pattern : str, d : int) -> List[int]:
    """
        Finds  the starting positions of the approximate occurrences of a pattern 
        in a genome string.

        Parameters
        ---
            genome : str
                The genome string to search in.
            pattern : str
                The pattern string to search for.
            d : int
                The maximum number of mismatches allowed for a pattern to be
                considered acceptable.

        Returns
        ---
            List[int] : Starting positions where pattern is a substring of genome
                with maximum d mismatches.
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
def aprox_pattern_count(genome : str, pattern : str, d : int) -> int:
    """
        Counts all the approximate ocurrences of a pattern in a genome string.

        Parameters
        ---
            genome : str
                The genome string to search in.
            pattern : str
                The pattern string to search for.
            d : int
                The maximum number of mismatches allowed for a pattern to be
                considered acceptable.

        Returns
        ---
            int : The amount of times the given pattern is present in the genome
                with maximum d mismatches.
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
    ecoli_file_path = os.path.join(os.path.dirname(__file__), "../decoded-genomes/escherichia_coli.txt")
    with open(ecoli_file_path, "r") as ecoli_file:
        ecoli_genome = ecoli_file.read()
    
    #this function takes forever to execute due
    #to intentional bad optimization
    #print(symbol_array(ecoli_genome, "C")) 
    
    print (improved_symbol_array(ecoli_genome, "C"))
    print (skew_array("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))
    print (minimum_skew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))
    print(aprox_pattern_matching("CCAAATCCCCTCATGGCATGCATTCCCGCAGTATTTAATCCTTTCATTCTGCATATAAGTAGTGAAGGTATAGAAACCCGTTCAAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGCGCCATAATCCAAACA", "AATCCTTTCA", 3))
    print (aprox_pattern_count("TTTAGAGCCTTCAGAGG", "GAGG", 2))
    

if __name__ == '__main__':
    test_functions()