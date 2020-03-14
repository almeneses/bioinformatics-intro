from typing import List, Dict
#3.3.1
def count_motifs(motifs : List[str]) -> Dict[str, int]:
    """
        Counts the instances per column of each nucleotide in the motifs list.

        Paramters
        ---
            motifs : List[str]
                List of k-mers strings.
            
        Returns
        ---
            Dict[str, int] : The count of each motif in each column of the given 
            motifs list of strings.
    """
    count_dict = {key : [0]*len(motifs[0]) for key in "ACGT"}

    for i in range(len(motifs)):
        for j in range(len(motifs[0])):
            nucleotide = motifs[i][j]
            count_dict[nucleotide][j] += 1
    
    return count_dict

#3.3.2
def profile(motifs : List[str]) -> Dict[str, List[float]]:
    """
        Calculates the profile matrix of the given motifs list of stirngs.

        Parameters
        ---
            motifs : List[str]
                List of k-mers strings.
            
        Returns
        ---
            Dict[str, List[float]] : The calculated profile of motifs, as a dictionary of lists.
    """
    count_dict  = count_motifs(motifs)
    motif_len = len(motifs)
    for key, value in count_dict.items():
        count_dict[key] = [x / motif_len for x in value]
    
    return count_dict

#3.3.3
def consensus(motifs : List[str]) -> str:
    """
        Calculates the consensus string of motifs.

        Parameters
        ---
            motifs : List[str]
                The list of k-mers motifs.

        Returns
        ---
            str : A consensus string of Motifs.
    """
    count_dict = count_motifs(motifs)
    j = len(motifs)
    concensus = ""
    count_dict_keys = count_dict.keys()

    for j in range(len(motifs[0])):
        frequent_nucleotide = ""
        quantity = 0
        for key in count_dict_keys:
            if count_dict[key][j] > quantity:
                quantity = count_dict[key][j]
                frequent_nucleotide = key
        concensus += frequent_nucleotide
    
    return concensus 

#3.3.4
def score(motifs : List[str]) -> int:
    """
        Calculates the total score of the given motifs.

        Parameters
        ---
            motifs : List[str]
                The list of k-mers motifs.

        Returns
        ---
            int : The score from the given list of string motifs.
    """
    most_frequent = consensus(motifs)
    count = count_motifs(motifs)
    elements_per_column = len(motifs)
    score = 0

    for j in range(len(most_frequent)):
        score += elements_per_column - count[most_frequent[j]][j]
        
    return score


# 3.4.1
def pr(text : str, profile : Dict[str, List[float]]) -> float:
    """
        Calculates the probability of the given profile to generate the given text.
        
        Parameters
        ---
            text : str
                The string to calculate the probability.
            profile : Dict[str, List[float]]
                The profile matrix needed for the probability calculation.
        
        Returns
        ---
            float : The probability that profile can generate the given text string.
    """
    p = 1
    for i in range(len(text)):
        p *= profile[text[i]][i]
    return p
    
# 3.4.2
def profile_most_probable_kmer(text : str, profile : Dict[str, List[float]], k : int) -> str:
    """
        Finds the most probable k-mer in a string with the given profile.

        Parameters
        ---
            text : str
                The string to search in the profile most probable k-mer.
            profile : Dict[str, List[float]]
                A profile matrix.
            k : int
                The length of the most probable k-mer.
        
        Returns
        ---
            str : The most probable k-mer in the given text based on the given profile.
    """
    probability = -1
    most_probable_text = ""
    
    for i in range(len(text) - k+1):

        current_pr = pr(text[i:i+k], profile)
        
        if current_pr > probability:
            probability = current_pr
            most_probable_text = text[i:i+k]

    return most_probable_text
    
#3.4.3
def greedy_motif_search(dna : List[str], k : int, t : int) -> List[str]:
    """
        Creates an array of motifs from the given dna string array 
        using a greedy algorithm approach.
        
        Parameters
        ---
            dna : List[str]
                The DNA string to search in.
            k : int
                The length of the resulting k-mers.
            t : int
                The amount of k-mers to return.

        Returns
        ---
            List[str] : Array of t motifs of length k from the
                given dna using a greedy algorithm approach.
    """
    best_motifs = [x[0:k] for x in dna]
    
    for i in range(len(dna[0]) - k + 1):
        motifs = [dna[0][i:i+k]]
        for j in range(1, t):
            prof = profile(motifs)
            motifs.append(profile_most_probable_kmer(dna[j], prof, k))

        if score(motifs) < score(best_motifs):
            best_motifs = motifs

    return best_motifs


def test_functions():

    count_motifs_input = [
    "AACGTA",
    "CCCGTT",
    "CACCTT",
    "GGATTA",
    "TTCCGG"
    ]
   
    score_input = [
    "GTACAACTGT",
    "CAACTATGAA",
    "TCCTACAGGA",
    "AAGCAAGGGT",
    "GCGTACGACC",
    "TCGTCAGCGT",
    "AACAAGGTCA",
    "CTCAGGCGTC",
    "GGATCCAGGT",
    "GGCAAGTACC"
    ]
    
    profile_most_probable_kmer_input = [
        "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT",
        {
        "A": [0.2, 0.2, 0.3, 0.2, 0.3],
        "C": [0.4, 0.3, 0.1, 0.5, 0.1],
        "G": [0.3, 0.3, 0.5, 0.2, 0.4],
        "T": [0.1, 0.2, 0.1, 0.1, 0.2]
        },
        5
    ]

    greedy_motif_search_input = [
        ["GAGGCGCACATCATTATCGATAACGATTCGCCGCATTGCC",
        "TCATCGAATCCGATAACTGACACCTGCTCTGGCACCGCTC",
        "TCGGCGGTATAGCCAGAAAGCGTAGTGCCAATAATTTCCT",
        "GAGTCGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG",
        "GACGGCAACTACGGTTACAACGCAGCAACCGAAGAATATT",
        "TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT",
        "AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG",
        "AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA"], 5, 8
    ]

    print (count_motifs(count_motifs_input))
    print(profile(count_motifs_input))
    print(consensus(count_motifs_input))
    print(score(score_input))
    print(profile_most_probable_kmer(profile_most_probable_kmer_input[0], 
        profile_most_probable_kmer_input[1], profile_most_probable_kmer_input[2]))
    print(greedy_motif_search(greedy_motif_search_input[0], greedy_motif_search_input[1], greedy_motif_search_input[2]))

if __name__ == '__main__':
    test_functions()