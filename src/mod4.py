from mod3 import *
#4.1.1
# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def count_with_pseudocounts(motifs:list):
    """
        Returns a dict with the count of each nucleotide in
        the motifs list using Laplace's rule of Succession
    """
    count_dict = {key : [1]*len(motifs[0]) for key in "ACGT"}
    matrix_len = len(motifs)
    motif_len = len(motifs[0])
    
    for i in range(matrix_len):
        for j in range(motif_len):
            nucleotide = motifs[i][j]
            count_dict[nucleotide][j] += 1
    
    return count_dict

#4.1.2
def profile_with_pseudocounts(motifs:list):
    """
        Takes a list of strings Motifs as input and returns the profile matrix of Motifs 
        with pseudocounts as a dictionary of lists.
        
        Parameters
        ---
            motifs : list<str>
                The list of k-mer motifs.
        
        Returns
        ---
            Dictionary with the pseudocount (adding 1 to each count) of each nucleotide 
            in the motifs list using Laplace's rule of Succession.

    """
    pseudocounts = count_with_pseudocounts(motifs)
    # we add 4 because we are adding 1 per nucleotide (ACTG), according
    # to Laplace's rule of succession 
    divider = len(motifs) + 4 
    for key, value in pseudocounts.items():
            pseudocounts[key] = [x / divider for x in value]

    return pseudocounts

#4.1.3
def greedy_motif_search_with_pseudocounts(dna, k, t):
    """
        Takes a list of strings dna followed by integers k and t and returns the result 
        of running GreedyMotifSearch where each profile matrix is generated with 
        pseudocounts.

        Parameters
        ---
            dna : list<str>
                List of k-mer motifs.
            k : int
                Length of the resulting motifs
            t : int
                Number of the resutling motifs

        Returns
        ---
            List of strings of length k containing the resulting t motifs. 
    """
    best_motifs = [x[0:k] for x in dna]
    
    for i in range(len(dna[0]) - k + 1):
        motifs = [dna[0][i:i+k]]
        for j in range(1, t):
            prof = profile_with_pseudocounts(motifs)
            motifs.append(profile_most_probable_kmer(dna[j], prof, k))

        if score(motifs) < score(best_motifs):
            best_motifs = motifs

    return best_motifs

def test_functions():
    count_with_pseudocounts_input = [
        "AACGTA",
        "CCCGTT",
        "CACCTT",
        "GGATTA",
        "TTCCGG"
    ]

    greedy_motif_search_pseudocounts_input =[ 
        [
        "GGCGTTCAGGCA",
        "AAGAATCAGTCA",
        "CAAGGAGTTCGC",
        "CACGTCAATCAC",
        "CAATAATATTCG"
        ], 3, 5]

    print(count_with_pseudocounts(count_with_pseudocounts_input))
    print(profile_with_pseudocounts(count_with_pseudocounts_input))
    print(greedy_motif_search_with_pseudocounts(
        greedy_motif_search_pseudocounts_input[0],
        greedy_motif_search_pseudocounts_input[1],
        greedy_motif_search_pseudocounts_input[2])
    )

test_functions()