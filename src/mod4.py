import random
from mod3 import *
#4.1.1
# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def count_with_pseudocounts(motifs: list) -> dict:
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
def profile_with_pseudocounts(motifs: list) -> dict:
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
def greedy_motif_search_with_pseudocounts(dna: list, k: int, t: int) -> list:
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
            list : Strings of length k containing the resulting t motifs. 
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

#4.2.1
def motifs(profile: dict, dna: list) -> list:
    """
        Takes a profile dictionary corresponding to a list of strings dna 
        and returns a list of the profile's most probable k-mers in each
        string from dna.

        Parameters
        ---
            profile : dict(list(str))
                Profile matrix.
            dna : list<str>
                List of k-mers.
        
        Returns
        ---
            list : The given profile's most probable k-mers in each string of dna.
    """
    result_motifs = []
    k = len(profile[next(iter(profile))])
    for i in dna:
        result_motifs.append(profile_most_probable_kmer(i, profile, k))
    
    return result_motifs

#4.2.2
def random_motifs(dna: list, k: int , t: int) -> list:
    """
        Generates a list of random t motifs of length k each from a list of string dna.

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
        list : randomly picked t motifs of lenght k.

    """
    result = []
    for i in dna:
        random_kmer = random.randint(0, len(dna[0]) - k)
        result.append(i[random_kmer:random_kmer+k])
    
    return result

#4.2.3
def randomized_motif_search(dna: list, k: int, t: int) -> list:
    """
        Generates a list of random t motifs of length k each from a list of string dna using
        a cotinuously improving best motifs search.

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
            List of best possible t motifs of lenght k.
    """

    mot = best_motifs = random_motifs(dna, k, t)
    while True:
        profile = profile_with_pseudocounts(mot)
        mot = motifs(profile, dna)
        if score(mot) < score(best_motifs):
            best_motifs = mot
        else:
            return best_motifs

#4.4.1
def normalize(probabilities : dict) -> dict:
    """
        Normalizes the probability of each k-mer in the probabilities dictionary. 
        
        Parameters
        ---
        probabilities : dict
            A dictionary of probabilities, where keys are k-mers and values 
            are the probabilities of these k-mers (which do not necessarily sum up to 1).
        
        Returns
        ---
            A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities.
    """
    total = sum(probabilities.values())
    return {key : value / total for key, value in probabilities.items()}


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

    motifs_input = [
        {"A":[0.8, 0.0, 0.0, 0.2],
        "C":[0.0, 0.6, 0.2, 0.0],
        "G":[0.2, 0.2, 0.8, 0.0],
        "T":[0.0, 0.2, 0.0, 0.8]},
        ["TTACCTTAAC",
        "GATGTCTGTC",
        "ACGGCGTTAG",
        "CCCTAACGAG",
        "CGTCAGAGGT"]
    ]

    random_motifs_input = [[

        "TTACCTTAAC",
        "GATGTCTGTC",
        "ACGGCGTTAG",
        "CCCTAACGAG",
        "CGTCAGAGGT"
    ], 3, 5]
    randomized_motif_search_input = [
        [
        "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
        "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
        "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
        "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
        "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
        ], 
        8, 
        5
    ]

    normalize_input = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}

    print(count_with_pseudocounts(count_with_pseudocounts_input))
    print(profile_with_pseudocounts(count_with_pseudocounts_input))
    print(greedy_motif_search_with_pseudocounts(
        greedy_motif_search_pseudocounts_input[0],
        greedy_motif_search_pseudocounts_input[1],
        greedy_motif_search_pseudocounts_input[2])
    )
    print(motifs(motifs_input[0], motifs_input[1]))
    print(random_motifs(random_motifs_input[0], random_motifs_input[1], random_motifs_input[2]))
    print(randomized_motif_search(
        randomized_motif_search_input[0], randomized_motif_search_input[1], randomized_motif_search_input[2]))
    
    print(normalize(normalize_input))

if __name__ == '__main__':
    test_functions()