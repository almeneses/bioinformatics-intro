import random
from typing import List, Dict
from mod3 import *

#4.1.1
def count_with_pseudocounts(motifs: List[str]) -> Dict[str, List[int]]:
    """
        Creates a dictionary with the count of each nucleotide in
        the motifs list using Laplace's rule of Succession.
        
        Parameters
        ---
            motifs : List[str]
                List of motifs.

        Returns
        ---
           Dict[str, List[int]]:  The count of each nucleotide.
    """
    motif_len = len(motifs[0])
    count_dict = {key : [1]*motif_len for key in "ACGT"}
    matrix_len = len(motifs)

    
    for i in range(matrix_len):
        for j in range(motif_len):
            nucleotide = motifs[i][j]
            count_dict[nucleotide][j] += 1
    
    return count_dict

#4.1.2
def profile_with_pseudocounts(motifs: List[str]) -> Dict[str, List[str]]:
    """
        Takes a list of strings Motifs as input and returns the profile matrix of Motifs 
        with pseudocounts as a dictionary of lists.
        
        Parameters
        ---
            motifs : List[str]
                The list of k-mer motifs.
        
        Returns
        ---
            Dict[str, List[str]] : Dictionary with the pseudocount (adding 1 to each count) of each nucleotide
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
def greedy_motif_search_with_pseudocounts(dna: List[str], k: int, t: int) -> List[str]:
    """
        Takes a list of strings dna followed by integers k and t and returns the result 
        of running GreedyMotifSearch where each profile matrix is generated with 
        pseudocounts.

        Parameters
        ---
            dna : List[str]
                List of k-mer motifs.
            k : int
                Length of the resulting motifs
            t : int
                Number of the resutling motifs

        Returns
        ---
            List[str] : Strings of length k containing the resulting t motifs. 
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
def motifs(profile: Dict[str, List[float]], dna: List[str]) -> List[str]:
    """
        Takes a profile dictionary corresponding to a list of strings dna 
        and returns a list of the profile's most probable k-mers in each
        string from dna.

        Parameters
        ---
            profile : Dict[str, List[float]]
                The profile matrix.

            dna : List[str]
                List of k-mers.
        
        Returns
        ---
            List[str] : The probable k-mers in each string of dna.
    """
    result_motifs = []
    k = len(profile[next(iter(profile))])
    for i in dna:
        result_motifs.append(profile_most_probable_kmer(i, profile, k))
    
    return result_motifs

#4.2.2
def random_motifs(dna: List[str], k: int , t: int) -> List[str]:
    """
        Generates a list of random t motifs of length k each from a list of string dna.

        Parameters
        ---
            dna : List[str]
                List of k-mer motifs.
            k : int
                Length of the resulting motifs
            t : int
                Number of the resutling motifs
        Returns
        ---
            List[str] : randomly picked t motifs of lenght k.

    """
    result = []
    for i in dna:
        random_kmer = random.randint(0, len(dna[0]) - k)
        result.append(i[random_kmer:random_kmer+k])
    
    return result

#4.2.3
def randomized_motif_search(dna: List[str], k: int, t: int) -> List[str]:
    """
        Generates a list of random t motifs of length k each from a list of string dna using
        a cotinuously improving best motifs search.

        Parameters
        ---
            dna : List[str]
                List of k-mer motifs.
            k : int
                Length of the resulting motifs
            t : int
                Number of the resutling motifs
        Returns
        ---
            List[str]: List of best possible t motifs of lenght k.
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
def normalize(probabilities : Dict[str, float]) -> Dict[str, float]:
    """
        Normalizes the probability of each k-mer in the probabilities dictionary. 
        
        Parameters
        ---
        probabilities : dict
            A dictionary of probabilities, where keys are k-mers and values 
            are the probabilities of these k-mers (which do not necessarily sum up to 1).
        
        Returns
        ---
            Dict[str, float]: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities.
    """
    total = sum(probabilities.values())
    return {key : value / total for key, value in probabilities.items()}

#4.4.2
def weighted_dice(probabilities : Dict[str, float]) -> str:
    """
        Randomly chooses a k-mer from the input.

        Parameters
        ---
            probabilities : Dict[str, float]
                A dictionary of probabilities, where keys are k-mers and values 
                are the probabilities of these k-mers (they have to sum up 1).

        Returns
        ---
            str : A randomly chosen k-mer with respect to the values in the probabilities input.
    """

    prob = random.uniform(0, 1)
    total = 0
    for key in probabilities:
        if total <= prob <= (total + probabilities[key]):
            return key
        total += probabilities[key]

#4.4.3
def profile_generated_string(text : str, profile : Dict[str, List[float]], k : int) -> Dict[str, List[float]]:
    """
        Randomly chooses a k-mer from a string 'text' based on a profile matrix 'profile'.

        Parameters
        ---
            text : str
                A string of nucleotides.

            profile : Dict[str, List[float]]
                The profile matrix.

            k : int
                Length of the chosen mer.

        Returns
        ---
            Dict[str, List[float]]: Randomly chosen str of size k.
            
    """
    n = len(text)
    probabilities = {}
    
    for i in range(n-k+1):
        probabilities[text[i:i+k]] = pr(text[i:i+k], profile)
    probabilities = normalize(probabilities)

    return weighted_dice(probabilities)
    
#4.4.4
def gibbs_sampler(dna : List[str], k : int, t : int, n : int = 100) -> List[str]:
    """
        Returns a list of best posible motifs using the Gibbs sampler method.

        Parameters
        ---
            dna : List[str]
                List of dna sequences.
            k : int
                Length of the resulting motifs.
            t : int
                Amount of generated motifs.
            n : int
                Number of iterations for finding the motifs (default 100).
        
        Returns
        ---
            Lis[str] : List of best possible motifs.
    """

    best_motifs = random_motifs(dna, k, t)
    
    for i in range(n):
        
        rand = random.randint(1, t)
        prof = profile_with_pseudocounts(best_motifs)
        motif = motifs(prof, dna)
        
        if score(motif) < score(best_motifs):
            best_motifs = motif
    
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
    weighted_dice_input = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25} 
    profile_generated_string_input = [
        "AAACCCAAACCC",
        {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]},
        2
    ]

    gibbs_sampler_input = [
        [
            "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
            "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
            "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
            "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
            "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
        ], 8, 5, 100]

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
    print(weighted_dice(weighted_dice_input))
    print(profile_generated_string(
        profile_generated_string_input[0], profile_generated_string_input[1], profile_generated_string_input[2]))
    print(gibbs_sampler(gibbs_sampler_input[0], gibbs_sampler_input[1], gibbs_sampler_input[2], gibbs_sampler_input[3]))

if __name__ == '__main__':
    test_functions()