#3.3.1
# Count
# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def count_motifs(motifs):

    count_dict = {key : [0]*len(motifs[0]) for key in "ACGT"}
    matrix_len = len(motifs)
    motif_len = len(motifs[0])
    
    for i in range(matrix_len):
        for j in range(motif_len):
            nucleotide = motifs[i][j]
            count_dict[nucleotide][j] += 1
    
    return count_dict

#3.3.2
# Profile
# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def profile(motifs):

    count_dict  = count_motifs(motifs)
    motif_len = len(motifs)
    for key, value in count_dict.items():
        count_dict[key] = [x / motif_len for x in value]
    
    return count_dict

#3.3.3
# Consensus
# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.

def consensus(motifs):
    
    count_dict = count_motifs(motifs)
    j = len(motifs)
    concensus = ""

    for j in range(len(motifs[0])):
        frequent_nucleotide = ""
        quantity = 0
        for key in count_dict.keys():
            if count_dict[key][j] > quantity:
                quantity = count_dict[key][j]
                frequent_nucleotide = key
        concensus += frequent_nucleotide
    
    return concensus        

def score(motifs):

    most_frequent = consensus(motifs)
    count = count_motifs(motifs)
    elements_per_column = len(motifs)
    score = 0

    for j in range(len(most_frequent)):
        score += elements_per_column - count[most_frequent[j]][j]
        
    return score


# 3.4.1
def pr(text, profile):
    """
        Returns the probability that profile
        generates the text string.
    """
    p = 1
    for i in range(len(text)):
        p *= profile[text[i]][i]
    return p
    
# 3.4.4
# Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
# Input: A string text, an integer k, and a 4 x k matrix profile.
# Output: A Profile-most probable k-mer in Text.
def profile_most_probable_kmer(text, profile, k):
    
    probability = -1
    most_probable_text = ""
    
    for i in range(len(text) - k+1):

        current_pr = pr(text[i:i+k], profile)
        
        if current_pr > probability:
            probability = current_pr
            most_probable_text = text[i:i+k]

    return most_probable_text
    

def test_functions():

    count_motifs_input = [
    "AACGTA",
    "CCCGTT",
    "CACCTT",
    "GGATTA",
    "TTCCGG"
    ]
    print (count_motifs(count_motifs_input))
    print(profile(count_motifs_input))
    print(consensus(count_motifs_input))


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
    print(score(score_input))

    #print(profile_most_probable_kmer("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", {
    #    "A": [0.2, 0.2, 0.3, 0.2, 0.3],
    #    "C": [0.4, 0.3, 0.1, 0.5, 0.1],
    #    "G": [0.3, 0.3, 0.5, 0.2, 0.4],
    #    "T": [0.1, 0.2, 0.1, 0.1, 0.2]}, 5)
    #)

    print(profile_most_probable_kmer("AACCGGTT",{
    "A":[1.0, 1.0, 1.0],
    "C":[0.0, 0.0, 0.0],
    "G":[0.0, 0.0, 0.0],
    "T":[0.0, 0.0, 0.0]}, 3))
test_functions()