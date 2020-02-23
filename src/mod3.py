#3.1
# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def count_motifs(motifs):

    count_dict = create_motif_count_dict(len(motifs[0]))
    matrix_len = len(motifs)
    motif_len = len(motifs[0])
    
    for i in range(matrix_len):
        for j in range(motif_len):
            nucleotide = motifs[i][j]
            count_dict[nucleotide][j] += 1
    
    return count_dict

#3.2
# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def profile(motifs):
    count_dict  = count_motifs(motifs)
    motif_len = len(motifs)
    for key, value in count_dict.items():
        count_dict[key] = [x / motif_len for x in value]
    
    return count_dict


# Helper function
# Input: the lenght of the strings in the motif matrix
# Output: Dictionary of list of zeroes, each one the size of length
def create_motif_count_dict(length):
    
    size = range(length)
    return {
        "A": [0 for x in size],
        "T": [0 for x in size],
        "C": [0 for x in size],
        "G": [0 for x in size],
    }



def test_functions():

    count_motifs_input = ["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]
    print (count_motifs(count_motifs_input))

    print(profile(count_motifs_input))

test_functions()