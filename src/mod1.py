from typing import Dict, List

#1.1
def pattern_count(text : str, pattern : str) -> int:
    """ 
		Counts the times the given pattern is present in text.
		
		Parameters
		---
			text : str
				The text to search in.
			pattern : str
				The pattern to search in text.
		
		Returns
		---
			int : How many times the pattern string is present in text.
	"""
    count = 0
    pattern_len = len(pattern)
    text_len = len(text)

    for i in range(text_len - pattern_len + 1):
        
        if(text[i : i + pattern_len] == pattern):
            count += 1
    
    return count

#1.2
def pattern_count_kmer(text : str, k : int) -> str:
	
	"""
		Returns the most frequent k-mer in text.
		
		Parameters
		---
			text : str
				Nucleotide equence to search in.
			k : int
				Length of the mer.

		Returns
		---
			int : Most frequent k-mer.
	"""
	result = ''
	count = 0
	for i in range(len(text) - k + 1):
		step = i + k
		mer = text[i:step]
		match_count = 1
		
		for j in range(step, len(text) - k +1):

			if ( text[j:j+k] == mer ):
				match_count += 1

		if( match_count > count):
			result = mer
			count = match_count
	
	return result


#1.3.1
def frequency_map(text : str, k : int) -> Dict[str, int]:
	"""
		Creates a frequency map whose key is a k-mer and its
		value is the times that mer is present in text.

		Parameters
		---
			text : str
				Nucleotide string to search for k-mers.
			k : int
				Length of the mers in the frequency map.
		
		Returns
		---
			Dict[str, int] : Dictionary of k-mers and the times they are present in text.
	"""

	text_len = len(text)
	freq_map = {}

	for i in range(text_len - k + 1):
		pattern = text[i:i+k]
		if (pattern in freq_map):
			freq_map[pattern] += 1
		else:
			freq_map[pattern] = 1

	return freq_map

#1.3.2
def frequent_words(text, k):
	"""
		Creates a frequency map whose key is a k-mer and its
		value is the times that mer is present in text.

		Parameters
		---
			text : str
				Nucleotide string to search for k-mers.
			k : int
				Length of the mers in the frequency map.
		
		Returns
		---
			Dict[str, int] : Dictionary of k-mers and the times they are present in text.
	"""
	freq_map = frequency_map(text, k)
	most_freq = max(freq_map.values())
	result = []
	
	for key in freq_map:
		if (freq_map[key] == most_freq):
			result.append(key)
	
	return result

#1.4.1
def reverse(pattern):
	"""
		Returns the reverse of the string pattern.

		Parameters
		---
			pattern : str
				String of nucleotides to be reversed.
		
		Returns
		---
			str : String of nucleotides reversed.
	"""
	return pattern[::-1]

#1.4.2
def complement(pattern):
	"""
		Returns the complement of the nucleotide(pattern).

		Parameters
		---
			pattern : str
				String of nucleotides.
		
		Returns
		---
			str : The complement of the given string of nucleotides.

	"""
	comp_map = {"T":"A", "A":"T", "C":"G", "G":"C"}
	result = ""

	for char in pattern:
		result += comp_map[char]
	
	return result

#1.4.3
def reverse_complement(pattern):
	"""
		Returns the reverse complement of a nucleotide string.

		Parameters
		---
			pattern : str
				String of nucleotides.

		Returns
		---
			str : The reverse complement of the given string of nucleotides.
	"""
	return reverse(complement(pattern))

#1.4.4
def pattern_match_positions(pattern : str, genome : str) -> List[int]:
	"""
		Returns a list with the positions in genome where
		pattern starts.

		Parameters
		---
			pattern : str
				The pattern to find.
			genome : str
				The genome string where pattern positions are going to be searched.
		
		Returns
		---
			List[int] : List of all positions where the given pattern starts 
				in the given genome string.
	"""
	positions = []
	genome_len = len(genome)
	pattern_len = len(pattern)
    
	for i in range(genome_len - pattern_len + 1):
		if(genome[i:i + pattern_len] == pattern):
			positions.append(i)
	
	return positions
