#Written in CamelCase because the code validation tool
#only accepts the function to be named like this
#1.1
def PatternCount(Text, Pattern):
    """ 
		Returns how many times the pattern string Pattern
		is present in the string Text
	"""
    count = 0
    pattern_len = len(Pattern)
    text_len = len(Text)

    for i in range(text_len - pattern_len + 1):
        
        if(Text[i : i + pattern_len] == Pattern):
            count += 1
    
    return count

#1.2
def pattern_count_kmer(text, kmer):
	
	"""
		Returns a tuple with the most frequent k-mer.
		A k-mer is a pattern or string of nucleotides 
		of lenght k
	"""

	result = ""
	count = 0

	for i in range(len(text) - kmer + 1):
		step = i + kmer
		mer = text[i:step]
		match_count = 1
		
		for j in range(step, len(text) - kmer +1):

			if ( text[j:j+kmer] == mer ):
				match_count += 1

		if( match_count > count):
			result = mer
			count = match_count
	
	return (result, count)


#1.3
def frequency_map(text, k):
	"""
		Returns a map whose key is a k-mer and its
		value is the times that mer is present in text.
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

#1.3(2)
def frequent_words(text, k):

	freq_map = frequency_map(text, k)
	most_freq = max(freq_map.values())
	result = []
	
	for key in freq_map:
		if (freq_map[key] == most_freq):
			result.append(key)
	
	return result

#1.4(1)
def reverse(pattern):
	"""
		Returns the reverse of the string pattern
	"""
	return pattern[::-1]

#1.4(2)
def complement(pattern):
	"""
		Returns the complement of the nucleotide(pattern)
	"""
	comp_map = {"T":"A", "A":"T", "C":"G", "G":"C"}
	result = ""

	for char in pattern:
		result += comp_map[char]
	
	return result

#1.4(3)
def reverse_complement(pattern):
	"""
		Returns the reverse complement of a nucleotide(pattern)
	"""
	return reverse(complement(pattern))

#1.4(4)
def pattern_match_positions(pattern, genome):
	"""
		Returns a list with the positions in genome where
		pattern starts.
	"""
	positions = []
	genome_len = len(genome)
	pattern_len = len(pattern)
    
	for i in range(genome_len - pattern_len + 1):
		if(genome[i:i + pattern_len] == pattern):
			positions.append(i)
	
	return positions
