#Written in CamelCase because the code validation tool
#only accepts the function to be named like this
def PatternCount(Text, Pattern):
    
    count = 0
    pattern_len = len(Pattern)
    text_len = len(Text)

    for i in range(text_len - pattern_len + 1):
        
        if(Text[i : i + pattern_len] == Pattern):
            count += 1
    
    return count


def pattern_find(text, kmer):

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
	
	return [result, count]

