import pattern.en as en
import csv
import sys
import codecs
import re
from array import array

# Adapted from: 
# https://github.com/alevchuk/pairwise-alignment-in-python/blob/master/alignment.py
def zeros(shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval

match_award      = 20
mismatch_penalty = -5
gap_penalty      = -5 
small_dist_award = 5


# Adapted from: 
# https://github.com/alevchuk/pairwise-alignment-in-python/blob/master/alignment.py
def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    elif levenshtein(alpha,beta) <= 2:
    	return small_dist_award
    else:
        return mismatch_penalty


# Adapted from: 
# https://github.com/alevchuk/pairwise-alignment-in-python/blob/master/alignment.py
def finalize(align1, align2):
    align1.reverse()    #reverse sequence 1
    align2.reverse()    #reverse sequence 2
    
    i,j = 0,0
    
    #calcuate identity, score and aligned sequeces

    found = 0
    score = 0
    identity = 0
    for i in range(0,len(align1)):
        if align1[i] == align2[i]:                

            identity = identity + 1
            score += match_score(align1[i], align2[i])
    
        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-': 
            score += match_score(align1[i], align2[i])
            found = 0
    
        #if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':          
            score += gap_penalty
    
    identity = float(identity) / len(align1)
    
    #print align1
    #print align2
    return align1, align2


# Adapted from: 
# https://github.com/alevchuk/pairwise-alignment-in-python/blob/master/alignment.py
def global_align(seq1, seq2):
    m, n = len(seq1), len(seq2)  # length of two sequences
    
    # Generate DP table and traceback path pointer matrix
    score = zeros((m+1, n+1))      # the DP table
   
    # Calculate DP table
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[i-1], seq2[j-1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    # Traceback and compute the alignment 
    align1, align2 = [], []
    i,j = m,n # start from the bottom right cell
    while i > 0 and j > 0: 
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]

        if score_current == score_diagonal + match_score(seq1[i-1], seq2[j-1]):
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1.append(seq1[i-1])
        align2.append('-')
        i -= 1
    while j > 0:
        align1.append('-')
        align2.append(seq2[j-1])
        j -= 1

    return finalize(align1, align2)

#checks if word is a function word
def is_function(word):
	function_pos = ["CC", "DT", "EX", "IN", "MD", "PDT", "PRP", 
		"PRP$", "TO", "WDT", "WP", "WP$"]
	if (check_pos(word) in function_pos):
		return True
	return False

#checks part of speech of word
def check_pos(word):
	for word, pos in en.tag(word):
		return pos

#sorts words into content and function buckets
def sort_words(sentence):
	buckets = []
	function_bucket = []
	content_bucket = []
	for word in sentence:
		if(is_function(word) or word == "is"):
			function_bucket.append(word.lower())
		else:
			content_bucket.append(word.lower())
	buckets.append(function_bucket)
	buckets.append(content_bucket)

	return buckets

#fixes special case compound words in our input
def special_cases(resp_words_to_align):

	if "dishtowel" in resp_words_to_align:
		index = resp_words_to_align.index("dishtowel")
		resp_words_to_align.remove("dishtowel")
		resp_words_to_align.insert(index, "dish")
		resp_words_to_align.insert(index+1, "towel")

	if "dishtowels" in resp_words_to_align:
		index = resp_words_to_align.index("dishtowels")
		resp_words_to_align.remove("dishtowels")
		resp_words_to_align.insert(index, "dish")
		resp_words_to_align.insert(index+1, "towels")
	if "3" in resp_words_to_align:
		index = resp_words_to_align.index("3")
		resp_words_to_align.remove("3")
		resp_words_to_align.insert(index, "three")
	if "2" in resp_words_to_align:
		index = resp_words_to_align.index("2")
		resp_words_to_align.remove("2")
		resp_words_to_align.insert(index, "two")
	if "layed" in resp_words_to_align:
		index = resp_words_to_align.index("layed")
		resp_words_to_align.remove("layed")
		resp_words_to_align.insert(index, "laid")
	resp_words_to_align = combine_compound("goal", "post", "goalpost", resp_words_to_align)
	resp_words_to_align = combine_compound("dish", "cloth", "dishcloth", resp_words_to_align)
	resp_words_to_align = combine_compound("mail", "man", "mailman", resp_words_to_align)
	resp_words_to_align = combine_compound("sauce", "pan", "saucepan", resp_words_to_align)
	resp_words_to_align = combine_compound("rain", "coat", "raincoat", resp_words_to_align)
	return resp_words_to_align


#combines compound words from resp 
def combine_compound(part1, part2, compound, resp_words_to_align):
	if part1 in resp_words_to_align and part2 in resp_words_to_align:
		index = resp_words_to_align.index(part1)
		resp_words_to_align.remove(part1)
		resp_words_to_align.remove(part2)
		resp_words_to_align.insert(index, compound)
	return resp_words_to_align


#finds errors in sentences
def find_errors(key_words, resp_words):
	all_errors = []
	all_alignments = []
	num_errors = 0
	#distractor 
	distractors = {
		frozenset(('daughter', 'set', 'table')):["coach", "racecar", "can", "go", "very"],
		frozenset(('train', 'stops', 'station')):["sister", "i", "are", "family", "team", "was", "trained"],
		frozenset(('driver', 'started', 'engine')):["breakfast", "drank", "some", "orange", "juice", "heard"],
		frozenset(('man', 'turning', 'faucet')):["birds", "build", "nests", "trees"],
		frozenset(('nice', 'people', 'are', 'coming')):["sister", "i", "are", "family", "team"],
		frozenset(('boy', 'hurried', 'school')):["plants", "are","full", "green", "leaves"],
		frozenset(('tiny', 'baby', 'was', 'pretty')):["juice", "heard", "ticking", "clock"],
		frozenset(('walked', 'grass')):["team", "was", "trained", "coach", "racecar", "can"],
		frozenset(('puppy', 'plays', 'ball')):["feet", "quarter", "is", "worth", "twenty-five"],
		frozenset(('football', 'hit', 'goalpost')):["mother", "father", "made", "bed", "clean"],
		frozenset(('met', 'friends')):["elephants", "are", "big", "animals"],
		frozenset(('sucking', 'thumb')):["first", "day", "week", "washed", "hands"],
		frozenset(('mother', 'tied', 'string')):["wrist", "elephants", "are", "big", "animals"],
		frozenset(('oven', 'too', 'hot')):["dessert", "had", "apple", "pie"],
		frozenset(('raincoat', 'hanging', 'up')):["necks", "parents", "sister", "i", "are"],
		frozenset(('bumped', 'head')):["go", "very", "fast", "people", "wear", "shoes"]}
	
	num_sentences = len(key_words)

	for sentence in range(num_sentences):

		errors = {"add_func":0, "omit_func":0, "sub_func":0, "add_cont":0, 
			"omit_cont":0, "sub_cont":0, "morph":0, "DNH_nothing":0, "DNH_incorrect":0, "perfect":0,
			"distractor_words":0, "total_cont_words":0}	
		alignments = []
		resp_words[sentence] = special_cases(resp_words[sentence])
		key_cont_bucket = sort_words(key_words[sentence])[1]
		resp_cont_bucket = sort_words(resp_words[sentence])[1]

		key_words_to_align = key_words[sentence]
		key_words_to_check = set()
		for word in key_words[sentence]:
			key_words_to_check.add(word)
		#print key_words_to_check

		i=0
		for word in key_words_to_align:
			if word != 'his':
				key_words_to_align[i] = en.lemma(word)
			i+=1

		
		resp_words_to_align = resp_words[sentence]
		resp_words_to_check = set()
		for word in resp_words[sentence]:
			resp_words_to_check.add(word)
		i=0
		for word in resp_words_to_align:
			suggestions = en.suggest(word)

			for suggestion in suggestions:

				if suggestion[0] in key_words_to_align:

					resp_words_to_align[i] = suggestion[0]

			if resp_words_to_align[i] != 'his':
				resp_words_to_align[i] = en.lemma(resp_words_to_align[i])

			i+=1

		alignment = global_align(key_words_to_align, resp_words_to_align)
		key_align = alignment[0]
		resp_align = alignment[1]
		print key_align
		print resp_align
		i=0
		#print key_words_to_check
		for aligned_word in key_align:
			for check_word in key_words_to_check:
				if en.lemma(check_word) == aligned_word and check_word != aligned_word:
					key_align[i] = check_word
			i+=1
		i=0
		for aligned_word in resp_align:
			for check_word in resp_words_to_check:
				#print aligned_word
				#print check_word
				if en.lemma(check_word) == aligned_word and check_word != aligned_word:
					resp_align[i] = check_word
			i+=1

		print key_align
		print resp_align
		alignments.append(alignment)

		#if any words in resp match target
		if any(word in key_cont_bucket for word in resp_cont_bucket): 
			
			if key_align == resp_align:
				errors["perfect"] += 1
			else:
				for i in range(len(key_align)):
					if key_align[i] != resp_align[i]:
						if(en.lemma(key_align[i]) == en.lemma(resp_align[i])):
							errors["morph"] += 1
						elif key_align[i] == '-':
							if is_function(resp_align[i]):
								errors["add_func"] += 1
							else:
								errors["add_cont"] += 1
						elif resp_align[i] == '-':
							if is_function(key_align[i]):
								errors["omit_func"] += 1
							else:
								errors["omit_cont"] += 1
						elif is_function(key_align[i]) and is_function(resp_align[i]):
							errors["sub_func"] += 1
						elif (not is_function(key_align[i])) and (not is_function(resp_align[i])):
							errors["sub_cont"] += 1
						elif is_function(key_align[i]) and (not is_function(resp_align[i])):
							errors["omit_func"] += 1
							errors["add_cont"] += 1
						elif (not is_function(key_align[i])) and is_function(resp_align[i]):
							errors["omit_cont"] += 1
							errors["add_func"] += 1

		#did not hear target sentence
		else:
			if len(resp_cont_bucket) == 0:
				errors["DNH_nothing"] = 1
			else:
				#look for distractor words in response
				kcb = set()
				for elem in key_cont_bucket:
					kcb.add(elem)
				key_with_distractor = frozenset((kcb))
				if key_with_distractor in distractors.keys():
					num=0
					for word in resp_cont_bucket:
						if word in distractors[key_with_distractor]:
							#print word
							num+=1

					errors["distractor_words"] +=num
					errors["total_cont_words"] +=len(resp_cont_bucket)


				errors["DNH_incorrect"] = 1

		all_alignments.append(alignment)
		all_errors.append(errors)
	ret = []
	ret.append(all_errors)
	ret.append(all_alignments)
	return ret


#removes formatting of words 
def parse_words(col_num, is_key, roots_only, file_name, stim_type_col):
	with open(file_name, 'rb') as f:
		reader = csv.reader(f, delimiter=',', quotechar = '|')
		all_words = []

		for line in reader:

			#check format
			try:
				sentence = line[col_num]
			except:
				col_num -= 3
				sentence = line[col_num]
			#check format
			if is_key and sentence[:2] == "A1":
				sentence = line[col_num + 1]
			sentence = replace_special_chars(sentence)  #replaces special characters

			sentence = sentence.split()  #splits sentence into words

			words = []
			for word in sentence:

				if roots_only:
					words.append(en.lemma(word.lower()))
				else:
					if (word.lower()[0]) != 'x':  #ignores words that are only x's
						words.append(word.lower())
			all_words.append(words)

		return all_words

#removes special from source file 
#change to regex??
def replace_special_chars(sentence):
	chars_to_replace = ['{TAB\}','{', '}', '?', ':', '(', ',', '"', ';' 
		'`', '-', '\'', 'LEFTBRACKET', 'LEFTARROW', '\INSERT', 'UPARROW',  
		 'CAPSLOCK', 'ESCAPE', 'ENTER', 'SHIFT', 'CONTROL', 'F1',
		'RIGHTARROW', 'RIGHTBRACKET']
	backmatch_pattern = re.compile(".{\BACKSPACE}")
	while backmatch_pattern.match(sentence):
		sentence = backmatch_pattern.sub('', sentence)

	for char in chars_to_replace:
		sentence = sentence.replace(char, '')

	sentence = sentence.replace('SPACE', ' ')
	sentence = sentence.replace('.', ' ')


	return sentence

#formats sentences (adds spaces in between words for output file)
def format_sentences(words):
	all_formatted_sentences = []
	for sentence in words:
		formatted_sentence = ""
		for word in sentence:
			formatted_sentence += word + " "
		all_formatted_sentences.append(formatted_sentence)
	return all_formatted_sentences

def main():

	subject_id_col = 1
	key_col = 12
	resp_col = 16
	stim_type_col = 27
	condition_col = 13
	
	#take in subject ID numbers from command line to run program
	for arg in sys.argv[1:]:
		file_name = 'sub%s_UTSPIN.csv'%arg
		key_words = parse_words(key_col, 1, 0, file_name, stim_type_col)
		resp_words = parse_words(resp_col, 0, 0, file_name, stim_type_col)
		
		find_errs = find_errors(key_words, resp_words)
		all_errors = find_errs[0]
		all_alignments = find_errs[1]

		all_key_sentences = format_sentences(key_words)
		all_resp_sentences = format_sentences(resp_words)


		with open('output_%s.csv' %arg, 'w') as csvfile:
			file_writer = csv.writer(csvfile, delimiter=',', quotechar=',', quoting=csv.QUOTE_MINIMAL)
			with open(file_name, "rb") as f:
				#reads in subject-level information for output file 
				subject_ids = []
				conditions = []
				reader = csv.reader(f,delimiter=',')
				for line in reader:

					subject_ids.append(line[subject_id_col])
					conditions.append(line[condition_col])


				#subject_ids.pop(0)
				#conditions.pop(0)


				temp = ["subject_id", "condition", "key_sentence", "resp_sentence"]
				for key,value in all_errors[0].items():
					temp.append(key)
				temp.append("total_errors")
				temp.append("total_morph_errors")
				temp.append("total_cont_errors")
				temp.append("total num words")
				file_writer.writerow(temp)

				for i in range(len(all_errors)):

					subject_id = subject_ids[i]
					line_to_write = []
					line_to_write.append(subject_ids[i])
					line_to_write.append(conditions[i])

					key_align_string = ""
					for word in all_alignments[i][0]:
						if word != '-':
							key_align_string += word
						else:
							key_align_string += '_'
						key_align_string += " "
					resp_align_string = ""
					for word in all_alignments[i][1]:
						if word != '-':
							resp_align_string += word
						else:
							resp_align_string += '_'
						resp_align_string += " "
					line_to_write.append(key_align_string)
					line_to_write.append(resp_align_string)

					total_errors = 0
					total_morph_errors = 0
					total_cont_errors = 0

					for key, value in all_errors[i].items():
						line_to_write.append(str(value))
						if (key == "sub_func" or key == "omit_func" or key== "add_func" or key== "morph"):
							total_morph_errors += value
						if (key == "sub_cont" or key == "omit_cont" or key== "add_cont"):
							total_cont_errors += value


						if not (key == "DNH_nothing" or key == "DNH_incorrect" or key == "perfect" or key == "distractor_words" or key == "total_cont_words"):
							total_errors += value
						else:
							if (key == "DNH_nothing" or key == "DNH_incorrect" or key == "distractor_words") and value > 0:
								total_errors = -1

					line_to_write.append(total_errors)
					line_to_write.append(total_morph_errors)
					line_to_write.append(total_cont_errors)
					i=0
					for token in key_align_string.split():
						if token != '_':
							i+=1
					line_to_write.append(i)
					if line_to_write[0] != 'Subject':
						file_writer.writerow(line_to_write)

# From:
# https://github.com/doukremt/distance/blob/master/distance/_levenshtein.py

def levenshtein(seq1, seq2, normalized=False, max_dist=-1):
	"""Compute the absolute Levenshtein distance between the two sequences
	`seq1` and `seq2`.
	
	The Levenshtein distance is the minimum number of edit operations necessary
	for transforming one sequence into the other. The edit operations allowed are:
	
		* deletion:     ABC -> BC, AC, AB
		* insertion:    ABC -> ABCD, EABC, AEBC..
		* substitution: ABC -> ABE, ADC, FBC..
	
	The `max_dist` parameter controls at which moment we should stop computing the
	distance between the provided sequences. If it is a negative integer, the
	distance will be computed until the sequences are exhausted; otherwise, the
	computation will stop at the moment the calculated distance is higher than
	`max_dist`, and then return -1. For example:
	
		>>> levenshtein("abc", "abcd", max_dist=1)  # dist = 1
		1
		>>> levenshtein("abc", "abcde", max_dist=1) # dist = 2
		-1
	
	This can be a time saver if you're not interested in the exact distance, but
	only need to check if the distance between the given sequences is below a
	given threshold.
	
	The `normalized` parameter is here for backward compatibility; providing
	it will result in a call to `nlevenshtein`, which should be used directly
	instead. 
	"""
	if normalized:
		return nlevenshtein(seq1, seq2, method=1)
		
	if seq1 == seq2:
		return 0
	
	len1, len2 = len(seq1), len(seq2)
	if max_dist >= 0 and abs(len1 - len2) > max_dist:
		return -1
	if len1 == 0:
		return len2
	if len2 == 0:
		return len1
	if len1 < len2:
		len1, len2 = len2, len1
		seq1, seq2 = seq2, seq1
	
	column = array('L', range(len2 + 1))
	
	for x in range(1, len1 + 1):
		column[0] = x
		last = x - 1
		for y in range(1, len2 + 1):
			old = column[y]
			cost = int(seq1[x - 1] != seq2[y - 1])
			column[y] = min(column[y] + 1, column[y - 1] + 1, last + cost)
			last = old
		if max_dist >= 0 and min(column) > max_dist:
			return -1
	
	if max_dist >= 0 and column[len2] > max_dist:
		# stay consistent, even if we have the exact distance
		return -1
	return column[len2]

# From:
# https://github.com/doukremt/distance/blob/master/distance/_levenshtein.py
def nlevenshtein(seq1, seq2, method=1):
	"""Compute the normalized Levenshtein distance between `seq1` and `seq2`.
	
	Two normalization methods are provided. For both of them, the normalized
	distance will be a float between 0 and 1, where 0 means equal and 1
	completely different. The computation obeys the following patterns:
	
		0.0                       if seq1 == seq2
		1.0                       if len(seq1) == 0 or len(seq2) == 0
		edit distance / factor    otherwise
	
	The `method` parameter specifies which normalization factor should be used.
	It can have the value 1 or 2, which correspond to the following:
	
		1: the length of the shortest alignment between the sequences
		   (that is, the length of the longest sequence)
		2: the length of the longest alignment between the sequences
	
	Which normalization factor should be chosen is a matter of taste. The first
	one is cheap to compute. The second one is more costly, but it accounts
	better than the first one for parallelisms of symbols between the sequences.
		
	For the rationale behind the use of the second method, see:
	Heeringa, "Measuring Dialect Pronunciation Differences using Levenshtein
	Distance", 2004, p. 130 sq, which is available online at:
	http://www.let.rug.nl/~heeringa/dialectology/thesis/thesis.pdf
	"""
	
	if seq1 == seq2:
		return 0.0
	len1, len2 = len(seq1), len(seq2)
	if len1 == 0 or len2 == 0:
		return 1.0
	if len1 < len2: # minimize the arrays size
		len1, len2 = len2, len1
		seq1, seq2 = seq2, seq1
	
	if method == 1:
		return levenshtein(seq1, seq2) / float(len1)
	if method != 2:
		raise ValueError("expected either 1 or 2 for `method` parameter")
	
	column = array('L', range(len2 + 1))
	length = array('L', range(len2 + 1))
	
	for x in range(1, len1 + 1):
	
		column[0] = length[0] = x
		last = llast = x - 1
		
		for y in range(1, len2 + 1):
		
			# dist
			old = column[y]
			ic = column[y - 1] + 1
			dc = column[y] + 1
			rc = last + (seq1[x - 1] != seq2[y - 1])
			column[y] = min(ic, dc, rc)
			last = old
			
			# length
			lold = length[y]
			lic = length[y - 1] + 1 if ic == column[y] else 0
			ldc = length[y] + 1 if dc == column[y] else 0
			lrc = llast + 1 if rc == column[y] else 0
			length[y] = max(ldc, lic, lrc)
			llast = lold
	
	return column[y] / float(length[y])
print en.suggest('fotball')
print levenshtein('hold', 'full')
print levenshtein('hold', 'of')

main()
