def remove_overlap(ranges):
	"""
	
	Merges overlapping ranges. Rnages are given in a format: [form, to, id]
	:param ranges: ranges to be merged
	:return: merged ranges
	"""
	
	result = []
	current_start = -float("inf")
	current_stop = -float("inf")

	for start, stop, Type in sorted(ranges, key=lambda i:i[:2]):
		Type = [Type, start, stop]
				
		if start > current_stop:
			# this segment starts after the last segment stops
			# just add a new segment
			result.append( (start, stop, [Type]) )
			current_start, current_stop = start, stop
		else:
			# segments overlap, replace
			oldType = result[-1][2]

			#result[-1] = (current_start, stop, list(set(oldType+[Type])))
			#print oldType
			#print Type
			result[-1] = (current_start, stop, oldType+[Type])
			# current_start already guaranteed to be lower
			current_stop = max(current_stop, stop)

	return result
    
def calc_sov(true, pred, label='I', calc_sigma=True):
	"""
	:param true: true assignment
	:param pred: predicted assignment
	:param label: label for which SOV is calculated
	:param calc_sigma: determines whether sigma factor is calculated
	:return: SOV value for protein
	"""

	assert len(true) == len(pred)
	if true.count(label) > 0:
		Ncc = 0
		true_segments = []
		c = 0
		for match in re.finditer(r"([%s]+)" % label, true):
			true_segments.append([])
			for i in range(match.start(), match.end()):
				true_segments[c].append(i)
			c += 1
		c = 0
		pred_segments = []
		for match in re.finditer(r"([%s]+)" % label, pred):
			pred_segments.append([])
			for i in range(match.start(), match.end()):
				pred_segments[c].append(i)
			c += 1
		summ = 0
		for segment in true_segments:
			matching_segment = False
			for segment2 in pred_segments:
				minov = len((set(segment) & set(segment2)))
				maxov = len((set(segment) | set(segment2)))
				if minov > 0:
					Ncc += len(segment)
					matching_segment = True
					len_s1_adj = int(0.5 * len(segment))
					len_s2_adj = int(0.5 * len(segment2))
					sigma = min(maxov - minov, minov, len_s1_adj, len_s2_adj)
					if not calc_sigma:
						sigma = 0
					summ += float(len(segment) * (minov + sigma) / maxov)
			if not matching_segment:
				Ncc += len(segment)
		return (float(summ / Ncc))