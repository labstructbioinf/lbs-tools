def remove_overlap(ranges):
	"""
	
	Merges overlapping ranges. Rnages are given in a format:
	[form, to, id]
	
	"""
    result = []
    current_start = -1
    current_stop = -1 

    for start, stop, Type in sorted(ranges, key=lambda i:i[:2]):
        if start > current_stop:
            # this segment starts after the last segment stops
            # just add a new segment
            result.append( (start, stop, [Type]) )
            current_start, current_stop = start, stop
        else:
            # segments overlap, replace
            oldType = result[-1][2]

            result[-1] = (current_start, stop, list(set(oldType+[Type])))
            # current_start already guaranteed to be lower
            current_stop = max(current_stop, stop)

    return result