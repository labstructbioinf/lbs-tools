def parse_pdb_ranges(pdb_ranges):
    
    for pdb_range in pdb_ranges.split(','):
        tmp = pdb_range.split(":")[1]

        if tmp[0] == '-':
            minus_start = True
            tmp = tmp[1:]
        else:
            minus_start = False

        if tmp.find('--')>0:
            tmp = tmp.replace('--', '-')
            minus_end = True
        else:
            minus_end = False
            
        ecod_start, ecod_end = tmp.split("-")

        if minus_start:
            ecod_start = '-' + ecod_start
            
        if minus_end:
            ecod_end = '-' + ecod_end
            
        yield ecod_start, ecod_end
        
        

