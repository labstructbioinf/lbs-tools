def parse_pcoils_toolkit(fn, window='14'):
    window2columns = {'14': 3, '21': 5, '28': 7}
    f = open(fn)
    data = f.readlines()
    assignment = []
    for i in range(4, len(data)):
        line = data[i].rstrip()
        results = line.split()
        if len(results) > 0:
            assignment.append(float(results[window2columns[window]]))
    return(assignment)

def parse_marcoil(fn):
    f = open(fn)
    data = f.readlines()
    for i in range(0, len(data)):
        line = data[i].rstrip()
        if 'cc-probability in percent and best heptad phase' in line:
            probs = []
            j = i+1
            while data[j].rstrip() != '':
                line2 = data[j].rstrip()
                data2 = line2.split()
                if len(data2) == 4:
                    probs.append(float(data2[2]))
                j += 1
            return(probs)

def parse_deepcoil(fn):
    f = open(fn)
    data = f.readlines()
    probs = []
    for a in data:
        a = a.rstrip()
        aa, prob = a.split()
        probs.append(float(prob))
    return probs

def parse_pcoils(fn):
    f = open(fn)
    data = f.readlines()
    probs = []
    seq = []
    for i in range(2, len(data)):
        line = data[i].rstrip()
        results = line.split()
        probs.append(float(results[4]))
    return probs

def parse_multicoil2(fn):
    f = open(fn)
    data = f.readlines()
    preds = []
    for i in range(3, len(data)):
        line = data[i].rstrip()
        preds.append(float(line.split()[3]))
    return preds
