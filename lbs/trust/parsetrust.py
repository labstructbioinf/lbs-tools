import tempfile
from Bio import SeqIO
import os
import subprocess


# usage example:

# all_reps = []
# 
# for fname in tqdm(glob.glob('./data_taa/trust/*.txt')):
#         
#     ident = os.path.splitext(os.path.basename(fname))[0]    
#     #if ident !='WP_080934149.1': continue
#     
#     res = parse_trust_result_file(fname, debug=False)
#     if len(res) > 0: 
#         for r in res:
#             all_reps.append([ident] + list(r))

def run_command(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode == 0:
        return result.stderr.strip()  # Return the output of the command
    else:
        raise ValueError
        
        
# not really used
charge = {'G':	abs(6.0 - 6.0),
'A':	abs(6.0 - 6.0),
'S':	abs(5.7 - 6.0),
'C':	abs(5.0 - 6.0),
'D':	abs(3.0 - 6.0),
'P':	abs(6.3 - 6.0),
'N':	abs(5.4 - 6.0),
'T':	abs(5.6 - 6.0),
'E':	abs(3.2 - 6.0),
'Q':	abs(5.7 - 6.0),
'H':	abs(7.6 - 6.0),
'M':	abs(5.7 - 6.0),
'I':	abs(6.0 - 6.0),
'L':	abs(6.0 - 6.0),
'K':	abs(9.7 - 6.0),
'R':	abs(10.8 - 6.0),
'F':	abs(5.5 - 6.0),
'Y':	abs(5.7 - 6.0),
'W':	abs(5.9 - 6.0),
'V':	abs(6.0 - 6.0)}


def main(commandline=False, data=None):
    if commandline:
        infile = open(sys.argv[1])
    else:
        res = ""
        infile = open(data)

    number = 1
    result = ""

    for fasta in SeqIO.parse(infile, 'fasta'):
        if commandline:
            print(str(number).ljust(13, " ") + str(fasta.seq).replace("~", "-").upper())
        else:
            result += str(number).ljust(13, " ") + str(fasta.seq).replace("~", "-").upper() + "\n"
        number += 1

    if not commandline:
        return result


    
def find_first_digit_position(string):
    for i, char in enumerate(string):
        if char.isdigit():
            return i
        
    print(string)
    raise ValueError
    

def runal2co(handle, weight=2, conservation=0, matrixtransf=0, debug=False, binary = '~/apps/al2co/al2co'):
    tr = tempfile.NamedTemporaryFile(delete=False, mode='wt')
    tr.close()
    
    # create input file
    tf = tempfile.NamedTemporaryFile(delete=False, mode='wt')
    tf.write(main(data=handle))
    tf.close()
    cmd = "%s -i %s -t %s -b 5000 -c %s -f %s -m %s" % (binary, tf.name, tr.name, conservation, weight, matrixtransf)
    
    txt_res = run_command(cmd)

    if txt_res == 'variance is zero':
        return []

    if debug:
        print(cmd)
    else:
        os.remove(tf.name)
        
    # get results
    res = open(tr.name)
    res.readline()
    res.readline()
    finalres = res.readline()
    res.close()    
    os.remove(tr.name)
    return list(finalres[find_first_digit_position(finalres):-1])



def parse_fasta(handle, scale):
    pos = None
    seqcount = 0
    slen = None

    for fasta in SeqIO.parse(handle, 'fasta'):
        seqcount += 1

        if pos is None:
            pos = [0] * len(fasta.seq)
            c = [0] * len(fasta.seq)
            slen = len(fasta.seq)

        if slen != len(fasta.seq):
            print("aln2stat: all sequences should be equal length!")
            print("first seq {}, other seq {}".format(slen, len(fasta.seq)))
            sys.exit(-1)

        for p in range(len(fasta.seq)):
            if fasta.seq[p] in scale:
                pos[p] += scale[fasta.seq[p]]
                c[p] += 1
            elif fasta.seq[p] not in ['-', '~']:
                pass
                # print("alntools: unknown character in sequence ->", fasta.sequence[p])
                # print(fasta)
                # sys.exit(-1)

    if seqcount == 0:
        print("aln2stat: there are no sequences in handle", handle)
        sys.exit(-1)

    return pos, c, seqcount



def parse_trust_result_file(filename, debug=False, MINREPLEN = 10, MAXREPLEN = 1500, MINREPEATNR = 3):

    
    """
    
        # MINREPLEN - Minimal length of repeat (aln columns)
        # MAXREPLEN - Maximal length of repeat
        # MINREPEATNR - Minimal number of repeat units

    """
    
    # detected repeats will be stored in
    res = []
    
    shortname = os.path.splitext(os.path.basename(filename))[0]
    
    sequence = False
    end_repeat = False
    
    infile = open(filename)
    l = infile.readline()

    first = True

    while l:
        
        #if debug:
        #    print(l)
        
        if 'REPEAT_TYPE' in l:
            tf = tempfile.NamedTemporaryFile(delete=False, mode='wt')
            REPEATNR = 0
            sequence = False

        elif l.find(">Repeat") != -1:
            REPEATNR += 1
            tf.write(">%s %s %s" % (shortname, REPEATNR, l[1:]))  # FASTA header
            sequence = True

        elif sequence:
            tf.write(l.upper())
            sequence = False
            end_repeat = True

        elif end_repeat:
            tf.close()

            cont = REPEATNR >= MINREPEATNR
            if cont:
                add = None
                avgcons = None
                try:
                    pos, _, seqcount = parse_fasta(tf.name, charge)

                    if len(pos) > MAXREPLEN or len(pos) < MINREPLEN:
                        wronglen = True
                    else:
                        wronglen = False
                        conservation = runal2co(tf.name, weight=0, conservation=0, matrixtransf=0, debug=debug)

                        if conservation == []:
                            raise IOError
                         
                        #print (pos)
                        #print (conservation)
                                
                        assert len(pos) == len(conservation)

                        cons = 0
                        count = 0

                        for p in range(len(pos)):
                            #gaps = c[p] / float(seqcount)

                            cons += int(conservation[p])
                            count += 1

                        if count > 0:
                            avgcons = round(cons * 1.0 / count, 1)

                except IOError:
                    avgcons = None

                if avgcons is not None and not wronglen:
                    res.append(
                        (round(avgcons, 1), REPEATNR, len(pos))
                    )

            os.remove(tf.name)

            end_repeat = False

        l = infile.readline()
        
    return res