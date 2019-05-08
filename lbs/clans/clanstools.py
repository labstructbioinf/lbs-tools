import clansparser, sys
import string


def valid_file_name(filename):
	valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
	return ''.join(c for c in filename if c in valid_chars)
    
def groups2files(clans_object):
	for c in clans_object.clusters.values():
		filename = valid_file_name(c.name)
		print(filename)

		f = open(filename + ".fas", 'w')
		
		for s in c.seqs:	
			f.write(">%s|%s\n%s\n" % (filename, s.title, s.sequence))
			
		f.close()
		
if __name__ == '__main__':
	c = clansparser.clansrun(open(sys.argv[1]), onlywarn=True, readHSP=False)
	groups2files(c)