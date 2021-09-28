import pickle, zlib, os
from functools import partial

class cache(object):

	def __init__(self, func, cachedir=""):
		self.func = func
		self.verb = False
		self.cachedir = cachedir
		
	def __call__(self, *args, **kw):		
		keyraw = list(args)
		keyraw += kw.values()
				
		key = self.func.__name__+"_"+str(zlib.crc32(str.encode(str(keyraw))))
		file_name = os.path.join(self.cachedir, key)

		if os.path.isfile(file_name):
			data = pickle.load(open(file_name, 'rb'))
			if self.verb: print (f'memoized: data for `{self.func.__name__}` read from cache')
			return data 
		else:
			value = self.func(*args, **kw)
			if self.verb: print (f'memoized: data for `{self.func.__name__}` wrote to cache')
			pickle.dump(value, open(file_name, 'wb'))
			return value
