import os
def search_deep(n_path,func,args=[]):
	try:
		for i in os.listdir(n_path):
			if os.path.isdir(n_path+'/'+i):
				search_deep(n_path+'/'+i,func,args=[])
			else:
				for j in args:
					if j in i:
						func(n_path+'/'+i)
	except OSError:
		func(n_path+'/'+i)
		#raise Exception('Not a directory ')