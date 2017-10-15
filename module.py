import os
def search_deep(n_path,func,args=[],ref=1):
        if 'Trash' in n_path:
		return
	try:
                for i in os.listdir(n_path):
                        if os.path.isdir(n_path+'/'+i):
                                search_deep(n_path+'/'+i,func,args,ref)
                        else:
				if ref:
                                	for j in args:
                                        	if j==i[-len(j):]:
                                                	func(n_path+'/'+i)
				else:
					for j in args:
                                                if j in i:
                                                        func(n_path+'/'+i)
        except OSError:
                 func(n_path)
                
                #raise Exception('Not a directory ')
