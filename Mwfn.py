import os
import sys
import threading
import multiprocessing
if len(sys.argv)>1:
	key = sys.argv[1]
else:
	key = '.wfx'
def gaussian(path,string):
	filename=path[:-4]
	os.system('mkdir '+filename)
	os.chdir(filename)
	os.system("cp ../"+path+' ./')
	os.system(" printf '"+string+"' | Multiwfn "+path)
	
processes=[]
for i in os.listdir('.'):
	if i[-4:]==key:
		p=multiprocessing.Process(target=gaussian, args=(i,'5\n1\n3\n2\n'))
                processes.append(p)
		q=multiprocessing.Process(target=gaussian, args=(i,'5\n12\n1\n2\n'))
                processes.append(q)
		#r=multiprocessing.Process(target=gaussian, args=(i,'12\n0\n'))
                #processes.append(r)

for p in processes:
	p.start()


print 'All calculation are done!!!'


