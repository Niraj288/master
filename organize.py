#!/usr/bin/python
import os
d={}
for file in os.listdir(os.getcwd()):
	f=file.split('.')[0]
	if len(f)==0:
		continue
	if f not in d:
		d[f]=[1,[[file,f+'/'+file]]]
		#os.mkdir(f)
	else:
		d[f][0]+=1
		d[f][1].append([file,f+'/'+file])
print d
for f in d:
	a,l=d[f]
	if a>1:
		os.mkdir(f)
		for b,c in l:
			try :
				os.rename(b,c)
			except OSError:
				pass

print "Done!!"


	
