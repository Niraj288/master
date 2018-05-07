#!/usr/bin/python
import os
os.system('squeue > test.txt')
file=open('test.txt','r')
data=file.readlines()
file.close()
d={}
for i in data:
	lis=i.split()
	name=lis[1]
	q=lis[4]
	if name in d:
		if q=='PD':
			d[name].append(q)
	else:
		d[name]=[]

for i in d:
	if i=='PARTITION':
		continue
	if len(d[i])==0:
		print i	
	
os.system('rm -rf test.txt')



