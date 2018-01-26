#!/usr/bin/python
import os
import module
import sys

def search_deep(n_path):
       	if 'Trash' in n_path:
		return
	try:	
		global li
		#print n_path
                for i in os.listdir(n_path):
                        if os.path.isdir(n_path+'/'+i) and ('.g16_' in i or '.g09_' in i):
				li.append(n_path+'/'+i)
                	search_deep(n_path+'/'+i)
	except OSError:
		pass


li=[]
path=raw_input("Enter path : ")
search_deep(path)

d={}
for i in li:
	f='_'.join(i.split('_')[:-1])
	if f in d:
		d[f][0]+=1
		d[f].append(i)
	else:
		d[f]=[1,i]
r_lis=[]
if '-all' in sys.argv:
	ref=1
else:
	ref=0

for i in d:
	if  d[i][0]>1:
		li=d[i][1:]
		li=sorted(li, key=lambda j_id: int(j_id.split('_')[-1]))
		print '*****************'
		for j in li:
			print j
		print '*****************'
		print ''
		if ref or raw_input('clean '+i+'(y/n) : ')=='y':
			r_lis+=li[:-1]
if len(r_lis)>0:
	os.system('mkdir Trash')
for item in  r_lis:
	os.system('mv '+item+' '+'Trash')

