#!/usr/bin/python
import os 
import module

def func(path):
	direc=os.getcwd()
	li=[]
	f=path.split('/')[-1][:-4]
	for i in os.listdir(direc):
		if f in i and os.path.isdir(i):
			print i,'i'
			li.append(i)
	li.sort()
	if len(li)==0:
		return
	d=li[-1]
	os.system('cp '+direc+'/'+d+'/'+f+'.chk ./')
	os.system('f16 '+f)
	os.system('lmode -b <'+f+'.alm> '+f+'.out')


path=raw_input("Enter path : ")
module.search_deep(path,func,['.alm'])

