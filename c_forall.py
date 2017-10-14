#!/usr/bin/python
import os 
import module
import sys

def func(path):
	global args
	filename=path.strip().split('/')[-1].split('.')[0]
	for arg in args:
		if '-p' in arg:
			s=arg.replace('-f',filename)
			s=s.replace('-p',path)
			print 'Performing :',s,'...'
			os.system(s)
			sys.stdout.flush()
		else:
			os.system(arg)
			sys.stdout.flush()

path=raw_input("Enter path : ")
lis='['+raw_input("Enter keyword for filename : ")+']'
print '-f is for filename'
print '-p is for path'
args=eval('['+raw_input("Enter comma separated Commands : ")+']')
module.search_deep(path,func,eval(lis))
