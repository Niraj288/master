#!/usr/bin/python
import os 
import module
import sys

def func(path):
	global args
	filename=path.replace('/','_')[1:]
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
lis=raw_input("Enter keyword for filename : ")
lis=lis.replace(',','","')
lis=eval('["'+lis+'"]')
print '-f is for filename'
print '-p is for path'
args=raw_input("Enter comma separated Commands : ")
args=args.replace(',','","')
args=eval('["'+args+'"]')
module.search_deep(path,func,lis)

