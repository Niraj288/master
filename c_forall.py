#!/usr/bin/python
import os 
import module
import sys
import multiprocessing

def threa(args1,path):
	os.chdir('/'.join(path.split('/')[:-1]))
	os.system(s)
        sys.stdout.flush()

def func(path):
	global args1
	filename=path.split('/')[-1].split('.')[0]
	if '-d' in args1:
        	os.chdir('/'.join(path.split('/')[:-1]))
                args1.remove('-d')
		for arg in args1:
                	if '-p' in arg or '-f' in arg:
                        	s=arg.replace('-f',filename)
                        	s=s.replace('-p',path)
                        	print 'Performing :',s,'...'
				multiprocessing.Process(target=threa, args=(s,path)).start()
                	else:
                	        multiprocessing.Process(target=threa, args=(arg,path)).start()
		return
	for arg in args1:
		if '-p' in arg or '-f' in arg:
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
print '-d to work in the same directory'
args1=raw_input("Enter comma separated Commands : ")
args1=args1.replace(',','","')
args1=eval('["'+args1+'"]')
module.search_deep(path,func,lis)

