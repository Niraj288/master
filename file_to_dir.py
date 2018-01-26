#!/usr/bin/python
import os 
import module

def func(path):
	global dire
	if path[:len(dire)]==dire:
		return 
	os.system("cp "+path+' '+dire)

def job(path,args,dire):
	os.chdir(path)
	os.system("mkdir "+dire)
	module.search_deep(path,func,args)

if __name__=='__main__':
	path=raw_input("Enter path : ")
	args=raw_input("Enter arguments ('.g16','.g09', ...) : ")
	args='["'+args.replace(',','","')+'"]'
	dire=path+'/'+raw_input("Enter directory name : ")+'/'
	os.chdir(path)
	os.system("mkdir "+dire)
	module.search_deep(path,func,eval(args))
